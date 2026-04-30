import argparse
import re
import requests
import pandas as pd
import time


SOURCE_PRESETS = {
    "wikidata": {
        "input": "data/wikidata/created/compounds_with_chebi_ids_homo_sapiens.tsv",
        "output": "data/wikidata/created/compounds_with_chebi_ids_homo_sapiens_updatedchebis.tsv",
        "smiles_columns": ["canonicalSmiles", "isomericSmiles"],
        "chebi_column": "chebi_id",
    },
    "hmdb": {
        "input": "data/hmdb_metabolites_extract_quantified_detected.tsv",
        "output": "data/hmdb_metabolites_extract_quantified_detected_updatedchebis.tsv",
        "smiles_columns": ["smiles"],
        "chebi_column": "chebi_id",
    },
}


def normalize_chebi_id(raw_value):
    if pd.isna(raw_value):
        return None
    value = str(raw_value).strip()
    if not value:
        return None

    normalized = []
    for part in re.split(r"[|;]", value):
        part = part.strip()
        if part.startswith("http://purl.obolibrary.org/obo/CHEBI_"):
            normalized.append(part.rsplit("/", 1)[-1])
        elif part.startswith("CHEBI:"):
            normalized.append(part.replace(":", "_", 1))
        elif part.startswith("CHEBI_"):
            normalized.append(part)
        elif part.isdigit():
            normalized.append(f"CHEBI_{part}")
        else:
            normalized.append(part)

    return "|".join(normalized)

def convert_smiles_to_chebi(smiles_string):

    chebi_ids = []
    was_resolved_directly = False
    parents_found = False

    # Get details from ChEBI lookup to check for a direct match to a ChEBI ID.
    response = requests.post(
        "https://chebifier.hastingslab.org/api/details",
        json={
            "type": "type",
            "smiles": smiles_string,
            "selectedModels": {
                "ChEBI Lookup": True,
            },
        },
    )

    lookup_infotext = response.json().get("models", {}).get("ChEBI Lookup", {}).get("highlights", [])

    # If the lookup highlights contain a ChEBI ID, use that directly.
    if lookup_infotext and "CHEBI:" in lookup_infotext[0][1]:
        chebi_id = lookup_infotext[0][1].split("CHEBI:")[1].split()[0].rstrip(".")
        chebi_ids.append(f"CHEBI:{chebi_id}")
        was_resolved_directly = True
        # print(f"Found ChEBI ID from lookup: CHEBI:{chebi_id} for SMILES {smiles_string}")
    else:
        # print(f"No direct ChEBI ID found from lookup for SMILES {smiles_string}, attempting classification...")
        response = requests.post(
            "https://chebifier.hastingslab.org/api/classify",
            json={
                "smiles": smiles_string,
                "ontology": False,
                "selectedModels": {
                    "ELECTRA (ChEBI50-3STAR)": True,
                },
            },
        )

        direct_parents = response.json().get("direct_parents")
        if direct_parents:
            for parent_list in direct_parents:
                if parent_list is not None:
                    parent_ids = [f"CHEBI:{parent[0]}" for parent in parent_list]
                    chebi_ids.extend(parent_ids)
                    parents_found = True
                    # print(f"Found direct parent ChEBI IDs from classification for SMILES {smiles_string}: {parent_ids}")
                else:
                    # print(f"No parents found in one of the classification results for SMILES {smiles_string}")
                    # print(f"Classification response content: {response.content}")
                    pass

            if chebi_ids:
                # print(f"Found {len(chebi_ids)} ChEBI IDs from classification for SMILES {smiles_string}")
                pass

    if not chebi_ids:
        # print(f"No ChEBI IDs found for SMILES {smiles_string} after lookup and classification.")
        pass

    return chebi_ids, was_resolved_directly, parents_found


def pick_smiles_candidates(raw_value):
    if pd.isna(raw_value):
        return []

    value = str(raw_value).strip()
    if not value:
        return []

    # Some rows contain multiple alternative SMILES joined by '|'.
    # Return all non-empty entries in order so we can try them one by one.
    candidates = []
    for part in value.split("|"):
        part = part.strip()
        if part:
            candidates.append(part)

    return candidates

def _choose_smiles_columns(df, requested_smiles_columns=None):
    """Pick SMILES columns that exist in the file.

    If `requested_smiles_columns` is provided, only those are considered and
    must exist. Otherwise, auto-detect common schemas.
    """

    if requested_smiles_columns:
        missing = [col for col in requested_smiles_columns if col not in df.columns]
        if missing:
            raise KeyError(f"Requested SMILES columns not found: {missing}")
        return list(requested_smiles_columns)

    auto_candidates = ["canonicalSmiles", "isomericSmiles", "smiles"]
    chosen = [col for col in auto_candidates if col in df.columns]
    if not chosen:
        raise KeyError(
            "Could not find a SMILES column. Expected at least one of "
            "['canonicalSmiles', 'isomericSmiles', 'smiles']."
        )
    return chosen


def find_missing_chebis(compounds_file, output_file_path=None, smiles_columns=None, chebi_column="chebi_id"):
    df = pd.read_csv(compounds_file, sep="\t")

    if chebi_column not in df.columns:
        raise KeyError(f"Missing required column '{chebi_column}' in {compounds_file}")

    selected_smiles_columns = _choose_smiles_columns(df, smiles_columns)
    print(f"Using SMILES columns: {selected_smiles_columns}")

    if "chebi_source" not in df.columns:
        df["chebi_source"] = None

    # Standardize existing IDs to CHEBI_##### format.
    df[chebi_column] = df[chebi_column].apply(normalize_chebi_id)

    # Treat both NaN and empty strings as missing ChEBI IDs.
    chebi_existing_mask = df[chebi_column].notna() & (df[chebi_column].astype(str).str.strip() != "")
    df.loc[chebi_existing_mask, "chebi_source"] = "already_existed"

    missing_indices = df.index[~chebi_existing_mask]
    print(f"Rows with missing ChEBI IDs: {len(missing_indices)}")

    for i, idx in enumerate(missing_indices, start=1):
        smiles_candidates = []
        for smiles_column in selected_smiles_columns:
            candidates = pick_smiles_candidates(df.at[idx, smiles_column])
            for candidate in candidates:
                if candidate not in smiles_candidates:
                    smiles_candidates.append(candidate)

        if not smiles_candidates:
            df.at[idx, "chebi_source"] = "unresolved"
            continue

        chebi_ids = []
        was_resolved_directly = False
        parents_found = False
        first_parent_ids = None

        # Prefer a direct match from any available SMILES candidate.
        # If none is found, fall back to the first parent-based match.
        for smiles_to_query in smiles_candidates:
            candidate_ids, candidate_direct, candidate_parents = convert_smiles_to_chebi(smiles_to_query)

            if candidate_direct and candidate_ids:
                chebi_ids = candidate_ids
                was_resolved_directly = True
                parents_found = False
                break

            if candidate_parents and candidate_ids and first_parent_ids is None:
                first_parent_ids = candidate_ids

        if not was_resolved_directly and first_parent_ids:
            chebi_ids = first_parent_ids
            parents_found = True

        if chebi_ids:
            df.at[idx, chebi_column] = normalize_chebi_id("|".join(chebi_ids))
            if was_resolved_directly:
                df.at[idx, "chebi_source"] = "found_directly"
            elif parents_found:
                df.at[idx, "chebi_source"] = "parents_compounds"
            else:
                df.at[idx, "chebi_source"] = "unresolved"
        else:
            df.at[idx, "chebi_source"] = "unresolved"

        if i % 50 == 0:
            print(f"Processed {i}/{len(missing_indices)} missing rows")

    save_path = output_file_path if output_file_path else compounds_file
    df.to_csv(save_path, sep="\t", index=False)
    print(f"Saved updated compounds file to {save_path}")

    direct_count = int((df["chebi_source"] == "found_directly").sum())
    parents_count = int((df["chebi_source"] == "parents_compounds").sum())
    still_missing_mask = df[chebi_column].isna() | (df[chebi_column].astype(str).str.strip() == "")
    unresolved_count = int(still_missing_mask.sum())

    print(f"Found directly: {direct_count}")
    print(f"Found using parents: {parents_count}")
    print(f"Still unresolved (no chebi_id): {unresolved_count}")

    return df

if __name__ == "__main__":
    start_time = time.time()

    # Quick source switch when editing the script directly.
    # Choices are "hmdb" and "wikidata".
    source = "wikidata"

    parser = argparse.ArgumentParser(
        description="Fill missing ChEBI IDs using SMILES for Wikidata/HMDB-style TSV files.",
    )
    parser.add_argument(
        "--source",
        choices=("wikidata", "hmdb"),
        default=source,
        help="Use built-in defaults for selected source (default: wikidata).",
    )
    parser.add_argument(
        "compounds_file",
        nargs="?",
        default=None,
        help="Input TSV file with at least a chebi_id column and one SMILES column.",
    )
    parser.add_argument(
        "output_file",
        nargs="?",
        default=None,
        help="Output TSV file path.",
    )
    parser.add_argument(
        "--smiles-columns",
        nargs="+",
        default=None,
        help="Optional explicit SMILES columns (e.g., --smiles-columns smiles or canonicalSmiles isomericSmiles).",
    )
    parser.add_argument(
        "--chebi-column",
        default=None,
        help="Name of the ChEBI ID column (default: chebi_id).",
    )

    args = parser.parse_args()

    preset = SOURCE_PRESETS[args.source]
    compounds_file = args.compounds_file if args.compounds_file else preset["input"]
    output_file = args.output_file if args.output_file else preset["output"]
    smiles_columns = args.smiles_columns if args.smiles_columns else preset["smiles_columns"]
    chebi_column = args.chebi_column if args.chebi_column else preset["chebi_column"]

    print(f"Source preset: {args.source}")
    print(f"Input file: {compounds_file}")
    print(f"Output file: {output_file}")

    find_missing_chebis(
        compounds_file,
        output_file,
        smiles_columns=smiles_columns,
        chebi_column=chebi_column,
    )

    elapsed_seconds = time.time() - start_time
    elapsed_minutes = elapsed_seconds / 60
    print(f"Total runtime: {elapsed_seconds:.2f} seconds ({elapsed_minutes:.2f} minutes)")


