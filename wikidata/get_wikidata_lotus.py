import pandas as pd
from pathlib import Path


def normalize_chebi_id(raw_value):
    if pd.isna(raw_value):
        return None
    value = str(raw_value).strip()
    if not value:
        return None

    values = []
    for part in value.split("|"):
        part = part.strip()
        if part.startswith("http://purl.obolibrary.org/obo/CHEBI_"):
            values.append(part.rsplit("/", 1)[-1])
        elif part.startswith("CHEBI:"):
            values.append(part.replace(":", "_", 1))
        else:
            values.append(part)

    return "|".join(values)


def connect_lotus_csv_to_chebi_ids(
    lotus_csv_path,
    output_tsv_path,
    inchikeys_csv="data/removed_leaf_classes_with_inchikeys.csv",
    smiles_csv="data/removed_leaf_classes_with_smiles.csv",
):
    """Match ChEBI IDs for a taxon-specific LOTUS explorer export.

    Unlike the old Zenodo Wikidata dump (see old_wikidata_scripts/), a LOTUS
    explorer CSV (e.g. data/lotus_homo_sapiens.csv) is already filtered to one
    taxon and already has compound/taxon metadata joined, so no separate
    taxon-tagging or filtering step is needed here.

    Matching is tried in order: InChIKey (compound_inchikey, always populated and exact),
    then connectivity SMILES (compound_smiles_conn), then isomeric SMILES (compound_smiles_iso).
    Anything still unmatched is left for find_missing_chebis (Chebifier API fallback).
    """

    lotus = pd.read_csv(lotus_csv_path)

    inchikey_lookup = pd.read_csv(inchikeys_csv)
    inchikey_to_chebi = dict(zip(inchikey_lookup["InChIKey"], inchikey_lookup["IRI"]))

    smiles_lookup = pd.read_csv(smiles_csv)
    smiles_lookup["SMILES"] = smiles_lookup["SMILES"].astype(str).str.strip('"').str.strip()
    smiles_to_chebi = dict(zip(smiles_lookup["SMILES"], smiles_lookup["IRI"]))

    lotus["chebi_id"] = None
    lotus["chebi_source"] = None

    for index, row in lotus.iterrows():
        inchikey = row.get("compound_inchikey")
        smiles_conn = row.get("compound_smiles_conn")
        smiles_iso = row.get("compound_smiles_iso")

        if pd.notna(inchikey) and inchikey in inchikey_to_chebi:
            lotus.at[index, "chebi_id"] = normalize_chebi_id(inchikey_to_chebi[inchikey])
            lotus.at[index, "chebi_source"] = "found_via_inchikey"
        elif pd.notna(smiles_conn) and smiles_conn in smiles_to_chebi:
            lotus.at[index, "chebi_id"] = normalize_chebi_id(smiles_to_chebi[smiles_conn])
            lotus.at[index, "chebi_source"] = "found_via_smiles_conn"
        elif pd.notna(smiles_iso) and smiles_iso in smiles_to_chebi:
            lotus.at[index, "chebi_id"] = normalize_chebi_id(smiles_to_chebi[smiles_iso])
            lotus.at[index, "chebi_source"] = "found_via_smiles_iso"

        if index % 1000 == 0:
            print(f"Processed {index} / {len(lotus)} rows")

    Path(output_tsv_path).parent.mkdir(parents=True, exist_ok=True)
    lotus.to_csv(output_tsv_path, index=False, sep="\t")

    matched_count = lotus["chebi_id"].notna().sum()
    print(f"Saved {len(lotus)} compounds with matched ChEBI IDs to {output_tsv_path}")
    print(f"Matched ChEBI IDs: {matched_count} out of {len(lotus)} compounds ({(matched_count/len(lotus))*100:.2f}%)")


if __name__ == "__main__":
    connect_lotus_csv_to_chebi_ids(
        "data/lotus_homo_sapiens.csv",
        "data/wikidata/created/lotus_homo_sapiens_with_chebi_ids.tsv",
    )
