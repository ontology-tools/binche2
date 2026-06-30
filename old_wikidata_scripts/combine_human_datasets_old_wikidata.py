"""Combine HMDB and Wikidata compound datasets into one harmonized TSV.

Output columns:
- chebi_id
- id
- smiles
- chebi_source
- source

Rules implemented:
- HMDB ChEBI IDs are normalized to CHEBI_<id> format.
- Wikidata prefers canonicalSmiles; falls back to isomericSmiles.
- Multiple ChEBI IDs are kept and written as ';'-separated values.
- Rows are dropped only when ChEBI ID(s) are missing.
- id is HMDB id, Wikidata id suffix, or both joined by '|'.
- source is 'hmdb', 'wikidata', or 'both'.
- At most one id per source is kept in merged rows.
"""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path


HMDB_INPUT = "data/hmdb_metabolites_extract_quantified_detected_updatedchebis.tsv"
WIKIDATA_INPUT = "data/wikidata/created/compounds_with_chebi_ids_homo_sapiens_updatedchebis.tsv"
OUTPUT_FILE = "data/combined_hmdb_wikidata.tsv"


def _clean(value: str | None) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _normalize_chebi_token(token: str) -> str:
    token = _clean(token)
    if not token:
        return ""

    if token.startswith("http://purl.obolibrary.org/obo/CHEBI_"):
        token = token.rsplit("/", 1)[-1]

    if token.startswith("CHEBI:"):
        return token.replace(":", "_", 1)

    if token.startswith("CHEBI_"):
        suffix = token.split("CHEBI_", 1)[1]
        if suffix.isdigit():
            return f"CHEBI_{int(suffix)}"
        return token

    # Handles HMDB-style values like "50599" or "50599.0"
    if re.fullmatch(r"\d+(?:\.0+)?", token):
        return f"CHEBI_{int(float(token))}"

    return token


def normalize_chebi_ids(raw_value: str | None) -> str:
    """Normalize ChEBI IDs and return ';'-joined unique values."""

    value = _clean(raw_value)
    if not value:
        return ""

    tokens = re.split(r"[|;]", value)
    normalized = []
    seen = set()
    for token in tokens:
        norm = _normalize_chebi_token(token)
        if norm and norm not in seen:
            normalized.append(norm)
            seen.add(norm)

    return ";".join(normalized)


def pick_preferred_smiles(row: dict[str, str]) -> str:
    """Prefer canonicalSmiles for Wikidata rows, then isomericSmiles."""

    canonical = _clean(row.get("canonicalSmiles"))
    isomeric = _clean(row.get("isomericSmiles"))

    if canonical:
        return canonical.split("|")[0].strip()
    if isomeric:
        return isomeric.split("|")[0].strip()
    return ""


def wikidata_id_suffix(wikidata_iri: str | None) -> str:
    iri = _clean(wikidata_iri)
    if not iri:
        return ""
    return iri.rsplit("/", 1)[-1]


def _merge_chebi_source(existing: str, new_value: str) -> str:
    parts = []
    seen = set()
    for text in [existing, new_value]:
        for token in _clean(text).split("|"):
            token = token.strip()
            if token and token not in seen:
                parts.append(token)
                seen.add(token)
    return "|".join(parts)


def _record_key(
    smiles: str,
    chebi_id: str,
    source: str,
    compound_id: str,
    row_index: int,
) -> tuple[str, str]:
    # Canonicalize ChEBI ordering for stable cross-source matching.
    chebis = [c.strip() for c in chebi_id.split(";") if c.strip()]
    chebi_key = ";".join(sorted(chebis))

    smiles_key = smiles.strip()
    if smiles_key:
        return (smiles_key, chebi_key)

    compound_key = compound_id.strip()
    if compound_key:
        # Avoid collapsing all empty-SMILES rows into one bucket.
        return (f"{source}:{compound_key}", chebi_key)

    # Last-resort stable key when both SMILES and id are missing.
    return (f"{source}:row{row_index}", chebi_key)


def combine_datasets(
    hmdb_path: str = HMDB_INPUT,
    wikidata_path: str = WIKIDATA_INPUT,
    output_path: str = OUTPUT_FILE,
) -> None:
    merged: dict[tuple[str, str], dict[str, str]] = {}
    drop_counts = {
        "hmdb": Counter(),
        "wikidata": Counter(),
    }
    kept_missing_smiles = {
        "hmdb": 0,
        "wikidata": 0,
    }
    total_counts = {
        "hmdb": 0,
        "wikidata": 0,
    }

    # HMDB rows
    with open(hmdb_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_index, row in enumerate(reader, start=1):
            total_counts["hmdb"] += 1
            chebi_id = normalize_chebi_ids(row.get("chebi_id"))
            smiles = _clean(row.get("smiles"))
            compound_id = _clean(row.get("hmdb_id"))
            chebi_source = _clean(row.get("chebi_source"))

            if not chebi_id:
                drop_counts["hmdb"]["chebi"] += 1
                continue
            if not smiles:
                kept_missing_smiles["hmdb"] += 1

            key = _record_key(smiles, chebi_id, "hmdb", compound_id, row_index)
            if key not in merged:
                merged[key] = {
                    "chebi_id": chebi_id,
                    "id": compound_id,
                    "smiles": smiles,
                    "chebi_source": chebi_source,
                    "source": "hmdb",
                    "hmdb_id": compound_id,
                    "wikidata_id": "",
                }
            else:
                rec = merged[key]
                if compound_id and not rec.get("hmdb_id"):
                    rec["hmdb_id"] = compound_id
                    rec["source"] = "both" if rec.get("wikidata_id") else "hmdb"
                rec["chebi_source"] = _merge_chebi_source(rec.get("chebi_source", ""), chebi_source)

    # Wikidata rows
    with open(wikidata_path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_index, row in enumerate(reader, start=1):
            total_counts["wikidata"] += 1
            chebi_id = normalize_chebi_ids(row.get("chebi_id"))
            smiles = pick_preferred_smiles(row)
            compound_id = wikidata_id_suffix(row.get("wikidataId"))
            chebi_source = _clean(row.get("chebi_source"))

            if not chebi_id:
                drop_counts["wikidata"]["chebi"] += 1
                continue
            if not smiles:
                kept_missing_smiles["wikidata"] += 1

            key = _record_key(smiles, chebi_id, "wikidata", compound_id, row_index)
            if key not in merged:
                merged[key] = {
                    "chebi_id": chebi_id,
                    "id": compound_id,
                    "smiles": smiles,
                    "chebi_source": chebi_source,
                    "source": "wikidata",
                    "hmdb_id": "",
                    "wikidata_id": compound_id,
                }
            else:
                rec = merged[key]
                if compound_id and not rec.get("wikidata_id"):
                    rec["wikidata_id"] = compound_id
                    rec["source"] = "both" if rec.get("hmdb_id") else "wikidata"
                rec["chebi_source"] = _merge_chebi_source(rec.get("chebi_source", ""), chebi_source)

    # Finalize id field and keep only requested output columns.
    rows = []
    for rec in merged.values():
        hmdb_id = _clean(rec.get("hmdb_id"))
        wikidata_id = _clean(rec.get("wikidata_id"))

        if rec["source"] == "both":
            rec["id"] = "|".join([value for value in [hmdb_id, wikidata_id] if value])
        elif rec["source"] == "hmdb":
            rec["id"] = hmdb_id
        else:
            rec["id"] = wikidata_id

        rows.append(
            {
                "chebi_id": rec["chebi_id"],
                "id": rec["id"],
                "smiles": rec["smiles"],
                "chebi_source": rec["chebi_source"],
                "source": rec["source"],
            }
        )

    rows.sort(key=lambda r: (r["source"], r["id"], r["smiles"]))

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["chebi_id", "id", "smiles", "chebi_source", "source"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Saved {len(rows)} combined compounds to {output_path}")

    for source in ["hmdb", "wikidata"]:
        dropped_total = sum(drop_counts[source].values())
        print(f"{source}: input={total_counts[source]}, dropped={dropped_total}, kept={total_counts[source] - dropped_total}")
        print(
            f"{source}: dropped_missing_chebi={drop_counts[source].get('chebi', 0)}, "
            f"kept_with_missing_smiles={kept_missing_smiles[source]}"
        )


if __name__ == "__main__":
    combine_datasets()
