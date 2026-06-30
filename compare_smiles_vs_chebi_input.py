"""Compare a SMILES-resolved study set against a directly-given ChEBI ID study set,
and compare their resulting enrichment CSVs.

Two independent comparisons are supported:

1. Study-set level (input level): paste the "Analysed leaves" text block the web app
   prints after a run (one for the SMILES run, one for the ChEBI ID run) and compare
   which ChEBI leaves they actually resolved to. This is the most direct way to answer
   "how many compounds match directly" between the two input modes, since it compares
   the resolved study sets *before* any graph pruning happens.

2. Result level (output level): compare the two final enrichment_results CSVs
   (Class, Raw p-value, Corrected p-value) for shared/unique classes and for how much
   the corrected p-values differ among classes both runs kept.

Usage:
    python compare_smiles_vs_chebi_input.py \
        --csv-a data/enrichment_results/enrichment_results_smiles_structure_prune.csv \
        --csv-b data/enrichment_results/enrichment_results_wine_structure_plain_prune.csv \
        --label-a smiles --label-b chebi_id
"""
import argparse
import re
import sys

import pandas as pd
from scipy.stats import spearmanr

CHEBI_RE = re.compile(r"([^(),]+?)\s*\((CHEBI_\d+)\)")


def parse_chebi_block(text: str) -> dict[str, str]:
    """Parse a comma-separated 'name (CHEBI_xxx), name (CHEBI_xxx), ...' block
    (as printed by the web app for 'Analysed leaves' / 'Removed nodes through
    pruning') into {chebi_id: name}.
    """
    return {chebi_id: name.strip() for name, chebi_id in CHEBI_RE.findall(text)}


def compare_study_sets(set_a: dict[str, str], set_b: dict[str, str], label_a: str, label_b: str) -> None:
    ids_a, ids_b = set(set_a), set(set_b)
    shared = ids_a & ids_b
    only_a = ids_a - ids_b
    only_b = ids_b - ids_a
    jaccard = len(shared) / len(ids_a | ids_b) if (ids_a | ids_b) else float("nan")

    print(f"\n=== Study set comparison: {label_a} vs {label_b} ===")
    print(f"{label_a}: {len(ids_a)} resolved leaves")
    print(f"{label_b}: {len(ids_b)} resolved leaves")
    print(f"Matching directly (same ChEBI leaf in both): {len(shared)}  (Jaccard = {jaccard:.3f})")

    if only_a:
        print(f"\nOnly resolved in {label_a} ({len(only_a)}):")
        for cid in sorted(only_a):
            print(f"  {set_a[cid]} ({cid})")
    if only_b:
        print(f"\nOnly resolved in {label_b} ({len(only_b)}):")
        for cid in sorted(only_b):
            print(f"  {set_b[cid]} ({cid})")


def load_enrichment_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["ChEBI_ID"] = df["Class"].str.extract(r"(CHEBI_\d+)")
    df["Name"] = df["Class"].str.replace(r"\s*\(CHEBI_\d+\)", "", regex=True)
    return df.set_index("ChEBI_ID")


def compare_enrichment_csvs(path_a: str, path_b: str, label_a: str, label_b: str, alpha: float = 0.05) -> pd.DataFrame:
    df_a = load_enrichment_csv(path_a)
    df_b = load_enrichment_csv(path_b)

    merged = df_a[["Name", "Raw p-value", "Corrected p-value"]].join(
        df_b[["Name", "Raw p-value", "Corrected p-value"]],
        how="outer",
        lsuffix=f"_{label_a}",
        rsuffix=f"_{label_b}",
    )
    merged["Name"] = merged[f"Name_{label_a}"].combine_first(merged[f"Name_{label_b}"])
    merged = merged.drop(columns=[f"Name_{label_a}", f"Name_{label_b}"])

    only_a = merged[merged[f"Corrected p-value_{label_b}"].isna()]
    only_b = merged[merged[f"Corrected p-value_{label_a}"].isna()]
    shared = merged.dropna(subset=[f"Corrected p-value_{label_a}", f"Corrected p-value_{label_b}"]).copy()

    shared["ratio"] = shared[f"Corrected p-value_{label_a}"] / shared[f"Corrected p-value_{label_b}"]
    shared[f"significant_{label_a}"] = shared[f"Corrected p-value_{label_a}"] < alpha
    shared[f"significant_{label_b}"] = shared[f"Corrected p-value_{label_b}"] < alpha
    shared["flipped"] = shared[f"significant_{label_a}"] != shared[f"significant_{label_b}"]

    rho, _ = spearmanr(shared[f"Corrected p-value_{label_a}"], shared[f"Corrected p-value_{label_b}"])

    print(f"\n=== Enrichment CSV comparison: {label_a} vs {label_b} (alpha={alpha}) ===")
    print(f"Classes in {label_a} only: {len(only_a)}")
    print(f"Classes in {label_b} only: {len(only_b)}")
    print(f"Classes in both: {len(shared)}")
    print(f"Spearman correlation of corrected p-values (shared classes): {rho:.3f}")
    print(f"Classes whose significance call flips at alpha={alpha}: {shared['flipped'].sum()}")

    print(f"\nTop 10 largest |log fold-change| in corrected p-value among shared classes:")
    top = shared.reindex(shared["ratio"].apply(lambda r: abs(__import__('math').log10(r))).sort_values(ascending=False).index)
    cols = ["Name", f"Corrected p-value_{label_a}", f"Corrected p-value_{label_b}", "ratio"]
    print(top[cols].head(10).to_string(index=False))

    return merged


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv-a", default="data/results/wine_enrichmentresults_prune_ids_new.csv")
    parser.add_argument("--csv-b", default="data/results/wine_enrichmentresults_prune_smiles_new.csv")
    parser.add_argument("--label-a", default="smiles")
    parser.add_argument("--label-b", default="chebi_id")
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()

    merged = compare_enrichment_csvs(args.csv_a, args.csv_b, args.label_a, args.label_b, args.alpha)
    merged.to_csv("data/results/comparison_smiles_vs_chebi.csv")
    print(f"\nFull merged comparison table written to data/enrichment_results/comparison_smiles_vs_chebi.csv")
