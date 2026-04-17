import requests
import pandas as pd
from pathlib import Path

"""
The scripts loads the latest Wikidata LOTUS dump from Zenodo, which contains pre-extracted compound-taxon pairs.

compounds.tsv - chemical structures metadata (wikidataId, canonicalSmiles, isomericSmiles, inchi, inchiKey)
references.tsv - bibliographical references metadata (wikidataId, pipe separated DOIs, titles)
taxa.tsv - biological organisms metadata (wikidataId, pipe separated names, taxa rank)
compound_reference_taxon.tsv - the documented structure-organism pairs

"""
def download_wikidata_lotus():
    SAVE_DIR = Path(__file__).resolve().parent.parent / "data" / "wikidata"
    SAVE_DIR.mkdir(parents=True, exist_ok=True)

    USER_AGENT = "MetabolomicsEnrichment/1.0 (miranda.carlsson@idiap.ch)"

    # Get latest version info from Zenodo API
    print("Fetching Zenodo metadata...")
    zenodo = requests.get(
        "https://zenodo.org/api/records/5668854",
        headers={"User-Agent": USER_AGENT}
    ).json()

    files = {f["key"]: f["links"]["self"] for f in zenodo["files"]}
    print("Available files:", list(files.keys()))

    # Download and save each file as-is
    dfs = {}
    for filename, url in files.items():
        save_path = SAVE_DIR / filename
        
        if save_path.exists():
            print(f"  {filename} already exists, loading from disk...")
            dfs[filename] = pd.read_csv(save_path, sep="\t")
            continue
        
        print(f"  Downloading {filename}...")
        r = requests.get(url, headers={"User-Agent": USER_AGENT}, stream=True)
        r.raise_for_status()
        
        with open(save_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print(f"  Saved to {save_path}")
        dfs[filename] = pd.read_csv(save_path, sep="\t")

    compounds = dfs["compounds.tsv"]
    taxa = dfs["taxa.tsv"]
    pairs = dfs["compound_reference_taxon.tsv"]

def print_summary():
    save_dir = Path(__file__).resolve().parent.parent / "data" / "wikidata"
    compounds_path = save_dir / "compounds.tsv"
    taxa_path = save_dir / "taxa.tsv"
    pairs_path = save_dir / "compound_reference_taxon.tsv"

    if not compounds_path.exists() or not taxa_path.exists() or not pairs_path.exists():
        print("Missing LOTUS files. Run download_wikidata_lotus first.")
        return

    compounds = pd.read_csv(compounds_path, sep="\t")
    taxa = pd.read_csv(taxa_path, sep="\t")
    pairs = pd.read_csv(pairs_path, sep="\t")

    print(f"\nTotal counts:")
    print(f"  Compounds:              {len(compounds)}")
    print(f"  Taxa:                   {len(taxa)}")
    print(f"  Structure-organism pairs: {len(pairs)}")

def count_human_entries():
    homo_sapiens_uri = "http://www.wikidata.org/entity/Q15978631"
    save_dir = Path(__file__).resolve().parent.parent / "data" / "wikidata"
    taxa_path = save_dir / "taxa.tsv"
    pairs_path = save_dir / "compound_reference_taxon.tsv"

    if not taxa_path.exists() or not pairs_path.exists():
        print("Missing LOTUS files. Run download_wikidata_lotus first.")
        return

    taxa = pd.read_csv(taxa_path, sep="\t")
    pairs = pd.read_csv(pairs_path, sep="\t")

    human_taxa = taxa[taxa["wikidataId"] == homo_sapiens_uri]

    if human_taxa.empty:
        print("Homo sapiens URI not found in taxa.tsv, falling back to name match.")
        human_taxa = taxa[taxa["names_pipe_separated"].str.contains("Homo sapiens", na=False)]

    human_pairs = pairs[pairs["taxon"].isin(human_taxa["wikidataId"])]

    print(f"Homo sapiens taxon entries:    {len(human_taxa)}") # From taxa.tsv
    print(f"Human structure-organism pairs: {len(human_pairs)}") # From compound_reference_taxon.tsv. Number of rows that have the taxon human sapiens
    print(f"Unique human compounds:         {human_pairs['compound'].nunique()}") # From compound_reference_taxon.tsv. Number of unique compounds associated with human taxa


def find_top_n_duplicate_human_compounds(top_n=10): # Just used for simple checks
    homo_sapiens_uri = "http://www.wikidata.org/entity/Q15978631"
    save_dir = Path(__file__).resolve().parent.parent / "data" / "wikidata"
    taxa_path = save_dir / "taxa.tsv"
    pairs_path = save_dir / "compound_reference_taxon.tsv"

    if not taxa_path.exists() or not pairs_path.exists():
        print("Missing LOTUS files. Run download_wikidata_lotus first.")
        return

    taxa = pd.read_csv(taxa_path, sep="\t")
    pairs = pd.read_csv(pairs_path, sep="\t")

    human_taxa = taxa[taxa["wikidataId"] == homo_sapiens_uri]
    if human_taxa.empty:
        print("Homo sapiens URI not found in taxa.tsv, falling back to name match.")
        human_taxa = taxa[taxa["names_pipe_separated"].str.contains("Homo sapiens", na=False)]

    human_pairs = pairs[pairs["taxon"].isin(human_taxa["wikidataId"])]
    if human_pairs.empty:
        print("No human pairs found.")
        return

    compound_counts = human_pairs["compound"].value_counts()
    repeated = compound_counts[compound_counts > 1]

    if repeated.empty:
        print("No repeated compounds found in human pairs.")
        return

    top_repeated = repeated.head(top_n)

    print(f"Top {len(top_repeated)} repeated compounds in human pairs:")
    print("compound\toccurrences")
    for compound_id, count in top_repeated.items():
        print(f"{compound_id}\t{int(count)}")

def connect_smiles_to_chebi_ids(new_file_path):
    chebi_smiles_file = 'data/removed_leaf_classes_with_smiles.csv'
    chebi_smiles = pd.read_csv(chebi_smiles_file)
    chebi_smiles["SMILES"] = chebi_smiles["SMILES"].astype(str).str.strip('"').str.strip()

    wiki_compunds_file = 'data/wikidata/compounds.tsv'
    wiki_compounds = pd.read_csv(wiki_compunds_file, sep="\t")


    chebi_smiles_dict = dict(zip(chebi_smiles['SMILES'], chebi_smiles['IRI']))

    # Add a new column to wiki_compounds for the matched ChEBI ID
    # For each row in wiki_compounds, check if the canonicalSmiles matches any SMILES in chebi_smiles_dict and if so, add the corresponding ChEBI ID to the new column
    # If theres is no match, try to match the isomericSmiles column instead. If still no match, leave the new column as None
    
    wiki_compounds['chebi_id'] = None

    for index, row in wiki_compounds.iterrows():
        canonical_smiles = row['canonicalSmiles']
        isomeric_smiles = row['isomericSmiles']

        if pd.notna(canonical_smiles) and canonical_smiles in chebi_smiles_dict:
            wiki_compounds.at[index, 'chebi_id'] = chebi_smiles_dict[canonical_smiles]
        elif pd.notna(isomeric_smiles) and isomeric_smiles in chebi_smiles_dict:
            wiki_compounds.at[index, 'chebi_id'] = chebi_smiles_dict[isomeric_smiles]
    
        # Print progress every 1000 rows
        if index % 1000 == 0:
            print(f"Processed {index} / {len(wiki_compounds)} rows")

    wiki_compounds.to_csv(new_file_path, index=False, sep="\t")

def count_matched_chebi_ids(file_path):
    df = pd.read_csv(file_path, sep="\t")
    matched_count = df['chebi_id'].notna().sum()
    total_count = len(df)
    print(f"Matched ChEBI IDs: {matched_count} out of {total_count} compounds ({(matched_count/total_count)*100:.2f}%)")


def add_taxon_and_taxon_names(input_file_path, output_file_path):
    compounds = pd.read_csv(input_file_path, sep="\t")
    pairs = pd.read_csv("data/wikidata/compound_reference_taxon.tsv", sep="\t")
    taxa = pd.read_csv("data/wikidata/taxa.tsv", sep="\t")

    compound_col = "wikidataId"
    taxon_lookup_col = "wikidataId"

    if "compound" not in pairs.columns or "taxon" not in pairs.columns:
        raise KeyError("compound_reference_taxon.tsv must contain 'compound' and 'taxon' columns")
    if "names_pipe_separated" not in taxa.columns:
        raise KeyError("taxa.tsv must contain 'names_pipe_separated' column")

    # A compound can map to multiple taxa; store all unique matches as a pipe-separated string.
    compound_to_taxa = (
        pairs.dropna(subset=["compound", "taxon"])
        .groupby("compound")["taxon"]
        .agg(lambda values: "|".join(sorted(set(values.astype(str)))))
        .to_dict()
    )

    taxa_to_names = (
        taxa.dropna(subset=[taxon_lookup_col, "names_pipe_separated"])
        .groupby(taxon_lookup_col)["names_pipe_separated"]
        .agg(lambda values: "|".join(sorted(set(values.astype(str)))))
        .to_dict()
    )

    compounds["taxon"] = compounds[compound_col].map(compound_to_taxa)

    def map_taxa_to_names(taxa_string):
        if pd.isna(taxa_string):
            return None
        names = []
        for taxon_id in str(taxa_string).split("|"):
            taxon_name = taxa_to_names.get(taxon_id)
            if taxon_name:
                names.append(taxon_name)
        if not names:
            return None
        return "|".join(sorted(set(names)))

    compounds["taxon_name"] = compounds["taxon"].apply(map_taxa_to_names)

    Path(output_file_path).parent.mkdir(parents=True, exist_ok=True)
    compounds.to_csv(output_file_path, index=False, sep="\t")

    matched_taxon_count = compounds["taxon"].notna().sum()
    matched_name_count = compounds["taxon_name"].notna().sum()
    print(f"Saved enriched file to {output_file_path}")
    print(
        f"Matched taxon for {matched_taxon_count} out of {len(compounds)} compounds "
        f"({(matched_taxon_count/len(compounds))*100:.2f}%)."
    )
    print(
        f"Matched taxon names for {matched_name_count} out of {len(compounds)} compounds "
        f"({(matched_name_count/len(compounds))*100:.2f}%)."
    )


def keep_homo_sapiens_compounds(input_file_path, output_file_path):
    homo_sapiens_uri = "http://www.wikidata.org/entity/Q15978631"
    compounds = pd.read_csv(input_file_path, sep="\t")

    if "taxon" not in compounds.columns:
        raise KeyError(
            "Input file must contain column 'taxon'. "
            "Run add_taxon_and_taxon_names first."
        )

    taxon_mask = compounds["taxon"].astype(str).str.contains(homo_sapiens_uri, na=False)
    homo_sapiens_compounds = compounds[taxon_mask].copy()

    Path(output_file_path).parent.mkdir(parents=True, exist_ok=True)
    homo_sapiens_compounds.to_csv(output_file_path, index=False, sep="\t")

    print(f"Saved Homo sapiens subset to {output_file_path}")
    print(
        f"Kept {len(homo_sapiens_compounds)} out of {len(compounds)} compounds "
        f"({(len(homo_sapiens_compounds)/len(compounds))*100:.2f}%)."
    )

def count_chebi_ids(input_file_path):
    compounds = pd.read_csv(input_file_path, sep="\t")
    matched_chebi_count = compounds["chebi_id"].notna().sum()
    print(f"Matched ChEBI IDs for {matched_chebi_count} out of {len(compounds)} compounds in {input_file_path}"
          f"({(matched_chebi_count/len(compounds))*100:.2f}%).")


if __name__ == "__main__":
    
    task = count_chebi_ids 
    # options: download_wikidata_lotus, count_human_entries, print_summary, find_top_n_duplicate_human_compounds, connect_smiles_to_chebi_ids, add_taxon_and_taxon_names, keep_homo_sapiens_compounds, count_chebi_ids

    if task == download_wikidata_lotus:
        download_wikidata_lotus()
        print_summary()
    elif task == count_human_entries:
        count_human_entries()
    elif task == print_summary:
        print_summary()
    elif task == find_top_n_duplicate_human_compounds:
        find_top_n_duplicate_human_compounds(top_n=10)
    elif task == connect_smiles_to_chebi_ids:
        new_file_path = 'data/wikidata/created/compounds_with_chebi_ids.tsv'
        # create the directory if it doesn't exist
        Path(new_file_path).parent.mkdir(parents=True, exist_ok=True)
        connect_smiles_to_chebi_ids(new_file_path)
        print(f"Saved compounds with ChEBI IDs to {new_file_path}")
        count_matched_chebi_ids(new_file_path)
    elif task == add_taxon_and_taxon_names:
        input_file_path = "data/wikidata/created/compounds_with_chebi_ids.tsv"
        output_file_path = input_file_path
        add_taxon_and_taxon_names(input_file_path, output_file_path)
    elif task == keep_homo_sapiens_compounds:
        input_file_path = "data/wikidata/created/compounds_with_chebi_ids.tsv"
        output_file_path = "data/wikidata/created/compounds_with_chebi_ids_homo_sapiens.tsv"
        keep_homo_sapiens_compounds(input_file_path, output_file_path)
    elif task == count_chebi_ids:
        input_file_path = "data/wikidata/created/compounds_with_chebi_ids_homo_sapiens.tsv"
        count_chebi_ids(input_file_path)
    else:
        print("Invalid task specified.")
        
    
