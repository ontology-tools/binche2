import pyhornedowl
import os
import urllib.request
import pandas as pd

"""
This script downloads the ChEBI ontology in OWL format and loads it using pyhornedowl.
"""

# URL and local save path
CHEBI_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.owl"
DATA_DIR = "data"
OWL_PATH = os.path.join(DATA_DIR, "chebi.owl")

def download_chebi():
    os.makedirs(DATA_DIR, exist_ok=True)

    # Download the OWL file if it is not already present
    if not os.path.exists(OWL_PATH):
        print("⬇ Downloading ChEBI ontology...")
        urllib.request.urlretrieve(CHEBI_URL, OWL_PATH)
        print("Download complete.")
    else:
        print("ChEBI ontology already exists locally.")

    return OWL_PATH

def load_chebi():
    owl_file = download_chebi()
    print("Loading ChEBI ontology...")
    ontology = pyhornedowl.open_ontology(owl_file)
    print(f"Loaded ontology with {len(ontology.get_classes())} classes "
          f"and {len(ontology.get_axioms())} axioms.")
    return ontology

def load_ontology(owl_file):
    print("Loading ontology...")
    ontology = pyhornedowl.open_ontology(owl_file)
    print(f"Loaded ontology from {owl_file} with {len(ontology.get_classes())} classes "
          f"and {len(ontology.get_axioms())} axioms.")
    return ontology

if __name__ == "__main__":
    # chebi_ontology = load_chebi()
    # load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl")

    csv = ("data/removed_leaf_classes_with_smiles.csv")
    # open and read the CSV file
    df = pd.read_csv(csv)
    len_df = len(df)
    print(f"Number of entries in CSV: {len_df}")
    
