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

def download_chebi(download_dir=None):
    """Download ChEBI ontology to specified directory (default: data/).
    
    Args:
        download_dir (str): Directory to download the OWL file to. If None, uses DATA_DIR ("data/").
    
    Returns:
        str: Path to the downloaded OWL file
    """
    if download_dir is None:
        download_dir = DATA_DIR
    
    os.makedirs(download_dir, exist_ok=True)
    owl_path = os.path.join(download_dir, "chebi.owl")

    # Download the OWL file if it is not already present
    if not os.path.exists(owl_path):
        print("⬇ Downloading ChEBI ontology...")
        urllib.request.urlretrieve(CHEBI_URL, owl_path)
        print("Download complete.")
    else:
        print("ChEBI ontology already exists locally.")

    return owl_path

def load_chebi(download_dir=None):
    """Load ChEBI ontology from specified directory (default: data/).
    
    Args:
        download_dir (str): Directory to download/load the OWL file from. If None, uses DATA_DIR ("data/").
                           Useful option: "data_new/" to keep the download separate from older versions.
    
    Returns:
        pyhornedowl ontology object
    """
    owl_file = download_chebi(download_dir=download_dir)
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
    chebi_ontology = load_chebi()

    # load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl")

    # csv = ("data/removed_leaf_classes_with_smiles.csv")
    # # open and read the CSV file
    # df = pd.read_csv(csv)
    # len_df = len(df)
    # print(f"Number of entries in CSV: {len_df}")
    
