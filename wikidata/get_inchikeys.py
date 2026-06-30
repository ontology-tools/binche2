print("This script converts SMILES strings to InChIKeys using RDKit.")

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
import pandas as pd

def smiles_to_inchikey(smiles, starcount=0):
    # print(f"Converting SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, starcount
    if "*" in smiles:
        starcount += 1
        return None, starcount
    inchi = MolToInchi(mol)
    return InchiToInchiKey(inchi), starcount

def canonical_smiles(smiles):
    """RDKit-canonical form of a SMILES string, or None if it can't be parsed.

    ChEBI's asserted SMILES aren't guaranteed to be in any particular canonical
    form, so baking the RDKit-canonical form into this file lets the website's
    exact-string lookup match incoming SMILES regardless of how they were written.
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol is not None else None

# Example
smiles = "[H]C(=Nc1ccccc1NS(=O)(=O)c1ccc(C)cc1)c1ccccc1NS(=O)(=O)c1ccc(C)cc1" 
inchikey = smiles_to_inchikey(smiles)
print(inchikey) 

def convert_smiles_file(input_file, output_file):

    starcount = 0 # count how many SMILES had a star in them, which means they are not valid for conversion

    # Open the input CSV file and read it
    df = pd.read_csv(input_file)
    print(f"Read {len(df)} rows from {input_file}")
    # Strip triple quotes from the SMILES column
    df['SMILES'] = df['SMILES'].str.strip('"')
    # Canonicalize so the website's exact-string lookup is toolkit-consistent
    df['SMILES'] = df['SMILES'].apply(lambda x: canonical_smiles(x) or x)
    # Convert SMILES to InChIKey and store in a new column
    df['InChIKey'] = df['SMILES'].apply(lambda x: smiles_to_inchikey(x, starcount)[0])
    starcount = df['SMILES'].apply(lambda x: smiles_to_inchikey(x, starcount)[1]).sum()
    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_file, index=False)

    # Print summary statistics
    total_rows = len(df)
    generated_keys = df['InChIKey'].notnull().sum()
    failed_conversions = total_rows - generated_keys
    print(f"Processed {total_rows} rows.")
    print(f"Generated {generated_keys} InChIKeys.")
    print(f"Failed conversions: {failed_conversions}.")
    print(f"SMILES with stars (not converted): {starcount}.")

def count_nans(csv_file):
    df = pd.read_csv(csv_file)
    # check how many of the rows that do not have an inchikey that have Classification 'strictural'
    column_name = 'InChIKey'
    na_count = df[column_name].isna().sum()
    structural_na_count = df[df['Classification'] == 'structural'][column_name].isna().sum()

    print(f"Total NaN in '{column_name}': {na_count}")
    print(f"NaN in '{column_name}' where Classification is 'structural': {structural_na_count}")




    
if __name__ == "__main__":
    input_file = "data/removed_leaf_classes_with_smiles.csv"
    output_file = "data/removed_leaf_classes_with_inchikeys.csv"

    convert_smiles_file(input_file, output_file)
    count_nans(output_file)