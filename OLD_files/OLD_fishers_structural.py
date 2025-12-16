"I need to count classes removed from a given class. This will only be based on the structural split."
"The total number of classes in the structural split will also be counted."
"OBS: Need to make sure to only count a removed class once."

from sqlite3 import Time
from symtable import Class
from load_chebi import load_chebi, load_ontology
from pruning_split_up_structure import get_descendants
import xml.etree.ElementTree as ET
import os
import pandas as pd
import time

def count_removed_classes_for_structural_class(class_iri, full_ontology, structural_ontology, removed_leaves_csv):
    """Count how many classes under the given class_iri were removed in the structural ontology."""

    structural_classes = structural_ontology.get_classes()
    connected_leaf_classes = set()

    # Check that class is in structural ontology
    if class_iri not in structural_classes:
        print(f"Class {class_iri} not found in structural ontology.")
        return 0, 0
    
    # Get all descendants of the class in the full ontology
    descendants = get_descendants(full_ontology, class_iri)

    # If no descendants, it already is a leaf
    if not descendants:
        print(f"⚠️ Class {class_iri} has no descendants in full ontology.")
        return 0, 0

    df = pd.read_csv(removed_leaves_csv)
    removed_classes = set(df["IRI"].tolist())
    # Check if descendants are in the CSV
    for desc in descendants:
        if desc in removed_classes:
            if desc not in connected_leaf_classes:
                connected_leaf_classes.add(desc)
                # Check that the column "Classification" is "structural
                # row = df[df["IRI"] == desc]
                # if not row.empty and row.iloc[0]["Classification"] != "structural":
                    # print(f"⚠️ Leaf class {desc} for {class_iri} found in removed CSV but not classified as structural.")
 
    # Count total leaf classes under the given class
    n_removed_leaf_classes = len(connected_leaf_classes)
    return n_removed_leaf_classes, connected_leaf_classes

def create_csv_for_counts(csv):
    # Create csv 
    df = pd.DataFrame(columns=["Class", "N_removed_leaves"])
    df.to_csv(csv, index=False)

def count_removed_classes_for_structural_class_using_saved_counts(class_iri, structural_ontology_file, removed_leaves_csv, counts_csv):
    """Count how many classes under the given class_iri were removed in the structural ontology."""
    # TODO: Is there a smart way to use descendants saved counts to speed this up?

    # Check if class has already been counted in the CSV
    df_counts = pd.read_csv(counts_csv)
    # Check if class_iri is in the "Class" column
    if class_iri in df_counts["Class"].values:
        n_removed = df_counts[df_counts["Class"] == class_iri]["N_removed_leaves"].values[0]
        print(f"Class {class_iri} already counted in CSV with {n_removed} removed leaves.")
        return n_removed, None
    
    full_ontology = load_chebi()
    structural_ontology = load_ontology(structural_ontology_file)
    
    structural_classes = structural_ontology.get_classes()
    connected_leaf_classes = set()

    # Check that class is in structural ontology
    if class_iri not in structural_classes:
        print(f"Class {class_iri} not found in structural ontology.")
        return 0, 0
    
    # Get all descendants of the class in the full ontology
    descendants = get_descendants(full_ontology, class_iri)

    # If no descendants, it already is a leaf
    if not descendants:
        print(f"⚠️ Class {class_iri} has no descendants in full ontology.")
        return 0, 0

    df_removed_classes = pd.read_csv(removed_leaves_csv)
    removed_classes = set(df_removed_classes["IRI"].tolist())
    # Check if descendants are in the CSV
    for desc in descendants:
        if desc in removed_classes:
            if desc not in connected_leaf_classes:
                connected_leaf_classes.add(desc)
                # Check that the column "Classification" is "structural
                # row = df_removed_classes[df_removed_classes["IRI"] == desc]
                # if not row.empty and row.iloc[0]["Classification"] != "structural":
                #     print(f"⚠️ Leaf class {desc} for {class_iri} found in removed CSV but not classified as structural.")
 
    # Count total leaf classes under the given class
    n_removed_leaf_classes = len(connected_leaf_classes)

    # Append to CSV
    new_row = pd.DataFrame({"Class": [class_iri], "N_removed_leaves": [n_removed_leaf_classes]})
    new_row.to_csv(counts_csv, mode='a', header=False, index=False)

    return n_removed_leaf_classes, connected_leaf_classes


def count_structured_leaves(csv):
    removed_classes = pd.read_csv(csv)
    # Count how many are that are set to "structural" as their Classification
    removed_structural = removed_classes[removed_classes["Classification"] == "structural"]
    return len(removed_structural)



if __name__ == "__main__":

    task = "count_removed_for_class_using_saved_counts"  # Options: "count_removed_for_class", "count_total_structured_leaves", "load_ontology", "count_removed_for_class_using_saved_counts"

    class_iri = "http://purl.obolibrary.org/obo/CHEBI_35505"

    if task == "count_removed_for_class":

        full_ontology = load_chebi()
        structural_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl"
        structural_ontology = load_ontology(structural_file)

        removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"

        # Timer start
        start_time = time.time()

        n_removed, removed_classes = count_removed_classes_for_structural_class(
            class_iri, full_ontology, structural_ontology, removed_leaves_csv
        )
        # Timer end
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Time taken to find classes: {elapsed_time:.2f} seconds")

        print(f"Class {class_iri} has {n_removed} removed leaf classes.")
        #if n_removed > 0:
            #print(f"Removed classes for {class_iri}:")
            # for rc in removed_classes:
                # print(f" - {rc}")

        #Time taken to find classes: 1.10 seconds (new version 1.01)
        #Class http://purl.obolibrary.org/obo/CHEBI_17234 had 13 removed leaf classes.

        # Time taken to find classes: 243.79 seconds (new version 225.40)
        # Class http://purl.obolibrary.org/obo/CHEBI_16646 had 5131 removed leaf classes.

    elif task == "count_removed_for_class_using_saved_counts":
        removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
        counts_csv = "data/leaf_classes_counts_structural.csv"

        if not os.path.exists(counts_csv):
            create_csv_for_counts(counts_csv)
            print(f"Created CSV file for counts: {counts_csv}")

        structural_ontology_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl"

        n_removed, removed_classes = count_removed_classes_for_structural_class_using_saved_counts(
            class_iri, structural_ontology_file, removed_leaves_csv, counts_csv
        )

        print(f"Class {class_iri} has {n_removed} removed leaf classes.")

        # full_ontology = load_chebi()
        # descendants = get_descendants(full_ontology, class_iri)
        # print(f"Class {class_iri} has {len(descendants)} descendants in full ontology.")
        # # print descendants
        # for d in descendants:
        #     print(f" - {d}")

    elif task == "count_total_structured_leaves":

        removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
        n_structured_leaves = count_structured_leaves(removed_leaves_csv)
        print(f"Total number of removed structural leaf classes: {n_structured_leaves}")

        # Output - Total number of removed structural leaf classes: 184 436

    elif task == "load_ontology": # Can be removed later, just to count classes
        structural_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl"
        structural_ontology = load_ontology(structural_file)

        functional_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl"
        functional_ontology = load_ontology(functional_file)

    else:
        print("No valid task selected. Please choose a valid task.")