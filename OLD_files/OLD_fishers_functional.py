from sqlite3 import Time
from symtable import Class
from load_chebi import load_chebi, load_ontology
from pruning_split_up_structure import get_descendants
import xml.etree.ElementTree as ET
import os
import pandas as pd
import time

def count_functional_leaves(csv):
    removed_classes = pd.read_csv(csv)
    # Count how many are that are set to "functional" as their Classification
    removed_functional = removed_classes[removed_classes["Classification"] == "functional"]
    return len(removed_functional)

def create_csv_for_counts(csv):
    # Create csv 
    df = pd.DataFrame(columns=["Class", "N_removed_leaves"])
    df.to_csv(csv, index=False)

def count_removed_classes_for_functional_class_using_saved_counts(class_iri, functional_ontology_file, removed_leaves_csv, counts_csv):
    """Count how many classes under the given class_iri were removed in the functional ontology."""

    # Check if class has already been counted in the CSV
    df_counts = pd.read_csv(counts_csv)
    # Check if class_iri is in the "Class" column
    if class_iri in df_counts["Class"].values:
        n_removed = df_counts[df_counts["Class"] == class_iri]["N_removed_leaves"].values[0]
        print(f"Class {class_iri} already counted in CSV with {n_removed} removed leaves.")
        return n_removed, None
    
    full_ontology = load_chebi()
    functional_ontology = load_ontology(functional_ontology_file)

    functional_classes = functional_ontology.get_classes()
    connected_leaf_classes = set()

    # Check that class is in functional ontology
    if class_iri not in functional_classes:
        print(f"Class {class_iri} not found in functional ontology.")
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
                
    # Count total leaf classes under the given class
    n_removed_leaf_classes = len(connected_leaf_classes)

    # Append to CSV
    new_row = pd.DataFrame({"Class": [class_iri], "N_removed_leaves": [n_removed_leaf_classes]})
    new_row.to_csv(counts_csv, mode='a', header=False, index=False)

    return n_removed_leaf_classes, connected_leaf_classes

if __name__ == "__main__":

    task = "count_removed_for_class_using_saved_counts"  # Options: "count_total_functional_leaves", "count_removed_for_class_using_saved_counts"
    class_iri = "http://purl.obolibrary.org/obo/CHEBI_232401"

    if task == "count_total_functional_leaves":
        removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
        n_functional_leaves = count_functional_leaves(removed_leaves_csv)
        print(f"Total number of removed functional leaf classes: {n_functional_leaves}")

        # Total number of removed functional leaf classes: 41

    elif task == "count_removed_for_class_using_saved_counts":
        removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
        counts_csv = "data/leaf_classes_counts_functional.csv"

        if not os.path.exists(counts_csv):
            create_csv_for_counts(counts_csv)
            print(f"Created CSV file for counts: {counts_csv}")

        functional_ontology_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl"

        n_removed, removed_classes = count_removed_classes_for_functional_class_using_saved_counts(
            class_iri, functional_ontology_file, removed_leaves_csv, counts_csv
        )

        print(f"Class {class_iri} has {n_removed} removed leaf classes.")

    else:
        print("Unknown task. Please set 'task' variable correctly.")