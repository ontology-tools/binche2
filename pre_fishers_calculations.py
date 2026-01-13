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
import json
from collections import defaultdict
from tqdm import tqdm
from collections import deque

def count_removed_leaves(removed_classes_csv, classification):
    removed_classes = pd.read_csv(removed_classes_csv)
    # Count how many classes that are set to "structural"/"functional" as their Classification
    if classification == "full":
        return len(removed_classes)
    elif classification not in ["structural", "functional"]:
        print("Classification must be either 'structural' or 'functional'.")
        return 0
    removed_structural = removed_classes[removed_classes["Classification"] == classification]
    return len(removed_structural)

def build_class_to_leaf_map(leaf_to_ancestors_file, class_to_leaf_output_file):
    """ Build a JSON map from each class IRI to ALL its leaf descendants using an existing leaf-to-ancestors map."""

    print(f"Loading leaf to ancestors map from {leaf_to_ancestors_file}...")
    with open(leaf_to_ancestors_file, 'r') as f:
        leaf_to_ancestors = json.load(f)
    print(f"Loaded leaf to ancestors map with {len(leaf_to_ancestors)} leaf classes.")

    class_to_leafs = defaultdict(set)
    print("Building class to leaf descendants map...")

    # Collect all leaf classes
    all_leaf_classes = set(leaf_to_ancestors.keys())

    # For each leaf, add it to all its ancestors
    for leaf, ancestors in leaf_to_ancestors.items():
        for ancestor in ancestors:
            class_to_leafs[ancestor].add(leaf)

    ### Include if I also want leaf classes to appear in the map. Now these will give the same output as something not in the ontology at all. 
    # # Ensure all leaf classes appear with empty lists
    # for leaf in all_leaf_classes:
    #     if leaf not in class_to_leafs:
    #         class_to_leafs[leaf] = set()
    
    # Convert sets to lists for JSON serialization
    class_to_leaf_json = {
        cls: list(leafs)
        for cls, leafs in class_to_leafs.items()} 
    
    # Save to JSON
    with open(class_to_leaf_output_file, 'w') as f:
        json.dump(class_to_leaf_json, f, indent=2)

    print(f"Saved class to leaf descendants map with {len(class_to_leaf_json)} classes to {class_to_leaf_output_file}.")
    
def count_removed_classes_for_class(class_iri, subclass_map, classification, check_leaf_classes = False, removed_classes_csv = None):
    """Count how many classes under the given class_iri were removed in the structural ontology."""

    if str(class_iri) not in subclass_map: 
        print(f"⚠️ Class {class_iri} is not found in map file.")
        return 0, 0
    elif len(subclass_map[str(class_iri)]) == 0: # Should not happen since then it would be a leaf
        print(f"⚠️ Class {class_iri} has no descendants in full ontology.")
        return 0, 0
    
    # Get string of sublclasses from the map
    subclasses = subclass_map[str(class_iri)]
    n_subclasses = len(subclasses)

    """ Checks that all found subclasses are of the expected type (structural or functional)""" # Probably not needed
    if check_leaf_classes:
        if classification not in ["structural", "functional"]:
            print("Classification must be either 'structural' or 'functional' to check leaf classes.")
            return subclasses, n_subclasses
        
        # Read csv with removed leaf classes
        all_leaf_classes_correct = True
        try:
            removed_classes_df = pd.read_csv(removed_classes_csv)
            print("Checking leaf classes")
        except FileNotFoundError:
            print("⚠️ Removed leaf classes CSV file not found. Cannot double check leaf classes.")
            return subclasses, n_subclasses
        # Check that all subclasses exist in removed classes and that the column "Classification" are of the correct type
        for sub in subclasses:
            if sub not in removed_classes_df["IRI"].values:
                print(f"⚠️ Subclass {sub} not found in removed classes csv.")
                all_leaf_classes_correct = False
            else:
                class_type = removed_classes_df[removed_classes_df["IRI"] == sub]["Classification"].values[0]
                if class_type != classification:
                    print(f"⚠️ Subclass {sub} is of type {class_type}, expected {classification}.")
                    all_leaf_classes_correct = False
        if all_leaf_classes_correct:
            print(f"All leaf classes are found and correctly identified as {classification}.")

    return subclasses, n_subclasses
    
if __name__ == "__main__":

    """ Select task to perform. To calculate the removed leaf classes for a given class, 
    the file with the class to leaf descendants map must be created first.
    This is done in "build_class_to_leaf_map" task. """

    task = "count_removed_classes_for_class" 
    # Options: "count_total_removed_leaves" "count_removed_classes_for_class" "build_class_to_leaf_map" "enrichment_analysis_plain"

    # Variables used in "count_removed_classes_for_class":
    # - classification (only if check_leaf_classes is True), check_leaf_classes, class_iri

    # Variables used in "count_total_removed_leaves":
    # - classification

    # Variables used in "enrichment_analysis_plain":
    # - classification, class_iri, (check_leaf_classes)

    # Relevant for task "count_removed_classes_for_class"
    classification = "structural" # "functional" or "structural" or "full"
    check_leaf_classes = True # Checks that the found the leaf classes are of the exppected type (Functional or Structural)
    class_iri = "http://purl.obolibrary.org/obo/CHEBI_83822"
    # Not found in map: "http://purl.obolibrary.org/obo/CHEBI_38870"
    #Children for "http://purl.obolibrary.org/obo/CHEBI_38867" are not found in the csv but are found in map

    # Files
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
    map_file = "data/class_to_leaf_descendants_map.json"
    leaf_to_ancestors_file = "data/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_output_file = "data/class_to_leaf_descendants_map.json"

    if task == "count_total_removed_leaves":

        n_removed_leaves = count_removed_leaves(removed_leaves_csv, classification)
        print(f"Total number of removed {classification} leaf classes: {n_removed_leaves}")

        # Output
        # - Total number of removed structural leaf classes: 184 436
        # - Total number of removed functional leaf classes: 41

    elif task == "count_removed_classes_for_class":

        subclasses, n_subclasses = count_removed_classes_for_class(class_iri, map_file, classification, check_leaf_classes, removed_leaves_csv)
        print(f"Class {class_iri} has {n_subclasses} removed subclasses in the ontology.")

        if subclasses and n_subclasses < 50: 
            for sub in subclasses:
                print(f" - {sub}")

    elif task == "build_class_to_leaf_map":

        build_class_to_leaf_map(leaf_to_ancestors_file, class_to_leaf_output_file)

    elif task == "enrichment_analysis_plain": # to be removed

        _, n_subclasses = count_removed_classes_for_class(class_iri, map_file, classification, check_leaf_classes, removed_leaves_csv)
        n_tot_removed_leaves = count_removed_leaves(removed_leaves_csv, classification)

        print(f"Calulated for {classification} ontology for {class_iri} (make sure classification is correct):")
        prob_of_success = n_subclasses / n_tot_removed_leaves

        print(f" The probability of success value for class {class_iri} is: {prob_of_success:.4g} ({n_subclasses} / {n_tot_removed_leaves})")

        # TODO?: Could check that the input class actually has the expected classification by checking which filtered entology it is found in



    else:
        print("No valid task selected.")

        # Compare length of map with removed classes csv
        # removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
        # removed_classes_df = pd.read_csv(removed_leaves_csv)
        # print(f"Removed classes CSV has {len(removed_classes_df)} entries.")
        # with open("data/removed_leaf_classes_to_parents.json", 'r') as f:
        #     removed_classes_map = json.load(f)
        # print(f"Removed classes map has {len(removed_classes_map)} entries.")

        # Output:
        # Removed classes CSV has 187470 entries.
        # Removed classes map has 187470 entries.

        # compare parent classes in filtered chebi and in json file
        # parent_map = "data/parent_to_leaf_class_map.json"
        # with open(parent_map, 'r') as f:
        #     parent_map_data = json.load(f)
        # print(f"Parent map has {len(parent_map_data)} entries.")
        # structural_file = "data/filtered_chebi_no_leaves_with_smiles_no_deprecated_new.owl"
        # structural_ontology = load_ontology(structural_file)

        # load_ontology("data/filtered_chebi_no_leaves_with_smiles_no_deprecated_new.owl")