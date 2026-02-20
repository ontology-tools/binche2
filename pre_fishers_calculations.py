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

def count_removed_leaves(removed_classes_csv):
    classification = 'structural' # Since 'functional' are errors in CHEBI ontology
    removed_classes = pd.read_csv(removed_classes_csv)
    # Count how many classes that are set to "structural"/"functional" as their Classification
    removed_set = removed_classes[removed_classes["Classification"] == classification]
    return len(removed_set)

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
    
def count_removed_classes_for_class(class_iri, leaf_descendants_map, classification, class_to_all_roles_map, roles_to_leaves_map):
    """
    Count how many leaf classes are associated with the given class_iri.
    
    Parameters:
        class_iri: The class to check
        leaf_descendants_map: Maps structural classes to their leaf descendants
        classification: "structural", "functional", or "full"
        roles_to_leaves_map: Maps role classes to their associated leaves
    
    Returns:
        leaves: Set of leaf class IRIs
        n_leaves: Count of leaves
    """

    # if classification not in ["structural", "functional", "full"]:
        # print("Classification must be either 'structural', 'functional' or 'full'.")
        # return 0, 0

    leaves = set()

    if classification in ["structural", "full"]:
        
        if str(class_iri) not in leaf_descendants_map: 
            print(f"⚠️ Class {class_iri} is not found in map file.")
            return 0, 0
        elif len(leaf_descendants_map[str(class_iri)]) == 0: # Should not happen since then it would be a leaf
            print(f"⚠️ Class {class_iri} has no descendants in full ontology.")
            return 0, 0

        # Get string of sublclasses from the map
        leaves_structure = leaf_descendants_map[str(class_iri)]
        leaves.update(leaves_structure)
    
    else:
        print(f"Classification {classification} is not supported for counting classes.")
        return 0, 0

    # if classification in ["functional", "full"]:
    #     # Get ALL roles for the class being tested (direct + inherited from ancestors)
    #     all_roles = class_to_all_roles_map.get(class_iri, [])

    #     if all_roles:
    #         # Collect all leaves associated with those roles
    #         leaves_role = set()
    #         for role in all_roles:
    #             leaves_role.update(roles_to_leaves_map.get(role, []))
    #             if role not in roles_to_leaves_map:
    #                 print(f"⚠️ Role {role} has no associated leaves in roles_to_leaves_map.")
    #             print(f"Role {role} has {len(roles_to_leaves_map.get(role, []))} associated leaves.")
    #         leaves.update(leaves_role)

    #         print(f"Class {class_iri} has {len(leaves_role)} functional leaf descendants from roles.")
    #         print(f"Class {class_iri} has {len(all_roles)} roles.")
    #     else:
    #         print(f"⚠️ Class {class_iri} has no associated roles in class_to_all_roles_map.")
        

    n_leaves = len(leaves)

    return leaves, n_leaves

def count_removed_classes_for_roles(class_iri, leaf_descendants_map, classification, roles_to_leaves_map):
    
    if classification not in ["functional", "full"]:
        print(f"Classification {classification} is not supported for counting roles.")
        return 0, 0

    leaves = set()
    leaves.update(roles_to_leaves_map.get(class_iri, []))

    if not leaves:
        print(f"⚠️ Role {class_iri} has no associated leaves in roles_to_leaves_map.")

    n_leaves = len(leaves)

    return leaves, n_leaves

if __name__ == "__main__":

    """ Select task to perform. To calculate the removed leaf classes for a given class, 
    the file with the class to leaf descendants map must be created first.
    This is done in "build_class_to_leaf_map" task. """

    task = "other" 
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
    class_to_all_roles_json = "data/class_to_all_roles_map.json"
    roles_to_leaves_map_json = "data/roles_to_leaves_map.json"

    with open(class_to_all_roles_json, "r") as f:
        class_to_all_roles_map = json.load(f)
    with open(roles_to_leaves_map_json, "r") as f:
        roles_to_leaves_map = json.load(f)

    if task == "count_total_removed_leaves":

        n_removed_leaves = count_removed_leaves(removed_leaves_csv)
        print(f"Total number of removed {classification} leaf classes: {n_removed_leaves}")

        # Output
        # - Total number of removed structural leaf classes: 184 436
        # - Total number of removed functional leaf classes: 41

    elif task == "count_removed_classes_for_class":

        subclasses, n_subclasses = count_removed_classes_for_class(
            class_iri,
            map_file,
            classification,
            class_to_all_roles_map,
            roles_to_leaves_map,
        )
        print(f"Class {class_iri} has {n_subclasses} removed subclasses in the ontology.")

        if subclasses and n_subclasses < 50: 
            for sub in subclasses:
                print(f" - {sub}")

    elif task == "build_class_to_leaf_map":

        build_class_to_leaf_map(leaf_to_ancestors_file, class_to_leaf_output_file)

    elif task == "enrichment_analysis_plain": # to be removed

        _, n_subclasses = count_removed_classes_for_class(
            class_iri,
            map_file,
            classification,
            class_to_all_roles_map,
            roles_to_leaves_map,
        )
        n_tot_removed_leaves = count_removed_leaves(removed_leaves_csv)

        print(f"Calulated for {classification} ontology for {class_iri} (make sure classification is correct):")
        prob_of_success = n_subclasses / n_tot_removed_leaves

        print(f" The probability of success value for class {class_iri} is: {prob_of_success:.4g} ({n_subclasses} / {n_tot_removed_leaves})")



    else:
        print("No valid task selected.")

        