# Fucntion to check subclasses of a given class through different json files
import pandas as pd
import json

def count_connected_classes(class_iri):
    chebi_subclass_json = 'data/chebi_subclass_map.json'
    chebi_leafdescendant_json = 'data/class_to_leaf_descendants_map.json'
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"

    with open(chebi_subclass_json, 'r') as f:
        subclass_map = json.load(f)

    with open(chebi_leafdescendant_json, 'r') as f:
        leaf_descendants_map = json.load(f)
        
    removed_classes = pd.read_csv(removed_leaves_csv)

    # find all subclasses of the given class
    subclasses = subclass_map.get(class_iri, [])
    leaf_descendants = leaf_descendants_map.get(class_iri, [])

    print(f"Class {class_iri} has {len(subclasses)} subclasses and {len(leaf_descendants)} leaf descendants.")

    subclasses_in_removed = 0
    leafs_in_removed = 0

    for this_class in subclasses:
        # check if this_class is in removed_classes
        if this_class in removed_classes['IRI'].values:
            subclasses_in_removed += 1

    for this_leaf in leaf_descendants:
        if this_leaf in removed_classes['IRI'].values:
            leafs_in_removed += 1

    print(f"Class {class_iri}")

    print("Subclasses:")
    for sub in subclasses:
        print(f" - {sub}")  
    
    print("Leaf Descendants:")
    for leaf in leaf_descendants:
        print(f" - {leaf}")

    print(f"Class {class_iri} has {len(subclasses)} subclasses and {len(leaf_descendants)} leaf descendants.")
    print(f"Out of these, {subclasses_in_removed} subclasses and {leafs_in_removed} leaf descendants are in the removed classes.")
    
    return subclasses, len(subclasses)


if __name__ == "__main__":

    class_iri = "http://purl.obolibrary.org/obo/CHEBI_33917"  # Example class IRI
    count_connected_classes(class_iri)