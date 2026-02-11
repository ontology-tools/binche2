# Fucntion to check subclasses of a given class through different json files
import pandas as pd
import json
from pre_fishers_calculations import count_removed_classes_for_class

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

# Check if class is in removed classes

def is_class_removed(class_iri):
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
    removed_classes = pd.read_csv(removed_leaves_csv)
    is_deleted = class_iri in removed_classes['IRI'].values
    print(f"Class {class_iri} removed: {is_deleted}")
    return is_deleted

# Check the leaves of a given class
def check_leaves_of_class(class_iri):
    leaf_descendants_json = 'data/class_to_leaf_descendants_map.json'
    with open(leaf_descendants_json, 'r') as f:
        leaf_descendants_map = json.load(f)
    leaf_descendants = leaf_descendants_map.get(class_iri, [])
    print(f"Class {class_iri} has {len(leaf_descendants)} leaf descendants:")
    for leaf in leaf_descendants:
        print(f" - {leaf}")
    return leaf_descendants

def check_direct_parents_of_class(class_iri):
    parent_map_json = 'data/chebi_parent_map.json'
    with open(parent_map_json, 'r') as f:
        parent_map = json.load(f)
    direct_parents = parent_map.get(class_iri, [])
    print(f"Class {class_iri} has {len(direct_parents)} direct parents:")
    for parent in direct_parents:
        print(f" - {parent}")
    return direct_parents

# Count how many of the removed leaf classes has classification "structure" and how many has "function"
def count_removed_leaves():
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
    removed_classes = pd.read_csv(removed_leaves_csv)
    structure_count = len(removed_classes[removed_classes['Classification'] == 'structural'])
    function_count = len(removed_classes[removed_classes['Classification'] == 'functional'])
    print(f"Total removed leaf classes: {len(removed_classes)}")
    print(f" - Structural (Struture): {structure_count}")
    print(f" - Functional (Role): {function_count}")
    # Print all leaf classes with a functional classification
    print("Functional leaf classes:")
    functional_classes = removed_classes[removed_classes['Classification'] == 'functional']
    for index, row in functional_classes.iterrows():
        print(f" - {row['IRI']}")
    return structure_count, function_count

def find_roles_for_leaf(leaf_iri):
    leaf_roles_json = 'data/removed_leaf_classes_to_ALL_roles_map.json'
    with open(leaf_roles_json, 'r') as f:
        leaf_roles_map = json.load(f)
    # Check if class exists in map
    if leaf_iri not in leaf_roles_map:
        print(f"Leaf class {leaf_iri} not found in roles map.")
        return []
    roles = leaf_roles_map.get(leaf_iri, [])
    print(f"Leaf class {leaf_iri} has {len(roles)} associated roles:")
    for role in roles:
        print(f" - {role}")
    return roles

def find_direct_roles_for_any_class(class_iri):
    roles_json = 'data/class_to_direct_roles_map.json'
    with open(roles_json, 'r') as f:
        roles_map = json.load(f)
    # Check if class exists in map
    if class_iri not in roles_map:
        print(f"Class {class_iri} not found in roles map.")
        return []
    roles = roles_map.get(class_iri, [])
    print(f"Class {class_iri} has {len(roles)} associated direct roles:")
    for role in roles:
        print(f" - {role}")
    return roles

def find_all_roles_for_any_class(class_iri):
    roles_json = 'data/class_to_all_roles_map.json'
    with open(roles_json, 'r') as f:
        roles_map = json.load(f)
    # Check if class exists in map
    if class_iri not in roles_map:
        print(f"Class {class_iri} not found in all roles map.")
        return []
    roles = roles_map.get(class_iri, [])
    print(f"Class {class_iri} has {len(roles)} associated direct and inherited roles:")
    return roles

if __name__ == "__main__":

    # class_iri = "http://purl.obolibrary.org/obo/CHEBI_33917"  # Example class IRI
    # count_connected_classes(class_iri)

    
    class_iri = "http://purl.obolibrary.org/obo/CHEBI_27732"  
    # is_class_removed(class_iri)
    # check_leaves_of_class(class_iri)

    # check_direct_parents_of_class(class_iri)
    # count_removed_leaves()

    # find_roles_for_leaf(class_iri)

    #  find_direct_roles_for_any_class(class_iri)

    class_iri = "http://purl.obolibrary.org/obo/CHEBI_25899"

    leaf_descendants_map_file = "data/class_to_leaf_descendants_map.json"
    with open(leaf_descendants_map_file, 'r') as f:
        leaf_descendants_map = json.load(f)

    classification = "functional" # "functional" or "structural" or "full"
    class_to_all_roles_json = "data/class_to_all_roles_map.json"
    with open(class_to_all_roles_json, "r") as f:
        class_to_all_roles_map = json.load(f)
    roles_to_leaves_map_json = "data/roles_to_leaves_map.json"
    with open(roles_to_leaves_map_json, "r") as f:
        roles_to_leaves_map = json.load(f)

    # count_removed_classes_for_class(class_iri, leaf_descendants_map, classification, class_to_all_roles_map, roles_to_leaves_map)

    find_all_roles_for_any_class(class_iri)