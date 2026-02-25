import os
from pathlib import Path
from load_chebi import load_chebi, load_ontology
from pruning_smiles import (
    find_leaf_classes_with_smiles_and_deprecated,
    save_filtered_owl,
    save_leaf_classes_with_smiles,
    build_parent_map,
    map_names_to_classes,
)
from prepare_role_calculations import (
    find_has_role_connections_from_owl,
    create_leaves_to_all_roles_map,
    create_roles_to_all_leaves_map,
    create_class_to_all_roles_map,
)
from pruning_split_up_structure import identify_structural_vs_functional, split_owl_by_type
from pre_fishers_calculations import build_class_to_leaf_map
import time

start_time = time.time()

"""
A script for creating all necessary files for the project.
"""

def rename_folder(old_name, new_name):
    """
    Rename a folder from old_name to new_name.
    
    Args:
        old_name (str): Current folder name
        new_name (str): New folder name
    
    Returns:
        bool: True if rename was successful, False otherwise
    """
    try:
        if os.path.exists(old_name):
            Path(old_name).rename(new_name)
            print(f"Folder renamed: '{old_name}' -> '{new_name}'")
            return True
        else:
            print(f"Folder '{old_name}' does not exist")
            return False
    except Exception as e:
        print(f"Error renaming folder: {e}")
        return False

def create_temp_data_folder():
    """
    Create a temporary 'data_new' folder for new files.
    """
    if os.path.exists("data_new"):
        print("Warning: 'data_new' already exists. Removing it.")
        import shutil
        shutil.rmtree("data_new")
    os.makedirs("data_new")
    print("Created temporary 'data_new' folder")

def finalize_folder_structure():
    """
    After all files are created, rename folders:
    1. data -> data_last_used_YYYY.MM.DD (or with counter if needed)
    2. data_new -> data
    """
    import datetime
    
    timestamp = datetime.datetime.now().strftime("%Y.%m.%d")
    
    # Rename old data folder if it exists
    if os.path.exists("data"):
        base_name = f"data_last_used_{timestamp}"
        old_data_name = base_name
        counter = 1
        while os.path.exists(old_data_name):
            old_data_name = f"{base_name}_{counter}"
            counter += 1
        rename_folder("data", old_data_name)
    
    # Rename data_new to data
    if os.path.exists("data_new"):
        rename_folder("data_new", "data")

if __name__ == "__main__":
    ### Create temporary data folder for new files
    print("Creating temporary data folder...")
    create_temp_data_folder()

    ### Define properties and file paths
    smiles_property = "https://w3id.org/chemrof/smiles_string"
    deprecated_property = "http://www.w3.org/2002/07/owl#deprecated"
    has_role_property = "http://purl.obolibrary.org/obo/RO_0000087"

    chebi_file = "data_new/chebi.owl"
    subclass_map_file = "data_new/chebi_subclass_map.json" 
    leaf_parents_map_file = "data_new/removed_leaf_classes_to_ALL_parents_map.json"
    removed_leaf_classes_file = "data_new/removed_leaf_classes_with_smiles.csv"

    roles_map_json = "data_new/class_to_direct_roles_map.json" # output file for roles map
    leaves_to_all_parents_json = "data_new/removed_leaf_classes_to_ALL_parents_map.json"
    leaves_to_all_roles_json = "data_new/removed_leaf_classes_to_ALL_roles_map.json" # output file for leaves to all roles map
    parent_map_json = "data_new/chebi_parent_map.json"
    roles_to_all_leaves_json = "data_new/roles_to_leaves_map.json" # output file for roles to all leaves map
    class_to_all_roles_json = "data_new/class_to_all_roles_map.json" # output file for class to all roles map
    url_to_id_map_json = "data_new/chebi_url_to_id_map.json"

    filtered_chebi_no_leaves_file = "data_new/filtered_chebi_no_leaves_with_smiles_no_deprecated.owl"
    structural_ontology_file = "data_new/filtered_chebi_no_leaves_with_smiles_no_deprecated_structural.owl"
    functional_ontology_file = "data_new/filtered_chebi_no_leaves_with_smiles_no_deprecated_functional.owl"

    leaf_to_ancestors_file = "data_new/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_output_file = "data_new/class_to_leaf_descendants_map.json"

    ### Download and load the ChEBI ontology
    print("Downloading and loading ChEBI ontology...")
    # Option: download directly to data_new/ to keep separate from old versions
    chebi_ontology = load_chebi(download_dir="data_new")

    ### Find and filter leaf classes with SMILES and deprecated classes, and save the filtered ontology
    print("Removing leaf classes with SMILES and deprecated classes...")
    classes_with_smiles, deprecated_classes = find_leaf_classes_with_smiles_and_deprecated(
        chebi_ontology,
        smiles_property,
        deprecated_property,
        subclass_map_file,
        leaf_parents_map_file,
        use_found_leaf_classes=False,
        removed_leaf_classes_file=removed_leaf_classes_file,
    )

    # Comment out if you do not want to save the filtered OWL in a new file
    save_filtered_owl(chebi_file, classes_with_smiles, deprecated_classes, filtered_chebi_no_leaves_file)

    ### Build map of all classes to their direct parents and save as JSON
    print("Building parent map...")
    build_parent_map(chebi_ontology, parent_map_json, deprecated_property)

    ### Build map with CHEBI urls to short CHEBI IDs and save as JSON
    print("Building URL to ID map...")
    map_names_to_classes(chebi_file, url_to_id_map_json)

    ### Map all classes to their direct "has_role" connections and save as JSON,
    ### then create a map of all leaf classes to all their "has_role" connections (not just direct ones) and save as JSON
    ### and the reverse map of all roles to all leaf classes that have that role (directly or indirectly) and save as JSON
    print("Building roles maps...")
    roles_map = find_has_role_connections_from_owl(chebi_file, has_role_property, deprecated_property, roles_map_json)
    create_leaves_to_all_roles_map(roles_map_json, leaves_to_all_parents_json, leaves_to_all_roles_json, parent_map_json)
    create_roles_to_all_leaves_map(leaves_to_all_roles_json, roles_to_all_leaves_json)

    ### Map each class to its direct roles, roles inherited from ancestor classes, and descendants of those roles in the role hierarchy.
    create_class_to_all_roles_map(roles_map_json, parent_map_json, class_to_all_roles_json)

    ### Split the filtered ontology into structural and functional classes and save as separate OWL files
    print("Splitting ontology into structural and functional classes...")
    structural_classes, functional_classes, unknown_classes = identify_structural_vs_functional(chebi_ontology)
    split_owl_by_type(structural_classes, functional_classes, unknown_classes, filtered_chebi_no_leaves_file, output_dir="data_new")

    ### Save a CSV file with the removed leaf classes 
    print("Saving CSV with removed leaf classes...")
    structural_ontology = load_ontology(structural_ontology_file)
    functional_ontology = load_ontology(functional_ontology_file)

    # Extract sets of class IRIs for quick membership checking
    structural_classes = set(structural_ontology.get_classes())
    functional_classes = set(functional_ontology.get_classes())

    classes_with_smiles, _ = find_leaf_classes_with_smiles_and_deprecated(
        chebi_ontology,
        smiles_property,
        deprecated_property,
        subclass_map_file,
        leaf_parents_map_file,
    )

    # Save CSV with classification
    save_leaf_classes_with_smiles(
        classes_with_smiles,
        chebi_ontology,
        smiles_property,
        removed_leaf_classes_file,
        structural_classes,
        functional_classes,
    )  

    ### Build a JSON map from each class IRI to ALL its leaf descendants
    print("Building class to leaf descendants map...")
    build_class_to_leaf_map(leaf_to_ancestors_file, class_to_leaf_output_file)

    ### Finalize folder structure: rename old data and move data_new to data
    print("Finalizing folder structure...")
    finalize_folder_structure()

    print("All files created successfully!")
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Total execution time: {elapsed_time:.2f} seconds or {elapsed_time/60:.2f} minutes or {elapsed_time/3600:.2f} hours")














