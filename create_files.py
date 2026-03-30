import os
from pathlib import Path
from load_chebi import load_chebi, load_ontology
from pruning_smiles import (
    find_leaf_classes_with_smiles_and_deprecated,
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
from pruning_split_up_structure import identify_structural_vs_functional
from pre_fishers_calculations import build_class_to_leaf_map
import time

start_time = time.time()

"""
A script for creating all necessary files for the project.
"""


def _run_stage(stage_name, stage_timings, func, *args, **kwargs):
    """Run a stage, record elapsed time, and print a compact timing line."""
    stage_start = time.time()
    result = func(*args, **kwargs)
    elapsed = time.time() - stage_start
    stage_timings.append((stage_name, elapsed))
    print(f"[TIMING] {stage_name}: {elapsed:.2f}s ({elapsed/60:.2f} min)")
    return result


def _print_timing_summary(stage_timings, total_elapsed):
    """Print stage timings and percentages of the total runtime."""
    print("\n" + "=" * 80)
    print("PIPELINE TIMING SUMMARY")
    print("=" * 80)

    if total_elapsed <= 0:
        total_elapsed = 1e-9

    sorted_timings = sorted(stage_timings, key=lambda x: x[1], reverse=True)
    for stage_name, elapsed in sorted_timings:
        pct = (elapsed / total_elapsed) * 100
        print(
            f"{stage_name}: {elapsed:.2f}s ({elapsed/60:.2f} min, {pct:.2f}% of total)"
        )

    measured_total = sum(elapsed for _, elapsed in stage_timings)
    overhead = total_elapsed - measured_total
    overhead_pct = (overhead / total_elapsed) * 100
    print("-" * 80)
    print(
        f"Untracked/overhead: {overhead:.2f}s ({overhead/60:.2f} min, {overhead_pct:.2f}% of total)"
    )
    print("=" * 80)

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
    stage_timings = []

    ### Create temporary data folder for new files
    print("Creating temporary data folder...")
    _run_stage("create_temp_data_folder", stage_timings, create_temp_data_folder)

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
    id_to_name_map_json = "data_new/chebi_id_to_name_map.json"

    leaf_to_ancestors_file = "data_new/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_output_file = "data_new/class_to_leaf_descendants_map.json"

    ### Download and load the ChEBI ontology
    print("Downloading and loading ChEBI ontology...")
    # Option: download directly to data_new/ to keep separate from old versions
    chebi_ontology = _run_stage(
        "download_and_load_chebi",
        stage_timings,
        load_chebi,
        download_dir="data_new",
    )

    ### Find and filter leaf classes with SMILES and deprecated classes, and save the filtered ontology
    print("Removing leaf classes with SMILES and deprecated classes...")
    classes_with_smiles, deprecated_classes = _run_stage(
        "find_leaf_classes_with_smiles_and_deprecated",
        stage_timings,
        find_leaf_classes_with_smiles_and_deprecated,
        chebi_ontology,
        smiles_property,
        deprecated_property,
        subclass_map_file,
        leaf_parents_map_file,
        use_found_leaf_classes=False,
        removed_leaf_classes_file=removed_leaf_classes_file,
    )

    ### Build map of all classes to their direct parents and save as JSON
    print("Building parent map...")
    _run_stage(
        "build_parent_map",
        stage_timings,
        build_parent_map,
        chebi_ontology,
        parent_map_json,
        deprecated_property,
        subclass_map_file=subclass_map_file,
        precomputed_deprecated_classes=deprecated_classes,
    )

    ### Build map with CHEBI short IDs to names and save as JSON
    print("Building ID to name map...")
    _run_stage(
        "map_names_to_classes",
        stage_timings,
        map_names_to_classes,
        chebi_file,
        id_to_name_map_json,
    )

    ### Map all classes to their direct "has_role" connections and save as JSON,
    ### then create a map of all leaf classes to all their "has_role" connections (not just direct ones) and save as JSON
    ### and the reverse map of all roles to all leaf classes that have that role (directly or indirectly) and save as JSON
    print("Building roles maps...")
    roles_map = _run_stage(
        "find_has_role_connections_from_owl",
        stage_timings,
        find_has_role_connections_from_owl,
        chebi_file,
        has_role_property,
        deprecated_property,
        roles_map_json,
    )
    _run_stage(
        "create_leaves_to_all_roles_map",
        stage_timings,
        create_leaves_to_all_roles_map,
        roles_map_json,
        leaves_to_all_parents_json,
        leaves_to_all_roles_json,
        parent_map_json,
    )
    _run_stage(
        "create_roles_to_all_leaves_map",
        stage_timings,
        create_roles_to_all_leaves_map,
        leaves_to_all_roles_json,
        roles_to_all_leaves_json,
    )

    ### Map each class to its direct roles, roles inherited from ancestor classes, and descendants of those roles in the role hierarchy.
    _run_stage(
        "create_class_to_all_roles_map",
        stage_timings,
        create_class_to_all_roles_map,
        roles_map_json,
        parent_map_json,
        class_to_all_roles_json,
    )

    ### Identify structural vs functional classes (in-memory, no OWL files needed)
    print("Identifying structural vs functional classes...")
    structural_classes, functional_classes, unknown_classes = _run_stage(
        "identify_structural_vs_functional",
        stage_timings,
        identify_structural_vs_functional,
        chebi_ontology,
    )
    
    print(f"Identified {len(structural_classes)} structural classes")
    print(f"Identified {len(functional_classes)} functional classes")
    print(f"Identified {len(unknown_classes)} unknown classes")

    ### Save a CSV file with the removed leaf classes (using in-memory class sets)
    print("Saving CSV with removed leaf classes...")
    _run_stage(
        "save_leaf_classes_with_smiles",
        stage_timings,
        save_leaf_classes_with_smiles,
        classes_with_smiles,
        chebi_ontology,
        smiles_property,
        removed_leaf_classes_file,
        structural_classes,
        functional_classes,
        owl_file=chebi_file,  # OPTIMIZATION: Pass OWL file for fast XML parsing of SMILES
        parent_map_file=parent_map_json,  # OPTIMIZATION: Use already-created parent map instead of ontology API
    )  

    ### Build a JSON map from each class IRI to ALL its leaf descendants
    print("Building class to leaf descendants map...")
    _run_stage(
        "build_class_to_leaf_map",
        stage_timings,
        build_class_to_leaf_map,
        leaf_to_ancestors_file,
        class_to_leaf_output_file,
    )

    ### Finalize folder structure: rename old data and move data_new to data
    print("Finalizing folder structure...")
    _run_stage("finalize_folder_structure", stage_timings, finalize_folder_structure)

    print("All files created successfully!")
    end_time = time.time()
    elapsed_time = end_time - start_time
    _print_timing_summary(stage_timings, elapsed_time)
    print(f"Total execution time: {elapsed_time:.2f} seconds or {elapsed_time/60:.2f} minutes or {elapsed_time/3600:.2f} hours")














