import os
import shutil
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

from wikidata.get_wikidata_lotus import create_wikidata_output_files, keep_taxon_compounds
from wikidata.get_inchikeys import convert_smiles_file
from hmdb.extract_hmdb import extract_hmdb_to_file
from hmdb.filter_hmdb_statuses import filter_hmdb_statuses_main
from wikidata.find_missing_chebis import run_find_missing_chebis
from wikidata.combine_human_datasets import combine_datasets
from wikidata.narrow_background_fishers import gather_narrow_leaves
from BiGG.get_model import download_model_json, gather_recon3d_leaves



import time
import json

start_time = time.time()

"""
A script for creating all necessary files for the project. 
The old data is moved to a timestamped folder.
If there are more than three folders with old data, the oldest ones are deleted to keep only the three most recent.
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

def cleanup_old_data_folders(max_folders=3):
    """
    Keep only the most recent 'max_folders' of old data folders and delete the rest.

    Args:
        max_folders (int): Maximum number of old data folders to keep
    """
    import re

    # Match folders named "data_last_used_YYYY.MM.DD" or "data_last_used_YYYY.MM.DD_N"
    pattern = re.compile(r"data_last_used_(\d{4}\.\d{2}\.\d{2})(?:_(\d+))?$")

    old_data_folders = []
    for f in os.listdir("."):
        if not os.path.isdir(f):
            continue
        match = pattern.match(f)
        if match:
            date_str, counter = match.groups()
            old_data_folders.append((f, date_str, int(counter) if counter else 0))

    # Sort by the date (and counter, for same-day folders) encoded in the name, newest first
    old_data_folders.sort(key=lambda item: (item[1], item[2]), reverse=True)

    # Keep only the most recent 'max_folders'
    folders_to_delete = old_data_folders[max_folders:]

    for folder, _, _ in folders_to_delete:
        print(f"Deleting old data folder: {folder}")
        shutil.rmtree(folder)

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

    ### Identify structural vs functional classes (fast path: traverse the precomputed
    ### subclass map instead of the slow per-node ontology API). The flattened map written by
    ### find_leaf_classes gives identical descendant sets for the (non-leaf) roots.
    print("Identifying structural vs functional classes...")
    with open(subclass_map_file, "r") as f:
        subclass_map_for_classification = json.load(f)
    structural_classes, functional_classes, unknown_classes = _run_stage(
        "identify_structural_vs_functional",
        stage_timings,
        identify_structural_vs_functional,
        chebi_ontology,
        subclass_map_for_classification,
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

    # Copy HMDB XML to data_new before folder rename
    print("Copying HMDB XML to temporary data folder...")
    hmdb_xml_src = "data/hmdb_metabolites.xml"
    hmdb_xml_dst = "data_new/hmdb_metabolites.xml"
    print(f"Copying HMDB XML to data_new...")
    if os.path.exists(hmdb_xml_src):
        shutil.copy2(hmdb_xml_src, hmdb_xml_dst)
        print(f"Copied {hmdb_xml_src} to {hmdb_xml_dst}")
    else:
        raise FileNotFoundError(f"{hmdb_xml_src} not found. Please download it before running.")
    

    ### Finalize folder structure: rename old data and move data_new to data
    print("Finalizing folder structure...")
    _run_stage("finalize_folder_structure", stage_timings, finalize_folder_structure)

    ### Convert SMILES to InChIKeys for the removed leaf classes (needed by find_missing_chebis and the website)
    print("Generating InChIKeys for removed leaf classes...")
    _run_stage(
        "convert_smiles_to_inchikeys",
        stage_timings,
        convert_smiles_file,
        "data/removed_leaf_classes_with_smiles.csv",
        "data/removed_leaf_classes_with_inchikeys.csv",
    )

    print("Creating Wikidata output files...")
    _run_stage(
        "create_wikidata_output_files",
        stage_timings,
        create_wikidata_output_files,
        include_download=True,
    )

    print("Filtering Wikidata compounds for Arabidopsis thaliana...")
    _run_stage(
        "filter_wikidata_arabidopsis_thaliana",
        stage_timings,
        keep_taxon_compounds,
        "data/wikidata/created/compounds_with_chebi_ids.tsv",
        "data/wikidata/created/compounds_with_chebi_ids_arabidopsis_thaliana.tsv",
        "arabidopsis_thaliana",
    )

    print("Creating HMDB output file...")
    _run_stage(
        "create_hmdb_output_file",
        stage_timings,
        extract_hmdb_to_file,
        "data/hmdb_metabolites.xml",
        "data/hmdb_metabolites_extract.tsv",
    )

    print("Filtering HMDB statuses...")
    _run_stage(
        "filter_hmdb_statuses",
        stage_timings,
        filter_hmdb_statuses_main,
    )

    print("Finding missing ChEBI IDs for Wikidata (Homo sapiens)...")
    _run_stage(
        "find_missing_chebis_wikidata_hs",
        stage_timings,
        run_find_missing_chebis,
        "wikidata_hs",
    )

    print("Finding missing ChEBI IDs for Wikidata (Arabidopsis thaliana)...")
    _run_stage(
        "find_missing_chebis_wikidata_at",
        stage_timings,
        run_find_missing_chebis,
        "wikidata_at",
    )

    print("Finding missing ChEBI IDs for HMDB...")
    _run_stage(
        "find_missing_chebis_hmdb",
        stage_timings,
        run_find_missing_chebis,
        "hmdb",
    )

    print("Combining human datasets...")
    _run_stage(
        "combine_datasets",
        stage_timings,
        combine_datasets,
    )

    print("Gathering narrow leaf classes (Homo sapiens)...")
    _run_stage(
        "gather_narrow_leaves_homo_sapiens",
        stage_timings,
        gather_narrow_leaves,
        compounds_tsv="data/combined_hmdb_wikidata.tsv",
        leaves_csv="data/removed_leaf_classes_with_smiles.csv",
        class_to_leaf_map="data/class_to_leaf_descendants_map.json",
        output_json="data/human_entities_leaves.json",
        taxon_label="homo_sapiens",
    )

    print("Gathering narrow leaf classes (Arabidopsis thaliana)...")
    _run_stage(
        "gather_narrow_leaves_arabidopsis_thaliana",
        stage_timings,
        gather_narrow_leaves,
        compounds_tsv="data/wikidata/created/compounds_with_chebi_ids_arabidopsis_thaliana_updatedchebis.tsv",
        leaves_csv="data/removed_leaf_classes_with_smiles.csv",
        class_to_leaf_map="data/class_to_leaf_descendants_map.json",
        output_json="data/arabidopsis_thaliana_leaves.json",
        taxon_label="arabidopsis_thaliana",
    )

    print("Downloading Recon3D model...")
    _run_stage(
        "download_recon3d_model",
        stage_timings,
        download_model_json,
        "Recon3D",
        "data/Recon3D.json",
    )

    print("Gathering narrow leaf classes (Endogenous human / Recon3D)...")
    _run_stage(
        "gather_narrow_leaves_endogenous_human",
        stage_timings,
        gather_recon3d_leaves,
        model_json_path="data/Recon3D.json",
        leaves_csv="data/removed_leaf_classes_with_smiles.csv",
        class_to_leaf_map_json="data/class_to_leaf_descendants_map.json",
        output_json="data/recon3d_leaves.json",
    )

    ### Clean up old data folders, keeping only the most recent ones
    print("Cleaning up old data folders...")
    _run_stage("cleanup_old_data_folders", stage_timings, cleanup_old_data_folders)

    print("All files created successfully!")
    end_time = time.time()
    elapsed_time = end_time - start_time
    _print_timing_summary(stage_timings, elapsed_time)
    print(f"Total execution time: {elapsed_time:.2f} seconds or {elapsed_time/60:.2f} minutes or {elapsed_time/3600:.2f} hours")














