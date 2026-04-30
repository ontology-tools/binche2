""" Enrichment script for narrow compounds """

import csv
import json
import os
import re
import sys
import pandas as pd
from pathlib import Path
import time

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from fishers_calculations import normalize_id, calculate_p_value, get_ancestors_for_inputs, get_n_ss_annotated, print_enrichment_results
from visualitations_and_pruning import root_children_pruner, linear_branch_collapser_pruner_remove_less, high_p_value_branch_pruner, zero_degree_pruner, create_graph_from_map, id_to_name, create_graph_with_roles_and_structures
from pre_fishers_calculations import count_removed_classes_for_class, count_removed_classes_for_roles
from multiple_test_corrections import bonferroni_correction, benjamini_hochberg_fdr_correction

# n_ss_leaves = total number of input classes in the study set (if they are all leaf classes, otherwise count corresponding leaf classes)
# n_ss_annotated = number of the input classes that are descendants of the given class
# n_bg_leaves = total number of leaf classes connected to the narrow set (human entities in this case)   
# n_bg_annotated = number of those leaf classes that are descendants of the given class)

# For the background, we want to include all classes that are connected to the narrow set (human entitites)
# If the CHEBI_Id is a leaf, we include it (and count is as 1 leaf), and all of its ancetsors.
# If the CHEBI_Id is not a leaf, we include all of its leaf descendants, and all of its ancestors.

def _calculate_depth_to_root(chebi_iri, parent_map, memo=None):
    """Recursively calculate depth (path length) from a node to the root.
    
    Uses memoization to avoid recalculating already-seen nodes.
    If multiple parents exist, returns the maximum depth among them.
    """
    if memo is None:
        memo = {}
    if chebi_iri in memo:
        return memo[chebi_iri]
    
    parents = parent_map.get(chebi_iri, [])
    if not parents:
        memo[chebi_iri] = 0
        return 0
    
    # depth = 1 + maximum depth among all parents
    depth = 1 + max(_calculate_depth_to_root(p, parent_map, memo) for p in parents)
    memo[chebi_iri] = depth
    return depth


def filter_chebifier_parents(parent_chebis, chebi_parent_map_json):
    """If Chebifier finds several parent classes for a given class (that does not have its own CHEBI ID),
    we want to only keep the parent(s) furthest down the hierarchy, i.e. the one(s) with the longest path to the root. 
    This is to avoid inflating the background with very high-level classes.

    Input: list of parent CHEBI IDs (strings)
    Output: list of parent CHEBI IDs (strings) that are furthest down the hierarchy. 
            should only be one parent unless there is a tie
    """

    # open json file with parent-child relationships in CHEBI
    with open(chebi_parent_map_json, "r", encoding="utf-8") as f:
        chebi_parent_map = json.load(f)

    # calculate path length to root for each parent CHEBI ID using memoization
    memo = {}
    parent_path_lengths = {p: _calculate_depth_to_root(p, chebi_parent_map, memo) for p in parent_chebis}

    # find the maximum path length
    max_path_length = max(parent_path_lengths.values())
    
    # choose the parent(s) with the maximum path length
    chosen_parents = [parent for parent, path_length in parent_path_lengths.items() if path_length == max_path_length]

    return chosen_parents

def gather_narrow_leaves(human_entities_tsv, leaves_csv, class_to_leaf_map, output_json, chebi_parent_map_json):
    """Build narrow background leaves from human entities TSV.

    Workflow:
    - Read all ChEBI IDs from `human_entities_tsv` (supports ';' or '|' separators).
    - For entities with multiple ChEBI IDs, filter to keep only the deepest ones in hierarchy.
    - If an ID is already a leaf (present in `leaves_csv` IRI column), keep it.
    - Otherwise, expand to all leaf descendants using `class_to_leaf_map`.
    - Save the resulting unique leaf IRIs to `output_json`.
    """

    # Reuse cached leaves when available.
    # if os.path.exists(output_json):
        # with open(output_json, "r", encoding="utf-8") as f:
          #  cached = json.load(f)
        # leaves = set(cached.get("narrow_leaves", []))
       # print(f"Loaded {len(leaves)} narrow leaves from {output_json}")
       # return leaves

    print(f"Building narrow leaves cache from {human_entities_tsv}...")
    print(f"Using leaf membership from {leaves_csv}")

    # Accept either a ready dict or a JSON path for class->leaf mapping.
    if isinstance(class_to_leaf_map, str):
        with open(class_to_leaf_map, "r", encoding="utf-8") as f:
            class_to_leaf_map = json.load(f)

    # Fast leaf lookup set from the existing removed-leaves CSV.
    leaf_iris = set()
    with open(leaves_csv, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if "IRI" not in (reader.fieldnames or []):
            raise KeyError(f"Expected 'IRI' column in {leaves_csv}")
        for row in reader:
            iri = (row.get("IRI") or "").strip()
            if iri:
                leaf_iris.add(iri)

    # Collect all human CHEBI seeds and expand non-leaf seeds to descendants.
    narrow_leaves = set()
    seeds_seen = 0

    with open(human_entities_tsv, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "chebi_id" not in (reader.fieldnames or []):
            raise KeyError(f"Expected 'chebi_id' column in {human_entities_tsv}")

        print(f"Found columns in human TSV: {reader.fieldnames}")

        for row in reader:
            raw_chebi = (row.get("chebi_id") or "").strip()
            if not raw_chebi:
                continue

            # Collect all ChEBI IDs for this entity and normalize them.
            tokens = [t.strip() for t in re.split(r"[|;]", raw_chebi) if t.strip()]
            seed_iris = [normalize_id(t) for t in tokens]

            # If multiple ChEBI IDs, filter to keep only the deepest ones in hierarchy.
            if len(seed_iris) > 1:
                seed_iris = filter_chebifier_parents(seed_iris, chebi_parent_map_json)
                print(f"Entity has {len(tokens)} ChEBI IDs; filtered to {len(seed_iris)} deepest: {seed_iris}")

            # Process each (filtered) ChEBI ID.
            for seed_iri in seed_iris:
                seeds_seen += 1

                if seed_iri in leaf_iris:
                    narrow_leaves.add(seed_iri)
                    continue

                descendants = class_to_leaf_map.get(seed_iri, [])

                # If a class has over 150 leaf decsendants, do not keep any of them, 
                # to avoid inflating the background with high-level classes.

                if len(descendants) > 150:
                    print(f"Class {seed_iri} has {len(descendants)} leaf descendants, which exceeds the threshold. Skipping expansion to avoid inflating background.")
                
                else:
                    narrow_leaves.update(descendants)

    print(f"Processed {seeds_seen} human ChEBI seed IDs")

    payload = {
        "human_entities_tsv": human_entities_tsv,
        "leaves_csv": leaves_csv,
        "n_seeds_processed": seeds_seen,
        "narrow_leaves": sorted(narrow_leaves),
    }

    os.makedirs(os.path.dirname(output_json) or ".", exist_ok=True)
    with open(output_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"Saved {len(narrow_leaves)} narrow leaves to {output_json}")
    return narrow_leaves

"""Enrichment Calcultaions"""

def get_studyset_leaves_narrow(studyset_list, human_entities_leaves_json, leaves_csv, class_to_leaf_map):
    """Get leaf descendants for the input classes and track leaves outside the narrow background.

    `human_entities_leaves_json` provides the narrow background leaf set.
    `leaves_csv` is used to recognize leaves in the broader CHEBI leaf universe.
    """

    # Normalize study set IDs to full IRIs
    studyset_list = [normalize_id(cls) for cls in studyset_list]
    if isinstance(class_to_leaf_map, str):
        with open(class_to_leaf_map, "r", encoding="utf-8") as f:
            class_to_leaf_map = json.load(f)

    with open(human_entities_leaves_json, "r", encoding="utf-8") as f:
        cached = json.load(f)
    human_leaves = set(cached.get("narrow_leaves", []))

    all_leaves_df = pd.read_csv(leaves_csv)
    all_leaves = set(all_leaves_df["IRI"].values)

    studyset_leaves = set()
    leaves_to_expand_background = set()
    parents_to_expand_background = set()

    for cls in studyset_list:
        print(f"Processing class {cls}...")

        if cls in human_leaves:
            print(f"Class {cls} is already in the narrow background.")
            studyset_leaves.add(cls)
        elif cls in all_leaves:
            print(f"Class {cls} is a leaf in the broader background, adding it.")
            studyset_leaves.add(cls)
            leaves_to_expand_background.add(cls)
        else:
            leaf_descendants = set(class_to_leaf_map.get(cls, []))
            print(f"Class {cls} is not a leaf, adding its {len(leaf_descendants)} leaf descendants.")
            studyset_leaves.update(leaf_descendants)

            for leaf in leaf_descendants:
                if leaf not in human_leaves:
                    print(f"Leaf descendant {leaf} of class {cls} is not in the narrow background; marking for expansion.")
                    leaves_to_expand_background.add(leaf)
                    parents_to_expand_background.add(cls)

    return list(studyset_leaves), leaves_to_expand_background, parents_to_expand_background

def count_narrow_leaves(narrow_leaves_json, leaves_to_expand_background):

    # Count the leaves in the narrow background
    with open(narrow_leaves_json, "r", encoding="utf-8") as f:
        cached = json.load(f)
    narrow_leaves = set(cached.get("narrow_leaves", []))

    # Add leaves outside the narrow background that are a part of the study set (these will be added to the background for enrichment calculations)
    if leaves_to_expand_background:
        print(f"Adding {len(leaves_to_expand_background)} leaves outside the narrow background to the background set for enrichment calculations.")
        narrow_leaves.update(leaves_to_expand_background)

    return len(narrow_leaves)

def count_narrow_leaves_for_class(class_iri, narrow_leaves_json, leaf_descendants_map):
    """Count how many narrow leaves are associated with the given class_iri."""
    with open(narrow_leaves_json, "r", encoding="utf-8") as f:
        cached = json.load(f)
    narrow_leaves = set(cached.get("narrow_leaves", []))

    if str(class_iri) not in leaf_descendants_map: 
        print(f"⚠️ Class {class_iri} is not found in map file.")
        return 0, 0
    elif len(leaf_descendants_map[str(class_iri)]) == 0: # Should not happen since then it would be a leaf            print(f"⚠️ Class {class_iri} has no descendants in full ontology.")
        return 0, 0
    
    # Get string of subclasses from the map
    leaves_structure = leaf_descendants_map[str(class_iri)]

    # Find the ones that are in the narrow background
    narrow_leaves_for_class = set(leaves_structure).intersection(narrow_leaves)

    return narrow_leaves_for_class, len(narrow_leaves_for_class)
 

def get_enrichment_values_narrow(narrow_leaves_csv, classification, studyset_leaves, studyset_ancestors, class_to_leaf_map, class_to_all_roles_map, roles_to_leaves_map, studyset_ancestors_roles, leaves_to_expand_background):

    # n_bg_leaves and n_ss_leaves will be the same for all classes
    n_bg_leaves = count_narrow_leaves(narrow_leaves_csv, leaves_to_expand_background)
    n_ss_leaves = len(studyset_leaves)

    results =  {} # dictionary to hold results


    if classification in ["structural", "full"]:

        # Calculate enrichment for structural ancestors
        for class_to_check in studyset_ancestors:

            # print(f"Calculating enrichment for class {class_to_check}...")

            _, n_bg_annotated = count_narrow_leaves_for_class(
                class_to_check,
                narrow_leaves_csv,
                class_to_leaf_map
            )
            n_ss_annotated = get_n_ss_annotated(
                studyset_leaves,
                class_to_check,
                class_to_leaf_map,
                classification,
                class_to_all_roles_map,
                roles_to_leaves_map,
            )

            odds, p_value = calculate_p_value(n_ss_annotated, n_ss_leaves, n_bg_annotated, n_bg_leaves)

            # Skip if calculation failed due to invalid counts
            if p_value is None:
                print(f"Skipping class {id_to_name(class_to_check)} due to invalid contingency table")
                continue

            results[class_to_check]={
                "class": id_to_name(class_to_check),
                # "class_id": strip_prefix(class_to_check),
                "n_ss_annotated": n_ss_annotated,
                "n_ss_leaves": n_ss_leaves,
                "n_bg_annotated": n_bg_annotated,
                "n_bg_leaves": n_bg_leaves,
                "odds_ratio": odds,
                "p_value": p_value
            }

    # Calculate enrichment for role classes
    if classification in ["functional", "full"] and studyset_ancestors_roles:

        print(f"Role-based enrichment is not implemented yet.")
        # TODO: IMPLEMENT LATER
        
        """
        print(f"Calculating enrichment for {len(studyset_ancestors_roles)} roles...")
        for role_to_check in studyset_ancestors_roles:
            _, n_bg_annotated = count_removed_classes_for_roles(role_to_check, class_to_leaf_map, classification, roles_to_leaves_map)
            n_ss_annotated = get_n_ss_annotated_for_roles(studyset_leaves, role_to_check, class_to_all_roles_map, roles_to_leaves_map)
           
            odds, p_value = calculate_p_value(n_ss_annotated, n_ss_leaves, n_bg_annotated, n_bg_leaves)
            
            # Skip if calculation failed due to invalid counts
            if p_value is None:
                print(f"Skipping role {id_to_name(role_to_check)} due to invalid contingency table")
                continue
            
            results[role_to_check]={
                "class": id_to_name(role_to_check),
                "n_ss_annotated": n_ss_annotated,
                "n_ss_leaves": n_ss_leaves,
                "n_bg_annotated": n_bg_annotated,
                "n_bg_leaves": n_bg_leaves,
                "odds_ratio": odds,
                "p_value": p_value
            }
            """

    return results


def run_narrow_background_enrichment_analysis(studyset_list,
                            bonferroni_correct=False,
                            benjamini_hochberg_correct=True,
                            root_children_prune=False,
                            levels=2,
                            linear_branch_prune=False,
                            n=2,
                            high_p_value_prune=False,
                            p_value_threshold=0.05,
                            zero_degree_prune=False,
                            classification="structural",
                            narrow_background_leaves_json="data/human_entities_leaves.json"):

    pruning_before_enrichment = root_children_prune or linear_branch_prune

    # Files
    removed_leaves_full_csv = "data/removed_leaf_classes_with_smiles.csv"
    leaf_to_ancestors_map_file = "data/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_map_file = "data/class_to_leaf_descendants_map.json"
    parent_map_file = "data/chebi_parent_map.json"
    class_to_all_roles_map_json = "data/class_to_all_roles_map.json"
    roles_to_leaves_map_json = "data/roles_to_leaves_map.json"

    with open(class_to_leaf_map_file, 'r') as f:
        class_to_leaf_map = json.load(f)
    with open(class_to_all_roles_map_json, 'r') as f:
        class_to_all_roles_map = json.load(f)
    with open(roles_to_leaves_map_json, 'r') as f:
        roles_to_leaves_map = json.load(f)
    with open(narrow_background_leaves_json, 'r') as f:
        cached = json.load(f)

    # Normalize study set IDs to full IRIs
    studyset_list = [normalize_id(cls) for cls in studyset_list]
   
    studyset_leaves, leaves_to_expand_background, parents_to_expand_background = get_studyset_leaves_narrow(
        studyset_list,
        narrow_background_leaves_json,
        removed_leaves_full_csv,
        class_to_leaf_map_file,
    )
    #print(f"Study set leaves: {studyset_leaves}")
    
    studyset_ancestors_all = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    
    # print(f"Study set ancestors: {studyset_ancestors_all}")
    print(f"Number of study set ancestors: {len(studyset_ancestors_all)}")

    """
    if classification in ["functional", "full"]:
        # Collect roles associated with structural ancestors AND leaves
        studyset_ancestors_roles = set()
        
        # Add roles from leaves
        for leaf in studyset_leaves:
            studyset_ancestors_roles.update(class_to_all_roles_map.get(leaf, []))
        
        # Add roles from structural ancestors
        for class_to_check in studyset_ancestors_all:
            studyset_ancestors_roles.update(class_to_all_roles_map.get(class_to_check, []))
        
        # Also include all ancestors of these roles (they will appear in the graph)
        with open(parent_map_file, 'r') as f:
            parent_map = json.load(f)
        
        roles_with_ancestors = set(studyset_ancestors_roles)
        to_process = list(studyset_ancestors_roles)
        
        while to_process:
            role = to_process.pop(0)
            for parent in parent_map.get(role, []):
                if parent not in roles_with_ancestors:
                    roles_with_ancestors.add(parent)
                    to_process.append(parent)
        
        studyset_ancestors_roles = roles_with_ancestors
        print(f"Number of roles (including ancestors): {len(studyset_ancestors_roles)}")
    else:
        studyset_ancestors_roles = set()
        """
    studyset_ancestors_roles = set()
    # print("Role class enrichment is not implemented yet, skipping role collection for now.")
  
    G = create_graph_with_roles_and_structures(studyset_leaves, studyset_ancestors_all, 
                            studyset_ancestors_roles, parent_map_file, 
                            class_to_all_roles_map, classification)

    
    pruned_G = G.copy()

    all_removed_nodes = set()

    if pruning_before_enrichment:

        if root_children_prune:
            print(f"studyset_leaves: {studyset_leaves}")
            print(f"Root children pruner activated, pruning {levels} levels from root")
            time_start_total = time.time()
            pruned_G, removed_nodes, execution_count = root_children_pruner(pruned_G, levels, allow_re_execution = False, execution_count = 0)
            time_end_total = time.time()
            print(f"Total time for root children pruning: {time_end_total - time_start_total} seconds")
            # print(f"Removed nodes by root children pruner: {removed_nodes}")
            all_removed_nodes.update(removed_nodes)

        if linear_branch_prune:
            print(f"Linear branch pruner activated, keeping only every {n}-th node in linear branches")
            
            pruned_G, removed_nodes = linear_branch_collapser_pruner_remove_less(pruned_G, n)
            # print(f"Removed nodes by linear branch pruner: {removed_nodes}")
            all_removed_nodes.update(removed_nodes)
        
        # Remove pruned nodes from studyset_ancestors_all
        studyset_ancestors = [cls for cls in studyset_ancestors_all if cls not in all_removed_nodes]
        print(f"Number of study set ancestors after before-enrichment pruning: {len(studyset_ancestors)}")
    
    else:
        studyset_ancestors = studyset_ancestors_all

    enrichment_results = get_enrichment_values_narrow(
        narrow_leaves_csv=narrow_background_leaves_json,
        classification=classification,
        studyset_leaves=studyset_leaves,
        studyset_ancestors=studyset_ancestors,
        class_to_leaf_map=class_to_leaf_map,
        class_to_all_roles_map=class_to_all_roles_map,
        roles_to_leaves_map=roles_to_leaves_map,
        studyset_ancestors_roles=studyset_ancestors_roles,
        leaves_to_expand_background=leaves_to_expand_background,
    )

    #print("Enrichment results:")
    #print_enrichment_results(enrichment_results)

    if bonferroni_correct:
        print("Applying Bonferroni correction to p-values...")
        enrichment_results, correction_map = bonferroni_correction(enrichment_results)
        # print("Enrichment results after Bonferroni correction:")
        # print_enrichment_results(enrichment_results)
    elif benjamini_hochberg_correct:
        print("Applying Benjamini-Hochberg FDR correction to p-values...")
        enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
        # print("Enrichment results after Benjamini-Hochberg correction:")
        # print_enrichment_results(enrichment_results)
    
    if high_p_value_prune: # Uses corrected p-values if correction was applied
        print(f"High p-value pruner activated, pruning nodes with p-value above {p_value_threshold}")

        time_start_total = time.time()
        pruned_G, removed_nodes = high_p_value_branch_pruner(pruned_G, enrichment_results, p_value_threshold)
        time_end_total = time.time()
        print(f"Total time for high p-value pruning: {time_end_total - time_start_total} seconds")
        # print(f"Removed nodes by high p-value pruner: {removed_nodes}")
        print(f"Number of pruned nodes: {len(removed_nodes)}")

        # update all_removed_nodes
        all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from enrichment results
        for cls in removed_nodes:
            if cls in enrichment_results:
                del enrichment_results[cls]

        # print("Final enrichment results after high p-value pruning:")
        # print_enrichment_results(enrichment_results)

    if zero_degree_prune:
        print("Applying zero-degree pruner to remove nodes with zero degree...")

        time_start_total = time.time()
        pruned_G, removed_nodes = zero_degree_pruner(pruned_G)
        time_end_total = time.time()
        print(f"Total time for zero-degree pruning: {time_end_total - time_start_total} seconds")
        # print(f"Removed nodes by zero-degree pruner: {removed_nodes}")
        print(f"Number of pruned nodes: {len(removed_nodes)}")

        # update all_removed_nodes
        all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from enrichment results
        for cls in removed_nodes:
            if cls in enrichment_results:
                del enrichment_results[cls]

        # print("Final enrichment results after zero-degree pruning:")
        # print_enrichment_results(enrichment_results)

    print("Final enrichment results:")
    print_enrichment_results(enrichment_results)
       
    print(f"Number of removed nodes in total: {len(all_removed_nodes)}")
    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {id_to_name(cls): vals for cls, vals in enrichment_results.items()}
    }
    return results, pruned_G, leaves_to_expand_background, parents_to_expand_background

def run_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy(
        studyset_list,
        levels=2,           # for root children pruner
        n=0,                # for linear branch pruner (0 = remove all intermediate nodes)
        p_value_threshold=0.05,
        classification="structural",
        narrow_background_leaves_json="data/human_entities_leaves.json"):

    # Files
    removed_leaves_full_csv = "data/removed_leaf_classes_with_smiles.csv"
    leaf_to_ancestors_map_file = "data/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_map_file = "data/class_to_leaf_descendants_map.json"
    parent_map_file = "data/chebi_parent_map.json"
    class_to_all_roles_map_json = "data/class_to_all_roles_map.json"
    roles_to_leaves_map_json = "data/roles_to_leaves_map.json"

    with open(class_to_leaf_map_file, 'r') as f:
        class_to_leaf_map = json.load(f)
    with open(class_to_all_roles_map_json, 'r') as f:
        class_to_all_roles_map = json.load(f)
    with open(roles_to_leaves_map_json, 'r') as f:
        roles_to_leaves_map = json.load(f)

    studyset_list = [normalize_id(cls) for cls in studyset_list]

    # Narrow background: resolve leaves and track any that fall outside the human background
    studyset_leaves, leaves_to_expand_background, parents_to_expand_background = get_studyset_leaves_narrow(
        studyset_list,
        narrow_background_leaves_json,
        removed_leaves_full_csv,
        class_to_leaf_map_file,
    )
    print(f"Study set leaves: {studyset_leaves}")

    studyset_ancestors = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    print(f"Number of study set ancestors: {len(studyset_ancestors)}")

    all_removed_nodes = set()

    # Role classification is not yet implemented for narrow background.
    # TODO: when implementing, mirror the role collection block from
    # run_enrichment_analysis_plain_enrich_pruning_strategy and replace
    # get_enrichment_values_narrow's stub with real role enrichment logic.
    studyset_ancestors_roles = set()

    # Enrichment is calculated up front (before graph construction),
    # using the fixed leaves_to_expand_background from leaf resolution above.
    # This set does not change across loop iterations.
    enrichment_results = get_enrichment_values_narrow(
        narrow_leaves_csv=narrow_background_leaves_json,
        classification=classification,
        studyset_leaves=studyset_leaves,
        studyset_ancestors=studyset_ancestors,
        class_to_leaf_map=class_to_leaf_map,
        class_to_all_roles_map=class_to_all_roles_map,
        roles_to_leaves_map=roles_to_leaves_map,
        studyset_ancestors_roles=studyset_ancestors_roles,
        leaves_to_expand_background=leaves_to_expand_background,
    )

    print("Enrichment results:")
    print_enrichment_results(enrichment_results)

    enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
    print("Enrichment results after Benjamini-Hochberg correction:")
    print_enrichment_results(enrichment_results)

    G = create_graph_with_roles_and_structures(
        studyset_leaves,
        studyset_ancestors,
        studyset_ancestors_roles,
        parent_map_file,
        class_to_all_roles_map,
        classification
    )

    ## Pre-loop phase ##
    print("Starting pre-loop pruning phase.")

    G, removed_nodes = high_p_value_branch_pruner(G, enrichment_results, p_value_threshold)
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by high p-value pruner: {removed_nodes}")

    G, removed_nodes = linear_branch_collapser_pruner_remove_less(G, n)
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by linear branch pruner: {removed_nodes}")

    G, removed_nodes, execution_count = root_children_pruner(G, levels, allow_re_execution=False, execution_count=0)
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by root children pruner: {removed_nodes}")

    ## Loop phase ##
    print("Starting loop pruning phase.")
    size_before = G.number_of_nodes()
    size_after = size_before
    first_iteration = True
    iteration = 0

    while size_after < size_before or first_iteration:
        size_before = size_after
        iteration += 1
        print(f"Loop iteration {iteration}")

        # Recalculate BH-corrected p-values on surviving nodes only.
        # leaves_to_expand_background remains fixed — it reflects the original
        # study set composition and should not shrink as nodes are pruned.
        current_enrichment = {
            cls: vals for cls, vals in enrichment_results.items()
            if cls not in all_removed_nodes and G.has_node(cls)
        }
        current_enrichment = benjamini_hochberg_fdr_correction(current_enrichment)

        G, removed_nodes = high_p_value_branch_pruner(G, current_enrichment, p_value_threshold)
        all_removed_nodes.update(removed_nodes)
        print(f"Removed nodes by high p-value pruner: {removed_nodes}")

        G, removed_nodes = zero_degree_pruner(G)
        all_removed_nodes.update(removed_nodes)
        print(f"Removed nodes by zero-degree pruner: {removed_nodes}")

        size_after = G.number_of_nodes()
        first_iteration = False

    ## No final phase pruners ##

    final_enrichment = current_enrichment

    print(f"Number of removed nodes in total: {len(all_removed_nodes)}")
    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {id_to_name(cls): vals for cls, vals in final_enrichment.items()}
    }

    return results, G, leaves_to_expand_background, parents_to_expand_background


if __name__ == "__main__":

    task = "gather_leaves"  # "gather_leaves"

    # human_entities_tsv = "data/combined_hmdb_wikidata.tsv"
    human_entities_tsv = "data/combined_hmdb_wikidata.tsv"
    leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
    class_to_leaf_map = "data/class_to_leaf_descendants_map.json"
    human_leaves_json = "data/human_entities_leaves.json" # output of gather_narrow_leaves

    chebi_parent_map_json="data/chebi_parent_map.json"

    if task == "gather_leaves":
        gather_narrow_leaves(human_entities_tsv, leaves_csv, class_to_leaf_map, human_leaves_json, chebi_parent_map_json)
    
    elif task == "get_studyset_leaves":
        studyset_list = [
            "CHEBI:15377",  # water
            "CHEBI:16113",  # cholesterol
            "CHEBI:17234",  # glucose
            "CHEBI:33917",  # aldohexose
        ]
        studyset_leaves, leaves_to_expand, parents_to_expand = get_studyset_leaves_narrow(
            studyset_list,
            human_entities_leaves_json=human_leaves_json,
            leaves_csv=leaves_csv,
            class_to_leaf_map=class_to_leaf_map,
        )
        print(f"Study set leaves: {studyset_leaves}")
        print(f"Leaves to expand background: {leaves_to_expand}")
        print(f"Parents to expand background: {parents_to_expand}")
    else:
        print(f"Unknown task: {task}")

