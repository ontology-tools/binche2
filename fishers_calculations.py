from math import inf
from scipy.stats import fisher_exact
import pandas as pd
import json
import copy 
import math
from visualitations_and_pruning import root_children_pruner, linear_branch_collapser_pruner_remove_less, high_p_value_branch_pruner, zero_degree_pruner, create_graph_from_map, id_to_name, create_graph_with_roles_and_structures
from pre_fishers_calculations import count_removed_classes_for_class, count_removed_leaves, count_removed_classes_for_roles
from multiple_test_corrections import bonferroni_correction, benjamini_hochberg_fdr_correction
import time

"""
from goat-tools:
a = study_count         # n_ss_annotated
b = study_n - study_count  # n_ss_leaves - n_ss_annotated
c = pop_count - study_count  # n_bg_annotated - n_ss_annotated
d = pop_n - pop_count - b    # n_bg_leaves - n_bg_annotated - (n_ss_leaves - n_ss_annotated)
table = [[a, b], [c, d]]
"""


# n_ss_leaves = total number of input classes in the study set (if they are all leaf classes, otherwise count corresponding leaf classes)
# n_ss_annotated = number of the input classes that are descendants of the given class

# n_bg_leaves = total number of leaf classes in the background set (the ones filtered out as leaves from the ontology)   
# n_bg_annotated = number of those leaf classes that are descendants of the given class)


"""
                    In Class     Not in Class
Study Set              a             b
Background (rest)      c             d

a = study compounds annotated to this class
b = study compounds NOT annotated to this class
c = background compounds in this class (excluding study set)
d = background compounds not in this class (excluding study set)
"""

def calculate_p_value(n_ss_annotated, n_ss_leaves, n_bg_annotated, n_bg_leaves):
    a = n_ss_annotated
    b = n_ss_leaves - n_ss_annotated
    c = n_bg_annotated - n_ss_annotated
    d = n_bg_leaves - n_bg_annotated - b
    odds, p = fisher_exact([[a, b], [c, d]], alternative='greater')
    # ‘greater’: thget_ne odds ratio of the underlying population is greater than one
    return odds, p

def normalize_id(raw_id: str) -> str:
    value = raw_id.strip().replace('"', "")
    if value.startswith("http://") or value.startswith("https://"):
        return value
    # Convert CHEBI:ID to http://purl.obolibrary.org/obo/CHEBI_ID format
    if value.startswith("CHEBI:"):
        chebi_id = value.replace(":", "_")
        return f"http://purl.obolibrary.org/obo/{chebi_id}"
    if not value.startswith("http://"):
        value = value.replace(":", "_")
        return f"http://purl.obolibrary.org/obo/{value}"
    return value

def get_leaves(studyset_list, leaves_csv, class_to_leaf_map):
    studyset_leaves = set()

    leaves_df = pd.read_csv(leaves_csv)
    print(leaves_df['IRI'].iloc[0].encode())
    
    for cls in studyset_list:
        print(f"Processing class {cls}...")
        if cls in leaves_df['IRI'].values:
            # add to list of studyset leaves
            print(f"Class {cls} is already leaf.")
            studyset_leaves.add(cls)
        else:
            # get leaf descendants from map and add them to studyset leaves
            leaf_descendants = set(class_to_leaf_map.get(cls, []))
            print(f"Class {cls} is not a leaf, adding its {len(leaf_descendants)} leaf descendants.")
            studyset_leaves.update(leaf_descendants)

    return list(studyset_leaves)

def get_ancestors_for_inputs(studyset_leaves, leaf_to_all_parents_map_json):
    with open(leaf_to_all_parents_map_json, 'r') as f:
        leaf_to_all_parents_map = json.load(f)

    studyset_ancestors = set()
    for leaf in studyset_leaves:
        parents = leaf_to_all_parents_map.get(leaf, [])
        studyset_ancestors.update(parents)

    return studyset_ancestors


def get_n_ss_annotated(studyset_leaves, class_to_check, class_to_leaf_map, classification, class_to_all_roles_map, roles_to_leaves_map):
    """
    n_ss_annotated = number of input classes that are leaf descendants of the given class.

    studyset_leaves: list of class IDs that were provided by the user (the study set)
    class_to_check: the ontology class for which we want n_ss_annotated
    map_file: JSON file mapping each class to all its leaf descendants
    """

    leaves = set()

    if classification in ["structural", "full"]:
        # descendants of the class we are calculating enrichment for
        leaf_descendants = set(class_to_leaf_map.get(class_to_check, []))
        leaves.update(leaf_descendants)

        # count how many study classes appear in those leaf descendants
        n_ss_annotated = len(leaves.intersection(set(studyset_leaves)))
        return n_ss_annotated
    
    # if classification in ["functional", "full"]:
    #     # Get all roles (direct + inherited from ancestors + role ancestors)
    #     all_roles = class_to_all_roles_map.get(class_to_check, [])
    #     for role in all_roles:
    #         leaves.update(roles_to_leaves_map.get(role, []))

    else:
        print(f"Classification {classification} is not supported for counting classes.")
        return 0

def get_n_ss_annotated_for_roles(studyset_leaves, class_to_check, class_to_all_roles_map, roles_to_leaves_map):
    leaves = set()
    leaves.update(roles_to_leaves_map.get(class_to_check, []))
    n_ss_annotated = len(leaves.intersection(set(studyset_leaves)))
    return n_ss_annotated

def get_enrichment_values(removed_leaves_csv, classification, studyset_leaves, studyset_ancestors, class_to_leaf_map, class_to_all_roles_map, roles_to_leaves_map, studyset_ancestors_roles):

    # n_bg_leaves and n_ss_leaves will be the same for all classes
    n_bg_leaves = count_removed_leaves(removed_leaves_csv)
    n_ss_leaves = len(studyset_leaves)

    results =  {} # dictionary to hold results

    if classification in ["structural", "full"]:

        # Calculate enrichment for structural ancestors
        for class_to_check in studyset_ancestors:

            # print(f"Calculating enrichment for class {class_to_check}...")

            _, n_bg_annotated = count_removed_classes_for_class(
                class_to_check,
                class_to_leaf_map,
                classification,
                class_to_all_roles_map,
                roles_to_leaves_map,
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

        # if classification == "functional" or classification == "full":
        # # Update studyset_ancestors_roles with the roles associated (direct + inherited from ancestors) with the current class being checked
        # # These roles will be added to the graph
        #     studyset_ancestors_roles.update(class_to_all_roles_map.get(class_to_check, []))

    # Calculate enrichment for role classes
    if classification in ["functional", "full"] and studyset_ancestors_roles:
        print(f"Calculating enrichment for {len(studyset_ancestors_roles)} roles...")
        
        for role_to_check in studyset_ancestors_roles:
            _, n_bg_annotated = count_removed_classes_for_roles(role_to_check, class_to_leaf_map, classification, roles_to_leaves_map)
            n_ss_annotated = get_n_ss_annotated_for_roles(studyset_leaves, role_to_check, class_to_all_roles_map, roles_to_leaves_map)
           
            odds, p_value = calculate_p_value(n_ss_annotated, n_ss_leaves, n_bg_annotated, n_bg_leaves)
            
            results[role_to_check]={
                "class": id_to_name(role_to_check),
                "n_ss_annotated": n_ss_annotated,
                "n_ss_leaves": n_ss_leaves,
                "n_bg_annotated": n_bg_annotated,
                "n_bg_leaves": n_bg_leaves,
                "odds_ratio": odds,
                "p_value": p_value
            }

    return results

def print_enrichment_results(enrichment_results):
    # Include corrected p-values if they exist
    # print(f"{'Class:':45} {'p-value:':15} {'p-value (corrected):':20} {'n_ss_annotated':20} {'n_bg_annotated':20}" )
    print(f"{'Class:':45} {'raw p-value:':15} {'p-value (corrected):':20}" )
    print("-" * 200)

    for r in enrichment_results.values():
        corrected_p = r.get('p_value_corrected', None)
        corrected_p_str = f"{corrected_p:.4e}" if corrected_p is not None else "N/A"
        #print(f"{r['class']:45} {r['p_value']:.4e}      {corrected_p_str:20} {r['n_ss_annotated']:20} {r['n_bg_annotated']:20}")
        print(f"{r['class']:45} {r['p_value']:.4e}      {corrected_p_str:20}")

""" Graphing and pruning stategies """

def run_enrichment_analysis(studyset_list,
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
                            check_leaf_classes=False):

    pruning_before_enrichment = root_children_prune or linear_branch_prune # add other pruning strategies when implemented

    # Files
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
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

    # Normalize study set IDs to full IRIs
    studyset_list = [normalize_id(cls) for cls in studyset_list]

   
    studyset_leaves = get_leaves(studyset_list, removed_leaves_csv, class_to_leaf_map)
    #print(f"Study set leaves: {studyset_leaves}")
    
    studyset_ancestors_all = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    
    # print(f"Study set ancestors: {studyset_ancestors_all}")
    print(f"Number of study set ancestors: {len(studyset_ancestors_all)}")

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

    enrichment_results = get_enrichment_values(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        studyset_ancestors_roles
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
    return results, pruned_G

####################################
# Combine pruning strategies
####################################

# Plain Enrichment Pruning Strategy: For the pre-loop phase this strategy applies the High Value Branch Pruner (0.05), 
# the Linear Branch Collapser Pruner, and the Root Children Pruner (3 (change to 2) levels, without repetition). 
# During the loop phase,
# this strategy applies the Molecule Leaves Pruner, the High P-Value Branch Pruner (0.05), the Linear Branch Collapser Pruner,
# and the Zero Degree Vertex Pruner. No pruners are applied in the final phase post-loop.

def run_enrichment_analysis_plain_enrich_pruning_strategy(studyset_list,
                            levels=2, # for root children pruner
                            n=0, # for linear branch pruner
                            p_value_threshold=0.05, # for high p-value pruner
                            classification="structural",
                            check_leaf_classes=False):

    # Files
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
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

    studyset_leaves = get_leaves(studyset_list, removed_leaves_csv, class_to_leaf_map)
    print(f"Study set leaves: {studyset_leaves}")

    studyset_ancestors = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    print(f"Study set ancestors: {studyset_ancestors}")
    print(f"Number of study set ancestors: {len(studyset_ancestors)}")

    all_removed_nodes = set()

    if classification in ["functional", "full"]:
        # Collect roles associated with leaves AND structural ancestors
        studyset_ancestors_roles = set()
        
        # Add roles from leaves
        for leaf in studyset_leaves:
            studyset_ancestors_roles.update(class_to_all_roles_map.get(leaf, []))
        
        # Add roles from structural ancestors
        for class_to_check in studyset_ancestors:
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
        print(f"Number of roles (including all ancestors): {len(studyset_ancestors_roles)}")
    else:
        studyset_ancestors_roles = set()

    enrichment_results = get_enrichment_values(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        studyset_ancestors_roles
    )

    print("Enrichment results:")
    print_enrichment_results(enrichment_results)

    enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
    print("Enrichment results after Benjamini-Hochberg correction:")
    print_enrichment_results(enrichment_results)

    pre_pruned_G = create_graph_with_roles_and_structures(studyset_leaves, studyset_ancestors, 
                            studyset_ancestors_roles, parent_map_file, 
                            class_to_all_roles_map, classification)
    G = pre_pruned_G.copy()

    ## Pre-loop phase ##
    print("Starting pre-loop pruning phase.")
    G, removed_nodes = high_p_value_branch_pruner(G, enrichment_results, p_value_threshold)
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by high p-value pruner: {removed_nodes}")

    G, removed_nodes = linear_branch_collapser_pruner_remove_less(G, n)
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by linear branch pruner: {removed_nodes}")

    G, removed_nodes, execution_count = root_children_pruner(G, levels, allow_re_execution = False, execution_count = 0)
    all_removed_nodes.update(removed_nodes)
    print(f"Removed nodes by root children pruner: {removed_nodes}")

    ## Loop phase ##
    print("Starting loop pruning phase.")
    # Count the number of nodes in G so it can be compared after each iteration
    size_before = G.number_of_nodes()
    size_after = size_before
    first_iteration = True
    iteration = 0

    # while the size changes, keep applying the loop phase pruners
    while size_after < size_before or first_iteration:
        size_before = size_after
        iteration += 1
        print(f"Loop iteration {iteration}")

        ## Recalculate corrected p-values

        # Remove pruned nodes from enrichment results
        current_enrichment = {
            cls: vals for cls, vals in enrichment_results.items()
            if cls not in all_removed_nodes and G.has_node(cls)
        }

        current_enrichment = benjamini_hochberg_fdr_correction(current_enrichment)
        
        G, removed_nodes = high_p_value_branch_pruner(G, current_enrichment, p_value_threshold)
        all_removed_nodes.update(removed_nodes)
        print(f"Removed nodes by high p-value pruner: {removed_nodes}")
        # not including since we do not want to remove too many nodes
        # G, removed_nodes = linear_branch_collapser_pruner_remove_less(G, n)
        # all_removed_nodes.update(removed_nodes)

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

    return results, G



if __name__ == "__main__":

    bonferroni_correct = False
    benjamini_hochberg_correct = True # Used in Binche1

    root_children_prune = True
    levels = 2 # Number of levels to prune including root. 1 only prunes root, 2 prunes root and it's direct neighbour, and so on.
    # allow_re_execution = False  # Currently not necessary. Whether the pruner can be executed multiple times on a given graph.
    # execution_count = 0  # Currently not necessary. Counter for the number of executions

    linear_branch_prune = True
    n = 2 # Keep only every n-th node in linear branches

    high_p_value_prune = True
    p_value_threshold = 0.05

    zero_degree_prune = True


    classification = "functional" # "functional" or "structural" or "full"
    check_leaf_classes = False # Checks that the found the leaf classes are of the expected type (Functional or Structural). If the classification is correct, 
                                # this should never be a problem and can be set to False.

    # studyset_list = ["http://purl.obolibrary.org/obo/CHEBI_77030","http://purl.obolibrary.org/obo/CHEBI_79036"]
    # studyset_list = ["http://purl.obolibrary.org/obo/CHEBI_77030"] 
    # Problem "Warning: p-value for node http://purl.obolibrary.org/obo/CHEBI_36357 not found. Assuming high p-value." Not found because it was removed in root children pruner but whyyy is it still in the graph???" Inte kollat om fixat!!!!!!!!!!!!!!!!!!!
    #studyset_list =["http://purl.obolibrary.org/obo/CHEBI_17234"]
    # studyset_list =["http://purl.obolibrary.org/obo/CHEBI_37626"]
    studyset_list =["http://purl.obolibrary.org/obo/CHEBI_77030"]


    results = run_enrichment_analysis(studyset_list,
                            bonferroni_correct=bonferroni_correct,
                            benjamini_hochberg_correct=benjamini_hochberg_correct,
                            root_children_prune=root_children_prune,
                            levels=levels,
                            linear_branch_prune=linear_branch_prune,
                            n=n,
                            high_p_value_prune=high_p_value_prune,
                            p_value_threshold=p_value_threshold,
                            zero_degree_prune=zero_degree_prune,
                            classification=classification,
                            check_leaf_classes=check_leaf_classes
                            )

    # run_enrichment_analysis_plain_enrich_pruning_strategy(studyset_list,
    #                         levels=2, # for root children pruner
    #                         n=2, # for linear branch pruner
    #                         p_value_threshold=0.05, # for high p-value pruner
    #                         classification="structural",
    #                         check_leaf_classes=False)

    
    
    # print("Final results:")
    # print(results)

    