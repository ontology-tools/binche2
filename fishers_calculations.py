from math import inf
from scipy.stats import fisher_exact
import pandas as pd
import json
import copy 
import math
from visualitations_and_pruning import root_children_pruner, linear_branch_collapser_pruner_remove_less, high_p_value_branch_pruner, zero_degree_pruner, create_graph_from_map, id_to_name
from pre_fishers_calculations import count_removed_classes_for_class, count_removed_leaves
from multiple_test_corrections import bonferroni_correction, benjamini_hochberg_fdr_correction

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

def calculate_p_value(n_ss_annotated, n_ss_leaves, n_bg_annotated, n_bg_leaves):
    a = n_ss_annotated
    b = n_ss_leaves - n_ss_annotated
    c = n_bg_annotated - n_ss_annotated
    d = n_bg_leaves - n_bg_annotated - b
    odds, p = fisher_exact([[a, b], [c, d]], alternative='greater')
    # ‘greater’: the odds ratio of the underlying population is greater than one
    return odds, p

def get_leaves(studyset_list, leaves_csv, class_to_leaf_map_json):
    studyset_leaves = set()

    leaves_df = pd.read_csv(leaves_csv)
    print(leaves_df['IRI'].iloc[0].encode())
    with open(class_to_leaf_map_json, 'r') as f:
        class_to_leaf_map = json.load(f)
    
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


def get_n_ss_annotated(studyset_leaves, class_to_check, map_file):
    """
    n_ss_annotated = number of input classes that are leaf descendants of the given class.

    studyset_leaves: list of class IDs that were provided by the user (the study set)
    class_to_check: the ontology class for which we want n_ss_annotated
    map_file: JSON file mapping each class to all its leaf descendants
    """

    # load mapping: {class: [list_of_leaf_descendants]}
    with open(map_file, 'r') as f:
        class_to_leaf_map = json.load(f)

    # descendants of the class we are calculating enrichment for
    leaf_descendants = set(class_to_leaf_map.get(class_to_check, []))

    # count how many study classes appear in those leaf descendants
    n_ss_annotated = sum(cls in leaf_descendants for cls in studyset_leaves)

    return n_ss_annotated

def get_enrichment_values(removed_leaves_csv, classification, studyset_leaves, studyset_ancestors, class_to_leaf_map_file, check_leaf_classes = False):

    #n_bg_leaves and n_ss_leaves will be the same for all classes
    n_bg_leaves = count_removed_leaves(removed_leaves_csv, classification)
    n_ss_leaves = len(studyset_leaves)

    results =  {} # dictionary to hold results

    for class_to_check in studyset_ancestors:

        # print(f"Calculating enrichment for class {class_to_check}...")

        _, n_bg_annotated = count_removed_classes_for_class(class_to_check, class_to_leaf_map_file, classification, check_leaf_classes, removed_leaves_csv)
        n_ss_annotated = get_n_ss_annotated(studyset_leaves, class_to_check, class_to_leaf_map_file)

        odds, p_value = calculate_p_value(n_ss_annotated, n_ss_leaves, n_bg_annotated, n_bg_leaves)

        results[class_to_check]={
            "class": id_to_name(class_to_check),
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
# TODO: look over pruning strategies and loops as explained in article
# TODO: look over p value pruning

# TODO: Add pruning methods to not have to enrich as many classes
# TODO: Add multiple testing correction

# TODO: Maybe remove leaves from graph? Now they are there to be able to start the graph from them. Becomes a proble in high p-value pruning because they do not have p-values.


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
    color_map = ['#FFB6C1', "#F44280", "#AA83A7", "#83163A", "#E63FE6", '#FFA07A', '#FF69B4'] 

    # Files
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"
    leaf_to_ancestors_map_file = "data/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_map_file = "data/class_to_leaf_descendants_map.json"
    parent_map_file = "data/chebi_parent_map.json"

    # Currently the code is working with the whole chebi iri needed. Make sure the input study set is in the same format:
    if not studyset_list[0].startswith("http://purl.obolibrary.org/obo/"):
        studyset_list = [f"http://purl.obolibrary.org/obo/{cls}" for cls in studyset_list]

    studyset_leaves = get_leaves(studyset_list, removed_leaves_csv, class_to_leaf_map_file)
    print(f"Study set leaves: {studyset_leaves}")

    studyset_ancestors_all = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    print(f"Study set ancestors: {studyset_ancestors_all}")
    print(f"Number of study set ancestors: {len(studyset_ancestors_all)}")

    pruned_G = None
    all_removed_nodes = set()

    if pruning_before_enrichment:
        # Currently the whole ontology has to load for pruning to work. Very slow.
        # TODO: solve this issue.

        if root_children_prune:
            print(f"studyset_leaves: {studyset_leaves}")
            print(f"Root children pruner activated, pruning {levels} levels from root")
            G = create_graph_from_map(studyset_leaves, parent_map_file, color_map, max_n_leaf_classes=inf)
            # draw_graph(G, graphing_layout="default", title="Graph")
            pruned_G = G.copy()
            pruned_G, removed_nodes, execution_count = root_children_pruner(pruned_G, levels, allow_re_execution = False, execution_count = 0)
            print(f"Removed nodes by root children pruner: {removed_nodes}")
            all_removed_nodes.update(removed_nodes)

        if linear_branch_prune:
            print(f"Linear branch pruner activated, keeping only every {n}-th node in linear branches")
            if pruned_G is None:
                pruned_G = create_graph_from_map(studyset_leaves, parent_map_file, color_map, max_n_leaf_classes=inf)
            
            pruned_G, removed_nodes = linear_branch_collapser_pruner_remove_less(pruned_G, n)
            print(f"Removed nodes by linear branch pruner: {removed_nodes}")
            all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from studyset_ancestors_all
        studyset_ancestors = [cls for cls in studyset_ancestors_all if cls not in all_removed_nodes]
        print(f"Number of study set ancestors after before-enrichment pruning: {len(studyset_ancestors)}")
    
    else:
        studyset_ancestors = studyset_ancestors_all

    enrichment_results = get_enrichment_values(removed_leaves_csv, classification, studyset_leaves, studyset_ancestors, class_to_leaf_map_file, check_leaf_classes)

    print("Enrichment results:")
    print_enrichment_results(enrichment_results)

    if bonferroni_correct:
        print("Applying Bonferroni correction to p-values...")
        enrichment_results, correction_map = bonferroni_correction(enrichment_results)
        print("Enrichment results after Bonferroni correction:")
        print_enrichment_results(enrichment_results)
    elif benjamini_hochberg_correct:
        print("Applying Benjamini-Hochberg FDR correction to p-values...")
        enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
        print("Enrichment results after Benjamini-Hochberg correction:")
        print_enrichment_results(enrichment_results)

    if high_p_value_prune: # Uses corrected p-values if Bonferroni correction was applied
        print(f"High p-value pruner activated, pruning nodes with p-value above {p_value_threshold}")
        if pruned_G is None:
            print("Creating new graph for high p-value pruning (no prior pruning applied)")
            pruned_G = create_graph_from_map(studyset_leaves, parent_map_file, color_map, max_n_leaf_classes=inf)
        else:
            print("Using previously pruned graph for high p-value pruning")

        pruned_G, removed_nodes = high_p_value_branch_pruner(pruned_G, enrichment_results, p_value_threshold)
        print(f"Removed nodes by high p-value pruner: {removed_nodes}")
        print(f"Number of pruned nodes: {len(removed_nodes)}")

        # update all_removed_nodes
        all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from enrichment results
        for cls in removed_nodes:
            if cls in enrichment_results:
                del enrichment_results[cls]

        print("Final enrichment results after high p-value pruning:")
        print_enrichment_results(enrichment_results)

    if zero_degree_prune:
        print("Applying zero-degree pruner to remove nodes with zero degree...")
        if pruned_G is None:
            print("Creating new graph for zero-degree pruning (no prior pruning applied)")
            pruned_G = create_graph_from_map(studyset_leaves, parent_map_file, color_map, max_n_leaf_classes=inf)
        else:
            print("Using previously pruned graph for zero-degree pruning")

        pruned_G, removed_nodes = zero_degree_pruner(pruned_G)
        print(f"Removed nodes by zero-degree pruner: {removed_nodes}")
        print(f"Number of pruned nodes: {len(removed_nodes)}")

        # update all_removed_nodes
        all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from enrichment results
        for cls in removed_nodes:
            if cls in enrichment_results:
                del enrichment_results[cls]
        
    # Create graph if it does not exist yet
    if pruned_G is None:
        pruned_G = create_graph_from_map(studyset_leaves, parent_map_file, color_map, max_n_leaf_classes=inf)

        print("Final enrichment results after zero-degree pruning:")
        print_enrichment_results(enrichment_results)
    
    print(f"Number of removed nodes in total: {len(all_removed_nodes)}")
    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {id_to_name(cls): vals for cls, vals in enrichment_results.items()}
    }
    return results, pruned_G

# OBS: in Binche1, this would be the input (adapt code to follow this?):

# String ontologyFile = binchePrefs.get(BiNChEOntologyPrefs.RoleAndStructOntology.name(), null);
# // the input path points to a file where the list of ChEBI IDs (one per line, CHEBI:03432) are stored.
# String elementsForEnrichFile = inputPath;

if __name__ == "__main__":

    bonferroni_correct = False
    benjamini_hochberg_correct = False # Used in Binche1

    root_children_prune = True
    levels = 2 # Number of levels to prune from root. 1 only prunes root and it's direct neighbour, and so on.
    # allow_re_execution = False  # Currently not necessary. Whether the pruner can be executed multiple times on a given graph.
    # execution_count = 0  # Currently not necessary. Counter for the number of executions

    linear_branch_prune = False
    n = 2 # Keep only every n-th node in linear branches
    # TODO: implement

    high_p_value_prune = False
    p_value_threshold = 0.05

    zero_degree_prune = False


    classification = "structural" # "functional" or "structural" or "full"
    check_leaf_classes = False # Checks that the found the leaf classes are of the expected type (Functional or Structural). If the classification is correct, 
                                # this should never be a problem and can be set to False.

    # studyset_list = ["http://purl.obolibrary.org/obo/CHEBI_77030","http://purl.obolibrary.org/obo/CHEBI_79036"]
    # studyset_list = ["http://purl.obolibrary.org/obo/CHEBI_77030"] 
    # Problem "Warning: p-value for node http://purl.obolibrary.org/obo/CHEBI_36357 not found. Assuming high p-value." Not found because it was removed in root children pruner but whyyy is it still in the graph???" Inte kollat om fixat!!!!!!!!!!!!!!!!!!!
    # studyset_list =["http://purl.obolibrary.org/obo/CHEBI_17234"]
    studyset_list =["http://purl.obolibrary.org/obo/CHEBI_37626"]

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
                            check_leaf_classes=check_leaf_classes)
    
    print("Final results:")
    # print(results)