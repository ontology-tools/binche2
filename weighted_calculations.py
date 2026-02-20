import numpy as np
from scipy.special import erf
from scipy.optimize import brentq, newton
from fishers_calculations import get_n_ss_annotated
from pre_fishers_calculations import count_removed_classes_for_class
from math import inf
import json
from fishers_calculations import get_leaves, get_ancestors_for_inputs, count_removed_leaves, print_enrichment_results
from visualitations_and_pruning import (
    root_children_pruner,
    linear_branch_collapser_pruner_remove_less,
    high_p_value_branch_pruner,
    zero_degree_pruner,
    create_graph_from_map,
    create_graph_with_roles_and_structures,
    id_to_name,
    )
from multiple_test_corrections import bonferroni_correction, benjamini_hochberg_fdr_correction


"""Using SaddleSum for weighted enrichment calculations"""

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

def calc_cgf_derivatives(lambda_val, weights_dict, N_tot):
    """
    The Cumulant Generating Function (CGF) K(λ) describes the distribution of random sums.
    We need its first and second derivatives for the SaddleSum test.

    Parameters:
        lambda_val : float
            The parameter λ at which to evaluate derivatives
        weights_dict : dict
            {element_id: weight} for elements with known weights
        N_tot : int
            Total population size (including elements with weight 0)
    
    Returns:
        d1k : float
            First derivative K'(λ) - relates to mean
        d2k : float
            Second derivative K''(λ) - relates to variance
        rho_t : float
            Normalization constant
    """
    
    wmax = max(weights_dict.values())

    n_with_weights = len(weights_dict)
    n_without_weights = N_tot - n_with_weights
    
    rho_t = 0.0      # Denominator of K'(t)
    d1_rho_t = 0.0   # Numerator of K'(t)
    d2_rho_t = 0.0   # Numerator of K''(t)

    # Sum over elements with known weights (the study set)
    for w in weights_dict.values():
        exp_term = np.exp(lambda_val * (w - wmax))
        rho_t += exp_term
        d1_rho_t += exp_term * w
        d2_rho_t += exp_term * w * w

    # Add contribution from elements without weights (background, weight = 0)

    rho_t += n_without_weights * np.exp(lambda_val * (-wmax))
    
    # Calculate derivatives
    d1k = d1_rho_t / rho_t
    d2k = d2_rho_t / rho_t - d1k * d1k

    return d1k, d2k, rho_t


def saddlepoint_equation(lambda_val, m, sum_weights, weights_dict, N_tot):
    """
    The saddlepoint equation: m * K'(λ) - sum_weights = 0
    
    We're solving for λ such that the expected sum under the tilted distribution
    equals the observed sum.

    Parameters:
        lambda_val : float
            The parameter λ to evaluate
        m : int
            Number of elements in the category (from background)
        sum_weights : float
            Sum of observed weights for elements in both study set AND category
        weights_dict : dict
            {element_id: weight} for elements with known weights
        N_tot : int
            Total population size (including elements with weight 0)
    
    Returns:
        float : Value of equation (should be ~0 when λ is correct)
    """
    d1k, _, _ = calc_cgf_derivatives(lambda_val, weights_dict, N_tot)
    return m * d1k - sum_weights # S hat


def saddlepoint_equation_derivative(lambda_val, m, weights_dict, N_tot):
    """
    Derivative of the saddlepoint equation with respect to λ.
    Used for Newton-Raphson optimization.

    Parameters:
        lambda_val : float
            The parameter λ to evaluate
        m : int
            Number of elements in the category
        weights_dict : dict
            {element_id: weight} for elements with known weights
        N_tot : int
            Total population size (including elements with weight 0)
    
    Returns:
        float : Derivative value
    """
    _, d2k, _ = calc_cgf_derivatives(lambda_val, weights_dict, N_tot)
    return m * d2k


def find_lambda(m, sum_weights, weights_dict, N_tot):
    """
    Find the λ parameter that satisfies the saddlepoint equation: m * K'(λ) = sum_weights
    
    Uses two-stage optimization:
    1. Bisection method (brentq) - Robust but slower, finds initial interval
    2. Newton-Raphson method (newton) - Faster convergence, refines solution

    Parameters:
        m : int
            Number of elements in the category
        sum_weights : float
            Sum of observed weights
        weights_dict : dict
            {element_id: weight} for elements with known weights
        N_tot : int
            Total population size (including elements with weight 0)
    
    Returns:
        float : Optimal λ value
    """
    # Step 1: Bisection to find initial λ in range [0.00001, 5.0]
    a = 0.00001
    b = 5.0
    try:
        lambda_init = brentq(
            lambda l: saddlepoint_equation(l, m, sum_weights, weights_dict, N_tot),
            a,
            b,
            maxiter=100,
            xtol=0.05
        )
    except ValueError:
        # Bisection failed - try expanding the interval
        print(f"Warning: Bisection failed for m={m}, sum_weights={sum_weights}") # Most common warning
        f_a = saddlepoint_equation(a, m, sum_weights, weights_dict, N_tot)
        f_b = saddlepoint_equation(b, m, sum_weights, weights_dict, N_tot)
        if f_a == f_b:
            lambda_init = 0.0
        else:
            lambda_init = 0.0
            for _ in range(5):
                a = b
                b = b + 10.0
                try:
                    lambda_init = brentq(
                        lambda l: saddlepoint_equation(l, m, sum_weights, weights_dict, N_tot),
                        a,
                        b,
                        maxiter=100,
                        xtol=0.05
                    )
                    break
                except ValueError:
                    continue

    # Step 2: Newton-Raphson refinement for higher precision (only if needed)
    f_init = saddlepoint_equation(lambda_init, m, sum_weights, weights_dict, N_tot)
    if abs(f_init) < 1e-8:
        lambda_optimal = lambda_init
    else:
        try:
            lambda_optimal = newton(
                func=lambda l: saddlepoint_equation(l, m, sum_weights, weights_dict, N_tot),
                x0=lambda_init,
                fprime=lambda l: saddlepoint_equation_derivative(l, m, weights_dict, N_tot),
                tol=0.01,
                maxiter=50
            )
        except RuntimeError:
            # Newton failed to converge
            print(f"Warning: Newton-Raphson failed to converge for m={m}, sum_weights={sum_weights}")
            lambda_optimal = lambda_init  # Use bisection result

    return lambda_optimal


def lugannani_rice_pvalue(lambda_val, m, weights_dict, N_tot):
    """
    Calculate p-value using the Lugannani-Rice saddlepoint approximation.
    
    This approximates the tail probability P(Sum ≥ observed_sum) for a sum
    of random variables. It's more accurate than normal approximation.
    
    The formula combines:
    1. Normal approximation (ndtr)
    2. Two correction terms for skewness and kurtosis

    Parameters:
        lambda_val : float
            The saddlepoint parameter (solution to saddlepoint equation)
        m : int
            Number of observed study-set elements in the category
        weights_dict : dict
            {element_id: weight} for elements with known weights
        N_tot : int
            Total population size (including elements with weight 0)
    
    Returns:
        float : P-value (probability in [0, 1])
    """
    # Get CGF derivatives at the saddlepoint
    d1k, d2k, rho_t = calc_cgf_derivatives(lambda_val, weights_dict, N_tot)
    
    if len(weights_dict) == 0:
        return 1.0
    
    wmax = max(weights_dict.values())

    ## Calculate components for Lugannani-Rice formula

    # exp_h: Related to the probability mass function at the saddlepoint
    exp_h = rho_t * np.exp(lambda_val * (wmax - d1k)) / N_tot

    # C: Related to the curvature of the log probability at the saddlepoint
    C = 2 * lambda_val * np.sqrt(d2k)

    # D: Standardized distance from the mean to the observed sum
    D_argument = lambda_val * (d1k - wmax) - np.log(rho_t) + np.log(N_tot)
    
    if D_argument < 0:
        print(f"Warning: D_argument is negative ({D_argument}), setting p-value to 1.0")
        return 1.0  # Invalid configuration, return non-significant

    D = -np.sqrt(2) * np.sqrt(D_argument)

    # phi: Correction factor based on saddlepoint density
    phi = np.sqrt(2 / np.pi) * (exp_h ** m)

    # z: Standardized test statistic for normal CDF
    z = (D * np.sqrt(m)) / np.sqrt(2)

    # Calculate p-value only if D * sqrt(m) <= -1 to ensure we're in the enrichment tail
    if D * np.sqrt(m) <= -1:
        # Normal CDF using error function
        # erf(z) = 2/√π * integral from 0 to z of exp(-t²) dt. 
        # Standard normal CDF = (1 + erf(z/√2)) / 2
        ndtr = (1 + erf(z)) / 2 # Normal Distribution (cumulative) Tail Right
        
        # Lugannani-Rice correction terms
        # These adjust for skewness and kurtosis of the distribution
        correction1 = phi / C / np.sqrt(m)
        correction2 = phi / D / np.sqrt(m) / 2
        
        # Final p-valuendtr = Normal Distribution 
        p_value = ndtr + correction1 + correction2
        
        # Handle numerical issues
        if np.isnan(p_value) or np.isinf(p_value):
            print(f"Warning: p-value computation resulted in NaN or Inf, setting to 1.0")
            p_value = 1.0
    else:
        # Not in the enrichment tail (D * sqrt(m) > -1)
        print(f"D * sqrt(m) > -1 ({D * np.sqrt(m)}), setting p-value to 1.0")
        p_value = 1.0
    
    # # Clamp into valid probability range
    # if p_value < 0.0:
    #     p_value = 0.0
    # elif p_value > 1.0:
    #     p_value = 1.0

    print(f"\n--- Lugannani-Rice Debug ---")
    print(f"D: {D}")
    print(f"C: {C}")
    print(f"lambda*T - K: {lambda_val * (d1k - wmax) - np.log(rho_t) + np.log(N_tot)}")

    print("-----------------------------\n")

    return p_value


def calculate_weighted_pvalue(studyset_leaves, class_to_check, class_to_leaf_map,
                              weights_dict, N_total, classification,
                              class_to_all_roles_map, roles_to_leaves_map):
    """
    Main function to calculate weighted enrichment p-value for a single category.
    
    This integrates all the saddlepoint calculations to produce a final p-value
    similar to Fisher's exact test, but accounting for element weights.
    
    Parameters:
        studyset_leaves : Element IDs in the study set 
        class_to_check : Category/class IRI to test for enrichment
        class_to_leaf_map : Mapping from class IRI to its leaf descendants
        weights_dict : {element_id: weight} for DETECTED elements only
            Background elements without weights are implicitly weight=0
        N_total : int
            Total population size (all background elements)
        classification : str
            "structural", "functional", or "full"
        class_to_all_roles_map : dict
            Mapping from class to all its roles (direct + inherited + role descendants)
        roles_to_leaves_map : dict
            Mapping from role IRI to leaf classes with that role
    
    Returns:
        Weighted enrichment p-value
    """
    # Get all leaf descendants of this category based on classification
    leaves = set()
    
    if classification in ["structural", "full"]:
        leaf_descendants_structural = set(class_to_leaf_map.get(class_to_check, []))
        leaves.update(leaf_descendants_structural)
    
    if classification in ["functional", "full"]:
        # Get all roles (direct + inherited from ancestors + role descendants)
        all_roles = class_to_all_roles_map.get(class_to_check, [])
        for role in all_roles:
            leaves.update(roles_to_leaves_map.get(role, []))
    
    leaf_descendants = leaves
    
    if len(leaf_descendants) == 0:
        return None  # No descendants, no enrichment
    
    # Find overlap: elements in both study set AND this category
    studyset_set = set(studyset_leaves)
    overlap = studyset_set & leaf_descendants
    
    # m = number of background elements in this category
    m = len(leaf_descendants)
    
    # n_observed = number of study set elements in this category
    n_observed = len(overlap)
    
    if n_observed == 0:
        return   # No overlap, no enrichment. Shouldn't happen due to earlier check.
    
    # Calculate sum of weights for the overlap
    sum_observed_weights = sum(weights_dict.get(element_id, 0.0) for element_id in overlap)
    
    # Find the saddlepoint parameter λ
    try:
        lambda_optimal = find_lambda(m, sum_observed_weights, weights_dict, N_total)
    except Exception as e:
        print(f"Error finding lambda for {class_to_check}: {e}")
        return None

    # Debug diagnostics
    if sum_observed_weights > 0:
        d1k0, _, _ = calc_cgf_derivatives(0.0, weights_dict, N_total)
        expected_sum = m * d1k0

        print("\n--- DEBUG ---")
        print(f"class: {class_to_check}")
        print(f"m (category size): {m}")
        print(f"n_observed: {n_observed}")
        print(f"sum_observed_weights: {sum_observed_weights}")
        print(f"population mean weight (K'(0)): {d1k0}")
        print(f"expected_sum under null: {expected_sum}")
        print(f"lambda_optimal: {lambda_optimal}")
        print("----------------")
        
    # Calculate p-value using Lugannani-Rice formula
    # Use observed study-set count (n_observed) like the Java implementation
    try:
        p_value = lugannani_rice_pvalue(lambda_optimal, n_observed, weights_dict, N_total)
    except Exception as e:
        print(f"Error calculating p-value for {class_to_check}: {e}")
        return None
        
    return p_value

def get_enrichment_values_with_weights(removed_leaves_csv, classification, 
                                       studyset_leaves, studyset_ancestors, 
                                       class_to_leaf_map, weights_dict,
                                       class_to_all_roles_map, roles_to_leaves_map):
    """
    Calculate weighted p-values for all categories.
    Optionally applies multiple-testing correction using weighted p-values.
    Returns dict with weighted p-values plus p_value/p_value_corrected.
    """
    
    enrichment_results = {}
    
    n_bg_leaves = count_removed_leaves(removed_leaves_csv)
    n_ss_leaves = len(studyset_leaves)

    
    for class_iri in studyset_ancestors:
        # Calculate counts
        n_ss_annotated = get_n_ss_annotated( # Just used for reporting, not for p-value calculation
            studyset_leaves,
            class_iri,
            class_to_leaf_map,
            classification,
            class_to_all_roles_map,
            roles_to_leaves_map,
        )
        
        _, n_bg_annotated = count_removed_classes_for_class( # Just used for reporting, not for p-value calculation
            class_iri,
            class_to_leaf_map,
            classification,
            class_to_all_roles_map,
            roles_to_leaves_map,
        )
        
        # Weighted test
        p_weighted = calculate_weighted_pvalue(
            studyset_leaves, class_iri, class_to_leaf_map,
            weights_dict, n_bg_leaves, classification,
            class_to_all_roles_map, roles_to_leaves_map
        )
        
        if p_weighted is None:
            print(f"Warning: p_weighted is None for {id_to_name(class_iri)}, skipping")
            continue
        
        # Calculate sum of weights in overlap for reporting
        # Use same logic as calculate_weighted_pvalue to get leaves
        leaves_for_sum = set()
        if classification in ["structural", "full"]:
            leaves_for_sum.update(class_to_leaf_map.get(class_iri, []))
        if classification in ["functional", "full"]:
            all_roles = class_to_all_roles_map.get(class_iri, [])
            for role in all_roles:
                leaves_for_sum.update(roles_to_leaves_map.get(role, []))
        
        overlap = set(studyset_leaves) & leaves_for_sum
        sum_weights = sum(weights_dict.get(leaf_id, 0.0) for leaf_id in overlap)
        
        enrichment_results[class_iri] = {
            'n_ss_annotated': n_ss_annotated,
            'n_ss_leaves': n_ss_leaves,
            'n_bg_annotated': n_bg_annotated,
            'n_bg_leaves': n_bg_leaves,
            'p_value_weighted': p_weighted,
            'sum_weights': sum_weights
        }
    return enrichment_results


def run_weighted_enrichment_analysis_plain_enrich_pruning_strategy(weights_dict,
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

    if not weights_dict:
        return {"study_set": [], "removed_nodes": [], "enrichment_results": {}}, None

    weights_dict = {normalize_id(cls): weight for cls, weight in weights_dict.items()} # Normalize IDs in weights_dict to ensure they match the format used in maps
    studyset_list = list(weights_dict.keys())

    studyset_leaves = get_leaves(studyset_list, removed_leaves_csv, class_to_leaf_map)
    print(f"Study set leaves: {studyset_leaves}")

    studyset_ancestors = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    print(f"Study set ancestors: {studyset_ancestors}")
    print(f"Number of study set ancestors: {len(studyset_ancestors)}")


    all_removed_nodes = set()

    enrichment_results = get_enrichment_values_with_weights(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        weights_dict,
        class_to_all_roles_map,
        roles_to_leaves_map,
    )

    for cls, vals in enrichment_results.items():
        vals["p_value"] = vals.get("p_value_weighted")
        vals["class"] = id_to_name(cls)

    print("Enrichment results:")
    print_enrichment_results(enrichment_results)

    enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
    print("Enrichment results after Benjamini-Hochberg correction:")
    print_enrichment_results(enrichment_results)

    pre_pruned_G = create_graph_from_map(studyset_leaves, parent_map_file, max_n_leaf_classes=inf)
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
    
    # Initialize current_enrichment before the loop (in case the loop doesn't execute)
    current_enrichment = {
        cls: vals for cls, vals in enrichment_results.items()
        if cls not in all_removed_nodes and G.has_node(cls)
    }

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


def run_weighted_enrichment_analysis(weights_dict,
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
    """
    Weighted enrichment analysis with flexible pruning options.
    Similar to fishers_calculations.run_enrichment_analysis but for weighted inputs.
    """

    pruning_before_enrichment = root_children_prune or linear_branch_prune

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

    if not weights_dict:
        return {"study_set": [], "removed_nodes": [], "enrichment_results": {}}, None

    weights_dict = {normalize_id(cls): weight for cls, weight in weights_dict.items()}
    studyset_list = list(weights_dict.keys())

    studyset_leaves = get_leaves(studyset_list, removed_leaves_csv, class_to_leaf_map)
    print(f"Study set leaves: {studyset_leaves}")

    studyset_ancestors_all = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    print(f"Study set ancestors: {studyset_ancestors_all}")
    print(f"Number of study set ancestors: {len(studyset_ancestors_all)}")

    G = create_graph_from_map(studyset_leaves, parent_map_file, max_n_leaf_classes=inf)
    pruned_G = G.copy()

    all_removed_nodes = set()

    if pruning_before_enrichment:
        if root_children_prune:
            print(f"Root children pruner activated, pruning {levels} levels from root")
            pruned_G, removed_nodes, execution_count = root_children_pruner(pruned_G, levels, allow_re_execution = False, execution_count = 0)
            print(f"Removed nodes by root children pruner: {removed_nodes}")
            all_removed_nodes.update(removed_nodes)

        if linear_branch_prune:
            print(f"Linear branch pruner activated, keeping only every {n}-th node in linear branches")
            pruned_G, removed_nodes = linear_branch_collapser_pruner_remove_less(pruned_G, n)
            print(f"Removed nodes by linear branch pruner: {removed_nodes}")
            all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from studyset_ancestors_all
        studyset_ancestors = [cls for cls in studyset_ancestors_all if cls not in all_removed_nodes]
        print(f"Number of study set ancestors after before-enrichment pruning: {len(studyset_ancestors)}")
    
    else:
        studyset_ancestors = studyset_ancestors_all

    enrichment_results = get_enrichment_values_with_weights(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        weights_dict,
        class_to_all_roles_map,
        roles_to_leaves_map,
    )

    for cls, vals in enrichment_results.items():
        vals["p_value"] = vals.get("p_value_weighted")
        vals["class"] = id_to_name(cls)

    print("Enrichment results:")
    print_enrichment_results(enrichment_results)

    if bonferroni_correct:
        print("Applying Bonferroni correction to p-values...")
        enrichment_results, _ = bonferroni_correction(enrichment_results)
    elif benjamini_hochberg_correct:
        print("Applying Benjamini-Hochberg FDR correction to p-values...")
        enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)

    if high_p_value_prune:
        print(f"High p-value pruner activated, pruning nodes with p-value above {p_value_threshold}")
        pruned_G, removed_nodes = high_p_value_branch_pruner(pruned_G, enrichment_results, p_value_threshold)
        print(f"Removed nodes by high p-value pruner: {removed_nodes}")
        all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from enrichment results
        for cls in removed_nodes:
            if cls in enrichment_results:
                del enrichment_results[cls]

    if zero_degree_prune:
        print("Applying zero-degree pruner to remove nodes with zero degree...")
        pruned_G, removed_nodes = zero_degree_pruner(pruned_G)
        print(f"Removed nodes by zero-degree pruner: {removed_nodes}")
        all_removed_nodes.update(removed_nodes)

        # Remove pruned nodes from enrichment results
        for cls in removed_nodes:
            if cls in enrichment_results:
                del enrichment_results[cls]

    print("Final enrichment results:")
    print_enrichment_results(enrichment_results)

    print(f"Number of removed nodes in total: {len(all_removed_nodes)}")
    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {id_to_name(cls): vals for cls, vals in enrichment_results.items()}
    }
    return results, pruned_G

# TODO: Adjust so input wth weights works as intended. Maybe a change in website files?
# TODO: Make sure it functions as intended with weights input.
# TODO: Add function to website
# TODO: Check so that it works for all classifications


# Autoscaling weights to see if makes p-values more significant. Maybe not sensible but worth a try.
def auto_scale_weights(weights_dict, target_max=100):
    """
    Automatically scale weights so the maximum is target_max.
    
    Only scales up if current max is below target_max.
    This ensures the weighted test has enough signal to detect enrichment.
    """
    if not weights_dict:
        return {}
    
    current_max = max(weights_dict.values())
    
    # Only scale if current max is below target
    if current_max < target_max:
        scale_factor = target_max / current_max
        print(f"Auto-scaling weights by factor {scale_factor:.2f} to set max to {target_max}")
        return {k: v * scale_factor for k, v in weights_dict.items()}
    
    # If already >= target_max, return unchanged
    print(f"No scaling applied. Current max weight {current_max:.2f} is >= target max {target_max}.")
    return weights_dict




if __name__ == "__main__":


    classification = "structural"  # "structural", "functional", or "full"
    removed_leaves_csv = "data/removed_leaf_classes_with_smiles.csv"

    weights_dict = {
        "http://purl.obolibrary.org/obo/CHEBI_17079": 0.7665,
        "http://purl.obolibrary.org/obo/CHEBI_46816": 0.7465,
        "http://purl.obolibrary.org/obo/CHEBI_28658": 0.7465,
        "http://purl.obolibrary.org/obo/CHEBI_28611": 0.7465,
        "http://purl.obolibrary.org/obo/CHEBI_28594": 0.6915,
        "http://purl.obolibrary.org/obo/CHEBI_17048": 0.6915,
        "http://purl.obolibrary.org/obo/CHEBI_7852": 0.60575,
        "http://purl.obolibrary.org/obo/CHEBI_164200": 0.2342,
        "http://purl.obolibrary.org/obo/CHEBI_8489": 0.25321,
        "http://purl.obolibrary.org/obo/CHEBI_9630": 0.2543,
        "http://purl.obolibrary.org/obo/CHEBI_59477": 0.2335,
        "http://purl.obolibrary.org/obo/CHEBI_9495": 0.2433,
        "http://purl.obolibrary.org/obo/CHEBI_3540": 0.509,
    }

    # Optionally, auto-scale weights
    # weights_dict = auto_scale_weights(weights_dict, target_max=100)
    studyset_leaves = list(weights_dict.keys())
    
    # Load class to leaf descendants map
    class_to_leaf_map_file = "data/class_to_leaf_descendants_map.json"
    with open(class_to_leaf_map_file, "r") as f:
        class_to_leaf_map = json.load(f)

    # Load class to all roles map
    class_to_all_roles_json = "data/class_to_all_roles_map.json"
    with open(class_to_all_roles_json, "r") as f:
        class_to_all_roles_map = json.load(f)

    # Load roles to leaves map
    roles_to_leaves_map_json = "data/roles_to_leaves_map.json"
    with open(roles_to_leaves_map_json, "r") as f:
        roles_to_leaves_map = json.load(f)

    # Get all ancestors of the study set leaves
    leaf_to_ancestors_map_file = "data/removed_leaf_classes_to_ALL_parents_map.json"
    studyset_ancestors = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    
    print(f"Study set: {len(studyset_leaves)} leaves")
    print(f"Ancestors: {len(studyset_ancestors)} classes")
    
    # Calculate enrichment values
    enrichment_results = get_enrichment_values_with_weights(
        removed_leaves_csv, classification,
        studyset_leaves, studyset_ancestors,
        class_to_leaf_map, weights_dict,
        class_to_all_roles_map, roles_to_leaves_map
    )
    
    print(f"\nFound {len(enrichment_results)} enriched classes:")
    for class_iri, results in enrichment_results.items():
        print(f"\nClass: {class_iri}")
        if results['p_value_weighted'] is not None:
            print(f"  Weighted p-value: {results['p_value_weighted']:.4e}")
        else:
            print(f"  Weighted p-value: N/A (calculation failed)")
            # print(f"  Sum of weights in overlap: {results['sum_weights']:.2f}")
