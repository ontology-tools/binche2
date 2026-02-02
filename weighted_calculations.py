import numpy as np
from scipy.special import erf
from scipy.optimize import brentq, newton

"""Using SaddleSum for weighted enrichment calculations"""

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
    if len(weights_dict) == 0:
        # Handle edge case of no weights
        return 0.0, 1.0, 1.0
    
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
    if n_without_weights > 0:
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
    return m * d1k - sum_weights


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
    # Step 1: Bisection to find initial λ in range [0.00001, 15.0]
    try:
        lambda_init = brentq(
            lambda l: saddlepoint_equation(l, m, sum_weights, weights_dict, N_tot),
            0.00001,
            15.0,
            maxiter=100
        )
    except ValueError:
        # Bisection failed - equation may not have root in this range
        print(f"Warning: Bisection failed for m={m}, sum_weights={sum_weights}")
        lambda_init = 0.1  # Fallback initial guess

    # Step 2: Newton-Raphson refinement for higher precision
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
            Number of elements in the category (from background)
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
        return 1.0  # Invalid configuration, return non-significant

    D = -np.sqrt(2) * np.sqrt(D_argument)

    # phi: Correction factor based on saddlepoint density
    phi = np.sqrt(2 / np.pi) * (exp_h ** m)

    # z: Standardized test statistic for normal CDF
    z = (D * np.sqrt(m)) / np.sqrt(2)

    # Calculate p-value only if D * sqrt(m) <= -1 to ensure we're in the enrichment tail
    if D * np.sqrt(m) <= -1:
        # Normal CDF using error function
        # erf(z) = 2/√π * integral from 0 to z of exp(-t²) dt
        # Standard normal CDF = (1 + erf(z/√2)) / 2
        ndtr = (1 + erf(z)) / 2
        
        # Lugannani-Rice correction terms
        # These adjust for skewness and kurtosis of the distribution
        correction1 = phi / C / np.sqrt(m)
        correction2 = phi / D / np.sqrt(m) / 2
        
        # Final p-value
        p_value = ndtr + correction1 + correction2
        
        # Handle numerical issues
        if np.isnan(p_value) or np.isinf(p_value):
            p_value = 1.0
    else:
        # Not in the enrichment tail (D * sqrt(m) > -1)
        p_value = 1.0
    
    return p_value





####################### SE OVER HARIFRAN #############################
def calculate_weighted_pvalue(studyset_leaves, class_to_check, class_to_leaf_map,
                              weights_dict, N_total):
    """
    Main function to calculate weighted enrichment p-value for a single category.
    
    This integrates all the saddlepoint calculations to produce a final p-value
    similar to Fisher's exact test, but accounting for element weights.
    
    Parameters:
        studyset_leaves : list or set
            Element IDs in your study set (detected compounds)
        class_to_check : str
            Category/class IRI to test for enrichment
        class_to_leaf_map : dict
            Mapping from class IRI to its leaf descendants
        weights_dict : dict
            {element_id: weight} for DETECTED elements only
            Background elements without weights are implicitly weight=0
        N_total : int
            Total population size (all background elements)
    
    Returns:
        float : Weighted enrichment p-value
    """
    # Get all leaf descendants of this category
    leaf_descendants = set(class_to_leaf_map.get(class_to_check, []))
    
    if len(leaf_descendants) == 0:
        return 1.0  # No descendants, no enrichment
    
    # Find overlap: elements in both study set AND this category
    studyset_set = set(studyset_leaves)
    overlap = studyset_set & leaf_descendants
    
    # m = number of background elements in this category
    m = len(leaf_descendants)
    
    # n_observed = number of study set elements in this category
    n_observed = len(overlap)
    
    if n_observed == 0:
        return 1.0  # No overlap, no enrichment
    
    # Calculate sum of weights for the overlap
    sum_observed_weights = sum(weights_dict.get(element_id, 0.0) for element_id in overlap)
    
    if sum_observed_weights == 0:
        return 1.0  # No weight in overlap
    
    # Find the saddlepoint parameter λ
    try:
        lambda_optimal = find_lambda(m, sum_observed_weights, weights_dict, N_total)
    except Exception as e:
        print(f"Error finding lambda for {class_to_check}: {e}")
        return 1.0
    
    # Calculate p-value using Lugannani-Rice formula
    try:
        p_value = lugannani_rice_pvalue(lambda_optimal, m, weights_dict, N_total)
    except Exception as e:
        print(f"Error calculating p-value for {class_to_check}: {e}")
        return 1.0
    
    # Ensure p-value is in valid range
    return min(max(p_value, 0.0), 1.0)

def get_enrichment_values_with_weights(removed_leaves_csv, classification, 
                                       studyset_leaves, studyset_ancestors, 
                                       class_to_leaf_map, weights_dict,
                                       check_leaf_classes=False):
    """
    Calculate both Fisher's exact and weighted p-values for all categories.
    
    Returns dict with both p-values for comparison.
    """
    from fishers_calculations import calculate_p_value, count_removed_leaves
    
    enrichment_results = {}
    
    # Get background size
    n_bg_leaves = count_removed_leaves(removed_leaves_csv, classification)
    n_ss_leaves = len(studyset_leaves)
    
    for class_iri in studyset_ancestors:
        # Calculate counts
        n_ss_annotated = get_n_ss_annotated(studyset_leaves, class_iri, class_to_leaf_map)
        
        if n_ss_annotated == 0:
            continue
        
        n_bg_annotated, _ = count_removed_classes_for_class(
            class_iri, class_to_leaf_map, classification, 
            check_leaf_classes, removed_leaves_csv
        )
        
        if n_bg_annotated == 0:
            continue
        
        # Fisher's exact test
        odds, p_fisher = calculate_p_value(
            n_ss_annotated, n_ss_leaves, n_bg_annotated, n_bg_leaves
        )
        
        # Weighted test
        p_weighted = calculate_weighted_pvalue(
            studyset_leaves, class_iri, class_to_leaf_map,
            weights_dict, n_bg_leaves
        )
        
        # Calculate sum of weights in overlap for reporting
        leaf_descendants = set(class_to_leaf_map.get(class_iri, []))
        overlap = set(studyset_leaves) & leaf_descendants
        sum_weights = sum(weights_dict.get(leaf_id, 0.0) for leaf_id in overlap)
        
        enrichment_results[class_iri] = {
            'n_ss_annotated': n_ss_annotated,
            'n_ss_leaves': n_ss_leaves,
            'n_bg_annotated': n_bg_annotated,
            'n_bg_leaves': n_bg_leaves,
            'odds_ratio': odds,
            'p_value_fisher': p_fisher,
            'p_value_weighted': p_weighted,
            'sum_weights': sum_weights
        }
    
    return enrichment_results