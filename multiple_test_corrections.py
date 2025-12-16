from math import inf
import pandas as pd
import copy 



def bonferroni_correction(enrichment_results):

    """ Performs Bonferroni correction on the p-values in the enrichment_results dictionary to adjust for multiple hypothesis testing."""

    adjusted_results = copy.deepcopy(enrichment_results)

    # Extract raw p-values
    class_p_list = [(cls, res['p_value']) for cls, res in enrichment_results.items()]


    m = len(class_p_list)  # number of tests
    correction_map = {}

    for cls, raw_p in class_p_list:
        if raw_p is None:
            corrected_p = None
        else:
            corrected_p = min(raw_p * m, 1.0)  # Bonferroni correction

        adjusted_results[cls]['p_value_corrected'] = corrected_p 
        correction_map[cls] = corrected_p # Maybe not necessary to have both adjusted_results and correction_map

    return adjusted_results, correction_map


def benjamini_hochberg_fdr_correction(enrichment_results):
    """ Performs Benjamini-Hochberg FDR (false discovery rate) correction on the p-values in the enrichment_results dictionary to adjust for multiple hypothesis testing."""
    # Extract p-values into a sorted list
    ordered_p_values = [(cls, info["p_value"]) for cls, info in enrichment_results.items()]
    ordered_p_values.sort(key=lambda x: x[1])  # Sort by p-value ascending

    m = len(ordered_p_values)  # number of tests
    adj = [None] * m  # list for adjusted p-values

    # Apply Benjamini-Hochberg procedure
    running_min = 1.0 # p-value cannot be higher than 1
    for i in range(m, 0, -1): 
        raw_p = ordered_p_values[i - 1][1]
        bh_value = (raw_p * m) / i
        running_min = min(running_min, bh_value)
        adj[i - 1] = running_min

    # Map adjusted p-values back into the results dictionary
    for (cls, _), corrected_p in zip(ordered_p_values, adj):
        enrichment_results[cls]["p_value_corrected"] = corrected_p

    return enrichment_results
