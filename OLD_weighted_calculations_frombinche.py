"""
weighted_calculations.py

SaddleSum (Lugannani-Rice) enrichment for ChEBI ontology.
Mirrors fishers_calculations.py exactly — same function signatures,
same file/data dependencies, same pruning pipeline.

Reference: Stojmirovic & Yu (2010), ArXiv:1004.5088
"""

import json
import math
import time
import numpy as np
import pandas as pd
from scipy.special import ndtr

from visualitations_and_pruning import (
    root_children_pruner,
    linear_branch_collapser_pruner_remove_less,
    high_p_value_branch_pruner,
    zero_degree_pruner,
    create_graph_with_roles_and_structures,
    id_to_name,
)
from fishers_calculations import (
    get_leaves,
    get_ancestors_for_inputs,
    normalize_id,
    print_enrichment_results,
)
from pre_fishers_calculations import (
    count_removed_leaves,
    count_removed_classes_for_class,
    count_removed_classes_for_roles,
)
from multiple_test_corrections import (
    bonferroni_correction,
    benjamini_hochberg_fdr_correction,
)


# ---------------------------------------------------------------------------
# Saddlepoint engine  (translated directly from saddlesum.c)
# ---------------------------------------------------------------------------

class _SaddleSum:
    MAX_ITERS = 50
    TOLERANCE = 1.0e-11

    def __init__(self, background_weights):
        w = np.asarray(background_weights, dtype=float)
        self._weights = np.sort(w)
        self._N = len(w)
        self._mean = float(w.mean())
        self._wmax = float(self._weights[-1])
        self._cache = []

    def _compute_item(self, lmbd):
        w = self._weights
        tmp = np.exp(lmbd * (w - self._wmax))
        Nrho  = tmp.sum()
        Nrho1 = (tmp * w).sum()
        Nrho2 = (tmp * w * w).sum()

        if Nrho <= 0.0:
            return None

        D1K = Nrho1 / Nrho
        D2K = Nrho2 / Nrho - D1K ** 2

        if D2K <= 0.0:
            return None

        expH = Nrho * np.exp(lmbd * (self._wmax - D1K)) / self._N
        C = 2.0 * lmbd * math.sqrt(D2K)

        inner = lmbd * (D1K - self._wmax) - math.log(Nrho) + math.log(self._N)
        if inner < 0.0:
            return None

        D = -math.sqrt(2.0) * math.sqrt(inner)
        return dict(lambda_=lmbd, mean=D1K, D2K=D2K, expH=expH, C=C, D=D)

    @staticmethod
    def _item_pvalue(item, m):
        sqrtm = math.sqrt(m)
        if item['D'] * sqrtm > -1.0:
            return 1.0
        # Guard against C or D being zero (divide by zero)
        if item['C'] == 0.0 or item['D'] == 0.0:
            return 1.0
        phi = math.sqrt(2.0 / math.pi) * (item['expH'] ** m)
        result = (float(ndtr(item['D'] * sqrtm))
                  + phi / item['C'] / sqrtm
                  + phi / item['D'] / sqrtm / 2.0)
        # Guard against nan/inf from numerical edge cases
        if not math.isfinite(result):
            return 1.0
        return result

    def _bisect(self, x):
        lo, hi = 0, len(self._cache)
        while lo < hi:
            mid = (lo + hi) // 2
            if self._cache[mid]['mean'] < x:
                lo = mid + 1
            else:
                hi = mid
        return lo

    def _insert(self, item):
        self._cache.insert(self._bisect(item['mean']), item)

    def pvalue(self, score, num_hits):
        x = score / num_hits
        if x <= self._mean:
            return 1.0

        min_pval = (1.0 / self._N) ** num_hits
        i = self._bisect(x)
        ya = self._cache[i - 1]['lambda_'] if i > 0 else 0.0 

        if i < len(self._cache):
            item = self._cache[i]
            yb = item['lambda_']
            yc = 0.5 * (ya + yb)
            pval = self._item_pvalue(item, num_hits)
            if (yb - ya) < self.TOLERANCE:
                return min(max(pval, min_pval), 1.0)
        else:
            yb = 0.5
            while True:
                yb *= 2.0
                if yb > 1e6:
                    return min_pval
                item = self._compute_item(yb)
                if item is None:
                    return min_pval
                self._insert(item)
                if abs(item['mean'] - self._wmax) < self.TOLERANCE:
                    break
            yc = 0.5 * (ya + yb)

        pval = 1.0
        for _ in range(self.MAX_ITERS):
            item = self._compute_item(yc)
            if item is None:
                break
            pval = self._item_pvalue(item, num_hits)
            diff = item['mean'] - x
            if diff < 0.0:
                ya = yc
            else:
                yb = yc
            y = yc - diff / item['D2K']
            if y < ya or y > yb:
                y = 0.5 * (ya + yb)
            if abs(y - yc) < self.TOLERANCE or abs(diff) < self.TOLERANCE:
                break
            self._insert(item)
            yc = y

        return min(max(pval, min_pval), 1.0)


# ---------------------------------------------------------------------------
# Weight helpers
# ---------------------------------------------------------------------------

def auto_scale_weights(weights_dict, target_max=1000.0):
    """Scale weights so max absolute value equals target_max (if below it)."""
    if not weights_dict:
        return weights_dict
    max_abs = max(abs(v) for v in weights_dict.values())
    if max_abs == 0 or max_abs >= target_max:
        return weights_dict
    scale = target_max / max_abs
    return {k: v * scale for k, v in weights_dict.items()}


def _propagate_weights_to_leaves(weights_dict, class_to_leaf_map, removed_leaves_csv):
    """
    Build weights_with_leaves: expand any non-leaf submitted IDs so their
    leaf descendants inherit the weight. Directly submitted leaves take
    priority over inherited values.
    Returns weights_with_leaves dict keyed by leaf IRI.
    """
    leaves_df = pd.read_csv(removed_leaves_csv)
    leaf_set = set(leaves_df['IRI'].values)

    weights_with_leaves = {}

    # First pass: propagate non-leaf weights to descendants
    for cls, weight in weights_dict.items():
        if cls not in leaf_set:
            for leaf in class_to_leaf_map.get(cls, []):
                if leaf not in weights_with_leaves:
                    weights_with_leaves[leaf] = weight

    # Second pass: directly submitted leaves always overwrite inherited values
    for cls, weight in weights_dict.items():
        if cls in leaf_set:
            weights_with_leaves[cls] = weight

    return weights_with_leaves


def _build_background_weights(weights_with_leaves, removed_leaves_csv):
    """
    Build the full background weight array.
    Every leaf in the background gets its weight from weights_with_leaves
    (if measured) or 0.0 otherwise.

    Must be called with weights_with_leaves (post-propagation), not
    the raw weights_dict.
    """
    leaves_df = pd.read_csv(removed_leaves_csv)
    all_bg_leaves = list(leaves_df['IRI'].values)
    return [weights_with_leaves.get(leaf, 0.0) for leaf in all_bg_leaves]


# ---------------------------------------------------------------------------
# calculate_weighted_p_value  (mirrors calculate_p_value in fishers_calculations)
# ---------------------------------------------------------------------------

def calculate_weighted_p_value(saddler, term_leaves, studyset_leaves_set, weights_with_leaves):
    """
    Compute the SaddleSum p-value for a single term.
    Returns (score, n_ss_annotated, p_value).
    """
    annotated = [
        leaf for leaf in term_leaves
        if leaf in studyset_leaves_set and leaf in weights_with_leaves
    ]
    n_ss_annotated = len(annotated)
    if n_ss_annotated == 0:
        return 0.0, 0, 1.0

    score = sum(weights_with_leaves[leaf] for leaf in annotated)
    p_value = saddler.pvalue(score, n_ss_annotated)
    return score, n_ss_annotated, p_value


# ---------------------------------------------------------------------------
# get_weighted_enrichment_values  (mirrors get_enrichment_values)
# ---------------------------------------------------------------------------

def get_weighted_enrichment_values(
    removed_leaves_csv,
    classification,
    studyset_leaves,
    studyset_ancestors,
    class_to_leaf_map,
    class_to_all_roles_map,
    roles_to_leaves_map,
    studyset_ancestors_roles,
    weights_with_leaves,
):
    """
    Compute SaddleSum enrichment for every ancestor class (and role class).
    Stores ALL results — BH/Bonferroni correction filters them afterwards,
    exactly as get_enrichment_values does for Fisher.
    """
    background_weights = _build_background_weights(weights_with_leaves, removed_leaves_csv)
    saddler = _SaddleSum(background_weights)

    studyset_leaves_set = set(studyset_leaves)
    n_bg_leaves = count_removed_leaves(removed_leaves_csv)
    n_ss_leaves = len(studyset_leaves)

    results = {}

    if classification in ["structural", "full"]:
        for cls in studyset_ancestors:
            term_leaves = set(class_to_leaf_map.get(cls, []))
            score, n_ss_annotated, p_value = calculate_weighted_p_value(
                saddler, term_leaves, studyset_leaves_set, weights_with_leaves
            )
            _, n_bg_annotated = count_removed_classes_for_class(
                cls, class_to_leaf_map, classification,
                class_to_all_roles_map, roles_to_leaves_map,
            )
            results[cls] = {
                "class": id_to_name(cls),
                "score": score,
                "n_ss_annotated": n_ss_annotated,
                "n_ss_leaves": n_ss_leaves,
                "n_bg_annotated": n_bg_annotated,
                "n_bg_leaves": n_bg_leaves,
                "odds_ratio": float('inf') if n_ss_annotated > 0 else 0.0,
                "p_value": p_value,
            }

    if classification in ["functional", "full"] and studyset_ancestors_roles:
        print(f"Calculating weighted enrichment for {len(studyset_ancestors_roles)} roles...")
        for role in studyset_ancestors_roles:
            term_leaves = set(roles_to_leaves_map.get(role, []))
            score, n_ss_annotated, p_value = calculate_weighted_p_value(
                saddler, term_leaves, studyset_leaves_set, weights_with_leaves
            )
            _, n_bg_annotated = count_removed_classes_for_roles(
                role, class_to_leaf_map, classification, roles_to_leaves_map
            )
            results[role] = {
                "class": id_to_name(role),
                "score": score,
                "n_ss_annotated": n_ss_annotated,
                "n_ss_leaves": n_ss_leaves,
                "n_bg_annotated": n_bg_annotated,
                "n_bg_leaves": n_bg_leaves,
                "odds_ratio": float('inf') if n_ss_annotated > 0 else 0.0,
                "p_value": p_value,
            }

    return results


# ---------------------------------------------------------------------------
# Shared setup logic
# ---------------------------------------------------------------------------

def _setup_weighted_analysis(weights_dict, classification,
                             removed_leaves_csv, leaf_to_ancestors_map_file,
                             class_to_leaf_map, class_to_all_roles_map,
                             roles_to_leaves_map, parent_map_file):
    """Extract leaves, propagate weights, find ancestors and roles."""
    normalized_weights_dict = {
        normalize_id(cls): weight for cls, weight in weights_dict.items()
    }

    studyset_list = list(normalized_weights_dict.keys())
    studyset_leaves = get_leaves(studyset_list, removed_leaves_csv, class_to_leaf_map)
    print(f"Study set leaves: {len(studyset_leaves)}")

    weights_with_leaves = _propagate_weights_to_leaves(
        normalized_weights_dict, class_to_leaf_map, removed_leaves_csv
    )
    print(f"Weights with leaf propagation: {len(weights_with_leaves)} leaves have weights")

    studyset_ancestors_all = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    print(f"Number of study set ancestors: {len(studyset_ancestors_all)}")

    if classification in ["functional", "full"]:
        studyset_ancestors_roles = set()
        for leaf in studyset_leaves:
            studyset_ancestors_roles.update(class_to_all_roles_map.get(leaf, []))
        for cls in studyset_ancestors_all:
            studyset_ancestors_roles.update(class_to_all_roles_map.get(cls, []))

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

    return studyset_leaves, weights_with_leaves, studyset_ancestors_all, studyset_ancestors_roles


# ---------------------------------------------------------------------------
# run_weighted_enrichment_analysis  (mirrors run_enrichment_analysis)
# ---------------------------------------------------------------------------

def run_weighted_enrichment_analysis(
    weights_dict,
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
):
    removed_leaves_csv          = "data/removed_leaf_classes_with_smiles.csv"
    leaf_to_ancestors_map_file  = "data/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_map_file      = "data/class_to_leaf_descendants_map.json"
    parent_map_file             = "data/chebi_parent_map.json"
    class_to_all_roles_map_json = "data/class_to_all_roles_map.json"
    roles_to_leaves_map_json    = "data/roles_to_leaves_map.json"

    with open(class_to_leaf_map_file, 'r') as f:
        class_to_leaf_map = json.load(f)
    with open(class_to_all_roles_map_json, 'r') as f:
        class_to_all_roles_map = json.load(f)
    with open(roles_to_leaves_map_json, 'r') as f:
        roles_to_leaves_map = json.load(f)

    studyset_leaves, weights_with_leaves, studyset_ancestors_all, studyset_ancestors_roles = (
        _setup_weighted_analysis(
            weights_dict, classification,
            removed_leaves_csv, leaf_to_ancestors_map_file,
            class_to_leaf_map, class_to_all_roles_map,
            roles_to_leaves_map, parent_map_file,
        )
    )

    G = create_graph_with_roles_and_structures(
        studyset_leaves, studyset_ancestors_all,
        studyset_ancestors_roles, parent_map_file,
        class_to_all_roles_map, classification,
    )
    pruned_G = G.copy()
    all_removed_nodes = set()

    pruning_before_enrichment = root_children_prune or linear_branch_prune
    if pruning_before_enrichment:
        if root_children_prune:
            print(f"Root children pruner activated, pruning {levels} levels from root")
            t0 = time.time()
            pruned_G, removed_nodes, _ = root_children_pruner(
                pruned_G, levels, allow_re_execution=False, execution_count=0
            )
            print(f"Root children pruning: {time.time()-t0:.2f}s")
            all_removed_nodes.update(removed_nodes)

        if linear_branch_prune:
            print(f"Linear branch pruner activated, n={n}")
            pruned_G, removed_nodes = linear_branch_collapser_pruner_remove_less(pruned_G, n)
            all_removed_nodes.update(removed_nodes)

        studyset_ancestors = [c for c in studyset_ancestors_all if c not in all_removed_nodes]
        print(f"Ancestors after pre-enrichment pruning: {len(studyset_ancestors)}")
    else:
        studyset_ancestors = studyset_ancestors_all

    enrichment_results = get_weighted_enrichment_values(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        studyset_ancestors_roles,
        weights_with_leaves,
    )

    if bonferroni_correct:
        print("Applying Bonferroni correction...")
        enrichment_results, _ = bonferroni_correction(enrichment_results)
    elif benjamini_hochberg_correct:
        print("Applying Benjamini-Hochberg FDR correction...")
        enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)

    if high_p_value_prune:
        print(f"High p-value pruner, threshold={p_value_threshold}")
        t0 = time.time()
        pruned_G, removed_nodes = high_p_value_branch_pruner(
            pruned_G, enrichment_results, p_value_threshold
        )
        print(f"High p-value pruning: {time.time()-t0:.2f}s, removed {len(removed_nodes)}")
        all_removed_nodes.update(removed_nodes)
        for cls in removed_nodes:
            enrichment_results.pop(cls, None)

    if zero_degree_prune:
        print("Zero-degree pruner activated...")
        t0 = time.time()
        pruned_G, removed_nodes = zero_degree_pruner(pruned_G)
        print(f"Zero-degree pruning: {time.time()-t0:.2f}s, removed {len(removed_nodes)}")
        all_removed_nodes.update(removed_nodes)
        for cls in removed_nodes:
            enrichment_results.pop(cls, None)

    print("Final weighted enrichment results:")
    print_enrichment_results(enrichment_results)
    print(f"Total removed nodes: {len(all_removed_nodes)}")

    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {
            id_to_name(cls): vals for cls, vals in enrichment_results.items()
        },
    }
    return results, pruned_G


# ---------------------------------------------------------------------------
# run_weighted_enrichment_analysis_plain_enrich_pruning_strategy
# ---------------------------------------------------------------------------

def run_weighted_enrichment_analysis_plain_enrich_pruning_strategy(
    weights_dict,
    levels=2,
    n=0,
    p_value_threshold=0.05,
    classification="structural",
):
    removed_leaves_csv          = "data/removed_leaf_classes_with_smiles.csv"
    leaf_to_ancestors_map_file  = "data/removed_leaf_classes_to_ALL_parents_map.json"
    class_to_leaf_map_file      = "data/class_to_leaf_descendants_map.json"
    parent_map_file             = "data/chebi_parent_map.json"
    class_to_all_roles_map_json = "data/class_to_all_roles_map.json"
    roles_to_leaves_map_json    = "data/roles_to_leaves_map.json"

    with open(class_to_leaf_map_file, 'r') as f:
        class_to_leaf_map = json.load(f)
    with open(class_to_all_roles_map_json, 'r') as f:
        class_to_all_roles_map = json.load(f)
    with open(roles_to_leaves_map_json, 'r') as f:
        roles_to_leaves_map = json.load(f)

    studyset_leaves, weights_with_leaves, studyset_ancestors, studyset_ancestors_roles = (
        _setup_weighted_analysis(
            weights_dict, classification,
            removed_leaves_csv, leaf_to_ancestors_map_file,
            class_to_leaf_map, class_to_all_roles_map,
            roles_to_leaves_map, parent_map_file,
        )
    )

    all_removed_nodes = set()

    enrichment_results = get_weighted_enrichment_values(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        studyset_ancestors_roles,
        weights_with_leaves,
    )

    print("Enrichment results (raw):")
    print_enrichment_results(enrichment_results)

    enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
    print("After BH correction:")
    print_enrichment_results(enrichment_results)

    G = create_graph_with_roles_and_structures(
        studyset_leaves, studyset_ancestors,
        studyset_ancestors_roles, parent_map_file,
        class_to_all_roles_map, classification,
    )

    print("Starting pre-loop pruning phase.")
    G, removed_nodes = high_p_value_branch_pruner(G, enrichment_results, p_value_threshold)
    all_removed_nodes.update(removed_nodes)
    print(f"High p-value pruner removed: {len(removed_nodes)}")

    G, removed_nodes = linear_branch_collapser_pruner_remove_less(G, n)
    all_removed_nodes.update(removed_nodes)
    print(f"Linear branch pruner removed: {len(removed_nodes)}")

    G, removed_nodes, _ = root_children_pruner(
        G, levels, allow_re_execution=False, execution_count=0
    )
    all_removed_nodes.update(removed_nodes)
    print(f"Root children pruner removed: {len(removed_nodes)}")

    print("Starting loop pruning phase.")
    size_before = G.number_of_nodes()
    size_after = size_before
    first_iteration = True
    iteration = 0
    current_enrichment = enrichment_results

    while size_after < size_before or first_iteration:
        size_before = size_after
        iteration += 1
        print(f"Loop iteration {iteration}")

        current_enrichment = {
            cls: vals for cls, vals in enrichment_results.items()
            if cls not in all_removed_nodes and G.has_node(cls)
        }
        current_enrichment = benjamini_hochberg_fdr_correction(current_enrichment)

        G, removed_nodes = high_p_value_branch_pruner(G, current_enrichment, p_value_threshold)
        all_removed_nodes.update(removed_nodes)
        print(f"High p-value pruner removed: {len(removed_nodes)}")

        G, removed_nodes = zero_degree_pruner(G)
        all_removed_nodes.update(removed_nodes)
        print(f"Zero-degree pruner removed: {len(removed_nodes)}")

        size_after = G.number_of_nodes()
        first_iteration = False

    print(f"Total removed nodes: {len(all_removed_nodes)}")
    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {
            id_to_name(cls): vals for cls, vals in current_enrichment.items()
        },
    }
    return results, G