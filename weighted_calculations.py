"""
weighted_calculations.py

SaddleSum (Lugannani-Rice) enrichment for ChEBI ontology.
Mirrors fishers_calculations.py exactly — same function signatures,
same file/data dependencies, same pruning pipeline.

Reference: Stojmirovic & Yu (2010), ArXiv:1004.5088
"""

import json
import logging
import math
import time
from collections import Counter
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
    get_structural_leaf_ids,
)
from multiple_test_corrections import (
    bonferroni_correction,
    benjamini_hochberg_fdr_correction,
)
from wikidata.narrow_background_fishers import (
    get_studyset_leaves_narrow,
    count_narrow_leaves,
    count_narrow_leaves_for_class,
    count_narrow_leaves_for_role,
)


# ---------------------------------------------------------------------------
# Saddlepoint engine  (main code translated from saddlesum.c, public domain code by Stojmirovic)

# A. Stojmirovic and Y.-K. Yu. Robust and accurate data enrichment statistics via distribution function of sum of weights. _Bioinformatics_, **26** (21):2752-2759, 2010.
# ---------------------------------------------------------------------------

logger = logging.getLogger(__name__)

class _SaddleSum:
    """
    Precomputes background distribution, answers pvalue queries via
    the Lugannani-Rice saddlepoint approximation.
    """
    MAX_ITERS = 50
    TOLERANCE = 1.0e-11
    MAX_WARNING_EXAMPLES_PER_REASON = 5

    def __init__(self, background_weights):
        w = np.asarray(background_weights, dtype=float)
        self._N = len(w)                   # total number of leaves in background
        self._mean = float(w.mean())       # mean over ALL weights including zeros
        self._wmax = float(w.max())        # max over ALL weights
        self._cache = []
        self._fallback_counts = Counter()  # Track numerical fallback reasons

        # Split into measured (non-zero) and unmeasured (zero) — after computing
        # mean and wmax which must reflect the full background population
        w_sorted = np.sort(w)
        nonzero_mask = w_sorted != 0.0
        self._weights = w_sorted[nonzero_mask]         # sorted non-zero weights only
        self._n_zeros = int((~nonzero_mask).sum())     # count of zero-weight leaves

    def _debug_weight_summary(self): # For logging numerical fallback contexts. Addition.
        if self._weights.size == 0:
            return {
                "N": self._N,
                "nonzero_count": 0,
                "zero_count": self._n_zeros,
                "mean_all": self._mean,
                "wmax": self._wmax,
            }

        return {
            "N": self._N,
            "nonzero_count": int(self._weights.size),
            "zero_count": self._n_zeros,
            "mean_all": self._mean,
            "wmax": self._wmax,
            "nonzero_min": float(self._weights[0]),
            "nonzero_max": float(self._weights[-1]),
            "nonzero_mean": float(np.mean(self._weights)),
            "nonzero_std": float(np.std(self._weights)),
        }

    def _log_numerical_fallback(self, reason, score, num_hits, lambda_value=None, extra=None):
        self._fallback_counts[reason] += 1
        occurrence = self._fallback_counts[reason]

        payload = {
            "reason": reason,
            "score": float(score),
            "num_hits": int(num_hits),
            "x": float(score / num_hits) if num_hits else None,
            "lambda": float(lambda_value) if lambda_value is not None else None,
            "occurrence": occurrence,
            "background": self._debug_weight_summary(),
        }
        if extra is not None:
            payload["extra"] = extra
        if occurrence <= self.MAX_WARNING_EXAMPLES_PER_REASON:
            logger.warning("Weighted SaddleSum numerical fallback: %s", payload)
        elif occurrence == self.MAX_WARNING_EXAMPLES_PER_REASON + 1:
            logger.warning(
                "Weighted SaddleSum numerical fallback: further '%s' messages suppressed after %d examples",
                reason,
                self.MAX_WARNING_EXAMPLES_PER_REASON,
            )

    def log_fallback_summary(self):
        if not self._fallback_counts:
            return
        logger.warning(
            "Weighted SaddleSum numerical fallback summary (counts by reason): %s",
            dict(self._fallback_counts),
        )

    def _compute_item(self, lmbd):
        w = self._weights
        tmp = np.exp(lmbd * (w - self._wmax)) # Shift by wmax for numerical stability (largest term becomes exp(0)=1, others <= 1)
        Nrho  = tmp.sum()
        Nrho1 = (tmp * w).sum()
        Nrho2 = (tmp * w * w).sum()

        # Add contribution of all zero-weight leaves as a single term
        # Each contributes exp(lmbd * (0 - wmax)) to Nrho, and 0 to Nrho1/Nrho2
        zero_contribution = self._n_zeros * math.exp(lmbd * (-self._wmax))
        Nrho += zero_contribution
        # Nrho1 and Nrho2 unchanged — zero weights contribute 0 * anything = 0

        if Nrho <= 0.0: 
            return None

        D1K = Nrho1 / Nrho # The mean of the tilted distribution, i.e. the saddlepoint mean at lambda=lmbd. First derivative of the cumulant generating function K'(lambda).
        D2K = Nrho2 / Nrho - D1K ** 2 # The variance of the tilted distribution, i.e. the second derivative of the cumulant generating function K''(lambda). Must be positive for a valid saddlepoint approximation.
        if D2K <= 0.0:
            return None

        expH = Nrho * np.exp(lmbd * (self._wmax - D1K)) / self._N # The main exponential term in the Lugannani-Rice formula
        C = 2.0 * lmbd * math.sqrt(D2K)

        inner = lmbd * (D1K - self._wmax) - math.log(Nrho) + math.log(self._N) # The argument inside the square root for D — must be non-negative for a valid saddlepoint approximation. 
        if inner < 0.0:
            return None

        D = -math.sqrt(2.0) * math.sqrt(inner)
        return dict(lambda_=lmbd, mean=D1K, D2K=D2K, expH=expH, C=C, D=D)

    @staticmethod
    def _item_pvalue(item, m):
        sqrtm = math.sqrt(m)
        if item['D'] * sqrtm > -1.0:
            return 1.0
        if item['C'] == 0.0 or item['D'] == 0.0:
            return 1.0
        phi = math.sqrt(2.0 / math.pi) * (item['expH'] ** m)
        result = (float(ndtr(item['D'] * sqrtm))
                  + phi / item['C'] / sqrtm
                  + phi / item['D'] / sqrtm / 2.0)
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

    def _exact_singleton_right_tail(self, score): # Addition. To handle the num_hits=1 case exactly, since the saddlepoint approximation is not accurate for small m.
        """
        Exact right-tail p-value for num_hits == 1:
        P(W >= score) under the empirical background distribution.
        """
        nonzero_ge = int(self._weights.size - np.searchsorted(self._weights, score, side='left'))
        if score <= 0.0:
            count_ge = nonzero_ge + self._n_zeros
        else:
            count_ge = nonzero_ge

        p = count_ge / self._N
        min_pval = 1.0 / self._N
        return min(max(p, min_pval), 1.0)

    def pvalue(self, score, num_hits):
        """
        Right-tail p-value for observing sum=score across num_hits members.

        No cutoff_pvalue argument — we compute the full p-value for every
        term and let BH correction decide significance, exactly as the
        Fisher pipeline does.
        """
        if num_hits == 1:
            return self._exact_singleton_right_tail(score)

        x = score / num_hits
        if x <= self._mean:
            return 1.0

        min_pval = (1.0 / self._N) ** num_hits # Minimum possible p-value if all hits had the maximum weight (wmax), which is the most extreme case. 
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
                    self._log_numerical_fallback(
                        reason="lambda_upper_bound_reached",
                        score=score,
                        num_hits=num_hits,
                        lambda_value=yb,
                    )
                    return min_pval
                item = self._compute_item(yb)
                if item is None:
                    self._log_numerical_fallback(
                        reason="invalid_item_during_bracketing",
                        score=score,
                        num_hits=num_hits,
                        lambda_value=yb,
                    )
                    return min_pval
                self._insert(item)
                if abs(item['mean'] - self._wmax) < self.TOLERANCE:
                    break
            yc = 0.5 * (ya + yb)

        # Newton's method with bisection fallback
        pval = 1.0
        for _ in range(self.MAX_ITERS):
            item = self._compute_item(yc)
            if item is None:
                self._log_numerical_fallback(
                    reason="invalid_item_during_newton_bisection",
                    score=score,
                    num_hits=num_hits,
                    lambda_value=yc,
                    extra={"ya": ya, "yb": yb},
                )
                break
            pval = self._item_pvalue(item, num_hits)
            diff = item['mean'] - x # trying to drive this to zero — if it is, then we are at the saddlepoint lambda that gives the desired mean, and thus the correct p-value.
            if diff < 0.0: 
                ya = yc
            else:
                yb = yc
            y = yc - diff / item['D2K'] 
            if y < ya or y > yb:
                y = 0.5 * (ya + yb) # Bisection fallback if Newton step goes out of bounds
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


def _propagate_weights_to_leaves(weights_dict, class_to_leaf_map, removed_leaves_csv, structural_leaf_ids=None):
    """
    Build weights_with_leaves: expand any non-leaf submitted IDs so their
    leaf descendants inherit the weight. Directly submitted leaves take
    priority over inherited values.

    structural_leaf_ids: if given, leaves mislabeled with a Classification
    other than 'structural' (a ChEBI ontology error) are excluded, mirroring
    get_leaves() in fishers_calculations.py.
    Returns weights_with_leaves dict keyed by leaf IRI.
    """
    leaves_df = pd.read_csv(removed_leaves_csv)
    leaf_set = set(leaves_df['IRI'].values)

    weights_with_leaves = {}

    # First pass: propagate non-leaf weights to descendants (take maximum if conflict)
    for cls, weight in weights_dict.items():
        if cls not in leaf_set:
            descendants = class_to_leaf_map.get(cls, [])
            if structural_leaf_ids is not None:
                descendants = [leaf for leaf in descendants if leaf in structural_leaf_ids]
            for leaf in descendants:
                current = weights_with_leaves.get(leaf)
                if current is None or weight > current:
                    weights_with_leaves[leaf] = weight

    # Second pass: directly submitted leaves always overwrite inherited values
    for cls, weight in weights_dict.items():
        if cls in leaf_set:
            if structural_leaf_ids is not None and cls not in structural_leaf_ids:
                print(f"Excluding class {cls}: not classified as 'structural' in ChEBI (likely a mislabeled leaf).")
                continue
            weights_with_leaves[cls] = weight

    return weights_with_leaves


def _build_background_weights(weights_with_leaves, removed_leaves_csv, background_leaf_ids=None):
    """
    Build the full background weight array.
    Every leaf in the background gets its weight from weights_with_leaves
    (if measured) or 0.0 otherwise.

    background_leaf_ids: if given, restricts the background population to
    this leaf set — either the genuine 'structural'-classified leaves (full
    background) or a narrow background's leaves (plus any expansion), so the
    population size used inside the SaddleSum statistic matches n_bg_leaves
    reported elsewhere.

    Must be called with weights_with_leaves (post-propagation), not
    the raw weights_dict.
    """
    if background_leaf_ids is not None:
        all_bg_leaves = list(background_leaf_ids)
    else:
        leaves_df = pd.read_csv(removed_leaves_csv)
        all_bg_leaves = list(leaves_df['IRI'].values)
    return [weights_with_leaves.get(leaf, 0.0) for leaf in all_bg_leaves]


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

def _print_non_finite_pvalue_diagnostics(enrichment_results, stage):
    """
    Check for non-finite p-values (inf, -inf, nan) and print diagnostics.
    This is a diagnostic tool to track when p-value computation fails.
    """
    non_finite_classes = {}
    for cls, vals in enrichment_results.items():
        p_value = vals.get("p_value")
        if not math.isfinite(p_value):
            non_finite_classes[cls] = {
                "p_value": p_value,
                "score": vals.get("score"),
                "n_ss_annotated": vals.get("n_ss_annotated"),
                "class_name": vals.get("class"),
            }
    
    if non_finite_classes:
        logger.warning(
            "Non-finite p-values detected at stage '%s': %d classes",
            stage,
            len(non_finite_classes),
        )
        for cls, info in list(non_finite_classes.items())[:5]:  # Log first 5 examples
            logger.debug("  Class %s: %s", cls, info)
        if len(non_finite_classes) > 5:
            logger.debug("  ... and %d more", len(non_finite_classes) - 5)


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
    structural_leaf_ids=None,
):
    """
    Compute SaddleSum enrichment for every ancestor class (and role class).
    Stores ALL results — BH/Bonferroni correction filters them afterwards,
    exactly as get_enrichment_values does for Fisher.
    """
    background_weights = _build_background_weights(weights_with_leaves, removed_leaves_csv, structural_leaf_ids)
    saddler = _SaddleSum(background_weights)

    studyset_leaves_set = set(studyset_leaves)
    n_bg_leaves = count_removed_leaves(removed_leaves_csv)
    n_ss_leaves = len(studyset_leaves)

    results = {}

    if classification in ["structural", "full"]:
        for cls in studyset_ancestors:
            if cls in studyset_leaves_set:
                continue

            term_leaves = set(class_to_leaf_map.get(cls, []))

            # Calculate p-value if leaves exist; otherwise use safe defaults
            if term_leaves:
                score, n_ss_annotated, p_value = calculate_weighted_p_value(
                    saddler, term_leaves, studyset_leaves_set, weights_with_leaves
                )
            else:
                score, n_ss_annotated, p_value = 0.0, 0, 1.0

            # Background annotation count: use mapped counts when available,
            # otherwise fall back to the computed leaf set size.
            if cls in class_to_leaf_map:
                _, n_bg_annotated = count_removed_classes_for_class(
                    cls, class_to_leaf_map, classification,
                    class_to_all_roles_map, roles_to_leaves_map,
                    structural_leaf_ids,
                )
            else:
                n_bg_annotated = len(term_leaves)

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

    saddler.log_fallback_summary()

    return results


# ---------------------------------------------------------------------------
# get_weighted_enrichment_values_narrow  (mirrors get_enrichment_values_narrow
# in wikidata/narrow_background_fishers.py)
# ---------------------------------------------------------------------------

def get_weighted_enrichment_values_narrow(
    narrow_background_leaves_json,
    classification,
    studyset_leaves,
    studyset_ancestors,
    class_to_leaf_map,
    class_to_all_roles_map,
    roles_to_leaves_map,
    studyset_ancestors_roles,
    weights_with_leaves,
    leaves_to_expand_background,
    expand_background=True,
):
    """
    Compute SaddleSum enrichment for every ancestor class (and role class)
    against a narrow (organism-specific) background, instead of the full
    ChEBI background. Mirrors get_weighted_enrichment_values, but the
    background population is the narrow leaf set (plus any study-set leaves
    outside it, if expand_background is True).
    """
    with open(narrow_background_leaves_json, 'r', encoding='utf-8') as f:
        cached = json.load(f)
    background_leaf_ids = set(cached.get("narrow_leaves", []))
    if expand_background and leaves_to_expand_background:
        background_leaf_ids.update(leaves_to_expand_background)

    background_weights = _build_background_weights(weights_with_leaves, None, background_leaf_ids)
    saddler = _SaddleSum(background_weights)

    studyset_leaves_set = set(studyset_leaves)
    n_bg_leaves = count_narrow_leaves(narrow_background_leaves_json, leaves_to_expand_background, expand_background)
    n_ss_leaves = len(studyset_leaves)

    results = {}

    if classification in ["structural", "full"]:
        for cls in studyset_ancestors:
            if cls in studyset_leaves_set:
                continue

            term_leaves = set(class_to_leaf_map.get(cls, []))

            if term_leaves:
                score, n_ss_annotated, p_value = calculate_weighted_p_value(
                    saddler, term_leaves, studyset_leaves_set, weights_with_leaves
                )
            else:
                score, n_ss_annotated, p_value = 0.0, 0, 1.0

            _, n_bg_annotated = count_narrow_leaves_for_class(
                cls, narrow_background_leaves_json, class_to_leaf_map,
                leaves_to_expand_background, expand_background,
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
            _, n_bg_annotated = count_narrow_leaves_for_role(
                role, narrow_background_leaves_json, roles_to_leaves_map,
                leaves_to_expand_background, expand_background,
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

    saddler.log_fallback_summary()

    return results


def _split_graph_nodes_for_enrichment(graph_nodes, role_nodes):
    """
    Split graph nodes into structural and role nodes for enrichment.
    """
    graph_nodes_set = set(graph_nodes)
    role_nodes_set = set(role_nodes)
    structural_nodes = [node for node in graph_nodes_set if node not in role_nodes_set]
    role_nodes_in_graph = [node for node in graph_nodes_set if node in role_nodes_set]
    return structural_nodes, role_nodes_in_graph



# ---------------------------------------------------------------------------
# Shared setup logic
# ---------------------------------------------------------------------------

def _setup_weighted_analysis(weights_dict, classification,
                             removed_leaves_csv, leaf_to_ancestors_map_file,
                             class_to_leaf_map, class_to_all_roles_map,
                             roles_to_leaves_map, parent_map_file,
                             structural_leaf_ids=None):
    """Extract leaves, propagate weights, find ancestors and roles."""
    normalized_weights_dict = {
        normalize_id(cls): weight for cls, weight in weights_dict.items()
    }

    studyset_list = list(normalized_weights_dict.keys())
    studyset_leaves = get_leaves(studyset_list, removed_leaves_csv, class_to_leaf_map, structural_leaf_ids)
    print(f"Study set leaves: {len(studyset_leaves)}")

    weights_with_leaves = _propagate_weights_to_leaves(
        normalized_weights_dict, class_to_leaf_map, removed_leaves_csv, structural_leaf_ids
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


def _setup_weighted_narrow_analysis(weights_dict, classification,
                                    narrow_background_leaves_json,
                                    removed_leaves_csv, leaf_to_ancestors_map_file,
                                    class_to_leaf_map, class_to_all_roles_map,
                                    roles_to_leaves_map, parent_map_file,
                                    expand_background=True):
    """Extract leaves, propagate weights, find ancestors and roles, against a
    narrow (organism-specific) background. Mirrors _setup_weighted_analysis,
    using get_studyset_leaves_narrow (the same leaf resolution the Fisher
    narrow-background path uses) instead of get_leaves.
    """
    normalized_weights_dict = {
        normalize_id(cls): weight for cls, weight in weights_dict.items()
    }
    studyset_list = list(normalized_weights_dict.keys())

    studyset_leaves, leaves_to_expand_background, parents_to_expand_background = get_studyset_leaves_narrow(
        studyset_list, narrow_background_leaves_json, removed_leaves_csv, class_to_leaf_map
    )
    print(f"Study set leaves: {len(studyset_leaves)}")

    # No structural_leaf_ids filter here: the Fisher narrow-background path
    # (get_studyset_leaves_narrow) doesn't apply it either, so weight
    # propagation stays consistent with how narrow-background leaves are resolved.
    weights_with_leaves = _propagate_weights_to_leaves(
        normalized_weights_dict, class_to_leaf_map, removed_leaves_csv
    )
    print(f"Weights with leaf propagation: {len(weights_with_leaves)} leaves have weights")

    studyset_ancestors_all = get_ancestors_for_inputs(studyset_leaves, leaf_to_ancestors_map_file)
    print(f"Number of study set ancestors: {len(studyset_ancestors_all)}")

    if classification in ["functional", "full"]:
        # Build expanded narrow leaf set (narrow background + any study-set
        # leaves outside it), mirroring run_narrow_background_enrichment_analysis.
        with open(narrow_background_leaves_json, 'r', encoding='utf-8') as f:
            cached_narrow = json.load(f)
        expanded_narrow_leaves = set(cached_narrow.get("narrow_leaves", []))
        if expand_background and leaves_to_expand_background:
            expanded_narrow_leaves.update(leaves_to_expand_background)

        studyset_ancestors_roles = set()
        for leaf in studyset_leaves:
            if leaf in expanded_narrow_leaves:
                studyset_ancestors_roles.update(class_to_all_roles_map.get(leaf, []))
        for cls in studyset_ancestors_all:
            leaf_descs = set(class_to_leaf_map.get(cls, []))
            if leaf_descs and leaf_descs.intersection(expanded_narrow_leaves):
                studyset_ancestors_roles.update(class_to_all_roles_map.get(cls, []))

        with open(parent_map_file, 'r') as f:
            parent_map = json.load(f)
        roles_with_ancestors = set(studyset_ancestors_roles)
        to_process = list(studyset_ancestors_roles)
        while to_process:
            role = to_process.pop(0)
            for parent in parent_map.get(role, []):
                if parent in roles_with_ancestors:
                    continue
                parent_leaves = set(roles_to_leaves_map.get(parent, []))
                if parent_leaves and parent_leaves.intersection(expanded_narrow_leaves):
                    roles_with_ancestors.add(parent)
                    to_process.append(parent)
        studyset_ancestors_roles = roles_with_ancestors
        print(f"Number of roles (including ancestors): {len(studyset_ancestors_roles)}")
    else:
        studyset_ancestors_roles = set()

    return (studyset_leaves, weights_with_leaves, studyset_ancestors_all, studyset_ancestors_roles,
            leaves_to_expand_background, parents_to_expand_background)


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

    structural_leaf_ids = get_structural_leaf_ids(removed_leaves_csv)

    studyset_leaves, weights_with_leaves, studyset_ancestors_all, studyset_ancestors_roles = (
        _setup_weighted_analysis(
            weights_dict, classification,
            removed_leaves_csv, leaf_to_ancestors_map_file,
            class_to_leaf_map, class_to_all_roles_map,
            roles_to_leaves_map, parent_map_file,
            structural_leaf_ids,
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

        structural_nodes_for_enrichment, role_nodes_for_enrichment = _split_graph_nodes_for_enrichment(
            pruned_G.nodes(), studyset_ancestors_roles
        )
        studyset_ancestors = structural_nodes_for_enrichment
        studyset_ancestors_roles_for_enrichment = role_nodes_for_enrichment
        print(f"Structural nodes for enrichment (from graph): {len(studyset_ancestors)}")
        print(f"Role nodes for enrichment (from graph): {len(studyset_ancestors_roles_for_enrichment)}")
    else:
        structural_nodes_for_enrichment, role_nodes_for_enrichment = _split_graph_nodes_for_enrichment(
            pruned_G.nodes(), studyset_ancestors_roles
        )
        studyset_ancestors = structural_nodes_for_enrichment
        studyset_ancestors_roles_for_enrichment = role_nodes_for_enrichment
        print(f"Structural nodes for enrichment (from graph): {len(studyset_ancestors)}")
        print(f"Role nodes for enrichment (from graph): {len(studyset_ancestors_roles_for_enrichment)}")

    enrichment_results = get_weighted_enrichment_values(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        studyset_ancestors_roles_for_enrichment,
        weights_with_leaves,
        structural_leaf_ids,
    )
    _print_non_finite_pvalue_diagnostics(enrichment_results, "post-raw-enrichment")

    if bonferroni_correct:
        print("Applying Bonferroni correction...")
        enrichment_results, _ = bonferroni_correction(enrichment_results)
        _print_non_finite_pvalue_diagnostics(enrichment_results, "post-bonferroni")
    elif benjamini_hochberg_correct:
        print("Applying Benjamini-Hochberg FDR correction...")
        enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
        _print_non_finite_pvalue_diagnostics(enrichment_results, "post-bh")

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
    _print_non_finite_pvalue_diagnostics(enrichment_results, "final-post-pruning")

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

    structural_leaf_ids = get_structural_leaf_ids(removed_leaves_csv)

    studyset_leaves, weights_with_leaves, studyset_ancestors, studyset_ancestors_roles = (
        _setup_weighted_analysis(
            weights_dict, classification,
            removed_leaves_csv, leaf_to_ancestors_map_file,
            class_to_leaf_map, class_to_all_roles_map,
            roles_to_leaves_map, parent_map_file,
            structural_leaf_ids,
        )
    )

    all_removed_nodes = set()

    G = create_graph_with_roles_and_structures(
        studyset_leaves, studyset_ancestors,
        studyset_ancestors_roles, parent_map_file,
        class_to_all_roles_map, classification,
    )

    structural_nodes_for_enrichment, role_nodes_for_enrichment = _split_graph_nodes_for_enrichment(
        G.nodes(), studyset_ancestors_roles
    )
    print(f"Structural nodes for enrichment (from graph): {len(structural_nodes_for_enrichment)}")
    print(f"Role nodes for enrichment (from graph): {len(role_nodes_for_enrichment)}")

    enrichment_results = get_weighted_enrichment_values(
        removed_leaves_csv,
        classification,
        studyset_leaves,
        structural_nodes_for_enrichment,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        role_nodes_for_enrichment,
        weights_with_leaves,
        structural_leaf_ids,
    )

    print("Enrichment results (raw):")
    print_enrichment_results(enrichment_results)

    enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
    print("After BH correction:")
    print_enrichment_results(enrichment_results)

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


# ---------------------------------------------------------------------------
# run_weighted_narrow_background_enrichment_analysis
# (mirrors run_narrow_background_enrichment_analysis in
#  wikidata/narrow_background_fishers.py)
# ---------------------------------------------------------------------------

def run_weighted_narrow_background_enrichment_analysis(
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
    narrow_background_leaves_json="data/human_entities_leaves.json",
    expand_background=True,
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

    (studyset_leaves, weights_with_leaves, studyset_ancestors_all, studyset_ancestors_roles,
     leaves_to_expand_background, parents_to_expand_background) = _setup_weighted_narrow_analysis(
        weights_dict, classification,
        narrow_background_leaves_json,
        removed_leaves_csv, leaf_to_ancestors_map_file,
        class_to_leaf_map, class_to_all_roles_map,
        roles_to_leaves_map, parent_map_file,
        expand_background,
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

    structural_nodes_for_enrichment, role_nodes_for_enrichment = _split_graph_nodes_for_enrichment(
        pruned_G.nodes(), studyset_ancestors_roles
    )
    studyset_ancestors = structural_nodes_for_enrichment
    studyset_ancestors_roles_for_enrichment = role_nodes_for_enrichment
    print(f"Structural nodes for enrichment (from graph): {len(studyset_ancestors)}")
    print(f"Role nodes for enrichment (from graph): {len(studyset_ancestors_roles_for_enrichment)}")

    enrichment_results = get_weighted_enrichment_values_narrow(
        narrow_background_leaves_json,
        classification,
        studyset_leaves,
        studyset_ancestors,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        studyset_ancestors_roles_for_enrichment,
        weights_with_leaves,
        leaves_to_expand_background,
        expand_background,
    )
    _print_non_finite_pvalue_diagnostics(enrichment_results, "post-raw-enrichment")

    if bonferroni_correct:
        print("Applying Bonferroni correction...")
        enrichment_results, _ = bonferroni_correction(enrichment_results)
        _print_non_finite_pvalue_diagnostics(enrichment_results, "post-bonferroni")
    elif benjamini_hochberg_correct:
        print("Applying Benjamini-Hochberg FDR correction...")
        enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
        _print_non_finite_pvalue_diagnostics(enrichment_results, "post-bh")

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
    _print_non_finite_pvalue_diagnostics(enrichment_results, "final-post-pruning")

    results = {
        "study_set": [id_to_name(c) for c in studyset_leaves],
        "removed_nodes": [id_to_name(c) for c in all_removed_nodes],
        "enrichment_results": {
            id_to_name(cls): vals for cls, vals in enrichment_results.items()
        },
    }

    return results, pruned_G, leaves_to_expand_background, parents_to_expand_background


# ---------------------------------------------------------------------------
# run_weighted_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy
# ---------------------------------------------------------------------------

def run_weighted_narrow_background_enrichment_analysis_plain_enrich_pruning_strategy(
    weights_dict,
    levels=2,
    n=0,
    p_value_threshold=0.05,
    classification="structural",
    narrow_background_leaves_json="data/human_entities_leaves.json",
    expand_background=True,
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

    (studyset_leaves, weights_with_leaves, studyset_ancestors, studyset_ancestors_roles,
     leaves_to_expand_background, parents_to_expand_background) = _setup_weighted_narrow_analysis(
        weights_dict, classification,
        narrow_background_leaves_json,
        removed_leaves_csv, leaf_to_ancestors_map_file,
        class_to_leaf_map, class_to_all_roles_map,
        roles_to_leaves_map, parent_map_file,
        expand_background,
    )

    all_removed_nodes = set()

    G = create_graph_with_roles_and_structures(
        studyset_leaves, studyset_ancestors,
        studyset_ancestors_roles, parent_map_file,
        class_to_all_roles_map, classification,
    )

    structural_nodes_for_enrichment, role_nodes_for_enrichment = _split_graph_nodes_for_enrichment(
        G.nodes(), studyset_ancestors_roles
    )
    print(f"Structural nodes for enrichment (from graph): {len(structural_nodes_for_enrichment)}")
    print(f"Role nodes for enrichment (from graph): {len(role_nodes_for_enrichment)}")

    enrichment_results = get_weighted_enrichment_values_narrow(
        narrow_background_leaves_json,
        classification,
        studyset_leaves,
        structural_nodes_for_enrichment,
        class_to_leaf_map,
        class_to_all_roles_map,
        roles_to_leaves_map,
        role_nodes_for_enrichment,
        weights_with_leaves,
        leaves_to_expand_background,
        expand_background,
    )

    print("Enrichment results (raw):")
    print_enrichment_results(enrichment_results)

    enrichment_results = benjamini_hochberg_fdr_correction(enrichment_results)
    print("After BH correction:")
    print_enrichment_results(enrichment_results)

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

    return results, G, leaves_to_expand_background, parents_to_expand_background