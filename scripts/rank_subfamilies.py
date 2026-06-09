#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rank_subfamilies.py

Ranks subfamilies of genus 2 hyperelliptic curves using Magma rank results.

Primary input:
    magma_rank_results_all_appearances.csv

Secondary inputs (for additional family metadata: q_family, H_family, rs, heights):
    families_construct_even_automorphism_two_stage_integral_models/all_stage2_integral_models_summary.csv
    families_construct_general_two_stage_integral_models/all_stage2_integral_models_summary.csv

Grouping key: (construction_kind, execution_relative_dir, family_index)

Within each subfamily, curves are already ordered by curve_index (scanner sorted best first).
Only the first TOP_K curves per subfamily are considered.

Outputs:
    ranked_subfamilies.csv
    ranked_subfamilies_detailed_curves.csv
    selected_subfamilies_report.txt
    ranked_subfamilies_latex.tex  (optional)
"""

import csv
import json
import math
import statistics
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ============================================================
# CONFIGURATION
# ============================================================

TOP_K = 5
RANK_THRESHOLDS = [6, 8, 10]   # R1, R2, R3

ALL_APPEARANCES_CSV = Path(
    "magma_rank_results_integral_models/magma_rank_results_all_appearances.csv"
)
# If the above does not exist yet, point to the latest saved iteration, e.g.:
# ALL_APPEARANCES_CSV = Path("../Dizertatie_SavedResults/magma_rank_results_integral_models_iteration3/magma_rank_results_all_appearances.csv")

INTEGRAL_MODEL_CSVS = [
    Path("families_construct_even_automorphism_two_stage_integral_models/all_stage2_integral_models_summary.csv"),
    Path("families_construct_general_two_stage_integral_models/all_stage2_integral_models_summary.csv"),
]

OUTPUT_DIR = Path("subfamily_ranking_results")

WRITE_LATEX = True
WRITE_JSON  = True

# Lexicographic ranking: tuple components, each maximised unless noted.
# Negative sign = minimise that component.
# The tuple is built in rank_key() below.

# ============================================================
# HELPERS
# ============================================================

def safe_int(v, default=None):
    try:
        return int(v)
    except (TypeError, ValueError):
        return default

def safe_float(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default

def safe_median(lst):
    return statistics.median(lst) if lst else None

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def write_csv(path: Path, rows: List[Dict], fieldnames: List[str]):
    ensure_dir(path.parent)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})

# ============================================================
# LOAD DATA
# ============================================================

def load_all_appearances(path: Path) -> List[Dict]:
    if not path.exists():
        raise FileNotFoundError(f"Primary input not found: {path}")
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_integral_model_metadata(paths: List[Path]) -> Dict[Tuple, Dict]:
    """
    Returns a dict keyed by (construction_kind, execution_relative_dir, family_index, curve_index)
    with extra metadata: q_family, H_family, Pq_family, rs, forced_xs,
    height_F_original, height_F_integral.
    """
    meta = {}
    for path in paths:
        if not path.exists():
            print(f"WARNING: integral model summary not found: {path}")
            continue
        with open(path, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                key = (
                    row.get("construction_kind", ""),
                    row.get("execution_relative_dir", ""),
                    row.get("family_index", ""),
                    row.get("curve_index", ""),
                )
                meta[key] = {
                    "q_family":         row.get("q_family", ""),
                    "H_family":         row.get("H_family", ""),
                    "Pq_family":        row.get("Pq_family", ""),
                    "rs":               row.get("rs", ""),
                    "forced_xs":        row.get("forced_xs", ""),
                    "height_F_original": row.get("height_F_original", ""),
                    "height_F_integral": row.get("height_F_integral", ""),
                    "m":                row.get("m", ""),
                    "q_params_family":  row.get("q_params_family", ""),
                    "top_k_min_affine_family": row.get("top_k_min_affine", ""),
                    "top_k_sum_affine_family": row.get("top_k_sum_affine", ""),
                    "best_affine_point_count_family": row.get("best_affine_point_count", ""),
                }
    return meta

# ============================================================
# DETERMINE RANK STATUS CATEGORIES
# ============================================================

DETERMINED_STATUSES = {
    "determined_by_RankBounds",
    "determined_by_MordellWeilGroupGenus2",
}

TIMEOUT_STATUSES = {
    "magma_time_limit",
    "transient_failure_after_retries",
}

MEMORY_STATUSES = {
    "magma_memory_limit",
}

REQUEST_ERROR_STATUSES = {
    "request_error",
}

UNRESOLVED_STATUSES = {
    "not_closed_after_MordellWeilGroupGenus2",
    "only_RankBounds_available",
    "magma_output_unparsed",
    "magma_runtime_error_unparsed",
    "magma_user_error_unparsed",
    "magma_calculator_offline",
    "empty_results",
    "xml_parse_error",
    "unknown",
}

def is_determined(status: str) -> bool:
    return status in DETERMINED_STATUSES

def is_successful_request(request_status: str) -> bool:
    return request_status in ("ok", "cached")

# ============================================================
# GROUP INTO SUBFAMILIES
# ============================================================

def group_by_subfamily(rows: List[Dict]) -> Dict[Tuple, List[Dict]]:
    """Groups rows by (construction_kind, execution_relative_dir, family_index)."""
    groups: Dict[Tuple, List[Dict]] = {}
    for row in rows:
        key = (
            row.get("construction_kind", ""),
            row.get("execution_relative_dir", ""),
            row.get("family_index", ""),
        )
        groups.setdefault(key, []).append(row)
    return groups


def sort_subfamily_curves(curves: List[Dict]) -> List[Dict]:
    """Sort curves within a subfamily by curve_index (ascending = best first per scanner)."""
    def curve_sort_key(r):
        return safe_int(r.get("curve_index"), 9999)
    return sorted(curves, key=curve_sort_key)

# ============================================================
# COMPUTE SUBFAMILY METRICS
# ============================================================

def compute_subfamily_metrics(
    key: Tuple,
    curves: List[Dict],
    meta: Dict[Tuple, Dict],
    top_k: int,
    r1: int,
    r2: int,
    r3: int,
) -> Dict:
    construction_kind, execution_relative_dir, family_index = key

    sorted_curves = sort_subfamily_curves(curves)
    top_curves = sorted_curves[:top_k]

    # --- family-level metadata from integral model summary (use first curve) ---
    first_meta = meta.get((
        construction_kind,
        execution_relative_dir,
        family_index,
        top_curves[0].get("curve_index", "") if top_curves else "",
    ), {})

    family_label   = top_curves[0].get("family_label", "") if top_curves else ""
    q_family       = first_meta.get("q_family", "")
    H_family       = first_meta.get("H_family", "")
    Pq_family      = first_meta.get("Pq_family", "")
    rs             = first_meta.get("rs", "")
    forced_xs      = first_meta.get("forced_xs", "")
    m              = first_meta.get("m", "")
    q_params       = first_meta.get("q_params_family", "")

    # --- basic counts ---
    top_k_count = len(top_curves)
    n_successful = sum(1 for c in top_curves if is_successful_request(c.get("request_status", "")))
    n_determined = sum(1 for c in top_curves if is_determined(c.get("rank_status", "")))
    n_only_rankbounds = sum(1 for c in top_curves if c.get("rank_status") == "only_RankBounds_available")
    n_not_closed = sum(1 for c in top_curves if c.get("rank_status") == "not_closed_after_MordellWeilGroupGenus2")
    n_timeout    = sum(1 for c in top_curves if c.get("rank_status", "") in TIMEOUT_STATUSES)
    n_memory     = sum(1 for c in top_curves if c.get("rank_status", "") in MEMORY_STATUSES)
    n_req_error  = sum(1 for c in top_curves if c.get("rank_status", "") in REQUEST_ERROR_STATUSES)
    n_unresolved = sum(1 for c in top_curves if c.get("rank_status", "") in UNRESOLVED_STATUSES)
    n_failed     = n_timeout + n_memory + n_req_error + n_unresolved

    # --- rank values ---
    determined_ranks = []
    for c in top_curves:
        if is_determined(c.get("rank_status", "")):
            r = safe_int(c.get("rank_Jacobian"))
            if r is not None:
                determined_ranks.append(r)

    n_rank_ge_r1 = sum(1 for r in determined_ranks if r >= r1)
    n_rank_ge_r2 = sum(1 for r in determined_ranks if r >= r2)
    n_rank_ge_r3 = sum(1 for r in determined_ranks if r >= r3)

    min_rank  = min(determined_ranks) if determined_ranks else None
    max_rank  = max(determined_ranks) if determined_ranks else None
    sum_rank  = sum(determined_ranks) if determined_ranks else None
    avg_rank  = sum(determined_ranks) / len(determined_ranks) if determined_ranks else None
    med_rank  = safe_median(determined_ranks)

    # rank source breakdown
    n_by_rankbounds   = sum(1 for c in top_curves if c.get("rank_status") == "determined_by_RankBounds")
    n_by_mordellweil  = sum(1 for c in top_curves if c.get("rank_status") == "determined_by_MordellWeilGroupGenus2")
    n_finite_index    = sum(1 for c in top_curves if c.get("finiteIndex", "").lower() == "true")
    n_proved          = sum(1 for c in top_curves if c.get("proved", "").lower() == "true")
    n_grh             = sum(1 for c in top_curves if c.get("assumption", "").upper() == "GRH")

    # magma timing
    magma_times = [safe_float(c.get("magma_time")) for c in top_curves]
    magma_times = [t for t in magma_times if t is not None]
    avg_magma_time = sum(magma_times) / len(magma_times) if magma_times else None
    max_magma_time = max(magma_times) if magma_times else None
    total_magma_time = sum(magma_times) if magma_times else None

    # --- point-search metrics ---
    affine_vals   = [safe_int(c.get("affine_point_count")) for c in top_curves]
    affine_vals   = [v for v in affine_vals if v is not None]
    proj_vals     = [safe_int(c.get("projective_lower_bound")) for c in top_curves]
    proj_vals     = [v for v in proj_vals if v is not None]
    extra_x_vals  = [safe_int(c.get("extra_x_count")) for c in top_curves]
    extra_x_vals  = [v for v in extra_x_vals if v is not None]

    top_k_min_affine  = min(affine_vals)  if affine_vals  else None
    top_k_max_affine  = max(affine_vals)  if affine_vals  else None
    top_k_sum_affine  = sum(affine_vals)  if affine_vals  else None
    top_k_avg_affine  = sum(affine_vals) / len(affine_vals) if affine_vals else None
    top_k_min_proj    = min(proj_vals)    if proj_vals    else None
    top_k_sum_proj    = sum(proj_vals)    if proj_vals    else None
    top_k_min_extra_x = min(extra_x_vals) if extra_x_vals else None
    top_k_sum_extra_x = sum(extra_x_vals) if extra_x_vals else None

    best_affine    = max(affine_vals) if affine_vals else None
    best_proj      = max(proj_vals)   if proj_vals   else None

    # heights from meta (per curve joined)
    heights = []
    for c in top_curves:
        cmeta = meta.get((construction_kind, execution_relative_dir, family_index, c.get("curve_index", "")), {})
        h = safe_float(cmeta.get("height_F_integral")) or safe_float(cmeta.get("height_F_original"))
        if h is not None:
            heights.append(h)
    best_height   = min(heights) if heights else None
    avg_height    = sum(heights) / len(heights) if heights else None
    max_height    = max(heights) if heights else None

    return {
        # identity
        "construction_kind":       construction_kind,
        "execution_relative_dir":  execution_relative_dir,
        "family_index":            family_index,
        "family_label":            family_label,
        "m":                       m,
        "rs":                      rs,
        "forced_xs":               forced_xs,
        "q_params":                q_params,
        "q_family":                q_family,
        "H_family":                H_family,
        "Pq_family":               Pq_family,
        # counts
        "top_k_count":             top_k_count,
        "n_successful":            n_successful,
        "n_determined":            n_determined,
        "n_only_rankbounds":       n_only_rankbounds,
        "n_not_closed":            n_not_closed,
        "n_timeout":               n_timeout,
        "n_memory":                n_memory,
        "n_req_error":             n_req_error,
        "n_unresolved":            n_unresolved,
        "n_failed":                n_failed,
        # rank thresholds
        f"n_rank_ge_{r1}":         n_rank_ge_r1,
        f"n_rank_ge_{r2}":         n_rank_ge_r2,
        f"n_rank_ge_{r3}":         n_rank_ge_r3,
        # rank stats
        "min_rank":                min_rank,
        "max_rank":                max_rank,
        "sum_rank":                sum_rank,
        "avg_rank":                avg_rank,
        "med_rank":                med_rank,
        "determined_ranks":        str(sorted(determined_ranks)),
        # rank source
        "n_by_rankbounds":         n_by_rankbounds,
        "n_by_mordellweil":        n_by_mordellweil,
        "n_finite_index":          n_finite_index,
        "n_proved":                n_proved,
        "n_grh":                   n_grh,
        # timing
        "avg_magma_time":          avg_magma_time,
        "max_magma_time":          max_magma_time,
        "total_magma_time":        total_magma_time,
        # point metrics
        "top_k_min_affine":        top_k_min_affine,
        "top_k_max_affine":        top_k_max_affine,
        "top_k_sum_affine":        top_k_sum_affine,
        "top_k_avg_affine":        top_k_avg_affine,
        "top_k_min_proj":          top_k_min_proj,
        "top_k_sum_proj":          top_k_sum_proj,
        "top_k_min_extra_x":       top_k_min_extra_x,
        "top_k_sum_extra_x":       top_k_sum_extra_x,
        "best_affine":             best_affine,
        "best_proj":               best_proj,
        # height metrics
        "best_height":             best_height,
        "avg_height":              avg_height,
        "max_height":              max_height,
    }

# ============================================================
# RANKING KEY
# ============================================================

def rank_key(m: Dict, r1: int, r2: int, r3: int) -> Tuple:
    """
    Lexicographic ranking tuple. Python sorts ascending, so negate to maximise.

    The philosophy: a subfamily is strong if ALL its top-k curves have high rank,
    not just one exceptional outlier. This reflects the heuristic that a high floor
    (min_rank) across all top-k curves is evidence that the forced rational points
    generate independent divisors on the Jacobian, systematically contributing rank
    regardless of which Q-coefficients are chosen.

    Components (in priority order):
      1. maximize n_determined          — ranks must be known
      2. maximize min_rank              — high floor = systematic forced-point contribution
      3. maximize n_rank_ge_R3          — how many hit the exceptional threshold
      4. maximize n_rank_ge_R2
      5. maximize n_rank_ge_R1
      6. maximize sum_rank              — total rank mass across top-k
      7. maximize avg_rank
      8. maximize n_successful          — prefer fewer Magma failures
      9. minimize n_failed
      10. maximize top_k_min_affine     — point-search quality tiebreakers
      11. maximize top_k_sum_affine
      12. maximize top_k_min_extra_x
      13. minimize avg_height
    """
    def neg(v, default=0):
        return -(v if v is not None else default)

    return (
        neg(m.get("n_determined", 0)),
        neg(m.get("min_rank"), 0),             # floor first: high floor = systematic contribution
        neg(m.get(f"n_rank_ge_{r3}", 0)),
        neg(m.get(f"n_rank_ge_{r2}", 0)),
        neg(m.get(f"n_rank_ge_{r1}", 0)),
        neg(m.get("sum_rank"), 0),
        neg(m.get("avg_rank"), 0.0),
        neg(m.get("n_successful", 0)),
        m.get("n_failed", 0),                  # minimise → no negation
        neg(m.get("top_k_min_affine"), 0),
        neg(m.get("top_k_sum_affine"), 0),
        neg(m.get("top_k_min_extra_x"), 0),
        m.get("avg_height") if m.get("avg_height") is not None else float("inf"),  # minimise
    )

# ============================================================
# DETAILED CURVES OUTPUT
# ============================================================

def build_detailed_curve_rows(
    ranked_subfamilies: List[Dict],
    groups: Dict[Tuple, List[Dict]],
    meta: Dict[Tuple, Dict],
    top_k: int,
) -> List[Dict]:
    detailed = []
    for sf in ranked_subfamilies:
        key = (sf["construction_kind"], sf["execution_relative_dir"], sf["family_index"])
        curves = sort_subfamily_curves(groups[key])[:top_k]
        for c in curves:
            cmeta = meta.get((
                sf["construction_kind"],
                sf["execution_relative_dir"],
                sf["family_index"],
                c.get("curve_index", ""),
            ), {})
            row = {
                "subfamily_rank":          sf["subfamily_rank"],
                "construction_kind":       sf["construction_kind"],
                "execution_relative_dir":  sf["execution_relative_dir"],
                "family_index":            sf["family_index"],
                "family_label":            sf["family_label"],
                "curve_index":             c.get("curve_index", ""),
                "curve_label":             c.get("curve_label", ""),
                "unique_id":               c.get("unique_id", ""),
                "F_integral":              c.get("F_integral", ""),
                "affine_point_count":      c.get("affine_point_count", ""),
                "projective_lower_bound":  c.get("projective_lower_bound", ""),
                "extra_x_count":           c.get("extra_x_count", ""),
                "height_F_integral":       cmeta.get("height_F_integral", ""),
                "rank_Jacobian":           c.get("rank_Jacobian", ""),
                "rank_status":             c.get("rank_status", ""),
                "rank_source":             c.get("rank_source", ""),
                "rankbounds_lb":           c.get("rankbounds_lb", ""),
                "rankbounds_ub":           c.get("rankbounds_ub", ""),
                "finiteIndex":             c.get("finiteIndex", ""),
                "proved":                  c.get("proved", ""),
                "assumption":              c.get("assumption", ""),
                "magma_time":              c.get("magma_time", ""),
                "request_status":          c.get("request_status", ""),
                "q_family":                cmeta.get("q_family", ""),
                "H_family":                cmeta.get("H_family", ""),
            }
            detailed.append(row)
    return detailed

# ============================================================
# REPORT
# ============================================================

def write_report(
    path: Path,
    ranked: List[Dict],
    top_n: int,
    r1: int, r2: int, r3: int,
    top_k: int,
):
    ensure_dir(path.parent)
    lines = []
    lines.append("=" * 100)
    lines.append("SUBFAMILY RANKING REPORT")
    lines.append(f"TOP_K = {top_k}  |  RANK_THRESHOLDS = R1={r1}, R2={r2}, R3={r3}")
    lines.append(f"Showing top {top_n} subfamilies out of {len(ranked)} total")
    lines.append("=" * 100)

    for sf in ranked[:top_n]:
        lines.append("")
        lines.append(f"RANK #{sf['subfamily_rank']}  |  {sf['construction_kind']}  |  {sf['family_label']}  |  family_index={sf['family_index']}")
        lines.append(f"  execution_dir : {sf['execution_relative_dir']}")
        if sf.get("rs"):
            lines.append(f"  rs            : {sf['rs']}")
        if sf.get("forced_xs"):
            lines.append(f"  forced_xs     : {sf['forced_xs']}")
        if sf.get("q_family"):
            lines.append(f"  q(x)          : {sf['q_family']}")
        if sf.get("H_family"):
            lines.append(f"  H(x)          : {sf['H_family']}")
        lines.append(f"  top_k curves  : {sf['top_k_count']}")
        lines.append(f"  determined    : {sf['n_determined']} / {sf['top_k_count']}")
        lines.append(f"  ranks         : {sf['determined_ranks']}")
        lines.append(f"  min/max/sum   : {sf['min_rank']} / {sf['max_rank']} / {sf['sum_rank']}")
        lines.append(f"  avg / median  : {sf['avg_rank']} / {sf['med_rank']}")
        lines.append(f"  >= R1={r1}       : {sf.get(f'n_rank_ge_{r1}', 0)}")
        lines.append(f"  >= R2={r2}       : {sf.get(f'n_rank_ge_{r2}', 0)}")
        lines.append(f"  >= R3={r3}       : {sf.get(f'n_rank_ge_{r3}', 0)}")
        lines.append(f"  n_failed      : {sf['n_failed']} (timeout={sf['n_timeout']}, mem={sf['n_memory']}, req_err={sf['n_req_error']})")
        lines.append(f"  affine counts : min={sf['top_k_min_affine']}  sum={sf['top_k_sum_affine']}  best={sf['best_affine']}")
        lines.append(f"  extra_x       : min={sf['top_k_min_extra_x']}  sum={sf['top_k_sum_extra_x']}")
        if sf.get("avg_height"):
            lines.append(f"  height        : best={sf['best_height']}  avg={sf['avg_height']:.1f}")
        lines.append("-" * 100)

    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Report written to {path}")

# ============================================================
# LATEX OUTPUT
# ============================================================

def write_latex(path: Path, ranked: List[Dict], top_n: int, r1: int, r2: int, r3: int):
    ensure_dir(path.parent)
    rows_tex = []
    for sf in ranked[:top_n]:
        kind_short = "even" if sf["construction_kind"] == "even_automorphism" else "gen"
        rows_tex.append(
            f"  {sf['subfamily_rank']} & {kind_short} & "
            f"{sf['family_label'].replace('#', r'\#')} & "
            f"{sf['n_determined']} & "
            f"{sf.get(f'n_rank_ge_{r3}', 0)} & "
            f"{sf['min_rank'] if sf['min_rank'] is not None else '--'} & "
            f"{sf['max_rank'] if sf['max_rank'] is not None else '--'} & "
            f"{sf['avg_rank']:.2f} if {sf['avg_rank'] is not None} else '--' & "
            f"{sf['top_k_min_affine']} \\\\"
        )

    content = (
        r"\begin{table}[ht]" + "\n"
        r"\centering" + "\n"
        r"\small" + "\n"
        r"\begin{tabular}{rlllrrrrl}" + "\n"
        r"\toprule" + "\n"
        rf"Rank & Kind & Family & Det. & $\geq${r3} & Min & Max & Avg & MinAff \\" + "\n"
        r"\midrule" + "\n"
        + "\n".join(rows_tex) + "\n"
        + r"\bottomrule" + "\n"
        r"\end{tabular}" + "\n"
        r"\caption{Top subfamilies ranked by Mordell--Weil rank performance.}" + "\n"
        r"\label{tab:subfamily_ranking}" + "\n"
        r"\end{table}" + "\n"
    )
    path.write_text(content, encoding="utf-8")
    print(f"LaTeX table written to {path}")

# ============================================================
# MAIN
# ============================================================

def main():
    ensure_dir(OUTPUT_DIR)

    r1, r2, r3 = RANK_THRESHOLDS

    print("Loading all_appearances CSV...")
    rows = load_all_appearances(ALL_APPEARANCES_CSV)
    print(f"  {len(rows)} rows loaded")

    print("Loading integral model metadata...")
    meta = load_integral_model_metadata(INTEGRAL_MODEL_CSVS)
    print(f"  {len(meta)} curve-level metadata entries loaded")

    print("Grouping into subfamilies...")
    groups = group_by_subfamily(rows)
    print(f"  {len(groups)} subfamilies found")

    print(f"Computing metrics (TOP_K={TOP_K})...")
    subfamily_metrics = []
    for key, curves in groups.items():
        m = compute_subfamily_metrics(key, curves, meta, TOP_K, r1, r2, r3)
        subfamily_metrics.append(m)

    print("Ranking subfamilies...")
    subfamily_metrics.sort(key=lambda m: rank_key(m, r1, r2, r3))
    for i, m in enumerate(subfamily_metrics, start=1):
        m["subfamily_rank"] = i

    # --- ranked_subfamilies.csv ---
    sf_fields = [
        "subfamily_rank", "construction_kind", "family_label", "family_index",
        "execution_relative_dir", "m", "rs", "forced_xs", "q_params",
        "q_family", "H_family", "Pq_family",
        "top_k_count", "n_successful", "n_determined",
        "n_only_rankbounds", "n_not_closed",
        "n_timeout", "n_memory", "n_req_error", "n_unresolved", "n_failed",
        f"n_rank_ge_{r1}", f"n_rank_ge_{r2}", f"n_rank_ge_{r3}",
        "min_rank", "max_rank", "sum_rank", "avg_rank", "med_rank",
        "determined_ranks",
        "n_by_rankbounds", "n_by_mordellweil", "n_finite_index", "n_proved", "n_grh",
        "avg_magma_time", "max_magma_time", "total_magma_time",
        "top_k_min_affine", "top_k_max_affine", "top_k_sum_affine", "top_k_avg_affine",
        "top_k_min_proj", "top_k_sum_proj",
        "top_k_min_extra_x", "top_k_sum_extra_x",
        "best_affine", "best_proj",
        "best_height", "avg_height", "max_height",
    ]
    write_csv(OUTPUT_DIR / "ranked_subfamilies.csv", subfamily_metrics, sf_fields)
    print(f"Written: {OUTPUT_DIR / 'ranked_subfamilies.csv'}")

    # --- ranked_subfamilies_detailed_curves.csv ---
    detailed = build_detailed_curve_rows(subfamily_metrics, groups, meta, TOP_K)
    detail_fields = [
        "subfamily_rank", "construction_kind", "family_label", "family_index",
        "execution_relative_dir", "curve_index", "curve_label", "unique_id",
        "F_integral", "affine_point_count", "projective_lower_bound",
        "extra_x_count", "height_F_integral",
        "rank_Jacobian", "rank_status", "rank_source",
        "rankbounds_lb", "rankbounds_ub", "finiteIndex", "proved", "assumption",
        "magma_time", "request_status",
        "q_family", "H_family",
    ]
    write_csv(OUTPUT_DIR / "ranked_subfamilies_detailed_curves.csv", detailed, detail_fields)
    print(f"Written: {OUTPUT_DIR / 'ranked_subfamilies_detailed_curves.csv'}")

    # --- report ---
    write_report(
        OUTPUT_DIR / "selected_subfamilies_report.txt",
        subfamily_metrics,
        top_n=50,
        r1=r1, r2=r2, r3=r3,
        top_k=TOP_K,
    )

    # --- latex ---
    if WRITE_LATEX:
        write_latex(OUTPUT_DIR / "ranked_subfamilies_latex.tex", subfamily_metrics, top_n=30, r1=r1, r2=r2, r3=r3)

    # --- json ---
    if WRITE_JSON:
        json_path = OUTPUT_DIR / "ranked_subfamilies.json"
        # convert None to null-safe
        def default_serial(obj):
            if isinstance(obj, float) and math.isnan(obj):
                return None
            raise TypeError
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(subfamily_metrics, f, indent=2, default=str)
        print(f"Written: {json_path}")

    # --- summary ---
    total = len(subfamily_metrics)
    n_any_determined = sum(1 for m in subfamily_metrics if m["n_determined"] > 0)
    n_full = sum(1 for m in subfamily_metrics if m["n_determined"] == TOP_K)
    n_ge_r3_any = sum(1 for m in subfamily_metrics if m.get(f"n_rank_ge_{r3}", 0) > 0)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print(f"  Total subfamilies ranked    : {total}")
    print(f"  With at least 1 determined  : {n_any_determined}")
    print(f"  With all {TOP_K} determined       : {n_full}")
    print(f"  With at least 1 rank >= {r3}  : {n_ge_r3_any}")
    top1 = subfamily_metrics[0]
    print(f"  Best subfamily              : {top1['family_label']} ({top1['construction_kind']})")
    print(f"    ranks = {top1['determined_ranks']}")
    print("=" * 70)


if __name__ == "__main__":
    main()
