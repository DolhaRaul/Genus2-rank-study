#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rank_subfamilies_by_kind.py

Same ranking logic as rank_subfamilies.py, but restricted to a single
construction kind. Useful for per-type analysis and for the D4/reciprocal
family once its Magma results are available.

Change FILTER_CONSTRUCTION_KIND to one of:
    "even_automorphism"
    "general"
    "d4"          (or whatever label the D4 scanner produces)
"""

# ============================================================
# CONFIGURATION — change these two lines as needed
# ============================================================

FILTER_CONSTRUCTION_KIND = "even_automorphism"   # or "general", "d4", etc.

TOP_K = 5
RANK_THRESHOLDS = [6, 8, 10]   # R1, R2, R3

ALL_APPEARANCES_CSV_CANDIDATES = [
    # current run (if available)
    "magma_rank_results_integral_models/magma_rank_results_all_appearances.csv",
    # fallback to latest saved iteration
    "../Dizertatie_SavedResults/magma_rank_results_integral_models_iteration3/magma_rank_results_all_appearances.csv",
]

INTEGRAL_MODEL_CSVS = [
    "families_construct_even_automorphism_two_stage_integral_models/all_stage2_integral_models_summary.csv",
    "families_construct_general_two_stage_integral_models/all_stage2_integral_models_summary.csv",
    # add D4 summary CSV here when available:
    # "families_construct_d4_two_stage_integral_models/all_stage2_integral_models_summary.csv",
]

WRITE_LATEX = True
WRITE_JSON  = True

# ============================================================
# IMPORT SHARED LOGIC FROM rank_subfamilies.py
# ============================================================

import csv, json, math, statistics
from pathlib import Path
from typing import Dict, List, Tuple

# resolve paths
from rank_subfamilies import (
    load_integral_model_metadata,
    group_by_subfamily,
    compute_subfamily_metrics,
    rank_key,
    build_detailed_curve_rows,
    write_report,
    write_latex,
    write_csv,
    ensure_dir,
)

# ============================================================
# MAIN
# ============================================================

def main():
    kind = FILTER_CONSTRUCTION_KIND
    output_dir = Path(f"subfamily_ranking_results_{kind}")
    ensure_dir(output_dir)

    r1, r2, r3 = RANK_THRESHOLDS

    # resolve input CSV
    appearances_path = None
    for candidate in ALL_APPEARANCES_CSV_CANDIDATES:
        p = Path(candidate)
        if p.exists():
            appearances_path = p
            break
    if appearances_path is None:
        raise FileNotFoundError("No all_appearances CSV found. Check ALL_APPEARANCES_CSV_CANDIDATES.")

    print(f"Using: {appearances_path}")
    with open(appearances_path, newline="", encoding="utf-8") as f:
        all_rows = list(csv.DictReader(f))

    # filter to requested kind
    rows = [r for r in all_rows if r.get("construction_kind", "") == kind]
    print(f"Rows for '{kind}': {len(rows)} / {len(all_rows)} total")
    if not rows:
        raise ValueError(f"No rows found for construction_kind='{kind}'. Check the value.")

    meta = load_integral_model_metadata([Path(p) for p in INTEGRAL_MODEL_CSVS])

    groups = group_by_subfamily(rows)
    print(f"Subfamilies of kind '{kind}': {len(groups)}")

    subfamily_metrics = []
    for key, curves in groups.items():
        m = compute_subfamily_metrics(key, curves, meta, TOP_K, r1, r2, r3)
        subfamily_metrics.append(m)

    subfamily_metrics.sort(key=lambda m: rank_key(m, r1, r2, r3))
    for i, m in enumerate(subfamily_metrics, start=1):
        m["subfamily_rank"] = i

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
    write_csv(output_dir / "ranked_subfamilies.csv", subfamily_metrics, sf_fields)

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
    write_csv(output_dir / "ranked_subfamilies_detailed_curves.csv", detailed, detail_fields)

    write_report(
        output_dir / "selected_subfamilies_report.txt",
        subfamily_metrics, top_n=50, r1=r1, r2=r2, r3=r3, top_k=TOP_K,
    )

    if WRITE_LATEX:
        write_latex(output_dir / "ranked_subfamilies_latex.tex", subfamily_metrics, top_n=30, r1=r1, r2=r2, r3=r3)

    if WRITE_JSON:
        with open(output_dir / "ranked_subfamilies.json", "w", encoding="utf-8") as f:
            json.dump(subfamily_metrics, f, indent=2, default=str)

    top1 = subfamily_metrics[0]
    n_ge_r3 = sum(1 for m in subfamily_metrics if m.get(f"n_rank_ge_{r3}", 0) > 0)
    print(f"\nResults written to {output_dir}/")
    print(f"  {len(subfamily_metrics)} subfamilies ranked")
    print(f"  With at least 1 rank >= {r3}: {n_ge_r3}")
    print(f"  Best: {top1['family_label']}  ranks={top1['determined_ranks']}")


if __name__ == "__main__":
    main()
