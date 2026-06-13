# Genus 2 Hyperelliptic Curves — Mordell–Weil Rank Study

Computational pipeline and results for the dissertation *"Explicit families of genus 2 hyperelliptic curves with many rational points and high-rank Jacobians"*. The project constructs families of genus 2 hyperelliptic curves over ℚ with many rational points, then computes the Mordell–Weil ranks of their Jacobians.

## Overview

The central heuristic: subfamilies of curves with many rational points may yield Jacobians of large Mordell–Weil rank, because the rational points produce divisor classes that generate a large-rank subgroup of J(C)(ℚ).

Three families are studied:

| Family | Model | Symmetry | Free parameters |
|--------|-------|----------|-----------------|
| **General interpolation** | F = P_q + H·Q, Q = Ax²+Bx+C, A∈(ℚ×)² | None | 2–3 (B,C or A,B,C) |
| **Even interpolation** | Same, with F(-x)=F(x) | (x,y)→(-x,y) | 2 (A, B in Q = Ax²+B) |
| **D₄ reciprocal** | y² = a²x⁶+bx⁴+bx²+a² | Aut_ℚ ≅ D₄ | 1 (parameter t = b/a²) |

## Pipeline

```
1. Point Search        2. Integral Models      3. Magma Requests       4. Ranking
(SageMath + PARI)  →  (clear denominators)  →  (async HTTP)         →  (subfamily-level)
```

### Stage 1: Subfamily Construction and Rational Point Search

**Scripts:**
- `scripts/scan_general_common_forced_families_two_stage_parallel_balanced.py` — General family (monic Q)
- `scripts/scan_general_non_monic_Q_two_stage_parallel.py` — General family (non-monic Q = Ax²+Bx+C, A∈(ℚ×)²)
- `scripts/scan_even_automorphism_two_stage_parallel_rational_AB.py` — Even family
- `scripts/scan_d4_reciprocal_family.py` — D₄ family

Each script enumerates subfamilies defined by common forced rational points, then searches for additional rational points using PARI's `hyperellratpoints`. A two-stage strategy is used:
- **Stage 1:** broad screening with moderate PARI height
- **Stage 2:** refinement of only the best subfamilies with much larger height

The non-monic variant extends the general family by allowing Q(x) = Ax² + Bx + C with A a rational square, increasing the parameter space from 2 to 3 free coefficients per subfamily while preserving rationality of the points at infinity.

Configuration parameters (edit at top of each script):
- `M_VALUES` — number of forced abscissae (4, 5, or 6)
- `A_BASE_POOL` — leading coefficient values for Q (non-monic script)
- `PARI_HEIGHT_STAGE1`, `PARI_HEIGHT_STAGE2` — search height bounds
- `MAX_FORCED_FAMILIES` — subfamilies to screen in Stage 1
- `TOP_FAMILIES_TO_REFINE` — subfamilies promoted to Stage 2

**Requires:** SageMath ≥ 10.x (via conda). Run: `conda activate sage && python <script>.py`

### Stage 2: Integral Model Construction

**Script:** `scripts/create_integral_models_from_stage2_details.py`

Converts rational-coefficient models y² = F(x) to integral models Y² = F_integral(X) via the coordinate change (X, Y) → (X, L·Y), where L = lcm of denominators of F's coefficients. The script scans all configured family directories for `stage2_refined_details.txt` files and processes them uniformly. Configure `SOURCE_ROOTS` at the top to add/remove family directories.

### Stage 3: Asynchronous Magma Rank Computation

**Script:** `scripts/request_magma_ranks_from_integral_models_async_robust_v3_global_rate_limiter_FIXED.py`

A single script handles all families. Configure `INPUT_CSVS` and `OUTPUT_DIR` at the top for each run.

For each curve, a Magma script is submitted to the [Magma Online Calculator](https://magma.maths.usyd.edu.au/calc/) that:
1. Computes `RankBounds(J)` → [lower, upper]
2. If equal, rank is determined directly
3. Otherwise runs `MordellWeilGroupGenus2(J : RankOnly := true)`
4. If `finiteIndex = true` and rank matches upper bound, rank is certified

Features:
- Async HTTP with `aiohttp`, configurable global rate limiter
- XML response caching and resume logic
- Offline detection with automatic sleep/probe
- Transient failure retries
- All results computed under GRH assumption

**Requires:** Python 3.10+, `aiohttp` (`pip install aiohttp`)

### Stage 4: Subfamily Ranking

**Scripts:**
- `scripts/rank_subfamilies.py` — Rank all subfamilies combined
- `scripts/rank_subfamilies_by_kind.py` — Rank within a single construction kind

The ranking key prioritises consistency over outliers:
1. Number of top-k curves with determined rank
2. **Minimum determined rank** (high floor = systematic contribution from forced points)
3. Number with rank ≥ R3, R2, R1 (configurable thresholds, default 10, 8, 6)
4. Sum and average of determined ranks
5. Point-search quality tiebreakers

## Key Results

| Family | Subfamilies | Best subfamily ranks | Max rank |
|--------|-------------|---------------------|----------|
| General (monic + non-monic) | 525 | [9, 9, 10, 10, 11] | 11 |
| Even interpolation | 200 | [10, 11, 12, 12, 12] | 13 |
| D₄ reciprocal | 250 curves | [6, 8, 8, 8, 8] | 8 |

The even family significantly outperforms the general family — the imposed involution (x,y) → (-x,y) systematically produces higher Mordell–Weil ranks. The non-monic enhancement of the general family improved the best result from [9,9,9,9,11] to [9,9,10,10,11].

## Repository Structure

```
├── scripts/                              # Python/SageMath scripts (8 files)
├── data/
│   ├── general_family/                   # Integral model summary (2225 curves: monic + non-monic)
│   ├── even_family/                      # Integral model summary (1000 curves)
│   └── d4_family/                        # Top curves summary (250 curves)
├── results/
│   ├── general_and_even_families/        # Magma rank results (combined run)
│   ├── d4_family/                        # Magma rank results for D₄
│   └── subfamily_rankings/              # Ranked CSVs + human-readable reports
├── dissertation_pipeline_context.txt     # Full pipeline specification
└── README.md
```

## Data Format

### Magma Results CSV (`magma_rank_results_unique.csv`)

Key columns:
- `unique_id` — Unique curve identifier
- `F_integral` — Integral defining polynomial
- `rank_Jacobian` — Certified rank (empty if undetermined)
- `rank_status` — One of: `determined_by_RankBounds`, `determined_by_MordellWeilGroupGenus2`, `magma_memory_limit`, `request_error`
- `rankbounds_lb`, `rankbounds_ub` — Bounds from RankBounds
- `finiteIndex`, `proved` — Certification flags from MordellWeilGroupGenus2
- `assumption` — `GRH` if Generalized Riemann Hypothesis was assumed

### Subfamily Rankings CSV (`ranked_subfamilies_*.csv`)

Key columns:
- `subfamily_rank` — Ranking position
- `construction_kind` — `general`, `even_automorphism`, or `d4_reciprocal`
- `family_label`, `family_index`, `execution_relative_dir` — Subfamily identification
- `n_determined` — How many top-k curves have certified rank
- `min_rank`, `max_rank`, `sum_rank`, `avg_rank` — Rank statistics
- `determined_ranks` — List of determined ranks (e.g., `[10, 11, 12, 12, 12]`)

## Reproducing Results

```bash
# Install dependencies
conda create -n sage sage python=3.12
conda activate sage
pip install aiohttp

# 1. Point search (example: D₄ family, ~10 min)
python scripts/scan_d4_reciprocal_family.py

# 2. Integral models (< 1 min)
python scripts/create_integral_models_from_stage2_details.py

# 3. Magma rank requests (hours; edit INPUT_CSVS and OUTPUT_DIR first)
python scripts/request_magma_ranks_from_integral_models_async_robust_v3_global_rate_limiter_FIXED.py

# 4. Ranking (run from scripts/ directory)
cd scripts
python rank_subfamilies.py
```

**Note:** The `data/` and `results/` folders contain aggregated outputs from multiple execution runs with different hyperparameter configurations. The scripts document the methodology; the CSVs are the final merged results. To reproduce from scratch, run the scanner scripts with your chosen parameters, then feed outputs through the subsequent pipeline stages.

## Requirements

- **Point search:** SageMath 10.x (via conda/miniforge)
- **Magma requests:** Python 3.10+, `aiohttp`
- **Ranking scripts:** Python 3.10+ (standard library only)

## Citation

If you use this code or data, please cite the dissertation:
> Dolha, R. (2026). *Explicit families of genus 2 hyperelliptic curves with many rational points and high-rank Jacobians*. Master's dissertation.

## License

MIT
