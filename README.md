# Genus 2 Hyperelliptic Curves — Mordell–Weil Rank Study

This repository contains the computational pipeline and results for the dissertation *"Explicit families of genus 2 hyperelliptic curves with many rational points and high-rank Jacobians"*. The project constructs families of genus 2 hyperelliptic curves over ℚ with many rational points, then computes (or bounds) the Mordell–Weil ranks of their Jacobians.

## Overview

The central heuristic is that subfamilies of curves with many rational points may yield Jacobians of large Mordell–Weil rank, because the rational points produce divisor classes that can generate a large-rank subgroup of J(C)(ℚ).

Three families of curves are studied:

| Family | Model | Symmetry | Parameters |
|--------|-------|----------|------------|
| **General interpolation** | F(x) = P_q(x) + H(x)·Q(x) | None imposed | 2 (coefficients of Q) |
| **Even interpolation** | Same, with F(-x) = F(x) | (x,y) → (-x,y) | 2 (A, B in Q(x) = Ax² + B) |
| **D₄ reciprocal** | y² = a²x⁶ + bx⁴ + bx² + a² | Aut_ℚ ≅ D₄ | 1 (parameter t = b/a²) |

## Pipeline

The computational pipeline has four stages:

```
┌─────────────────────┐     ┌────────────────────────┐     ┌─────────────────────┐     ┌──────────────────┐
│  1. Point Search    │────▶│  2. Integral Models    │────▶│  3. Magma Requests  │────▶│  4. Ranking      │
│  (SageMath + PARI)  │     │  (clear denominators)  │     │  (async HTTP)       │     │  (subfamily-level)│
└─────────────────────┘     └────────────────────────┘     └─────────────────────┘     └──────────────────┘
```

### Stage 1: Subfamily Construction and Rational Point Search

**Scripts:**
- `scripts/scan_general_common_forced_families_two_stage_parallel_balanced.py` — General family
- `scripts/scan_even_automorphism_two_stage_parallel_rational_AB.py` — Even family
- `scripts/scan_d4_reciprocal_family.py` — D₄ family

These scripts enumerate subfamilies defined by common forced rational points, then search for additional rational points using PARI's `hyperellratpoints`. A two-stage strategy is used: Stage 1 screens broadly with moderate PARI height; Stage 2 refines only the best subfamilies with much larger height.

Within each subfamily, the top k = 5 curves are kept, ranked by:
1. Affine rational point count (descending)
2. Height of F (ascending)
3. Number of extra rational abscissae (descending)

**Requires:** SageMath (with PARI interface). Run with `conda activate sage && python <script>.py`

### Stage 2: Integral Model Construction

**Script:** `scripts/create_integral_models_from_stage2_details.py`

Converts rational-coefficient models y² = F(x) to integral models Y² = F_integral(X) by clearing denominators. The D₄ family already produces integral models directly.

### Stage 3: Asynchronous Magma Rank Computation

**Scripts:**
- `scripts/request_magma_ranks_from_integral_models_async_robust_v3_global_rate_limiter_FIXED.py` — General + Even families
- `scripts/request_magma_ranks_d4.py` — D₄ family

For each curve, a Magma script is generated and submitted to the [Magma Online Calculator](https://magma.maths.usyd.edu.au/calc/). The Magma script:
1. Computes `RankBounds(J)` → [lower, upper]
2. If lower = upper, rank is determined directly
3. Otherwise runs `MordellWeilGroupGenus2(J : RankOnly := true)`
4. If `finiteIndex = true` and rank matches upper bound, rank is certified

Features:
- Async HTTP with `aiohttp`, global rate limiter (configurable)
- XML response caching and resume logic
- Offline detection with automatic sleep/probe
- Transient failure retries
- All results under GRH assumption

**Requires:** Python 3.10+, `aiohttp` (`pip install aiohttp`)

### Stage 4: Subfamily Ranking

**Scripts:**
- `scripts/rank_subfamilies.py` — Ranks all subfamilies together
- `scripts/rank_subfamilies_by_kind.py` — Ranks within a single construction kind

The ranking key prioritises consistency over outliers:
1. Number of top-k curves with determined rank
2. **Minimum determined rank** (high floor = systematic contribution from forced points)
3. Number with rank ≥ R3 (default 10)
4. Number with rank ≥ R2 (default 8)
5. Number with rank ≥ R1 (default 6)
6. Sum of determined ranks
7. Average rank
8. Point-search quality metrics (affine count, extra abscissae, height)

## Key Results

| Family | Subfamilies | Best Subfamily Ranks | Max Single Rank |
|--------|-------------|---------------------|-----------------|
| General interpolation | 450 | [9, 9, 9, 9, 11] | 11 |
| Even interpolation | 200 | [10, 11, 12, 12, 12] | 13 |
| D₄ reciprocal | 250 curves | rank 8 (best) | 8 |

The even family significantly outperforms the general family — the imposed involution (x,y) → (-x,y) systematically produces higher Mordell–Weil ranks, suggesting the additional symmetry contributes independent divisor classes.

## Repository Structure

```
├── scripts/                           # All Python/SageMath scripts
├── data/
│   ├── general_family/               # Integral model summary CSV (1850 curves)
│   ├── even_family/                  # Integral model summary CSV (1000 curves)
│   └── d4_family/                    # Top curves CSV (250 curves)
├── results/
│   ├── general_family/               # Magma rank results (shared with even family)
│   ├── d4_family/                    # Magma rank results for D₄
│   └── subfamily_rankings/           # Ranked subfamily CSVs and reports
├── dissertation_pipeline_context.txt  # Full pipeline specification
└── README.md
```

## Data Format

### Magma Results CSV (`magma_rank_results_unique.csv`)

Key columns:
- `unique_id` — Unique curve identifier
- `F_integral` — Integral defining polynomial
- `rank_Jacobian` — Certified rank (empty if undetermined)
- `rank_status` — One of: `determined_by_RankBounds`, `determined_by_MordellWeilGroupGenus2`, `magma_memory_limit`, `request_error`
- `rankbounds_lb`, `rankbounds_ub` — Lower/upper bounds from RankBounds
- `finiteIndex`, `proved` — Certification flags from MordellWeilGroupGenus2
- `assumption` — `GRH` if Generalized Riemann Hypothesis was assumed

### Subfamily Rankings CSV (`ranked_subfamilies.csv`)

Key columns:
- `subfamily_rank` — Overall ranking position
- `construction_kind` — `general`, `even_automorphism`, or `d4_reciprocal`
- `family_label`, `rs`, `q_family`, `H_family` — Subfamily identification
- `n_determined` — How many top-k curves have certified rank
- `min_rank`, `max_rank`, `sum_rank`, `avg_rank` — Rank statistics
- `determined_ranks` — List of all determined ranks (e.g., `[10, 11, 12, 12, 12]`)

## Requirements

- **Point search:** SageMath 10.x (via conda/miniforge)
- **Magma requests:** Python 3.10+, `aiohttp`
- **Ranking scripts:** Python 3.10+ (standard library only)

## Reproducing Results

```bash
# 1. Install dependencies
conda create -n sage sage python=3.12
conda activate sage
pip install aiohttp

# 2. Run point search (example: D₄ family, ~10 min)
python scripts/scan_d4_reciprocal_family.py

# 3. Run Magma requests (requires internet; hours due to rate limiting)
python scripts/request_magma_ranks_d4.py

# 4. Run ranking
python scripts/rank_subfamilies.py
```

## Citation

If you use this code or data, please cite the dissertation:
> Dolha, R. (2026). *Explicit families of genus 2 hyperelliptic curves with many rational points and high-rank Jacobians*. Master's dissertation.

## License

MIT
