#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

"""
Two-stage parallel scanner for EVEN genus-2 hyperelliptic families with automorphism
    (x, y) -> (-x, y).

STAGE 1:
    Broad, cheap screening with a relatively small PARI height.

STAGE 2:
    Refine only the best STAGE1 families with a much larger PARI height.

This version allows rational values for the base a in A=a^2 and rational values for B.

Construction:
    C : y^2 = F(x),
    F(x) = P_q(x) + H(x) Q(x),

where:
    q(x) = u*x^2 + w,
    H(x) = prod_j (x^2 - r_j^2),
    Q(x) = A*x^2 + B,

and A is a rational square, so F has square leading coefficient.
All polynomials are even, hence the curve has the rational automorphism:
    (x, y) -> (-x, y).

Run with:
    sage -python scan_even_automorphism_two_stage_parallel.py
"""

from sage.all import *
from itertools import combinations, product
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import os
import time


# ============================================================
# RATIONAL PARAMETER POOLS
# ============================================================

def rational_pool_signed(num_bound, den_bound, include_zero=True):
    """
    Returns reduced rationals a/b with |a| <= num_bound and 1 <= b <= den_bound.
    Example: num_bound=3, den_bound=2 gives values such as -3, -3/2, -1/2, 0, 1/2, 3/2, 3.
    """
    vals = set()
    for b in range(1, den_bound + 1):
        for a in range(-num_bound, num_bound + 1):
            if a == 0 and not include_zero:
                continue
            if gcd(ZZ(a), ZZ(b)) == 1:
                vals.add(QQ(a) / QQ(b))
    return sorted(vals)


def rational_pool_positive(num_bound, den_bound):
    """
    Returns positive reduced rationals a/b with 1 <= a <= num_bound and 1 <= b <= den_bound.
    This is intended for the base a of A=a^2, not for A itself.
    """
    vals = set()
    for b in range(1, den_bound + 1):
        for a in range(1, num_bound + 1):
            if gcd(ZZ(a), ZZ(b)) == 1:
                vals.add(QQ(a) / QQ(b))
    return sorted(vals)


# ============================================================
# GLOBAL FAMILY CONFIGURATION
# ============================================================

N_FORCED_PAIRS = 2

RATIONAL_R_NUM_BOUND = 4
RATIONAL_R_DEN_BOUND = 4
MAX_R_TUPLES = 300

Q_COEFF_POOL = list(range(-3, 4))
ALLOW_ZERO_Q = False
MAX_ZERO_FORCED_ABS_ROOTS = 0

MAX_FORCED_FAMILIES = 2000

DEDUP_FAMILIES_BY_H_AND_PQ = True
DEDUP_CURVES_WITHIN_FAMILY = True


# ============================================================
# STAGE 1: BROAD SCREENING
# ============================================================

PARI_HEIGHT_STAGE1 = 3000

# A = a^2. Here a is allowed to be rational.
# Example: rational_pool_positive(10, 3) gives a in {1/3, 1/2, 2/3, ..., 10} after reduction.
A_BASE_POOL_STAGE1 = rational_pool_positive(num_bound=10, den_bound=5)

# B is allowed to be rational, independently of A.
# Example: rational_pool_signed(15, 3) gives B in [-15,15] with denominators <= 3.
B_PARAM_POOL_STAGE1 = rational_pool_signed(num_bound=15, den_bound=3, include_zero=True)

TOP_FAMILIES_TO_REFINE = 100


# ============================================================
# STAGE 2: EXPENSIVE REFINEMENT
# ============================================================

PARI_HEIGHT_STAGE2 = 80000

# Usually keep Stage 2 on the same A/B grid as Stage 1.
# You can make these larger, but Stage 2 is already expensive.
A_BASE_POOL_STAGE2 = A_BASE_POOL_STAGE1
B_PARAM_POOL_STAGE2 = B_PARAM_POOL_STAGE1

TOP_FINAL_FAMILIES = 50


# ============================================================
# SCORING CONFIGURATION
# ============================================================

TOP_CURVES_PER_FAMILY = 5
GOOD_EXTRA_X_THRESHOLD = 26
TARGET_AFFINE_1 = 55
TARGET_AFFINE_2 = 60


# ============================================================
# MULTIPROCESSING CONFIGURATION
# ============================================================

USE_PARALLEL = True
# MAX_WORKERS = 4
# Alternative:
MAX_WORKERS = max(1, (os.cpu_count() or 2) - 5)


# ============================================================
# OUTPUT CONFIGURATION
# ============================================================

OUTPUT_BASE_DIR = Path("families_construct_even_automorphism_two_stage")


# ============================================================
# SAGE SETUP
# ============================================================

R = PolynomialRing(QQ, "x")
x = R.gen()
OUTPUT_DIR = None


# ============================================================
# OUTPUT DIRECTORY AND RUN CONFIG
# ============================================================

def rational_to_name(q):
    """
    Transformă un număr rațional într-un string sigur pentru nume de folder.

    Exemple:
        3      -> "3"
        -3     -> "-3"
        1/5    -> "1over5"
        -2/3   -> "-2over3"
    """
    q = QQ(q)
    n = q.numerator()
    d = q.denominator()

    if d == 1:
        return str(n)

    return f"{n}over{d}"


def safe_path_component(s):
    """
    Elimină/caracterizează caracterele problematice pentru path-uri Windows/Linux.
    """
    s = str(s)

    replacements = {
        "/": "over",
        "\\": "over",
        ":": "-",
        "*": "",
        "?": "",
        '"': "",
        "<": "",
        ">": "",
        "|": "",
        " ": "",
    }

    for bad, good in replacements.items():
        s = s.replace(bad, good)

    return s


def pool_range_name(pool):
    """
    Returnează un nume sigur pentru intervalul unui pool.
    Funcționează și pentru valori raționale Sage QQ.
    """
    pool = list(pool)

    if not pool:
        return "empty"

    left = rational_to_name(min(pool))
    right = rational_to_name(max(pool))

    return safe_path_component(f"{left}to{right}")


def make_run_name():
    """
    Nume detaliat, dar safe, pentru folderul rulării.
    Nu conține '/', deci nu mai creează subdirectoare accidentale.
    """
    run_name = (
        f"even_pairs{N_FORCED_PAIRS}"
        f"_H1_{PARI_HEIGHT_STAGE1}"
        f"_H2_{PARI_HEIGHT_STAGE2}"
        f"_topRef{TOP_FAMILIES_TO_REFINE}"
        f"_goodX{GOOD_EXTRA_X_THRESHOLD}"
        f"_rN{RATIONAL_R_NUM_BOUND}_rD{RATIONAL_R_DEN_BOUND}"
        f"_maxR{MAX_R_TUPLES}"
        f"_maxFam{MAX_FORCED_FAMILIES}"
        f"_q{pool_range_name(Q_COEFF_POOL)}"
        f"_A1{pool_range_name(A_BASE_POOL_STAGE1)}sq"
        f"_nA{len(set(QQ(a)**2 for a in A_BASE_POOL_STAGE1))}"
        f"_B1{pool_range_name(B_PARAM_POOL_STAGE1)}"
        f"_nB{len(set(QQ(b) for b in B_PARAM_POOL_STAGE1))}"
        f"_workers{MAX_WORKERS if USE_PARALLEL else 1}"
    )

    return safe_path_component(run_name)

def init_output_dir():
    global OUTPUT_DIR
    OUTPUT_DIR = OUTPUT_BASE_DIR / make_run_name()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def save_run_config():
    path = Path(OUTPUT_DIR) / "run_config.txt"
    with open(path, "w", encoding="utf-8") as out:
        out.write("CONFIGURATION FOR TWO-STAGE EVEN-AUTOMORPHISM FAMILY SCAN\n")
        out.write("=" * 100 + "\n")
        out.write(f"N_FORCED_PAIRS = {N_FORCED_PAIRS}\n")
        out.write(f"RATIONAL_R_NUM_BOUND = {RATIONAL_R_NUM_BOUND}\n")
        out.write(f"RATIONAL_R_DEN_BOUND = {RATIONAL_R_DEN_BOUND}\n")
        out.write(f"MAX_R_TUPLES = {MAX_R_TUPLES}\n")
        out.write(f"Q_COEFF_POOL = {Q_COEFF_POOL}\n")
        out.write(f"ALLOW_ZERO_Q = {ALLOW_ZERO_Q}\n")
        out.write(f"MAX_ZERO_FORCED_ABS_ROOTS = {MAX_ZERO_FORCED_ABS_ROOTS}\n")
        out.write(f"MAX_FORCED_FAMILIES = {MAX_FORCED_FAMILIES}\n")
        out.write(f"DEDUP_FAMILIES_BY_H_AND_PQ = {DEDUP_FAMILIES_BY_H_AND_PQ}\n")
        out.write(f"DEDUP_CURVES_WITHIN_FAMILY = {DEDUP_CURVES_WITHIN_FAMILY}\n")
        out.write("\nSTAGE 1\n")
        out.write("-" * 100 + "\n")
        out.write(f"PARI_HEIGHT_STAGE1 = {PARI_HEIGHT_STAGE1}\n")
        out.write(f"A_BASE_POOL_STAGE1 = {A_BASE_POOL_STAGE1}\n")
        out.write(f"A_SQUARE_POOL_STAGE1 = {[QQ(a)**2 for a in A_BASE_POOL_STAGE1]}\n")
        out.write(f"B_PARAM_POOL_STAGE1 = {B_PARAM_POOL_STAGE1}\n")
        out.write(f"TOP_FAMILIES_TO_REFINE = {TOP_FAMILIES_TO_REFINE}\n")
        out.write("\nSTAGE 2\n")
        out.write("-" * 100 + "\n")
        out.write(f"PARI_HEIGHT_STAGE2 = {PARI_HEIGHT_STAGE2}\n")
        out.write(f"A_BASE_POOL_STAGE2 = {A_BASE_POOL_STAGE2}\n")
        out.write(f"A_SQUARE_POOL_STAGE2 = {[QQ(a)**2 for a in A_BASE_POOL_STAGE2]}\n")
        out.write(f"B_PARAM_POOL_STAGE2 = {B_PARAM_POOL_STAGE2}\n")
        out.write(f"TOP_FINAL_FAMILIES = {TOP_FINAL_FAMILIES}\n")
        out.write("\nSCORING\n")
        out.write("-" * 100 + "\n")
        out.write(f"TOP_CURVES_PER_FAMILY = {TOP_CURVES_PER_FAMILY}\n")
        out.write(f"GOOD_EXTRA_X_THRESHOLD = {GOOD_EXTRA_X_THRESHOLD}\n")
        out.write(f"TARGET_AFFINE_1 = {TARGET_AFFINE_1}\n")
        out.write(f"TARGET_AFFINE_2 = {TARGET_AFFINE_2}\n")
        out.write("\nMULTIPROCESSING\n")
        out.write("-" * 100 + "\n")
        out.write(f"USE_PARALLEL = {USE_PARALLEL}\n")
        out.write(f"MAX_WORKERS = {MAX_WORKERS}\n")
        out.write(f"os.cpu_count() = {os.cpu_count()}\n")
        out.write(f"OUTPUT_DIR = {OUTPUT_DIR}\n")


# ============================================================
# BASIC HELPERS
# ============================================================

def rational_height(a):
    a = QQ(a)
    return max(abs(a.numerator()), a.denominator())


def poly_height(F):
    F = R(F)
    return max([rational_height(c) for c in F.coefficients()] + [0])


def is_squarefree_sextic(F):
    F = R(F)
    return F.degree() == 6 and F.gcd(F.derivative()).degree() == 0


def is_even_polynomial(F):
    F = R(F)
    return all(F.monomial_coefficient(x**k) == 0 for k in range(1, F.degree() + 1, 2))


def leading_coeff_is_square_Q(F):
    F = R(F)
    lc = QQ(F.leading_coefficient())
    return lc != 0 and lc.is_square()


def polynomial_key(F):
    F = R(F)
    coeffs = list(F.list())
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return tuple(QQ(c) for c in coeffs)


def rational_positive_pool(num_bound, den_bound):
    vals = set()
    for b in range(1, den_bound + 1):
        for a in range(1, num_bound + 1):
            if gcd(ZZ(a), ZZ(b)) == 1:
                vals.add(QQ(a) / QQ(b))
    return sorted(vals)


def forced_abscissae_from_abs_roots(abs_roots):
    xs = []
    for r in abs_roots:
        r = QQ(r)
        xs.append(-r)
        xs.append(r)
    return tuple(sorted(set(xs)))


def forced_points(q, abs_roots):
    pts = set()
    for xx in forced_abscissae_from_abs_roots(abs_roots):
        yy = QQ(q(xx))
        pts.add((QQ(xx), yy))
        pts.add((QQ(xx), -yy))
    return sorted(pts, key=lambda P: (P[0], P[1]))


# ============================================================
# EVEN FAMILY CONSTRUCTION
# ============================================================

def build_even_family_member(abs_roots, q, A, B):
    abs_roots = tuple(QQ(r) for r in abs_roots)

    if any(r == 0 for r in abs_roots):
        raise ValueError("Absolute forced roots must be non-zero.")
    if len(set(abs_roots)) != len(abs_roots):
        raise ValueError("Absolute forced roots must be distinct.")

    q = R(q)
    A = QQ(A)
    B = QQ(B)

    H = R(prod(x**2 - r**2 for r in abs_roots))
    Pq = R((q**2) % H)
    Qpoly = R(A*x**2 + B)
    F = R(Pq + H*Qpoly)

    return {
        "abs_roots": abs_roots,
        "forced_xs": forced_abscissae_from_abs_roots(abs_roots),
        "q": q,
        "H": H,
        "Pq": Pq,
        "Q": Qpoly,
        "A": A,
        "B": B,
        "F": F,
    }


def A_square_pool_from_base(A_base_pool):
    vals = sorted(set(QQ(a)**2 for a in A_base_pool))
    vals = [a for a in vals if a != 0]
    return vals


def AB_grid(A_base_pool, B_param_pool):
    for A in A_square_pool_from_base(A_base_pool):
        for B in B_param_pool:
            yield QQ(A), QQ(B)


# ============================================================
# PARI POINT SEARCH
# ============================================================

def pari_height_argument(height):
    return ZZ(height)


def pari_points_to_python(pts):
    if hasattr(pts, "sage"):
        pts = pts.sage()
    elif hasattr(pts, "python"):
        pts = pts.python()
    return pts


def normalize_affine_point(pt):
    return (QQ(pt[0]), QQ(pt[1]))


def search_with_pari_on_poly(F, height):
    F = R(F)
    pts = pari(F).hyperellratpoints(pari_height_argument(height), 0)
    pts = pari_points_to_python(pts)

    out = set()
    for pt in pts:
        P = normalize_affine_point(pt)
        if P[1]**2 == F(P[0]):
            out.add(P)

    return sorted(out, key=lambda P: (P[0], P[1]))


def split_forced_and_extra(points, abs_roots):
    forced_xs = set(forced_abscissae_from_abs_roots(abs_roots))
    forced = [P for P in points if P[0] in forced_xs]
    extra = [P for P in points if P[0] not in forced_xs]
    return forced, extra


# ============================================================
# CURVE AND FAMILY SCORING
# ============================================================

def curve_record(data, pari_height):
    F = data["F"]
    abs_roots = data["abs_roots"]
    q = data["q"]

    if not is_squarefree_sextic(F):
        return None
    if not is_even_polynomial(F):
        return None
    if not leading_coeff_is_square_Q(F):
        return None

    try:
        pts = search_with_pari_on_poly(F, pari_height)
    except Exception:
        return None

    forced_found, extra = split_forced_and_extra(pts, abs_roots)
    extra_xs = sorted(set(P[0] for P in extra))

    rec = dict(data)
    rec.update({
        "height_F": poly_height(F),
        "forced_expected": forced_points(q, abs_roots),
        "forced_found": forced_found,
        "extra_points": extra,
        "extra_xs": extra_xs,
        "extra_x_count": len(extra_xs),
        "extra_point_count": len(extra),
        "affine_point_count": len(pts),
        "projective_lower_bound": len(pts) + 2,
        "pari_height": pari_height,
    })
    return rec


def curve_sort_key(rec):
    return (
        -rec["affine_point_count"],
        rec["height_F"],
        -rec["extra_x_count"],
    )


def summarize_family(abs_roots, q, q_params, q_coeffs, curve_records, stage_name):
    if not curve_records:
        return None

    records = sorted(curve_records, key=curve_sort_key)
    top = records[:TOP_CURVES_PER_FAMILY]
    best = top[0]

    top_affine_counts = [r["affine_point_count"] for r in top]
    top_projective_counts = [r["projective_lower_bound"] for r in top]
    top_extra_x_counts = [r["extra_x_count"] for r in top]

    good_curve_count = sum(
        1 for r in records
        if r["extra_x_count"] >= GOOD_EXTRA_X_THRESHOLD
    )

    return {
        "stage": stage_name,
        "abs_roots": tuple(QQ(r) for r in abs_roots),
        "forced_xs": forced_abscissae_from_abs_roots(abs_roots),
        "q": R(q),
        "q_params": dict(q_params),
        "q_coeffs": tuple(q_coeffs),
        "forced_points": forced_points(q, abs_roots),

        "num_curves_tested_good": len(records),
        "good_curve_count": good_curve_count,

        "best_extra_x_count": best["extra_x_count"],
        "best_extra_point_count": best["extra_point_count"],
        "best_affine_point_count": best["affine_point_count"],
        "best_projective_lower_bound": best["projective_lower_bound"],
        "best_height_F": best["height_F"],

        "top_k_count": len(top),
        "top_k_min_affine": min(top_affine_counts),
        "top_k_max_affine": max(top_affine_counts),
        "top_k_sum_affine": sum(top_affine_counts),
        "top_k_avg_affine": QQ(sum(top_affine_counts)) / QQ(len(top_affine_counts)),

        "top_k_min_projective": min(top_projective_counts),
        "top_k_sum_projective": sum(top_projective_counts),

        "top_k_min_extra_x": min(top_extra_x_counts),
        "top_k_sum_extra_x": sum(top_extra_x_counts),

        "top_k_affine_ge_target1": sum(1 for v in top_affine_counts if v >= TARGET_AFFINE_1),
        "top_k_affine_ge_target2": sum(1 for v in top_affine_counts if v >= TARGET_AFFINE_2),

        "top_curves": top,
    }


def family_sort_key(summary):
    return (
        -summary["top_k_affine_ge_target2"],
        -summary["top_k_affine_ge_target1"],
        -summary["top_k_min_affine"],
        -summary["top_k_sum_affine"],
        -summary["best_affine_point_count"],
        -summary["best_projective_lower_bound"],
        summary["best_height_F"],
    )


# ============================================================
# FORCED FAMILY GENERATION
# ============================================================

def generate_abs_root_tuples():
    pool = rational_positive_pool(RATIONAL_R_NUM_BOUND, RATIONAL_R_DEN_BOUND)
    tuples = list(combinations(pool, N_FORCED_PAIRS))

    tuples.sort(key=lambda rs: (
        max(rational_height(r) for r in rs),
        sum(rational_height(r) for r in rs),
        sum(abs(QQ(r)) for r in rs),
    ))

    if MAX_R_TUPLES is not None and len(tuples) > MAX_R_TUPLES:
        tuples = tuples[:MAX_R_TUPLES]

    return tuples


def generate_even_q_polys():
    qs = []

    for u, w in product(Q_COEFF_POOL, Q_COEFF_POOL):
        if not ALLOW_ZERO_Q and u == 0 and w == 0:
            continue

        q = R(QQ(u)*x**2 + QQ(w))
        q_params = {"u": u, "w": w}
        q_coeffs = (u, w)
        qs.append((q, q_params, q_coeffs))

    qs.sort(key=lambda item: (
        poly_height(item[0]),
        str(item[0]),
    ))

    return qs


def forced_family_candidates():
    root_tuples = generate_abs_root_tuples()
    q_items = generate_even_q_polys()

    count = 0
    seen_family_keys = set()

    for abs_roots in root_tuples:
        for q, q_params, q_coeffs in q_items:
            zero_abs_count = sum(1 for r in abs_roots if q(QQ(r)) == 0)
            if zero_abs_count > MAX_ZERO_FORCED_ABS_ROOTS:
                continue

            H = R(prod(x**2 - QQ(r)**2 for r in abs_roots))
            Pq = R((R(q)**2) % H)

            fam_key = (
                tuple(QQ(r) for r in abs_roots),
                polynomial_key(H),
                polynomial_key(Pq),
            )

            if DEDUP_FAMILIES_BY_H_AND_PQ and fam_key in seen_family_keys:
                continue

            seen_family_keys.add(fam_key)

            yield abs_roots, q, q_params, q_coeffs

            count += 1
            if MAX_FORCED_FAMILIES is not None and count >= MAX_FORCED_FAMILIES:
                return


def candidate_from_summary(summary):
    abs_roots = tuple(QQ(r) for r in summary["abs_roots"])
    q_coeffs = tuple(summary["q_coeffs"])
    q_params = dict(summary["q_params"])
    u, w = q_coeffs
    q = R(QQ(u)*x**2 + QQ(w))
    return abs_roots, q, q_params, q_coeffs


# ============================================================
# SCANNING
# ============================================================

def scan_one_even_family(abs_roots, q, q_params, q_coeffs,
                         pari_height, A_base_pool, B_param_pool, stage_name):
    curve_records = []
    seen_F = set()

    for A, B in AB_grid(A_base_pool, B_param_pool):
        data = build_even_family_member(abs_roots=abs_roots, q=q, A=A, B=B)

        if DEDUP_CURVES_WITHIN_FAMILY:
            key = polynomial_key(data["F"])
            if key in seen_F:
                continue
            seen_F.add(key)

        rec = curve_record(data, pari_height=pari_height)
        if rec is None:
            continue

        rec["q_params"] = dict(q_params)
        rec["q_coeffs"] = tuple(q_coeffs)
        curve_records.append(rec)

    return summarize_family(
        abs_roots=abs_roots,
        q=q,
        q_params=q_params,
        q_coeffs=q_coeffs,
        curve_records=curve_records,
        stage_name=stage_name,
    )


def scan_one_even_family_worker(args):
    (
        abs_roots_raw,
        q_coeffs,
        q_params,
        pari_height,
        A_base_raw,
        B_raw,
        stage_name,
    ) = args

    abs_roots = tuple(QQ(r) for r in abs_roots_raw)
    u, w = q_coeffs
    q = R(QQ(u)*x**2 + QQ(w))

    A_base_pool = [QQ(a) for a in A_base_raw]
    B_param_pool = [QQ(b) for b in B_raw]

    try:
        return scan_one_even_family(
            abs_roots=abs_roots,
            q=q,
            q_params=q_params,
            q_coeffs=q_coeffs,
            pari_height=pari_height,
            A_base_pool=A_base_pool,
            B_param_pool=B_param_pool,
            stage_name=stage_name,
        )
    except Exception:
        return None


def scan_candidate_list(candidates, stage_name, pari_height, A_base_pool, B_param_pool,
                        retain_top=None):
    print("=" * 100)
    print(f"{stage_name}: scanning {len(candidates)} forced families")
    print(f"PARI height = {pari_height}")
    print(f"A square pool = {[QQ(a)**2 for a in A_base_pool]}")
    print(f"B pool = {B_param_pool}")
    print("=" * 100)

    start = time.time()
    summaries = []

    if not USE_PARALLEL or MAX_WORKERS <= 1:
        completed = 0
        for abs_roots, q, q_params, q_coeffs in candidates:
            completed += 1

            summary = scan_one_even_family(
                abs_roots=abs_roots,
                q=q,
                q_params=q_params,
                q_coeffs=q_coeffs,
                pari_height=pari_height,
                A_base_pool=A_base_pool,
                B_param_pool=B_param_pool,
                stage_name=stage_name,
            )

            if summary is not None:
                summaries.append(summary)

            if completed % 25 == 0 or completed == len(candidates):
                elapsed = time.time() - start
                print(
                    f"{stage_name}: completed={completed}/{len(candidates)}, "
                    f"nonempty summaries={len(summaries)}, elapsed={elapsed:.1f}s"
                )

    else:
        tasks = []
        A_base_raw = [str(QQ(a)) for a in A_base_pool]
        B_raw = [str(QQ(b)) for b in B_param_pool]

        for abs_roots, q, q_params, q_coeffs in candidates:
            abs_roots_raw = tuple(str(QQ(r)) for r in abs_roots)
            tasks.append((
                abs_roots_raw,
                q_coeffs,
                q_params,
                pari_height,
                A_base_raw,
                B_raw,
                stage_name,
            ))

        completed = 0

        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
            future_to_task = {
                executor.submit(scan_one_even_family_worker, task): task
                for task in tasks
            }

            for future in as_completed(future_to_task):
                completed += 1

                try:
                    summary = future.result()
                except Exception:
                    summary = None

                if summary is not None:
                    summaries.append(summary)

                if completed % 25 == 0 or completed == len(candidates):
                    elapsed = time.time() - start
                    print(
                        f"{stage_name}: completed={completed}/{len(candidates)}, "
                        f"nonempty summaries={len(summaries)}, elapsed={elapsed:.1f}s"
                    )

    summaries.sort(key=family_sort_key)

    if retain_top is not None:
        return summaries[:retain_top]

    return summaries


# ============================================================
# OUTPUT
# ============================================================

def write_curve_record(out, rec, index):
    out.write("-" * 90 + "\n")
    out.write(f"CURBA #{index}\n")
    out.write(f"q_params = {rec['q_params']}\n")
    out.write(f"A = {rec['A']}\n")
    out.write(f"B = {rec['B']}\n")
    out.write(f"Q(x) = {rec['Q']}\n")
    out.write(f"F(x) = {rec['F']}\n")
    out.write(f"Height(F) = {rec['height_F']}\n")
    out.write(f"PARI height = {rec['pari_height']}\n")
    out.write(f"Nr. x-uri suplimentare distincte = {rec['extra_x_count']}\n")
    out.write(f"Nr. puncte suplimentare = {rec['extra_point_count']}\n")
    out.write(f"Nr. puncte afine gasite = {rec['affine_point_count']}\n")
    out.write(f"Lower bound projectiv = {rec['projective_lower_bound']}\n")
    out.write(f"x-uri suplimentare = {rec['extra_xs']}\n")
    out.write(f"Puncte suplimentare = {rec['extra_points']}\n")


def save_family_summaries(summaries, filename, title):
    path = Path(OUTPUT_DIR) / filename

    with open(path, "w", encoding="utf-8") as out:
        out.write("=" * 100 + "\n")
        out.write(title + "\n")
        out.write("EVEN FAMILIES OF GENUS-2 HYPERELLIPTIC CURVES\n")
        out.write("Automorphism: (x,y) -> (-x,y)\n")
        out.write("=" * 100 + "\n\n")

        for i, fam in enumerate(summaries, start=1):
            best = fam["top_curves"][0]

            out.write("=" * 100 + "\n")
            out.write(f"FAMILIA #{i}\n")
            out.write(f"stage = {fam['stage']}\n")
            out.write(f"abs_roots = {fam['abs_roots']}\n")
            out.write(f"forced_xs = {fam['forced_xs']}\n")
            out.write(f"q_params = {fam['q_params']}\n")
            out.write(f"q(x) = {fam['q']}\n")
            out.write(f"H(x) = {best['H']}\n")
            out.write(f"Pq(x) = {best['Pq']}\n")
            out.write(f"Puncte comune fortate = {fam['forced_points']}\n")
            out.write(f"num_curves_tested_good = {fam['num_curves_tested_good']}\n")
            out.write(f"good_curve_count(extra_x >= {GOOD_EXTRA_X_THRESHOLD}) = {fam['good_curve_count']}\n")
            out.write(f"top_k_count = {fam['top_k_count']}\n")
            out.write(f"top_k_min_affine = {fam['top_k_min_affine']}\n")
            out.write(f"top_k_max_affine = {fam['top_k_max_affine']}\n")
            out.write(f"top_k_sum_affine = {fam['top_k_sum_affine']}\n")
            out.write(f"top_k_avg_affine = {fam['top_k_avg_affine']}\n")
            out.write(f"top_k_affine_ge_{TARGET_AFFINE_1} = {fam['top_k_affine_ge_target1']}\n")
            out.write(f"top_k_affine_ge_{TARGET_AFFINE_2} = {fam['top_k_affine_ge_target2']}\n")
            out.write(f"top_k_min_extra_x = {fam['top_k_min_extra_x']}\n")
            out.write(f"top_k_sum_extra_x = {fam['top_k_sum_extra_x']}\n")
            out.write(f"best_extra_x_count = {fam['best_extra_x_count']}\n")
            out.write(f"best_affine_point_count = {fam['best_affine_point_count']}\n")
            out.write(f"best_projective_lower_bound = {fam['best_projective_lower_bound']}\n")
            out.write(f"best_height_F = {fam['best_height_F']}\n")
            out.write("\n")

            out.write("TOP CURVES IN THIS FAMILY\n")
            for j, rec in enumerate(fam["top_curves"], start=1):
                write_curve_record(out, rec, j)

            out.write("\n")

    print(f"Saved: {path}")


def save_family_csv(summaries, filename):
    path = Path(OUTPUT_DIR) / filename

    header = [
        "family_index",
        "stage",
        "abs_roots",
        "forced_xs",
        "q_params",
        "q",
        "H",
        "Pq",
        "forced_points",
        "num_curves_tested_good",
        "good_curve_count",
        "top_k_min_affine",
        "top_k_max_affine",
        "top_k_sum_affine",
        "top_k_avg_affine",
        f"top_k_affine_ge_{TARGET_AFFINE_1}",
        f"top_k_affine_ge_{TARGET_AFFINE_2}",
        "top_k_min_extra_x",
        "top_k_sum_extra_x",
        "best_extra_x_count",
        "best_affine_point_count",
        "best_projective_lower_bound",
        "best_height_F",
        "best_A",
        "best_B",
        "best_Q",
        "best_F",
    ]

    with open(path, "w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout)
        writer.writerow(header)

        for i, fam in enumerate(summaries, start=1):
            best = fam["top_curves"][0]

            row = [
                i,
                fam["stage"],
                str(fam["abs_roots"]),
                str(fam["forced_xs"]),
                str(fam["q_params"]),
                str(fam["q"]),
                str(best["H"]),
                str(best["Pq"]),
                str(fam["forced_points"]),
                fam["num_curves_tested_good"],
                fam["good_curve_count"],
                fam["top_k_min_affine"],
                fam["top_k_max_affine"],
                fam["top_k_sum_affine"],
                str(fam["top_k_avg_affine"]),
                fam["top_k_affine_ge_target1"],
                fam["top_k_affine_ge_target2"],
                fam["top_k_min_extra_x"],
                fam["top_k_sum_extra_x"],
                fam["best_extra_x_count"],
                fam["best_affine_point_count"],
                fam["best_projective_lower_bound"],
                fam["best_height_F"],
                str(best["A"]),
                str(best["B"]),
                str(best["Q"]),
                str(best["F"]),
            ]

            writer.writerow(row)

    print(f"Saved: {path}")


def save_selected_candidates(stage1_summaries):
    path = Path(OUTPUT_DIR) / "selected_families_for_stage2.csv"

    header = [
        "selected_index",
        "abs_roots",
        "q_params",
        "q",
        "stage1_best_affine_point_count",
        "stage1_best_extra_x_count",
        "stage1_top_k_min_affine",
        "stage1_top_k_sum_affine",
        "stage1_best_A",
        "stage1_best_B",
        "stage1_best_F",
    ]

    with open(path, "w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout)
        writer.writerow(header)

        for i, fam in enumerate(stage1_summaries, start=1):
            best = fam["top_curves"][0]
            writer.writerow([
                i,
                str(fam["abs_roots"]),
                str(fam["q_params"]),
                str(fam["q"]),
                fam["best_affine_point_count"],
                fam["best_extra_x_count"],
                fam["top_k_min_affine"],
                fam["top_k_sum_affine"],
                str(best["A"]),
                str(best["B"]),
                str(best["F"]),
            ])

    print(f"Saved: {path}")


# ============================================================
# MAIN
# ============================================================

def main():
    init_output_dir()
    save_run_config()

    print("=" * 100)
    print("STARTING TWO-STAGE EVEN-AUTOMORPHISM FAMILY SCAN")
    print("=" * 100)
    print(f"OUTPUT_DIR = {OUTPUT_DIR}")
    print(f"USE_PARALLEL = {USE_PARALLEL}")
    print(f"MAX_WORKERS = {MAX_WORKERS if USE_PARALLEL else 1}")
    print()

    total_start = time.time()

    all_candidates = list(forced_family_candidates())

    print("=" * 100)
    print(f"Generated {len(all_candidates)} forced families for Stage 1")
    print("=" * 100)

    stage1_summaries = scan_candidate_list(
        candidates=all_candidates,
        stage_name="STAGE_1_SCREENING",
        pari_height=PARI_HEIGHT_STAGE1,
        A_base_pool=A_BASE_POOL_STAGE1,
        B_param_pool=B_PARAM_POOL_STAGE1,
        retain_top=TOP_FAMILIES_TO_REFINE,
    )

    save_family_summaries(
        stage1_summaries,
        "stage1_screening_details.txt",
        "STAGE 1 SCREENING RESULTS",
    )
    save_family_csv(stage1_summaries, "stage1_screening_summary.csv")
    save_selected_candidates(stage1_summaries)

    stage2_candidates = [candidate_from_summary(s) for s in stage1_summaries]

    stage2_summaries = scan_candidate_list(
        candidates=stage2_candidates,
        stage_name="STAGE_2_REFINEMENT",
        pari_height=PARI_HEIGHT_STAGE2,
        A_base_pool=A_BASE_POOL_STAGE2,
        B_param_pool=B_PARAM_POOL_STAGE2,
        retain_top=TOP_FINAL_FAMILIES,
    )

    save_family_summaries(
        stage2_summaries,
        "stage2_refined_details.txt",
        "STAGE 2 REFINED RESULTS",
    )
    save_family_csv(stage2_summaries, "stage2_refined_summary.csv")

    elapsed = time.time() - total_start

    print("=" * 100)
    print("TWO-STAGE SCAN FINISHED")
    print(f"Total elapsed time: {elapsed:.1f}s")
    print(f"OUTPUT_DIR = {OUTPUT_DIR}")
    print("=" * 100)


if __name__ == "__main__":
    main()
