#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

"""
Two-stage parallel scanner for the GENERAL family of genus-2 hyperelliptic curves
with NON-MONIC Q(x) (leading coefficient A = rational square).

    C : y^2 = F(x),
    F(x) = P_q(x) + H(x) Q(x),

where:
    H(x)   = prod_i (x - r_i),
    P_q(x) = q(x)^2 mod H(x),
    Q(x)   = A*x^(6-m) + ...   with A = a^2 for some a in Q^*,

so every member of the family passes through the forced rational points:
    (r_i, ±q(r_i))

and the leading coefficient of F is A * leading_coeff(H), which is a rational
square (since A is a square and leading_coeff(H) = 1 for monic H). This ensures
the two points at infinity are rational.

This is the general case: no rational automorphism is imposed.

Main improvement over the previous version:
    - the generator is balanced using MAX_Q_PER_RS;
    - it avoids spending almost all MAX_FORCED_FAMILIES on the first few r-tuples;
    - coefficient pools for Q/L may be rational, not only integral.

STAGE 1:
    Broad screening with smaller PARI height and moderate Q/L coefficient pool.

STAGE 2:
    Refine only the best TOP_FAMILIES_TO_REFINE families using larger PARI height
    and, optionally, a larger rational Q/L coefficient pool.

Run with:
    sage -python scan_general_common_forced_families_two_stage_parallel_balanced.py

Recommended:
    Run as a .py file from WSL/Linux, not from Jupyter, because multiprocessing
    is more reliable when the __main__ module is importable.
"""

from sage.all import *
from itertools import combinations, product
from itertools import islice
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import csv
import os
import time


# ============================================================
# SAGE SETUP
# ============================================================

R = PolynomialRing(QQ, "x")
x = R.gen()

OUTPUT_DIR = None


# ============================================================
# POOL HELPERS
# ============================================================

def rational_pool_signed(num_bound, den_bound, include_zero=True):
    """
    Produces reduced rationals a/b with:
        |a| <= num_bound,
        1 <= b <= den_bound.

    Example:
        rational_pool_signed(3, 2) includes
        -3, -5/2, -2, -3/2, -1, -1/2, 0, 1/2, 1, ...
    """
    vals = set()

    for b in range(1, den_bound + 1):
        for a in range(-num_bound, num_bound + 1):
            if gcd(ZZ(a), ZZ(b)) == 1:
                q = QQ(a) / QQ(b)
                if include_zero or q != 0:
                    vals.add(q)

    return sorted(vals)


def rational_pool_unsigned(num_bound, den_bound, include_zero=False):
    """
    Produces nonnegative reduced rationals a/b with:
        0 <= a <= num_bound,
        1 <= b <= den_bound.
    """
    vals = set()

    for b in range(1, den_bound + 1):
        for a in range(0, num_bound + 1):
            if gcd(ZZ(a), ZZ(b)) == 1:
                q = QQ(a) / QQ(b)
                if include_zero or q != 0:
                    vals.add(q)

    return sorted(vals)


# ============================================================
# GLOBAL FAMILY CONFIGURATION
# ============================================================

# m = number of forced abscissae, so deg H = m.
# m=4: Q(x)=A*x^2+B*x+C  (A is a rational square for rational points at infinity)
# m=5: Q(x)=A*x+B         (A is a rational square)
# m=6: Q(x)=A             (A is a rational square)
M_VALUES = [4]

# Leading coefficient pool for Q(x): A = a^2 where a is a positive rational.
# 10 values covering small and moderate rationals.
A_BASE_POOL = [QQ(1)/QQ(3), QQ(1)/QQ(2), QQ(2)/QQ(3), QQ(1), QQ(3)/QQ(2), QQ(2), QQ(5)/QQ(2), QQ(3), QQ(4), QQ(5)]

# Pool for rational abscissae r_i = a/b, with |a| <= num_bound, 1 <= b <= den_bound.
RATIONAL_RS_NUM_BOUND = 8
RATIONAL_RS_DEN_BOUND = 8

# Number of rs-tuples to generate and sort for each m.
MAX_RS_TUPLES_PER_M = 50000
MAX_RS_POOL_SIZE = 30

# q(x)=u*x^2+v*x+w, with u,v,w from this pool.
# This may contain integers or Sage rationals.
Q_COEFF_POOL = list(range(-10, 11))

# Exclude q = 0.
ALLOW_ZERO_Q = False

# Avoid cases where q(r_i)=0 too often.
MAX_ZERO_FORCED_ABSCISSAE = 1

# Balanced generation:
#   for each fixed rs-tuple, test at most MAX_Q_PER_RS polynomials q.
# This prevents the Stage 1 budget from being spent almost entirely on the first few rs-tuples.
MAX_Q_PER_RS = 150

# Optional cap on the total number of forced families tested in Stage 1, across all m-values.
MAX_FORCED_FAMILIES = 2000

# Deduplicate q and variants producing the same P_q modulo H.
DEDUP_FAMILIES_BY_H_AND_PQ = True

# Deduplicate identical curves inside the same family.
DEDUP_CURVES_WITHIN_FAMILY = True


# ============================================================
# STAGE 1: BROAD SCREENING
# ============================================================

# Stage 1 should be relatively cheap.
PARI_HEIGHT_STAGE1 = 5000

# Tail coefficients for Q(x).
# For m=4, Q(x)=x^2+B*x+C, so two tail coefficients are drawn from this pool.
# For m=5, Q(x)=x+D, so one tail coefficient is drawn from this pool.
# For m=6, Q(x)=1, so there is no tail coefficient.
#
# Stage 1: moderate rational pool.
L_PARAM_POOL_STAGE1 = rational_pool_signed(num_bound=6, den_bound=2, include_zero=True)

# Retain this many families after Stage 1 and send them to Stage 2.
TOP_FAMILIES_TO_REFINE = 50


# ============================================================
# STAGE 2: EXPENSIVE REFINEMENT
# ============================================================

PARI_HEIGHT_STAGE2 = 200000

# Stage 2: denser rational pool, applied only to the selected families.
L_PARAM_POOL_STAGE2 = rational_pool_signed(num_bound=10, den_bound=3, include_zero=True)

# Number of final families retained.
TOP_FINAL_FAMILIES = 100


# ============================================================
# SCORING CONFIGURATION
# ============================================================

TOP_CURVES_PER_FAMILY = 5

# For m=4, affine_count is usually 8 + 2*extra_x_count.
# Thus extra_x_count=26 corresponds roughly to 60 affine points.
GOOD_EXTRA_X_THRESHOLD = 26

TARGET_AFFINE_1 = 55
TARGET_AFFINE_2 = 60


# ============================================================
# MULTIPROCESSING CONFIGURATION
# ============================================================

USE_PARALLEL = True

# Conservative default. Increase if your machine has enough CPU/RAM.
MAX_WORKERS = max(1, (os.cpu_count() or 2) - 4)
# Alternative:
# MAX_WORKERS = 4


# ============================================================
# OUTPUT CONFIGURATION
# ============================================================

OUTPUT_BASE_DIR = Path("families_construct_general_two_stage")

# If True, uses a short timestamped folder name and stores full config in run_config.txt.
# This avoids Windows MAX_PATH issues.
USE_SHORT_RUN_NAME = False


# ============================================================
# SAFE RUN NAME HELPERS
# ============================================================

def rational_to_name(q):
    q = QQ(q)
    n = q.numerator()
    d = q.denominator()
    if d == 1:
        return str(n)
    return f"{n}over{d}"


def safe_path_component(s):
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
    pool = list(pool)
    if not pool:
        return "empty"
    return safe_path_component(f"{rational_to_name(min(pool))}to{rational_to_name(max(pool))}")


def make_run_name():
    if USE_SHORT_RUN_NAME:
        run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        return safe_path_component(
            f"general_H1{PARI_HEIGHT_STAGE1}_H2{PARI_HEIGHT_STAGE2}_top{TOP_FAMILIES_TO_REFINE}_{run_id}"
        )

    run_name = (
        f"general_m{'-'.join(map(str, M_VALUES))}"
        f"_H1_{PARI_HEIGHT_STAGE1}"
        f"_H2_{PARI_HEIGHT_STAGE2}"
        f"_topRef{TOP_FAMILIES_TO_REFINE}"
        f"_goodX{GOOD_EXTRA_X_THRESHOLD}"
        f"_rsN{RATIONAL_RS_NUM_BOUND}_rsD{RATIONAL_RS_DEN_BOUND}"
        f"_maxRS{MAX_RS_TUPLES_PER_M}"
        f"_maxFam{MAX_FORCED_FAMILIES}"
        f"_maxQperRS{MAX_Q_PER_RS}"
        f"_q{pool_range_name(Q_COEFF_POOL)}"
        f"_L1{pool_range_name(L_PARAM_POOL_STAGE1)}"
        f"_nL1{len(set(QQ(a) for a in L_PARAM_POOL_STAGE1))}"
        f"_L2{pool_range_name(L_PARAM_POOL_STAGE2)}"
        f"_nL2{len(set(QQ(a) for a in L_PARAM_POOL_STAGE2))}"
        f"_nA{len(A_BASE_POOL)}"
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
        out.write("CONFIGURATION FOR TWO-STAGE GENERAL COMMON-FORCED-POINT FAMILY SCAN\n")
        out.write("=" * 100 + "\n")
        out.write(f"M_VALUES = {M_VALUES}\n")
        out.write(f"RATIONAL_RS_NUM_BOUND = {RATIONAL_RS_NUM_BOUND}\n")
        out.write(f"RATIONAL_RS_DEN_BOUND = {RATIONAL_RS_DEN_BOUND}\n")
        out.write(f"MAX_RS_TUPLES_PER_M = {MAX_RS_TUPLES_PER_M}\n")
        out.write(f"Q_COEFF_POOL = {Q_COEFF_POOL}\n")
        out.write(f"ALLOW_ZERO_Q = {ALLOW_ZERO_Q}\n")
        out.write(f"MAX_ZERO_FORCED_ABSCISSAE = {MAX_ZERO_FORCED_ABSCISSAE}\n")
        out.write(f"MAX_Q_PER_RS = {MAX_Q_PER_RS}\n")
        out.write(f"MAX_FORCED_FAMILIES = {MAX_FORCED_FAMILIES}\n")
        out.write(f"DEDUP_FAMILIES_BY_H_AND_PQ = {DEDUP_FAMILIES_BY_H_AND_PQ}\n")
        out.write(f"DEDUP_CURVES_WITHIN_FAMILY = {DEDUP_CURVES_WITHIN_FAMILY}\n")
        out.write("\n")
        out.write("STAGE 1\n")
        out.write("-" * 100 + "\n")
        out.write(f"PARI_HEIGHT_STAGE1 = {PARI_HEIGHT_STAGE1}\n")
        out.write(f"L_PARAM_POOL_STAGE1 = {L_PARAM_POOL_STAGE1}\n")
        out.write(f"len(L_PARAM_POOL_STAGE1) = {len(L_PARAM_POOL_STAGE1)}\n")
        out.write(f"TOP_FAMILIES_TO_REFINE = {TOP_FAMILIES_TO_REFINE}\n")
        out.write("\n")
        out.write("STAGE 2\n")
        out.write("-" * 100 + "\n")
        out.write(f"PARI_HEIGHT_STAGE2 = {PARI_HEIGHT_STAGE2}\n")
        out.write(f"L_PARAM_POOL_STAGE2 = {L_PARAM_POOL_STAGE2}\n")
        out.write(f"len(L_PARAM_POOL_STAGE2) = {len(L_PARAM_POOL_STAGE2)}\n")
        out.write(f"TOP_FINAL_FAMILIES = {TOP_FINAL_FAMILIES}\n")
        out.write("\n")
        out.write("SCORING\n")
        out.write("-" * 100 + "\n")
        out.write(f"TOP_CURVES_PER_FAMILY = {TOP_CURVES_PER_FAMILY}\n")
        out.write(f"GOOD_EXTRA_X_THRESHOLD = {GOOD_EXTRA_X_THRESHOLD}\n")
        out.write(f"TARGET_AFFINE_1 = {TARGET_AFFINE_1}\n")
        out.write(f"TARGET_AFFINE_2 = {TARGET_AFFINE_2}\n")
        out.write("\n")
        out.write("MULTIPROCESSING\n")
        out.write("-" * 100 + "\n")
        out.write(f"USE_PARALLEL = {USE_PARALLEL}\n")
        out.write(f"MAX_WORKERS = {MAX_WORKERS}\n")
        out.write(f"os.cpu_count() = {os.cpu_count()}\n")
        out.write("\n")
        out.write("OUTPUT\n")
        out.write("-" * 100 + "\n")
        out.write(f"OUTPUT_BASE_DIR = {OUTPUT_BASE_DIR}\n")
        out.write(f"OUTPUT_DIR = {OUTPUT_DIR}\n")
        out.write(f"USE_SHORT_RUN_NAME = {USE_SHORT_RUN_NAME}\n")


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


def rational_pool(num_bound, den_bound, include_zero=True):
    vals = set()
    for b in range(1, den_bound + 1):
        for a in range(-num_bound, num_bound + 1):
            if gcd(ZZ(a), ZZ(b)) == 1:
                q = QQ(a) / QQ(b)
                if include_zero or q != 0:
                    vals.add(q)
    return sorted(vals)


def forced_points(q, rs):
    pts = set()
    for r in rs:
        y = QQ(q(r))
        pts.add((QQ(r), y))
        pts.add((QQ(r), -y))
    return sorted(pts, key=lambda P: (P[0], P[1]))


# ============================================================
# GENERAL FAMILY CONSTRUCTION
# ============================================================

def make_poly_from_tail(tail_coeffs, deg, leading_coeff=QQ(1)):
    """
    Build polynomial of given degree with specified leading coefficient.
    deg=2, [B,C], A -> A*x^2 + B*x + C
    deg=1, [D], A   -> A*x + D
    deg=0, [], A    -> A
    """
    if len(tail_coeffs) != deg:
        raise ValueError(f"deg={deg}, but tail_coeffs has length {len(tail_coeffs)}")

    if deg == 0:
        return R(leading_coeff)

    L = QQ(leading_coeff) * x**deg
    for i, c in enumerate(tail_coeffs):
        L += QQ(c) * x**(deg - 1 - i)
    return R(L)


def build_family_general(rs, q, L_tail_coeffs, leading_A=QQ(1)):
    """
    For fixed rs, q and Q of degree 6-len(rs), construct:
        H  = prod (x-r_i)
        Pq = q^2 mod H
        Q  = A*x^deg + tail  (A is the leading coefficient, a rational square)
        F  = Pq + H*Q
    """
    rs = tuple(QQ(r) for r in rs)
    if len(set(rs)) != len(rs):
        raise ValueError("The r_i must be distinct.")

    m = len(rs)
    if m not in [4, 5, 6]:
        raise ValueError("This scanner is designed for m=4,5,6.")

    H = R(prod(x - r for r in rs))
    q = R(q)
    Pq = R((q**2) % H)

    degL = 6 - m
    L = make_poly_from_tail(L_tail_coeffs, degL, leading_coeff=leading_A)
    F = R(Pq + H * L)

    return {
        "m": m,
        "rs": rs,
        "q": q,
        "H": H,
        "Pq": Pq,
        "L": L,
        "L_tail": tuple(QQ(c) for c in L_tail_coeffs),
        "A": QQ(leading_A),
        "F": F,
    }


def L_tail_grid_for_m(m, L_param_pool):
    """
    Generate all (A, tail_coeffs) tuples for the given m.
    A is drawn from A_BASE_POOL (rational squares).
    tail_coeffs are drawn from L_param_pool.
    """
    degL = 6 - m
    pool = [QQ(a) for a in L_param_pool]
    A_pool = [a**2 for a in A_BASE_POOL]

    if degL == 2:
        return [(A, (B, C)) for A in A_pool for B, C in product(pool, pool)]
    if degL == 1:
        return [(A, (D,)) for A in A_pool for D in pool]
    if degL == 0:
        return [(A, ()) for A in A_pool]

    raise ValueError("Invalid m.")


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


def split_forced_and_extra(points, rs):
    rs_set = set(QQ(r) for r in rs)
    forced = [P for P in points if P[0] in rs_set]
    extra = [P for P in points if P[0] not in rs_set]
    return forced, extra


# ============================================================
# CURVE AND FAMILY SCORING
# ============================================================

def curve_record(data, pari_height):
    F = data["F"]
    rs = data["rs"]
    q = data["q"]

    if not is_squarefree_sextic(F):
        return None

    # In the monic construction, the leading coefficient is 1, hence a square.
    # This check is kept for safety if the user later modifies the leading coefficient.
    if not leading_coeff_is_square_Q(F):
        return None

    try:
        pts = search_with_pari_on_poly(F, pari_height)
    except Exception:
        return None

    forced_found, extra = split_forced_and_extra(pts, rs)
    extra_xs = sorted(set(P[0] for P in extra))

    rec = dict(data)
    rec.update({
        "height_F": poly_height(F),
        "forced_expected": forced_points(q, rs),
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
    """
    All curves here are sextics with square leading coefficient, so the projective
    lower bound is affine_point_count + 2. We sort mainly by affine count.
    """
    return (
        -rec["affine_point_count"],
        rec["height_F"],
        -rec["extra_x_count"],
    )


def summarize_family(m, rs, q, q_params, q_coeffs, curve_records, stage_name):
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
        "m": m,
        "rs": tuple(QQ(r) for r in rs),
        "q": R(q),
        "q_params": dict(q_params),
        "q_coeffs": tuple(QQ(c) for c in q_coeffs),
        "forced_points": forced_points(q, rs),

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
    """
    Sort families by the quality of their best TOP_CURVES_PER_FAMILY curves.
    """
    return (
        -summary["top_k_affine_ge_target2"],
        -summary["top_k_affine_ge_target1"],
        -summary["top_k_min_affine"],
        -summary["top_k_sum_affine"],
        -summary["best_affine_point_count"],
        -summary["best_projective_lower_bound"],
        summary["best_height_F"],
    )


def best_member_sort_key(summary):
    """
    Optional sort key: use it if you only care about the best single curve in a family.
    Not used by default.
    """
    return (
        -summary["best_affine_point_count"],
        -summary["best_projective_lower_bound"],
        -summary["top_k_sum_affine"],
        summary["best_height_F"],
    )


# ============================================================
# FORCED FAMILY GENERATION
# ============================================================

def generate_rs_tuples(m):
    """
    Memory-safe version for generating rs-tuples.

    Important for m=6, because list(combinations(pool, 6)) can be enormous.
    """

    pool = rational_pool(
        RATIONAL_RS_NUM_BOUND,
        RATIONAL_RS_DEN_BOUND,
        include_zero=True,
    )

    # Sortăm valorile individuale după înălțime și mărime.
    pool = sorted(pool, key=lambda r: (
        rational_height(r),
        abs(QQ(r)),
        QQ(r),
    ))

    # Limităm pool-ul înainte de a face combinații.
    # Pentru m=6, valori recomandate: 20-30.
    if MAX_RS_POOL_SIZE is not None and len(pool) > MAX_RS_POOL_SIZE:
        pool = pool[:MAX_RS_POOL_SIZE]

    tuples = list(combinations(pool, m))

    # Acum lista este controlată, deci putem sorta combinațiile.
    tuples.sort(key=lambda rs: (
        max(rational_height(r) for r in rs),
        sum(rational_height(r) for r in rs),
        sum(abs(QQ(r)) for r in rs),
    ))

    if MAX_RS_TUPLES_PER_M is not None and len(tuples) > MAX_RS_TUPLES_PER_M:
        tuples = tuples[:MAX_RS_TUPLES_PER_M]

    return tuples


def generate_q_polys():
    """
    q(x)=u*x^2+v*x+w.
    Returns (q, q_params, q_coeffs).
    """
    qs = []

    for u, v, w in product(Q_COEFF_POOL, Q_COEFF_POOL, Q_COEFF_POOL):
        if not ALLOW_ZERO_Q and QQ(u) == 0 and QQ(v) == 0 and QQ(w) == 0:
            continue

        q = R(QQ(u)*x**2 + QQ(v)*x + QQ(w))
        q_params = {"u": QQ(u), "v": QQ(v), "w": QQ(w)}
        q_coeffs = (QQ(u), QQ(v), QQ(w))
        qs.append((q, q_params, q_coeffs))

    qs.sort(key=lambda item: (
        poly_height(item[0]),
        str(item[0]),
    ))

    return qs


def forced_family_candidates():
    """
    Generates tuples:
        (m, rs, q, q_params, q_coeffs).

    Balanced version:
      - avoids spending almost all MAX_FORCED_FAMILIES on the first few rs-tuples;
      - tests at most MAX_Q_PER_RS q-polynomials for each rs-tuple.
    """
    q_items = generate_q_polys()
    count = 0
    seen_family_keys = set()

    for m in M_VALUES:
        rs_tuples = generate_rs_tuples(m)

        for rs in rs_tuples:
            q_count_for_this_rs = 0

            for q, q_params, q_coeffs in q_items:
                zero_y_count = sum(1 for r in rs if q(QQ(r)) == 0)
                if zero_y_count > MAX_ZERO_FORCED_ABSCISSAE:
                    continue

                H = R(prod(x - QQ(r) for r in rs))
                Pq = R((R(q)**2) % H)

                fam_key = (
                    m,
                    tuple(QQ(r) for r in rs),
                    polynomial_key(H),
                    polynomial_key(Pq),
                )

                if DEDUP_FAMILIES_BY_H_AND_PQ and fam_key in seen_family_keys:
                    continue

                seen_family_keys.add(fam_key)

                yield m, rs, q, q_params, q_coeffs

                count += 1
                q_count_for_this_rs += 1

                if q_count_for_this_rs >= MAX_Q_PER_RS:
                    break

                if MAX_FORCED_FAMILIES is not None and count >= MAX_FORCED_FAMILIES:
                    return


def candidate_from_summary(summary):
    """
    Reconstructs a candidate tuple from a summary.
    """
    m = int(summary["m"])
    rs = tuple(QQ(r) for r in summary["rs"])
    q_coeffs = tuple(QQ(c) for c in summary["q_coeffs"])
    u, v, w = q_coeffs
    q = R(QQ(u)*x**2 + QQ(v)*x + QQ(w))
    q_params = {"u": QQ(u), "v": QQ(v), "w": QQ(w)}
    return m, rs, q, q_params, q_coeffs


# ============================================================
# SCANNING
# ============================================================

def scan_one_general_family(m, rs, q, q_params, q_coeffs,
                            pari_height, L_param_pool, stage_name):
    curve_records = []
    seen_F = set()

    for leading_A, L_tail in L_tail_grid_for_m(m, L_param_pool):
        data = build_family_general(rs=rs, q=q, L_tail_coeffs=L_tail, leading_A=leading_A)

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
        m=m,
        rs=rs,
        q=q,
        q_params=q_params,
        q_coeffs=q_coeffs,
        curve_records=curve_records,
        stage_name=stage_name,
    )


def scan_one_general_family_worker(args):
    """
    Multiprocessing worker.
    """
    (
        m,
        rs_raw,
        q_coeffs_raw,
        q_params_raw,
        pari_height,
        L_param_raw,
        stage_name,
    ) = args

    rs = tuple(QQ(r) for r in rs_raw)
    q_coeffs = tuple(QQ(c) for c in q_coeffs_raw)
    u, v, w = q_coeffs
    q = R(QQ(u)*x**2 + QQ(v)*x + QQ(w))
    q_params = {"u": QQ(u), "v": QQ(v), "w": QQ(w)}
    L_param_pool = [QQ(a) for a in L_param_raw]

    try:
        return scan_one_general_family(
            m=m,
            rs=rs,
            q=q,
            q_params=q_params,
            q_coeffs=q_coeffs,
            pari_height=pari_height,
            L_param_pool=L_param_pool,
            stage_name=stage_name,
        )
    except Exception:
        return None


def scan_candidate_list(candidates, stage_name, pari_height, L_param_pool, retain_top=None):
    """
    Scans a given list of forced families with a given PARI height and Q/L coefficient grid.
    """
    print("=" * 100)
    print(f"{stage_name}: scanning {len(candidates)} forced families")
    print(f"PARI height = {pari_height}")
    print(f"len(L parameter pool) = {len(L_param_pool)}")
    print(f"For m=4, this means {len(L_param_pool)**2} curves per family before filtering.")
    print("=" * 100)

    start = time.time()
    summaries = []

    if not USE_PARALLEL or MAX_WORKERS <= 1:
        completed = 0

        for m, rs, q, q_params, q_coeffs in candidates:
            completed += 1

            summary = scan_one_general_family(
                m=m,
                rs=rs,
                q=q,
                q_params=q_params,
                q_coeffs=q_coeffs,
                pari_height=pari_height,
                L_param_pool=L_param_pool,
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
        L_param_raw = [str(QQ(a)) for a in L_param_pool]

        for m, rs, q, q_params, q_coeffs in candidates:
            rs_raw = tuple(str(QQ(r)) for r in rs)
            q_coeffs_raw = tuple(str(QQ(c)) for c in q_coeffs)
            q_params_raw = {k: str(QQ(v)) for k, v in q_params.items()}

            tasks.append((
                int(m),
                rs_raw,
                q_coeffs_raw,
                q_params_raw,
                pari_height,
                L_param_raw,
                stage_name,
            ))

        completed = 0

        with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
            future_to_task = {
                executor.submit(scan_one_general_family_worker, task): task
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
    out.write(f"L_tail = {rec['L_tail']}\n")
    out.write(f"L(x) = {rec['L']}\n")
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
        out.write("GENERAL FAMILIES OF GENUS-2 HYPERELLIPTIC CURVES\n")
        out.write("Construction: F(x)=P_q(x)+H(x)Q(x), with no rational automorphism imposed.\n")
        out.write("=" * 100 + "\n\n")

        for i, fam in enumerate(summaries, start=1):
            best = fam["top_curves"][0]

            out.write("=" * 100 + "\n")
            out.write(f"FAMILIA #{i}\n")
            out.write(f"stage = {fam['stage']}\n")
            out.write(f"m = {fam['m']}\n")
            out.write(f"rs = {fam['rs']}\n")
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
        "m",
        "rs",
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
        "best_L_tail",
        "best_A",
        "best_L",
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
                fam["m"],
                str(fam["rs"]),
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
                str(best["L_tail"]),
                str(best.get("A", 1)),
                str(best["L"]),
                str(best["F"]),
            ]

            writer.writerow(row)

    print(f"Saved: {path}")


def save_selected_candidates(stage1_summaries):
    path = Path(OUTPUT_DIR) / "selected_families_for_stage2.csv"

    header = [
        "selected_index",
        "m",
        "rs",
        "q_params",
        "q",
        "stage1_best_affine_point_count",
        "stage1_best_extra_x_count",
        "stage1_top_k_min_affine",
        "stage1_top_k_sum_affine",
        "stage1_best_L_tail",
        "stage1_best_F",
    ]

    with open(path, "w", newline="", encoding="utf-8") as fout:
        writer = csv.writer(fout)
        writer.writerow(header)

        for i, fam in enumerate(stage1_summaries, start=1):
            best = fam["top_curves"][0]
            writer.writerow([
                i,
                fam["m"],
                str(fam["rs"]),
                str(fam["q_params"]),
                str(fam["q"]),
                fam["best_affine_point_count"],
                fam["best_extra_x_count"],
                fam["top_k_min_affine"],
                fam["top_k_sum_affine"],
                str(best["L_tail"]),
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
    print("STARTING TWO-STAGE GENERAL COMMON-FORCED-POINT FAMILY SCAN")
    print("=" * 100)
    print(f"OUTPUT_DIR = {OUTPUT_DIR}")
    print(f"USE_PARALLEL = {USE_PARALLEL}")
    print(f"MAX_WORKERS = {MAX_WORKERS if USE_PARALLEL else 1}")
    print(f"MAX_Q_PER_RS = {MAX_Q_PER_RS}")
    print(f"len(L_PARAM_POOL_STAGE1) = {len(L_PARAM_POOL_STAGE1)}")
    print(f"len(L_PARAM_POOL_STAGE2) = {len(L_PARAM_POOL_STAGE2)}")
    print()

    total_start = time.time()

    # -------------------------------
    # Stage 1: broad cheap screening
    # -------------------------------
    all_candidates = list(forced_family_candidates())

    print("=" * 100)
    print(f"Generated {len(all_candidates)} forced families for Stage 1")
    print("=" * 100)

    stage1_summaries = scan_candidate_list(
        candidates=all_candidates,
        stage_name="STAGE_1_SCREENING",
        pari_height=PARI_HEIGHT_STAGE1,
        L_param_pool=L_PARAM_POOL_STAGE1,
        retain_top=TOP_FAMILIES_TO_REFINE,
    )

    save_family_summaries(
        stage1_summaries,
        "stage1_screening_details.txt",
        "STAGE 1 SCREENING RESULTS",
    )
    save_family_csv(stage1_summaries, "stage1_screening_summary.csv")
    save_selected_candidates(stage1_summaries)

    # -------------------------------
    # Stage 2: expensive refinement
    # -------------------------------
    stage2_candidates = [candidate_from_summary(s) for s in stage1_summaries]

    stage2_summaries = scan_candidate_list(
        candidates=stage2_candidates,
        stage_name="STAGE_2_REFINEMENT",
        pari_height=PARI_HEIGHT_STAGE2,
        L_param_pool=L_PARAM_POOL_STAGE2,
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
