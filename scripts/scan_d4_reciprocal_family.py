#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
scan_d4_reciprocal_family.py

Scanner for the D4-symmetric (reciprocal) genus 2 hyperelliptic family:

    C_{s,b}: y^2 = s*x^6 + b*x^4 + b*x^2 + s

with s = a^2 (a rational square), so the two points at infinity are rational.
The parameter t = b/s = b/a^2 determines the isomorphism class.

Rational automorphism group: Aut_Q(C_t) ≅ D4 for t ∉ {0, ±1, ±3, ±5, ±15}.

Since t determines the curve up to isomorphism, each valid (a,b) pair with
distinct t gives a distinct curve. We search over integer a and b to produce
integral models directly.

Run with:
    conda activate sage
    python scan_d4_reciprocal_family.py
"""

from sage.all import *
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import time

# ============================================================
# CONFIGURATION
# ============================================================

A_VALUES = list(range(1, 16))          # a = 1..15, so s = 1,4,9,...,225
B_RANGE = list(range(-200, 201))       # b = -200..200

EXCLUDED_T = frozenset([
    QQ(0), QQ(1), QQ(-1), QQ(3), QQ(-3), QQ(5), QQ(-5), QQ(15), QQ(-15)
])

PARI_HEIGHT = 200000
TOP_CURVES = 250
OUTPUT_DIR = Path("families_construct_d4_reciprocal")
NUM_WORKERS = 8

# ============================================================
# HELPERS
# ============================================================

def generate_ab_pairs():
    """Generate valid (a, b) pairs, excluding those where t = b/a^2 is forbidden."""
    seen_t = set()
    pairs = []
    for a in A_VALUES:
        s = a * a
        for b in B_RANGE:
            t = QQ(b) / QQ(s)
            if t in EXCLUDED_T:
                continue
            if t in seen_t:
                continue
            seen_t.add(t)
            pairs.append((a, b))
    return pairs


def d4_orbit(x0, y0):
    """Compute the D4 orbit of a rational affine point (x0, y0)."""
    points = set()
    if x0 == 0:
        points.add((x0, y0))
        points.add((x0, -y0))
        return points
    inv_x = QQ(1) / x0
    y_inv = y0 / x0**3
    for xv, yv in [(x0, y0), (-x0, y0), (inv_x, y_inv), (-inv_x, y_inv)]:
        points.add((xv, yv))
        points.add((xv, -yv))
    return points


def count_orbits(points_list):
    """Count distinct D4 orbits."""
    remaining = set(points_list)
    n_orbits = 0
    while remaining:
        p = next(iter(remaining))
        orbit = d4_orbit(p[0], p[1])
        remaining -= orbit
        n_orbits += 1
    return n_orbits


def process_one(ab_pair):
    """Process a single (a, b) pair."""
    a, b = ab_pair
    s = a * a

    Qx = PolynomialRing(QQ, 'x')
    x = Qx.gen()
    F = s * x**6 + b * x**4 + b * x**2 + s

    if not F.is_squarefree():
        return None

    try:
        pts_raw = pari(F).hyperellratpoints(PARI_HEIGHT)
    except Exception:
        return None

    affine_points = []
    for pt in pts_raw:
        affine_points.append((QQ(pt[0]), QQ(pt[1])))

    affine_point_count = len(affine_points)
    projective_lower_bound = affine_point_count + 2
    n_orbits = count_orbits(affine_points) if affine_points else 0

    coeffs = [abs(c) for c in F.list()]
    height_F = int(max(coeffs)) if coeffs else 0

    t = QQ(b) / QQ(s)

    return {
        "a": a,
        "b": b,
        "s": s,
        "t": str(t),
        "F_integral": str(F),
        "height_F": height_F,
        "affine_point_count": affine_point_count,
        "projective_lower_bound": projective_lower_bound,
        "n_orbits": n_orbits,
        "pari_height": PARI_HEIGHT,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pairs = generate_ab_pairs()

    print("=" * 80)
    print("D4 RECIPROCAL FAMILY SCANNER")
    print(f"  a range: 1..{max(A_VALUES)}")
    print(f"  b range: {min(B_RANGE)}..{max(B_RANGE)}")
    print(f"  Total (a,b) pairs (distinct t): {len(pairs)}")
    print(f"  PARI height: {PARI_HEIGHT}")
    print(f"  Workers: {NUM_WORKERS}")
    print("=" * 80)

    start = time.time()
    results = []

    with ProcessPoolExecutor(max_workers=NUM_WORKERS) as executor:
        futures = {executor.submit(process_one, p): p for p in pairs}
        done_count = 0
        for future in as_completed(futures):
            done_count += 1
            if done_count % 200 == 0:
                print(f"  [{done_count}/{len(pairs)}] processed...")
            result = future.result()
            if result is not None and result["affine_point_count"] > 0:
                results.append(result)

    elapsed = time.time() - start
    print(f"\nScan complete: {len(results)} curves with points in {elapsed:.1f}s")

    # Sort: maximize affine, then orbits, then minimize height
    results.sort(key=lambda r: (-r["affine_point_count"], -r["n_orbits"], r["height_F"]))
    top_results = results[:TOP_CURVES]

    for i, r in enumerate(top_results, start=1):
        r["curve_index"] = i
        r["curve_label"] = f"CURBA #{i}"
        r["family_label"] = f"D4_t={r['t']}"
        r["construction_kind"] = "d4_reciprocal"

    fieldnames = [
        "construction_kind", "curve_index", "curve_label", "family_label",
        "a", "b", "s", "t",
        "F_integral", "height_F",
        "affine_point_count", "projective_lower_bound", "n_orbits",
        "pari_height",
    ]

    csv_path = OUTPUT_DIR / "d4_top_curves_summary.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in top_results:
            writer.writerow({k: r.get(k, "") for k in fieldnames})

    print(f"\nTop {len(top_results)} curves written to: {csv_path}")

    print(f"\nTop 10:")
    for r in top_results[:10]:
        print(f"  t={r['t']:>10s}  a={r['a']:>2d}  b={r['b']:>4d}  "
              f"affine={r['affine_point_count']:>4d}  orbits={r['n_orbits']:>3d}  "
              f"proj={r['projective_lower_bound']:>4d}  height={r['height_F']}")

    if top_results:
        print(f"\n  Best F(x) = {top_results[0]['F_integral']}")

    print(f"\nElapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
