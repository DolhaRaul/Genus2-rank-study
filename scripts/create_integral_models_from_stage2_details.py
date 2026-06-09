#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-

"""
Create integral isomorphic models for the top curves saved in stage2_refined_details.txt.

Input roots expected by default:

    families_construct_even_automorphism_two_stage/
    families_construct_general_two_stage/

Each execution folder is expected to contain:

    run_config.txt
    selected_families_for_stage2.csv
    stage1_screening_details.txt
    stage1_screening_summary.csv
    stage2_refined_details.txt
    stage2_refined_summary.csv

The script reads stage2_refined_details.txt, extracts every curve F(x) appearing under
the CURBA #k blocks, and writes an integral isomorphic model.

Mathematical normalization used:

    C  : y^2 = F(x),     F in QQ[x]
    D  : lcm of denominators of all coefficients of F
    X* = X
    Y* = D Y

Then:

    (Y*)^2 = D^2 F(X*) = F_integral(X*)

and F_integral lies in ZZ[x].

Important:
    The coordinate scaling is Y* = D*Y, not D^2*Y.
    The RHS is multiplied by D^2.

Run with:

    sage -python create_integral_models_from_stage2_details.py
"""

from sage.all import QQ, ZZ, PolynomialRing, lcm
from sage.misc.sage_eval import sage_eval

from pathlib import Path
import csv
import re
import shutil
from datetime import datetime


# ============================================================
# CONFIGURATION
# ============================================================

SOURCE_ROOTS = [
    {
        "kind": "even_automorphism",
        "input_root": Path("families_construct_even_automorphism_two_stage"),
        "output_root": Path("families_construct_even_automorphism_two_stage_integral_models"),
    },
    {
        "kind": "general",
        "input_root": Path("families_construct_general_two_stage"),
        "output_root": Path("families_construct_general_two_stage_integral_models"),
    },
]

DETAILS_FILENAME = "stage2_refined_details.txt"

# If True, search recursively for stage2_refined_details.txt.
RECURSIVE = True

# If False, output execution folders mirror the original execution folder names.
# If True, output execution folders are run_0001, run_0002, ...
# and an execution_index.csv is written in the output root.
USE_SHORT_EXECUTION_FOLDER_NAMES = False

# Copy run_config.txt into the output execution folder for traceability.
COPY_RUN_CONFIG = True

# Keep only essential metadata in the details file; point lists are very large.
WRITE_LONG_POINT_LISTS = False

# Also write a Magma-friendly data file with integral polynomials.
WRITE_MAGMA_DATA_FILE = True

OUTPUT_DETAILS_FILENAME = "stage2_integral_models_details.txt"
OUTPUT_SUMMARY_FILENAME = "stage2_integral_models_summary.csv"
OUTPUT_MAGMA_DATA_FILENAME = "stage2_integral_models_magma_data.txt"
AGGREGATE_SUMMARY_FILENAME = "all_stage2_integral_models_summary.csv"
EXECUTION_INDEX_FILENAME = "execution_index.csv"


# ============================================================
# SAGE SETUP
# ============================================================

R = PolynomialRing(QQ, "x")
x = R.gen()

RZ = PolynomialRing(ZZ, "x")
xz = RZ.gen()


# ============================================================
# REGEX / PARSING HELPERS
# ============================================================

FAMILY_START_RE = re.compile(r"(?m)^\s*FAMILIA\s+#(\d+)\s*$")
CURVE_START_RE = re.compile(r"(?m)^\s*CURBA\s+#(\d+)\s*$")

def parse_sage_poly(expr):
    """
    Parse a Sage polynomial expression in x into QQ[x].
    """
    expr = expr.strip()
    return R(sage_eval(expr, locals={"x": x}))

# AICI regex este un obiect de tipu; Pattern ( obtinut prin re.compile(regex))
# IN loc sa refolosim acelasi regex de mai multe ori, CREEM un obiect Patterm dedicat
# Pattern.finditer() DE ne un iterator ce face yield pe obiecte de tipul Match ( OBIECTELE sunt create
# pe baza de continuturi ce dau match SI ACESTEA nu sunt overlapping!)
def split_blocks_by_regex(text, regex):
    """
    Split text into blocks that start at each regex match.
    Returns a list of (label_number, block_text).
    """
    matches = list(regex.finditer(text))
    blocks = []

    # NU putem folosi end = m.end() pentru a defini blocul, PENTRU ca asta ne ar da DOAR sfarsitul liniei "FAMILIA #n"
    # SFARSITUL blocului ESTE fix inainte de noua linie de tipul " FAMILIA #n+1", DECI e corect sa ne uitam la inceput
    # de nou match
    for i, m in enumerate(matches):
        start = m.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        number = int(m.group(1))
        block = text[start:end].strip()
        blocks.append((number, block))

    return blocks


def parse_key_value_lines(text):
    """
    Parse lines of the exact form:

        key = value

    using the delimiter ' = '. This is important because some keys contain '>=':
        good_curve_count(extra_x >= 26) = 9
    """
    fields = {}

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if " = " not in line:
            continue

        key, _, value = line.partition(" = ")
        key = key.strip()
        value = value.strip()

        if key:
            fields[key] = value

    return fields


def text_before_top_curves(family_block):
    """
    PENTRU fiecare familie ( bloc de familie cu ale ei curbe), SE returneaza FIX acea parte care caracterizeaza
    general familia ( DECI acea parte de dinainte de definitia primei curbe)
    """
    marker = "TOP CURVES IN THIS FAMILY"
    if marker in family_block:
        return family_block.split(marker, 1)[0]
    return family_block


def extract_family_records(details_text):
    """
    Extract family blocks and the curve blocks inside each family.

    Returns:
        [
          {
            "family_index": int,
            "family_label": "FAMILIA #i",
            "family_fields": {...},
            "curves": [
                {
                  "curve_index": int,
                  "curve_label": "CURBA #j",
                  "curve_fields": {...}
                },
                ...
            ]
          },
          ...
        ]
    """
    family_records = []

    for fam_no, fam_block in split_blocks_by_regex(details_text, FAMILY_START_RE):
        family_label = f"FAMILIA #{fam_no}"
        family_header_text = text_before_top_curves(fam_block)
        family_fields = parse_key_value_lines(family_header_text)

        curves = []
        curve_blocks = split_blocks_by_regex(fam_block, CURVE_START_RE)

        for curve_no, curve_block in curve_blocks:
            curve_label = f"CURBA #{curve_no}"
            curve_fields = parse_key_value_lines(curve_block)

            if "F(x)" not in curve_fields:
                # Robustness: skip malformed curve block.
                continue

            curves.append({
                "curve_index": curve_no,
                "curve_label": curve_label,
                "curve_fields": curve_fields,
            })

        family_records.append({
            "family_index": fam_no,
            "family_label": family_label,
            "family_fields": family_fields,
            "curves": curves,
        })

    return family_records


# ============================================================
# POLYNOMIAL NORMALIZATION
# ============================================================

def coeffs_from_poly_QQ(F):
    F = R(F)
    coeffs = list(F.list())
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return [QQ(c) for c in coeffs]


def polynomial_height_QQ(F):
    """
    Naive rational height: max over coefficients of max(abs(numerator), denominator).
    """
    coeffs = coeffs_from_poly_QQ(F)
    if not coeffs:
        return ZZ(0)
    return max(max(abs(c.numerator()), c.denominator()) for c in coeffs)


def polynomial_height_ZZ(FZ):
    coeffs = list(RZ(FZ).list())
    if not coeffs:
        return ZZ(0)
    return max(abs(ZZ(c)) for c in coeffs)


def is_squarefree_over_QQ(F):
    F = R(F)
    if F.degree() <= 0:
        return False
    return F.gcd(F.derivative()).degree() == 0


def is_even_polynomial(F):
    F = R(F)
    return all(F.monomial_coefficient(x**k) == 0 for k in range(1, F.degree() + 1, 2))


def integral_model_rhs(F):
    """
    For y^2 = F(x), F in QQ[x], set:

        D = lcm denominators(coefficients of F)
        Y* = D Y

    Then:

        (Y*)^2 = D^2 F(x) = F_integral(x)

    with F_integral in ZZ[x].
    """
    F = R(F)
    coeffs = coeffs_from_poly_QQ(F)

    if not coeffs:
        raise ValueError("Zero polynomial is not a valid curve model.")

    denoms = [ZZ(c.denominator()) for c in coeffs]
    D = ZZ(lcm(denoms)) if denoms else ZZ(1)

    F_scaled_QQ = R((D**2) * F)

    # Convert to ZZ[x]. If something is wrong, Sage will raise.
    F_integral = RZ(F_scaled_QQ)

    # Explicit verification.
    if any(c not in ZZ for c in F_integral.list()):
        raise ValueError("Internal error: scaled polynomial is not integral.")

    return D, F_integral


def magma_poly_string(FZ):
    """
    Sage's string representation is already close to Magma syntax:
        x^6 - 3*x^4 + 1
    We keep variable name x.
    """
    return str(FZ)


# ============================================================
# OUTPUT HELPERS
# ============================================================

def get_field(fields, key, default=""):
    return fields.get(key, default)


def get_first_existing_field(fields, keys, default=""):
    for key in keys:
        if key in fields:
            return fields[key]
    return default


def relevant_family_fields_for_output(kind, family_fields):
    """
    Keep compact family-level metadata.
    """
    if kind == "even_automorphism":
        keys = [
            "stage",
            "abs_roots",
            "forced_xs",
            "q_params",
            "q(x)",
            "H(x)",
            "Pq(x)",
            "num_curves_tested_good",
            "good_curve_count(extra_x >= 26)",
            "top_k_count",
            "top_k_min_affine",
            "top_k_max_affine",
            "top_k_sum_affine",
            "top_k_avg_affine",
            "top_k_affine_ge_55",
            "top_k_affine_ge_60",
            "top_k_min_extra_x",
            "top_k_sum_extra_x",
            "best_extra_x_count",
            "best_affine_point_count",
            "best_projective_lower_bound",
            "best_height_F",
        ]
    else:
        keys = [
            "stage",
            "m",
            "rs",
            "q_params",
            "q(x)",
            "H(x)",
            "Pq(x)",
            "Puncte comune fortate",
            "num_curves_tested_good",
            "good_curve_count(extra_x >= 26)",
            "top_k_count",
            "top_k_min_affine",
            "top_k_max_affine",
            "top_k_sum_affine",
            "top_k_avg_affine",
            "top_k_affine_ge_55",
            "top_k_affine_ge_60",
            "top_k_min_extra_x",
            "top_k_sum_extra_x",
            "best_extra_x_count",
            "best_affine_point_count",
            "best_projective_lower_bound",
            "best_height_F",
        ]

    return {k: family_fields[k] for k in keys if k in family_fields}


def relevant_curve_fields_for_output(kind, curve_fields):
    """
    Keep compact curve-level metadata.
    """
    common = [
        "q_params",
        "Height(F)",
        "PARI height",
        "Nr. x-uri suplimentare distincte",
        "Nr. puncte suplimentare",
        "Nr. puncte afine gasite",
        "Lower bound projectiv",
    ]

    if WRITE_LONG_POINT_LISTS:
        common += [
            "x-uri suplimentare",
            "Puncte suplimentare",
        ]

    if kind == "even_automorphism":
        keys = [
            "A",
            "B",
            "Q(x)",
        ] + common
    else:
        keys = [
            "L_tail",
            "L(x)",
        ] + common

    return {k: curve_fields[k] for k in keys if k in curve_fields}


def write_details_file(out_path, source_details_path, kind, execution_rel_dir, normalized_records):
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("=" * 100 + "\n")
        out.write("INTEGRAL ISOMORPHIC MODELS FOR STAGE 2 REFINED CURVES\n")
        out.write("=" * 100 + "\n")
        out.write(f"construction_kind = {kind}\n")
        out.write(f"source_details_path = {source_details_path}\n")
        out.write(f"execution_relative_dir = {execution_rel_dir}\n")
        out.write(f"generated_at = {datetime.now().isoformat(timespec='seconds')}\n")
        out.write("\n")
        out.write("Normalization convention:\n")
        out.write("  Original curve: y^2 = F(X), F in QQ[X]\n")
        out.write("  D = lcm of denominators of coefficients of F\n")
        out.write("  Coordinate change: X* = X, Y* = D * Y\n")
        out.write("  Integral model: (Y*)^2 = D^2 * F(X*) = F_integral(X*)\n")
        out.write("=" * 100 + "\n\n")

        current_family = None

        for rec in normalized_records:
            fam_key = rec["family_index"]

            if fam_key != current_family:
                current_family = fam_key
                out.write("=" * 100 + "\n")
                out.write(f"{rec['family_label']}\n")
                out.write("-" * 100 + "\n")
                for k, v in rec["family_fields_compact"].items():
                    out.write(f"{k} = {v}\n")
                out.write("\n")
                out.write("TOP CURVES IN THIS FAMILY, WRITTEN AS INTEGRAL ISOMORPHIC MODELS\n")

            out.write("-" * 90 + "\n")
            out.write(f"{rec['curve_label']}\n")

            for k, v in rec["curve_fields_compact"].items():
                out.write(f"{k} = {v}\n")

            out.write(f"F_original(x) = {rec['F_original']}\n")
            out.write(f"Height(F_original) = {rec['height_F_original']}\n")
            out.write(f"degree = {rec['degree']}\n")
            out.write(f"leading_coefficient_original = {rec['leading_coefficient_original']}\n")
            out.write(f"is_squarefree_original = {rec['is_squarefree_original']}\n")
            out.write(f"is_even_original = {rec['is_even_original']}\n")
            out.write(f"D_lcm_denominators = {rec['D_lcm_denominators']}\n")
            out.write("Schimbare de coordonate:\n")
            out.write("  X* = X\n")
            out.write(f"  Y* = {rec['D_lcm_denominators']} * Y\n")
            out.write("Model integral izomorf:\n")
            out.write("  (Y*)^2 = F_integral(X*)\n")
            out.write(f"F_integral(x) = {rec['F_integral']}\n")
            out.write(f"Height(F_integral) = {rec['height_F_integral']}\n")
            out.write(f"leading_coefficient_integral = {rec['leading_coefficient_integral']}\n")
            out.write(f"is_squarefree_integral_over_Q = {rec['is_squarefree_integral_over_Q']}\n")
            out.write("\n")


def write_summary_csv(csv_path, rows):
    fieldnames = [
        "construction_kind",
        "execution_relative_dir",
        "source_details_path",
        "family_index",
        "curve_index",
        "family_label",
        "curve_label",

        # family metadata
        "stage",
        "m",
        "rs",
        "abs_roots",
        "forced_xs",
        "q_params_family",
        "q_family",
        "H_family",
        "Pq_family",
        "top_k_min_affine",
        "top_k_max_affine",
        "top_k_sum_affine",
        "best_affine_point_count",
        "best_projective_lower_bound",

        # curve metadata
        "A",
        "B",
        "Q_curve",
        "L_tail",
        "L_curve",
        "pari_height",
        "extra_x_count",
        "extra_point_count",
        "affine_point_count",
        "projective_lower_bound",

        # model data
        "F_original",
        "height_F_original",
        "degree",
        "leading_coefficient_original",
        "is_squarefree_original",
        "is_even_original",
        "D_lcm_denominators",
        "F_integral",
        "height_F_integral",
        "leading_coefficient_integral",
        "is_squarefree_integral_over_Q",
    ]

    with open(csv_path, "w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def quote_magma_string(s):
    s = str(s)
    s = s.replace("\\", "\\\\")
    s = s.replace('"', '\\"')
    return f'"{s}"'


def write_magma_data_file(path, rows):
    """
    Writes a Magma-friendly list:

        Qx<x> := PolynomialRing(Rationals());
        data := [
            <"label", x^6 + ...>,
            ...
        ];

    This is useful as a starting point for rank-computation scripts.
    """
    with open(path, "w", encoding="utf-8") as out:
        out.write("Qx<x> := PolynomialRing(Rationals());\n")
        out.write("\n")
        out.write("data := [\n")

        for i, row in enumerate(rows):
            comma = "," if i + 1 < len(rows) else ""
            label = (
                f"{row['construction_kind']} | "
                f"{row['execution_relative_dir']} | "
                f"{row['family_label']} | "
                f"{row['curve_label']}"
            )
            out.write(
                f"    <{quote_magma_string(label)}, {row['F_integral']}>{comma}\n"
            )

        out.write("];\n")


# ============================================================
# PROCESSING
# ============================================================

def process_details_file(details_path, input_root, output_root, kind, short_name=None):
    details_path = Path(details_path)
    execution_dir = details_path.parent
    execution_rel_dir = execution_dir.relative_to(input_root)

    if short_name is None:
        output_exec_dir = output_root / execution_rel_dir
    else:
        output_exec_dir = output_root / short_name

    output_exec_dir.mkdir(parents=True, exist_ok=True)

    if COPY_RUN_CONFIG:
        run_config_src = execution_dir / "run_config.txt"
        if run_config_src.exists():
            shutil.copy2(run_config_src, output_exec_dir / "source_run_config.txt")

    text = details_path.read_text(encoding="utf-8", errors="replace")
    family_records = extract_family_records(text)

    normalized_records = []
    summary_rows = []

    for family in family_records:
        family_fields = family["family_fields"]
        family_fields_compact = relevant_family_fields_for_output(kind, family_fields)

        for curve in family["curves"]:
            curve_fields = curve["curve_fields"]

            F_original = parse_sage_poly(curve_fields["F(x)"])
            D, F_integral = integral_model_rhs(F_original)

            row = {
                "construction_kind": kind,
                "execution_relative_dir": str(execution_rel_dir),
                "source_details_path": str(details_path),
                "family_index": family["family_index"],
                "curve_index": curve["curve_index"],
                "family_label": family["family_label"],
                "curve_label": curve["curve_label"],

                "stage": get_field(family_fields, "stage"),
                "m": get_field(family_fields, "m"),
                "rs": get_field(family_fields, "rs"),
                "abs_roots": get_field(family_fields, "abs_roots"),
                "forced_xs": get_field(family_fields, "forced_xs"),
                "q_params_family": get_field(family_fields, "q_params"),
                "q_family": get_field(family_fields, "q(x)"),
                "H_family": get_field(family_fields, "H(x)"),
                "Pq_family": get_field(family_fields, "Pq(x)"),
                "top_k_min_affine": get_field(family_fields, "top_k_min_affine"),
                "top_k_max_affine": get_field(family_fields, "top_k_max_affine"),
                "top_k_sum_affine": get_field(family_fields, "top_k_sum_affine"),
                "best_affine_point_count": get_field(family_fields, "best_affine_point_count"),
                "best_projective_lower_bound": get_field(family_fields, "best_projective_lower_bound"),

                "A": get_field(curve_fields, "A"),
                "B": get_field(curve_fields, "B"),
                "Q_curve": get_field(curve_fields, "Q(x)"),
                "L_tail": get_field(curve_fields, "L_tail"),
                "L_curve": get_field(curve_fields, "L(x)"),
                "pari_height": get_field(curve_fields, "PARI height"),
                "extra_x_count": get_field(curve_fields, "Nr. x-uri suplimentare distincte"),
                "extra_point_count": get_field(curve_fields, "Nr. puncte suplimentare"),
                "affine_point_count": get_field(curve_fields, "Nr. puncte afine gasite"),
                "projective_lower_bound": get_field(curve_fields, "Lower bound projectiv"),

                "F_original": str(F_original),
                "height_F_original": str(polynomial_height_QQ(F_original)),
                "degree": str(F_original.degree()),
                "leading_coefficient_original": str(F_original.leading_coefficient()),
                "is_squarefree_original": str(is_squarefree_over_QQ(F_original)),
                "is_even_original": str(is_even_polynomial(F_original)),
                "D_lcm_denominators": str(D),
                "F_integral": str(F_integral),
                "height_F_integral": str(polynomial_height_ZZ(F_integral)),
                "leading_coefficient_integral": str(F_integral.leading_coefficient()),
                "is_squarefree_integral_over_Q": str(is_squarefree_over_QQ(R(F_integral))),
            }

            compact_curve_fields = relevant_curve_fields_for_output(kind, curve_fields)

            normalized_records.append({
                **row,
                "family_fields_compact": family_fields_compact,
                "curve_fields_compact": compact_curve_fields,
            })
            summary_rows.append(row)

    details_out = output_exec_dir / OUTPUT_DETAILS_FILENAME
    csv_out = output_exec_dir / OUTPUT_SUMMARY_FILENAME

    write_details_file(
        out_path=details_out,
        source_details_path=details_path,
        kind=kind,
        execution_rel_dir=execution_rel_dir,
        normalized_records=normalized_records,
    )
    write_summary_csv(csv_out, summary_rows)

    magma_out = None
    if WRITE_MAGMA_DATA_FILE:
        magma_out = output_exec_dir / OUTPUT_MAGMA_DATA_FILENAME
        write_magma_data_file(magma_out, summary_rows)

    return {
        "execution_relative_dir": str(execution_rel_dir),
        "details_path": str(details_path),
        "output_dir": str(output_exec_dir),
        "n_families": len(family_records),
        "n_curves": len(summary_rows),
        "details_out": str(details_out),
        "csv_out": str(csv_out),
        "magma_out": str(magma_out) if magma_out else "",
        "rows": summary_rows,
    }


def find_stage2_details_files(input_root):
    input_root = Path(input_root)
    if not input_root.exists():
        return []

    if RECURSIVE:
        return sorted(input_root.rglob(DETAILS_FILENAME))

    return sorted(input_root.glob(f"*/{DETAILS_FILENAME}"))


def process_source_root(spec):
    kind = spec["kind"]
    input_root = Path(spec["input_root"])
    output_root = Path(spec["output_root"])
    output_root.mkdir(parents=True, exist_ok=True)

    details_files = find_stage2_details_files(input_root)

    print("=" * 100)
    print(f"Processing source kind: {kind}")
    print(f"Input root : {input_root}")
    print(f"Output root: {output_root}")
    print(f"Found {len(details_files)} file(s) named {DETAILS_FILENAME}")
    print("=" * 100)

    aggregate_rows = []
    execution_index_rows = []

    for idx, details_path in enumerate(details_files, start=1):
        short_name = f"run_{idx:04d}" if USE_SHORT_EXECUTION_FOLDER_NAMES else None

        try:
            result = process_details_file(
                details_path=details_path,
                input_root=input_root,
                output_root=output_root,
                kind=kind,
                short_name=short_name,
            )
        except Exception as e:
            print(f"ERROR while processing {details_path}: {e}")
            continue

        aggregate_rows.extend(result["rows"])

        execution_index_rows.append({
            "index": idx,
            "short_name": short_name or "",
            "execution_relative_dir": result["execution_relative_dir"],
            "source_details_path": result["details_path"],
            "output_dir": result["output_dir"],
            "n_families": result["n_families"],
            "n_curves": result["n_curves"],
            "details_out": result["details_out"],
            "csv_out": result["csv_out"],
            "magma_out": result["magma_out"],
        })

        print(
            f"[{idx}/{len(details_files)}] {result['execution_relative_dir']} "
            f"-> families={result['n_families']}, curves={result['n_curves']}"
        )

    aggregate_csv = output_root / AGGREGATE_SUMMARY_FILENAME
    write_summary_csv(aggregate_csv, aggregate_rows)

    index_csv = output_root / EXECUTION_INDEX_FILENAME
    with open(index_csv, "w", newline="", encoding="utf-8") as fout:
        fieldnames = [
            "index",
            "short_name",
            "execution_relative_dir",
            "source_details_path",
            "output_dir",
            "n_families",
            "n_curves",
            "details_out",
            "csv_out",
            "magma_out",
        ]
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for row in execution_index_rows:
            writer.writerow(row)

    print(f"Aggregate CSV saved: {aggregate_csv}")
    print(f"Execution index saved: {index_csv}")
    print()

    return {
        "kind": kind,
        "input_root": str(input_root),
        "output_root": str(output_root),
        "n_details_files": len(details_files),
        "n_curves_total": len(aggregate_rows),
        "aggregate_csv": str(aggregate_csv),
        "execution_index_csv": str(index_csv),
    }


def main():
    print("=" * 100)
    print("CREATE INTEGRAL ISOMORPHIC MODELS FROM STAGE 2 DETAILS")
    print("=" * 100)

    results = []

    for spec in SOURCE_ROOTS:
        results.append(process_source_root(spec))

    print("=" * 100)
    print("DONE")
    print("=" * 100)
    for res in results:
        print(
            f"{res['kind']}: details_files={res['n_details_files']}, "
            f"curves={res['n_curves_total']}, output={res['output_root']}"
        )


if __name__ == "__main__":
    main()
