#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Robust v2 asynchronous Magma Calculator rank requests for integral genus-2 hyperelliptic models.

Input:
    CSV files produced by create_integral_models_from_stage2_details.py, especially:

        families_construct_even_automorphism_two_stage_integral_models/all_stage2_integral_models_summary.csv
        families_construct_general_two_stage_integral_models/all_stage2_integral_models_summary.csv

The script:
    1. Reads all input CSVs.
    2. Deduplicates curves by exact F_integral string.
    3. Sends one request per unique F_integral to Magma Online Calculator.
    4. Saves raw XML responses immediately.
    5. Parses RankBounds / MordellWeilGroupGenus2 output.
    6. Writes:
        - unique_curves_to_magma.csv
        - curve_appearances_mapping.csv
        - magma_rank_results_unique.csv
        - magma_rank_results_all_appearances.csv

Magma logic for each curve:
    - Compute RankBounds(J).
    - If lb = ub, rank_Jacobian = lb.
    - Else call MordellWeilGroupGenus2(J : RankOnly := true).
    - If finiteIndex=true and rank_G=ub2, rank_Jacobian = rank_G.
    - Otherwise the rank is marked as not fully determined by this request.

Important:
    Keep CONCURRENCY small. The Magma Calculator is a public service with an internal
    execution limit. Start with CONCURRENCY = 2.
"""

import asyncio
import csv
import hashlib
import html
import re
import time
import shutil
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import aiohttp
except ImportError as exc:
    raise SystemExit(
        "Missing dependency: aiohttp\n"
        "Install it with:\n"
        "    pip install aiohttp\n"
        "or, in conda:\n"
        "    conda install -c conda-forge aiohttp\n"
    ) from exc


# ============================================================
# CONFIGURATION
# ============================================================

MAGMA_URL = "https://magma.maths.usyd.edu.au/xml/calculator.xml"

INPUT_CSVS = [
    Path("families_construct_d4_reciprocal/d4_top_curves_summary.csv"),
]

OUTPUT_DIR = Path("magma_rank_results_d4_reciprocal")
RAW_XML_DIRNAME = "xml"
BAD_CACHE_XML_DIRNAME = "xml_bad_cache"
OFFLINE_XML_DIRNAME = "xml_offline"

# Balanced setting for the public Magma Calculator.
# This controls how many computations may be active simultaneously.
CONCURRENCY = 2

# Global request-start limiter shared by all coroutines.
# 10 seconds/request-start gives approximately 360 request starts/hour.
# This is NOT a per-coroutine delay; it spaces request starts globally.
GLOBAL_MIN_SECONDS_BETWEEN_REQUEST_STARTS = 10.0

# Local HTTP timeout. Magma has its own internal 120-second execution cap.
REQUEST_TIMEOUT_SECONDS = 180
MAX_HTTP_RETRIES = 3

# Resume mode: if XML exists and is a GOOD cache, do not re-submit the curve.
SKIP_EXISTING_XML = True

# Delete generated CSVs at script start so repeated runs do not accumulate append duplicates.
# This does NOT delete XML caches.
CLEAN_OUTPUT_CSVS_ON_START = True

# Existing XMLs with these statuses are treated as bad/transient and rerun.
# Memory limit is intentionally NOT included: it is usually a real computational failure.
RERUN_CACHE_STATUSES = {
    "unknown",
    "magma_calculator_offline",
    "invalid_calculator_response_no_headers",
    "invalid_calculator_response_no_time",
    "empty_results",
    "xml_parse_error",
    "magma_output_unparsed",
    "magma_user_error_unparsed",
    "magma_runtime_error_unparsed",
    "transient_failure_after_retries",
}
RERUN_BAD_CACHED_XML = True

# New responses with these statuses are considered transient and retried immediately
# for the SAME curve. HTTP/network retries are separate from these XML-level retries.
TRANSIENT_BAD_STATUSES = {
    "unknown",
    "invalid_calculator_response_no_headers",
    "invalid_calculator_response_no_time",
    "empty_results",
    "xml_parse_error",
    "magma_output_unparsed",
}

MAX_TRANSIENT_RERUNS_PER_CURVE = 3
TRANSIENT_RETRY_SLEEP_SECONDS = 90

# If Magma returns <offline>, pause and probe until it appears available again.
SLEEP_ON_OFFLINE = True
OFFLINE_SLEEP_SECONDS = 2 * 60
PROBE_INTERVAL_SECONDS = 60
OFFLINE_TRIGGER_THRESHOLD = 1
PROBE_MAGMA_SCRIPT = 'print "PING_OK";\n'

USE_GRH = True
USE_VERBOSE = False

# Useful for testing. Set to None for all curves.
PROCESS_ONLY_FIRST_N_UNIQUE_CURVES = None

UNIQUE_CURVES_CSV = "unique_curves_to_magma.csv"
APPEARANCES_CSV = "curve_appearances_mapping.csv"
RESULTS_UNIQUE_CSV = "magma_rank_results_unique.csv"
RESULTS_ALL_APPEARANCES_CSV = "magma_rank_results_all_appearances.csv"
FAILED_REQUESTS_CSV = "magma_rank_failed_requests.csv"


# ============================================================
# REGEX FOR PARSING MAGMA OUTPUT
# ============================================================

RANKBOUNDS_RE = re.compile(r"RankBounds(?:\(J\))?\s*[:=]\s*\[([^\],]+),\s*([^\]]+)\]")
DIRECT_RANK_RE = re.compile(r"RANK DETERMINAT DIRECT DIN RankBounds:\s*rank J\(Q\)\s*=\s*(-?\d+)")
RANK_G_RE = re.compile(r"rank\(G\)\s*=\s*(-?\d+)")
FINITE_RE = re.compile(
    r"finiteIndex\s*=\s*(true|false)\s*,\s*proved\s*=\s*(true|false)\s*,\s*ub2\s*=\s*(-?\d+)",
    re.IGNORECASE,
)
FINAL_RANK_RE = re.compile(r"RANK DETERMINAT:\s*rank J\(Q\)\s*=\s*(-?\d+)")
UNRESOLVED_RE = re.compile(r"RANG NEINCHIS", re.IGNORECASE)
ERROR_RE = re.compile(r"EROARE LA ACEASTA CURBA", re.IGNORECASE)


# ============================================================
# CSV HELPERS
# ============================================================

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_csv_rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        print(f"WARNING: input CSV does not exist: {path}")
        return []

    with open(path, "r", newline="", encoding="utf-8") as fin:
        reader = csv.DictReader(fin)
        rows = []
        for i, row in enumerate(reader, start=1):
            row = dict(row)
            row["_source_csv"] = str(path)
            row["_source_row_index"] = str(i)
            rows.append(row)
        return rows


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with open(path, "w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def append_csv_row(path: Path, row: Dict[str, object], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    exists = path.exists()

    with open(path, "a", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerow({k: row.get(k, "") for k in fieldnames})


def remove_file_if_exists(path: Path) -> None:
    try:
        path.unlink()
        print(f"[clean] removed {path}")
    except FileNotFoundError:
        pass


def clean_output_csvs(output_dir: Path) -> None:
    if not CLEAN_OUTPUT_CSVS_ON_START:
        return
    for name in [
        RESULTS_UNIQUE_CSV,
        RESULTS_ALL_APPEARANCES_CSV,
        FAILED_REQUESTS_CSV,
    ]:
        remove_file_if_exists(output_dir / name)


# ============================================================
# IDENTIFIERS / LABELS
# ============================================================

def normalize_poly_string(poly: str) -> str:
    """
    Exact-string deduplication with light whitespace normalization.

    This does NOT detect isomorphic curves or differently scaled integral models.
    It only deduplicates identical F_integral expressions.
    """
    return re.sub(r"\s+", " ", str(poly).strip())


def short_hash(s: str, n: int = 16) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:n]


def sanitize_filename(s: str, max_len: int = 120) -> str:
    s = str(s).strip()
    s = re.sub(r"[^\w\-\.]+", "_", s)
    s = re.sub(r"_+", "_", s)
    s = s.strip("_")
    if len(s) > max_len:
        s = s[:max_len].strip("_")
    return s or "curve"


def build_label(row: Dict[str, str], unique_id: str) -> str:
    """Short label printed by Magma. Full metadata remains in CSV files."""
    kind = row.get("first_construction_kind", "")
    if kind == "even_automorphism":
        kind_short = "even"
    elif kind == "general":
        kind_short = "general"
    else:
        kind_short = kind or "unknown"

    fam = row.get("first_family_label", "").replace("FAMILIA #", "F")
    cur = row.get("first_curve_label", "").replace("CURBA #", "C")
    return f"{unique_id} | {kind_short} | {fam} | {cur}"


# ============================================================
# DEDUPLICATION
# ============================================================

def load_and_deduplicate(input_csvs: List[Path]) -> Tuple[List[Dict[str, str]], List[Dict[str, str]]]:
    """
    Returns:
        unique_curves: one row per unique F_integral
        appearances: one row per original CSV row, mapped to unique_id
    """
    all_rows: List[Dict[str, str]] = []

    for path in input_csvs:
        all_rows.extend(read_csv_rows(path))

    unique_by_poly: Dict[str, Dict[str, str]] = {}
    appearances: List[Dict[str, str]] = []

    for row in all_rows:
        F_integral = normalize_poly_string(row.get("F_integral", ""))

        if not F_integral:
            print(f"WARNING: missing F_integral in {row.get('_source_csv')} row {row.get('_source_row_index')}")
            continue

        key = F_integral
        h = short_hash(key, 16)

        if key not in unique_by_poly:
            unique_id = f"curve_{len(unique_by_poly) + 1:06d}_{h}"
            unique_row = {
                "unique_id": unique_id,
                "F_integral": F_integral,
                "poly_hash": h,
                "first_source_csv": row.get("_source_csv", ""),
                "first_source_row_index": row.get("_source_row_index", ""),
                "first_construction_kind": row.get("construction_kind", ""),
                "first_execution_relative_dir": row.get("execution_relative_dir", ""),
                "first_family_index": row.get("family_index", ""),
                "first_curve_index": row.get("curve_index", ""),
                "first_family_label": row.get("family_label", ""),
                "first_curve_label": row.get("curve_label", ""),
                "first_affine_point_count": row.get("affine_point_count", ""),
                "first_projective_lower_bound": row.get("projective_lower_bound", ""),
                "height_F_integral": row.get("height_F_integral", ""),
                "degree": row.get("degree", ""),
                "is_squarefree_integral_over_Q": row.get("is_squarefree_integral_over_Q", ""),
                "appearance_count": "0",
            }
            unique_by_poly[key] = unique_row

        unique_row = unique_by_poly[key]
        unique_row["appearance_count"] = str(int(unique_row["appearance_count"]) + 1)

        app = {
            "unique_id": unique_row["unique_id"],
            "poly_hash": unique_row["poly_hash"],
            "_source_csv": row.get("_source_csv", ""),
            "_source_row_index": row.get("_source_row_index", ""),
            "construction_kind": row.get("construction_kind", ""),
            "execution_relative_dir": row.get("execution_relative_dir", ""),
            "family_index": row.get("family_index", ""),
            "curve_index": row.get("curve_index", ""),
            "family_label": row.get("family_label", ""),
            "curve_label": row.get("curve_label", ""),
            "affine_point_count": row.get("affine_point_count", ""),
            "projective_lower_bound": row.get("projective_lower_bound", ""),
            "extra_x_count": row.get("extra_x_count", ""),
            "F_integral": F_integral,
        }
        appearances.append(app)

    unique_curves = list(unique_by_poly.values())
    return unique_curves, appearances


# ============================================================
# MAGMA SCRIPT GENERATION
# ============================================================

def magma_quote_string(s: str) -> str:
    s = str(s)
    s = s.replace("\\", "\\\\")
    s = s.replace('"', '\\"')
    return f'"{s}"'


def build_single_curve_magma_script(label: str, poly: str) -> str:
    """
    Build a Magma script for one curve only.
    """
    parts: List[str] = []

    parts.append("Qx<x> := PolynomialRing(Rationals());")
    parts.append(f"f := {poly};")

    if USE_VERBOSE:
        parts.append('SetVerbose("MordellWeilGroup", 1);')

    if USE_GRH:
        parts.append('SetClassGroupBounds("GRH");')

    parts.append("")
    parts.append(f"label := {magma_quote_string(label)};")
    parts.append('print "============================================================";')
    parts.append("print label;")
    parts.append('print "f(x) =", f;')
    parts.append("")

    parts.append("try")
    parts.append("    C := HyperellipticCurve(f);")
    parts.append("    J := Jacobian(C);")
    parts.append("")
    parts.append("    lb, ub := RankBounds(J);")
    parts.append('    printf "RankBounds: [%o, %o]\\n", lb, ub;')
    parts.append("")
    parts.append("    if lb eq ub then")
    parts.append('        printf "RANK DETERMINAT DIRECT DIN RankBounds: rank J(Q) = %o\\n", lb;')
    parts.append("    else")
    parts.append("        G, phi, finiteIndex, proved, ub2 := MordellWeilGroupGenus2(J : RankOnly := true);")
    parts.append("        rg := TorsionFreeRank(G);")
    parts.append('        printf "rank(G) = %o\\n", rg;')
    parts.append('        printf "finiteIndex = %o, proved = %o, ub2 = %o\\n", finiteIndex, proved, ub2;')
    parts.append("")
    parts.append("        if finiteIndex and rg eq ub2 then")
    parts.append('            printf "RANK DETERMINAT: rank J(Q) = %o\\n", rg;')
    parts.append("        else")
    parts.append('            print "RANG NEINCHIS: lower bound intern =", rg, ", upper bound =", ub2;')
    parts.append("        end if;")
    parts.append("    end if;")
    parts.append("")
    parts.append("catch e")
    parts.append('    print "EROARE LA ACEASTA CURBA:";')
    parts.append("    print e`Object;")
    parts.append("end try;")
    parts.append("")

    return "\n".join(parts)


# ============================================================
# MAGMA XML PARSING
# ============================================================

def parse_magma_xml(xml_text: str) -> Dict[str, object]:
    out: Dict[str, object] = {
        "max_time": "",
        "max_input": "",
        "seed": "",
        "version": "",
        "time": "",
        "memory": "",
        "warning": "",
        "lines_text": "",
        "rankbounds_lb": "",
        "rankbounds_ub": "",
        "rank_direct": "",
        "rank_g": "",
        "finiteIndex": "",
        "proved": "",
        "ub2": "",
        "final_rank_printed": "",
        "rank_Jacobian": "",
        "rank_status": "",
        "rank_source": "",
        "assumption": "GRH" if USE_GRH else "",
        "magma_error": "",
        "unresolved": "",
    }

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        out["lines_text"] = xml_text
        out["rank_status"] = "xml_parse_error"
        out["magma_error"] = "could_not_parse_xml"
        return out

    offline_node = root.find("offline")
    if offline_node is not None:
        msg = offline_node.text.strip() if offline_node.text else "Magma calculator offline"
        out["lines_text"] = msg
        out["rank_status"] = "magma_calculator_offline"
        out["magma_error"] = "true"
        out["rank_source"] = "offline_response_from_calculator"
        return out

    headers = root.find("headers")
    if headers is None:
        out["lines_text"] = xml_text
        out["rank_status"] = "invalid_calculator_response_no_headers"
        out["magma_error"] = "true"
        out["rank_source"] = "no_headers_in_xml_response"
        return out

    for key in ["max_time", "max_input", "seed", "version", "time", "memory"]:
        node = headers.find(key)
        if node is not None and node.text is not None:
            out[key] = node.text.strip()

    warning_node = headers.find("warning")
    if warning_node is not None and warning_node.text is not None:
        out["warning"] = warning_node.text.strip()

    if out["time"] == "":
        out["lines_text"] = xml_text
        out["rank_status"] = "invalid_calculator_response_no_time"
        out["magma_error"] = "true"
        out["rank_source"] = "missing_time_header"
        return out

    results = root.find("results")
    lines = []

    if results is not None:
        for line in results.findall("line"):
            txt = line.text if line.text is not None else ""
            lines.append(html.unescape(txt))

    out["lines_text"] = "\n".join(lines)

    for line in lines:
        m = RANKBOUNDS_RE.search(line)
        if m:
            out["rankbounds_lb"] = m.group(1).strip()
            out["rankbounds_ub"] = m.group(2).strip()

        m = DIRECT_RANK_RE.search(line)
        if m:
            out["rank_direct"] = m.group(1).strip()

        m = RANK_G_RE.search(line)
        if m:
            out["rank_g"] = m.group(1).strip()

        m = FINITE_RE.search(line)
        if m:
            out["finiteIndex"] = m.group(1).lower()
            out["proved"] = m.group(2).lower()
            out["ub2"] = m.group(3).strip()

        m = FINAL_RANK_RE.search(line)
        if m:
            out["final_rank_printed"] = m.group(1).strip()

        if UNRESOLVED_RE.search(line):
            out["unresolved"] = "true"

        if ERROR_RE.search(line):
            out["magma_error"] = "true"

    combined = (out.get("warning", "") + "\n" + out.get("lines_text", "")).lower()

    if out["rank_direct"] != "":
        out["rank_Jacobian"] = out["rank_direct"]
        out["rank_status"] = "determined_by_RankBounds"
        out["rank_source"] = "RankBounds_lb_eq_ub"
        return out

    if out["rank_g"] != "" and out["finiteIndex"] == "true" and out["ub2"] != "" and out["rank_g"] == out["ub2"]:
        out["rank_Jacobian"] = out["rank_g"]
        out["rank_status"] = "determined_by_MordellWeilGroupGenus2"
        out["rank_source"] = "finiteIndex_true_and_rankG_eq_ub2"
        return out

    if out["rank_g"] != "":
        out["rank_status"] = "not_closed_after_MordellWeilGroupGenus2"
        out["rank_source"] = "rankG_available_but_not_certified_as_final"
        return out

    if out["rankbounds_lb"] != "" and out["rankbounds_ub"] != "":
        out["rank_status"] = "only_RankBounds_available"
        out["rank_source"] = "RankBounds_only"
        return out

    if "memory limit" in combined or "user memory limit" in combined:
        out["rank_status"] = "magma_memory_limit"
        out["magma_error"] = "true"
        out["rank_source"] = "magma_computation_exceeded_memory_limit"
        return out

    if "time limit" in combined or "timed out" in combined:
        out["rank_status"] = "magma_time_limit"
        out["magma_error"] = "true"
        out["rank_source"] = "magma_computation_exceeded_time_limit"
        return out

    if "runtime error" in combined:
        out["rank_status"] = "magma_runtime_error_unparsed"
        out["magma_error"] = "true"
        out["rank_source"] = "runtime_error_detected"
        return out

    if "user error" in combined:
        out["rank_status"] = "magma_user_error_unparsed"
        out["magma_error"] = "true"
        out["rank_source"] = "user_error_detected"
        return out

    if not out["lines_text"].strip():
        out["rank_status"] = "empty_results"
        out["rank_source"] = "xml_ok_but_no_result_lines"
        return out

    out["rank_status"] = "magma_output_unparsed"
    out["rank_source"] = "xml_ok_but_expected_rank_patterns_not_found"
    return out


def is_good_cache_status(rank_status: str) -> bool:
    return rank_status not in RERUN_CACHE_STATUSES


def is_offline_status(rank_status: str) -> bool:
    return rank_status == "magma_calculator_offline"


def is_transient_bad_status(rank_status: str) -> bool:
    return rank_status in TRANSIENT_BAD_STATUSES


def archive_xml_text(directory: Path, stem: str, reason: str, xml_text: str) -> Path:
    ensure_dir(directory)
    safe_reason = sanitize_filename(reason or "xml", max_len=60)
    target = directory / f"{stem}__{safe_reason}__{datetime.now().strftime('%Y%m%d_%H%M%S')}.xml"
    target.write_text(xml_text, encoding="utf-8")
    return target


def archive_bad_xml(xml_path: Path, bad_dir: Path, reason: str) -> None:
    ensure_dir(bad_dir)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_reason = sanitize_filename(reason or "bad_cache", max_len=60)
    target = bad_dir / f"{xml_path.stem}__{safe_reason}__{timestamp}{xml_path.suffix}"
    try:
        shutil.move(str(xml_path), str(target))
    except FileNotFoundError:
        return



# ============================================================
# GLOBAL RATE LIMITER
# ============================================================

class GlobalRateLimiter:
    """
    Global request-start limiter shared by all coroutines.

    CONCURRENCY controls how many computations may be active simultaneously.
    This limiter controls how close in time two POST requests may start globally.

    Example:
        CONCURRENCY = 2
        GLOBAL_MIN_SECONDS_BETWEEN_REQUEST_STARTS = 25.0

    Then up to two computations may overlap, but request starts are spaced:
        t=0, t=25, t=50, ...
    """

    def __init__(self, min_interval_seconds: float):
        self.min_interval_seconds = float(min_interval_seconds)
        self.lock = asyncio.Lock()
        self.last_start_time = 0.0

    async def wait_turn(self) -> None:
        async with self.lock:
            now = time.monotonic()
            elapsed = now - self.last_start_time
            wait_seconds = self.min_interval_seconds - elapsed

            if wait_seconds > 0:
                await asyncio.sleep(wait_seconds)

            self.last_start_time = time.monotonic()


# ============================================================
# ASYNC HTTP
# ============================================================

async def post_to_magma(session: aiohttp.ClientSession, script: str) -> str:
    last_error: Optional[Exception] = None

    for attempt in range(MAX_HTTP_RETRIES + 1):
        try:
            async with session.post(
                MAGMA_URL,
                data={"input": script},
                headers={
                    "User-Agent": "Mozilla/5.0",
                    "Referer": "https://magma.maths.usyd.edu.au/calc/",
                    "Content-Type": "application/x-www-form-urlencoded",
                },
            ) as response:
                response.raise_for_status()
                return await response.text()

        except Exception as exc:
            last_error = exc
            if attempt < MAX_HTTP_RETRIES:
                await asyncio.sleep(2.0 * (attempt + 1))

    raise RuntimeError(f"Magma POST failed after {MAX_HTTP_RETRIES + 1} attempts: {last_error}")


async def probe_magma_available(session: aiohttp.ClientSession, rate_limiter: GlobalRateLimiter) -> bool:
    try:
        await rate_limiter.wait_turn()
        xml_text = await post_to_magma(session, PROBE_MAGMA_SCRIPT)
    except Exception as exc:
        print(f"[probe] HTTP/request error: {exc}")
        return False

    parsed = parse_magma_xml(xml_text)

    if parsed.get("rank_status") == "magma_calculator_offline":
        print("[probe] Magma still offline.")
        return False

    if parsed.get("time", "") != "":
        print("[probe] Magma appears available again.")
        return True

    print(f"[probe] Unexpected status: {parsed.get('rank_status')}")
    return False


async def wait_until_magma_available(session: aiohttp.ClientSession, rate_limiter: GlobalRateLimiter) -> None:
    if not SLEEP_ON_OFFLINE:
        return

    print("=" * 100)
    print("Magma Calculator returned <offline>. Pausing requests.")
    print(f"Initial sleep: {OFFLINE_SLEEP_SECONDS} seconds.")
    print("=" * 100)
    await asyncio.sleep(OFFLINE_SLEEP_SECONDS)

    while True:
        if await probe_magma_available(session, rate_limiter):
            print("=" * 100)
            print("Magma Calculator appears available. Resuming.")
            print("=" * 100)
            return
        print(f"[probe] Sleeping {PROBE_INTERVAL_SECONDS} seconds before next probe.")
        await asyncio.sleep(PROBE_INTERVAL_SECONDS)


async def process_one_curve(
    session: aiohttp.ClientSession,
    semaphore: asyncio.Semaphore,
    unique_curve: Dict[str, str],
    result_fields: List[str],
    failed_fields: List[str],
    results_csv_path: Path,
    failed_csv_path: Path,
    xml_dir: Path,
    bad_cache_dir: Path,
    offline_dir: Path,
    completed_counter: Dict[str, int],
    total: int,
    offline_state: Dict[str, int],
    rate_limiter: GlobalRateLimiter,
) -> Dict[str, object]:
    async with semaphore:
        unique_id = unique_curve["unique_id"]
        F_integral = unique_curve["F_integral"]
        label = build_label(unique_curve, unique_id)

        xml_filename = sanitize_filename(f"{unique_id}_{unique_curve.get('poly_hash', '')}") + ".xml"
        xml_path = xml_dir / xml_filename
        xml_stem = xml_path.stem

        # 1. Existing cache: use only if it is a good, non-transient result.
        if SKIP_EXISTING_XML and xml_path.exists():
            xml_text = xml_path.read_text(encoding="utf-8", errors="replace")
            parsed = parse_magma_xml(xml_text)
            cached_status = str(parsed.get("rank_status", ""))

            if is_good_cache_status(cached_status):
                result_row = {
                    **unique_curve,
                    "xml_file": str(xml_path),
                    "request_status": "cached",
                    "elapsed_local_seconds": "0.000",
                    **parsed,
                }
                append_csv_row(results_csv_path, result_row, result_fields)

                completed_counter["done"] += 1
                print(
                    f"[{completed_counter['done']}/{total}] {unique_id} -> "
                    f"cached, rank_status={result_row.get('rank_status')}, "
                    f"rank={result_row.get('rank_Jacobian')}, magma_time={result_row.get('time')}"
                )
                return result_row

            if RERUN_BAD_CACHED_XML:
                archive_bad_xml(xml_path, bad_cache_dir, cached_status or "bad_cache")
                print(f"[cache] {unique_id}: existing XML status={cached_status}; moved to bad cache and will rerun.")
            else:
                result_row = {
                    **unique_curve,
                    "xml_file": str(xml_path),
                    "request_status": "cached_bad_not_rerun",
                    "elapsed_local_seconds": "0.000",
                    **parsed,
                }
                append_csv_row(results_csv_path, result_row, result_fields)
                completed_counter["done"] += 1
                return result_row

        # 2. Request this curve. Offline and transient responses are not accepted as final.
        script = build_single_curve_magma_script(label=label, poly=F_integral)
        start_time = time.time()
        transient_attempts = 0

        while True:
            await rate_limiter.wait_turn()

            try:
                xml_text = await post_to_magma(session, script)
            except Exception as exc:
                elapsed = time.time() - start_time
                failed_row = {
                    "unique_id": unique_id,
                    "poly_hash": unique_curve.get("poly_hash", ""),
                    "F_integral": F_integral,
                    "error": str(exc),
                    "elapsed_local_seconds": f"{elapsed:.3f}",
                }
                append_csv_row(failed_csv_path, failed_row, failed_fields)

                completed_counter["done"] += 1
                print(f"[{completed_counter['done']}/{total}] {unique_id} -> REQUEST ERROR: {exc}")
                return {
                    **unique_curve,
                    "xml_file": "",
                    "request_status": f"error: {exc}",
                    "elapsed_local_seconds": f"{elapsed:.3f}",
                    "rank_status": "request_error",
                }

            parsed = parse_magma_xml(xml_text)
            rank_status = str(parsed.get("rank_status", ""))

            # 2a. Global offline response: sleep/probe and retry the SAME curve.
            if is_offline_status(rank_status):
                offline_state["consecutive_offline"] += 1
                offline_path = archive_xml_text(offline_dir, xml_stem, "offline", xml_text)
                print(f"[offline] {unique_id}: Magma returned offline XML. Saved to {offline_path}")

                if offline_state["consecutive_offline"] >= OFFLINE_TRIGGER_THRESHOLD:
                    await wait_until_magma_available(session, rate_limiter)
                    offline_state["consecutive_offline"] = 0

                continue

            offline_state["consecutive_offline"] = 0

            # 2b. Transient XML-level bad response: retry the SAME curve a few times.
            if is_transient_bad_status(rank_status):
                transient_attempts += 1
                bad_path = archive_xml_text(
                    bad_cache_dir,
                    xml_stem,
                    f"{rank_status}_attempt{transient_attempts}",
                    xml_text,
                )
                print(
                    f"[transient] {unique_id}: status={rank_status}; saved to {bad_path}; "
                    f"attempt {transient_attempts}/{MAX_TRANSIENT_RERUNS_PER_CURVE}"
                )

                if transient_attempts <= MAX_TRANSIENT_RERUNS_PER_CURVE:
                    await asyncio.sleep(TRANSIENT_RETRY_SLEEP_SECONDS)
                    continue

                # Too many transient failures: save the last response as a final failure marker
                # so the run can progress. It is included in RERUN_CACHE_STATUSES, so a future
                # run can retry it again if desired.
                parsed["rank_status"] = "transient_failure_after_retries"
                parsed["rank_source"] = f"last_transient_status={rank_status}"
                parsed["magma_error"] = "true"
                xml_path.write_text(xml_text, encoding="utf-8")

                elapsed = time.time() - start_time
                result_row = {
                    **unique_curve,
                    "xml_file": str(xml_path),
                    "request_status": "transient_failure_after_retries",
                    "elapsed_local_seconds": f"{elapsed:.3f}",
                    **parsed,
                }
                append_csv_row(results_csv_path, result_row, result_fields)

                completed_counter["done"] += 1
                print(
                    f"[{completed_counter['done']}/{total}] {unique_id} -> "
                    f"transient_failure_after_retries, last_status={rank_status}, magma_time={result_row.get('time')}"
                )
                return result_row

            # 2c. Good/final/nontransient status. Save in main XML cache and complete.
            xml_path.write_text(xml_text, encoding="utf-8")

            elapsed = time.time() - start_time
            result_row = {
                **unique_curve,
                "xml_file": str(xml_path),
                "request_status": "ok",
                "elapsed_local_seconds": f"{elapsed:.3f}",
                **parsed,
            }
            append_csv_row(results_csv_path, result_row, result_fields)

            completed_counter["done"] += 1
            print(
                f"[{completed_counter['done']}/{total}] {unique_id} -> "
                f"ok, rank_status={result_row.get('rank_status')}, "
                f"rank={result_row.get('rank_Jacobian')}, magma_time={result_row.get('time')}"
            )
            return result_row


async def run_requests(unique_curves: List[Dict[str, str]], output_dir: Path) -> List[Dict[str, object]]:
    ensure_dir(output_dir)
    xml_dir = output_dir / RAW_XML_DIRNAME
    bad_cache_dir = output_dir / BAD_CACHE_XML_DIRNAME
    offline_dir = output_dir / OFFLINE_XML_DIRNAME
    ensure_dir(xml_dir)
    ensure_dir(bad_cache_dir)
    ensure_dir(offline_dir)

    results_csv_path = output_dir / RESULTS_UNIQUE_CSV
    failed_csv_path = output_dir / FAILED_REQUESTS_CSV

    result_fields = [
        "unique_id", "poly_hash", "F_integral",
        "first_source_csv", "first_source_row_index",
        "first_construction_kind", "first_execution_relative_dir",
        "first_family_index", "first_curve_index", "first_family_label", "first_curve_label",
        "first_affine_point_count", "first_projective_lower_bound",
        "height_F_integral", "degree", "is_squarefree_integral_over_Q", "appearance_count",
        "xml_file", "request_status", "elapsed_local_seconds",
        "max_time", "max_input", "seed", "version", "time", "memory", "warning",
        "rankbounds_lb", "rankbounds_ub", "rank_direct", "rank_g",
        "finiteIndex", "proved", "ub2", "final_rank_printed",
        "rank_Jacobian", "rank_status", "rank_source", "assumption",
        "magma_error", "unresolved", "lines_text",
    ]

    failed_fields = [
        "unique_id", "poly_hash", "F_integral", "error", "elapsed_local_seconds",
    ]

    timeout = aiohttp.ClientTimeout(total=REQUEST_TIMEOUT_SECONDS)
    connector = aiohttp.TCPConnector(limit_per_host=CONCURRENCY, limit=CONCURRENCY)

    semaphore = asyncio.Semaphore(CONCURRENCY)
    rate_limiter = GlobalRateLimiter(GLOBAL_MIN_SECONDS_BETWEEN_REQUEST_STARTS)
    completed_counter = {"done": 0}
    offline_state = {"consecutive_offline": 0}
    total = len(unique_curves)

    async with aiohttp.ClientSession(timeout=timeout, connector=connector) as session:
        tasks = [
            process_one_curve(
                session=session,
                semaphore=semaphore,
                unique_curve=curve,
                result_fields=result_fields,
                failed_fields=failed_fields,
                results_csv_path=results_csv_path,
                failed_csv_path=failed_csv_path,
                xml_dir=xml_dir,
                bad_cache_dir=bad_cache_dir,
                offline_dir=offline_dir,
                completed_counter=completed_counter,
                total=total,
                offline_state=offline_state,
                rate_limiter=rate_limiter,
            )
            for curve in unique_curves
        ]

        results = await asyncio.gather(*tasks)

    return list(results)


# ============================================================
# JOIN RESULTS BACK TO ALL APPEARANCES
# ============================================================

def results_by_unique_id(results: List[Dict[str, object]]) -> Dict[str, Dict[str, object]]:
    return {str(r.get("unique_id", "")): r for r in results if r.get("unique_id", "")}


def write_all_appearances_results(
    output_dir: Path,
    appearances: List[Dict[str, str]],
    unique_results: List[Dict[str, object]],
) -> None:
    by_id = results_by_unique_id(unique_results)

    rows = []

    for app in appearances:
        uid = app["unique_id"]
        result = by_id.get(uid, {})

        row = {
            **app,
            "rank_Jacobian": result.get("rank_Jacobian", ""),
            "rank_status": result.get("rank_status", ""),
            "rank_source": result.get("rank_source", ""),
            "rankbounds_lb": result.get("rankbounds_lb", ""),
            "rankbounds_ub": result.get("rankbounds_ub", ""),
            "rank_direct": result.get("rank_direct", ""),
            "rank_g": result.get("rank_g", ""),
            "finiteIndex": result.get("finiteIndex", ""),
            "proved": result.get("proved", ""),
            "ub2": result.get("ub2", ""),
            "assumption": result.get("assumption", ""),
            "magma_time": result.get("time", ""),
            "memory": result.get("memory", ""),
            "warning": result.get("warning", ""),
            "xml_file": result.get("xml_file", ""),
            "request_status": result.get("request_status", ""),
        }
        rows.append(row)

    fieldnames = [
        "unique_id", "poly_hash", "_source_csv", "_source_row_index",
        "construction_kind", "execution_relative_dir",
        "family_index", "curve_index", "family_label", "curve_label",
        "affine_point_count", "projective_lower_bound", "extra_x_count",
        "F_integral",
        "rank_Jacobian", "rank_status", "rank_source",
        "rankbounds_lb", "rankbounds_ub", "rank_direct", "rank_g",
        "finiteIndex", "proved", "ub2", "assumption",
        "magma_time", "memory", "warning", "xml_file", "request_status",
    ]

    write_csv(output_dir / RESULTS_ALL_APPEARANCES_CSV, rows, fieldnames)


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    ensure_dir(OUTPUT_DIR)
    clean_output_csvs(OUTPUT_DIR)

    print("=" * 100)
    print("ROBUST V3 MAGMA RANK REQUESTS FROM INTEGRAL MODELS")
    print("=" * 100)
    print(f"Output dir: {OUTPUT_DIR}")
    print(f"Concurrency: {CONCURRENCY}")
    print(f"Global min seconds between request starts: {GLOBAL_MIN_SECONDS_BETWEEN_REQUEST_STARTS}s")
    print(f"Transient retry sleep: {TRANSIENT_RETRY_SLEEP_SECONDS}s")
    print(f"Max transient reruns per curve: {MAX_TRANSIENT_RERUNS_PER_CURVE}")
    print(f"Clean CSVs on start: {CLEAN_OUTPUT_CSVS_ON_START}")
    print(f"Use GRH: {USE_GRH}")
    print(f"Rerun cache statuses: {sorted(RERUN_CACHE_STATUSES)}")
    print(f"Transient bad statuses: {sorted(TRANSIENT_BAD_STATUSES)}")
    print("Input CSVs:")
    for p in INPUT_CSVS:
        print(f"  {p}")
    print("=" * 100)

    unique_curves, appearances = load_and_deduplicate(INPUT_CSVS)

    if PROCESS_ONLY_FIRST_N_UNIQUE_CURVES is not None:
        unique_curves = unique_curves[:PROCESS_ONLY_FIRST_N_UNIQUE_CURVES]

    print(f"Total appearances: {len(appearances)}")
    print(f"Unique curves by exact F_integral: {len(unique_curves)}")

    unique_fields = [
        "unique_id", "poly_hash", "F_integral",
        "first_source_csv", "first_source_row_index",
        "first_construction_kind", "first_execution_relative_dir",
        "first_family_index", "first_curve_index", "first_family_label", "first_curve_label",
        "first_affine_point_count", "first_projective_lower_bound",
        "height_F_integral", "degree", "is_squarefree_integral_over_Q", "appearance_count",
    ]

    appearance_fields = [
        "unique_id", "poly_hash", "_source_csv", "_source_row_index",
        "construction_kind", "execution_relative_dir",
        "family_index", "curve_index", "family_label", "curve_label",
        "affine_point_count", "projective_lower_bound", "extra_x_count", "F_integral",
    ]

    write_csv(OUTPUT_DIR / UNIQUE_CURVES_CSV, unique_curves, unique_fields)
    write_csv(OUTPUT_DIR / APPEARANCES_CSV, appearances, appearance_fields)

    start = time.time()
    results = asyncio.run(run_requests(unique_curves, OUTPUT_DIR))
    elapsed = time.time() - start

    write_all_appearances_results(OUTPUT_DIR, appearances, results)

    print("=" * 100)
    print("DONE")
    print(f"Total unique curves processed: {len(results)}")
    print(f"Elapsed local time: {elapsed:.1f}s")
    print(f"Unique results CSV: {OUTPUT_DIR / RESULTS_UNIQUE_CSV}")
    print(f"All appearances CSV: {OUTPUT_DIR / RESULTS_ALL_APPEARANCES_CSV}")
    print("=" * 100)


if __name__ == "__main__":
    main()
