#!/usr/bin/env python3
"""
Parse gnina log files into a single tabular dataset.

Output schema (default):
identifier, pose, affinity_kcal_mol, intramol_kcal_mol, cnn_pose_score, cnn_affinity, log_path

Usage:
  python parse_gnina_logs.py --logdir dock_logs --out docking_results.parquet
  python parse_gnina_logs.py --logdir dock_logs --out docking_results.csv --format csv
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import pandas as pd


# --- regex helpers ---

# Commandline: gnina ... --out dock_out/CNP0415838.1_out.sdf --log ...
RE_OUT = re.compile(r"(?:^|\s)--out\s+(?P<out>\S+)")
RE_LIG = re.compile(r"(?:^|\s)-l\s+(?P<lig>\S+)")
RE_CMDLINE = re.compile(r"^Commandline:\s*(?P<cmd>.*)$", re.MULTILINE)

# A table row looks like:
# 1       -4.84       -0.69       0.0920      5.869
RE_ROW = re.compile(
    r"^\s*(?P<mode>\d+)\s+"
    r"(?P<affinity>-?\d+(?:\.\d+)?)\s+"
    r"(?P<intramol>-?\d+(?:\.\d+)?)\s+"
    r"(?P<cnn_pose>-?\d+(?:\.\d+)?)\s+"
    r"(?P<cnn_aff>-?\d+(?:\.\d+)?)\s*$"
)

# table starts after the dashed separator line
RE_SEPARATOR = re.compile(r"^-{3,}\+?-{3,}.*$")


def derive_identifier_from_out(out_path: str) -> str:
    """
    Given an --out path like dock_out/CNP0415838.1_out.sdf,
    return identifier like CNP0415838.1 (best-effort).
    """
    base = Path(out_path).name  # CNP0415838.1_out.sdf
    # common pattern: <id>_out.sdf
    if base.endswith("_out.sdf"):
        return base[: -len("_out.sdf")]
    # other: <id>.sdf or <id>_something.sdf -> strip extension first
    stem = Path(base).stem
    # if it ends with _out, remove that
    if stem.endswith("_out"):
        stem = stem[: -len("_out")]
    return stem


def derive_identifier_from_lig(lig_path: str) -> str:
    """
    Given a ligand path like ligands_pdbqt/CNP0415838.1.pdbqt
    return identifier like CNP0415838.1.
    """
    return Path(lig_path).stem


def extract_identifier(text: str) -> Optional[str]:
    """
    Prefer identifier from --out; fallback to -l ligand.
    """
    m_cmd = RE_CMDLINE.search(text)
    cmd = m_cmd.group("cmd") if m_cmd else text  # fallback: search whole text

    m_out = RE_OUT.search(cmd)
    if m_out:
        return derive_identifier_from_out(m_out.group("out"))

    m_lig = RE_LIG.search(cmd)
    if m_lig:
        return derive_identifier_from_lig(m_lig.group("lig"))

    return None


def parse_gnina_table_rows(text: str) -> List[Dict]:
    """
    Locate the results table and parse row lines into dicts.
    Returns list of rows with mode + scores.
    """
    lines = text.splitlines()

    # find separator line, then parse subsequent numeric rows
    start_idx = None
    for i, line in enumerate(lines):
        if RE_SEPARATOR.match(line.strip()):
            start_idx = i + 1
            break

    if start_idx is None:
        return []

    rows: List[Dict] = []
    for line in lines[start_idx:]:
        line = line.rstrip("\n")
        m = RE_ROW.match(line)
        if not m:
            # stop when table ends (blank line or non-row content)
            if rows and (not line.strip() or "Writing" in line or "Done" in line):
                break
            continue

        rows.append(
            {
                "pose": int(m.group("mode")),
                "affinity_kcal_mol": float(m.group("affinity")),
                "intramol_kcal_mol": float(m.group("intramol")),
                "cnn_pose_score": float(m.group("cnn_pose")),
                "cnn_affinity": float(m.group("cnn_aff")),
            }
        )

    return rows


def iter_log_files(logdir: Path, pattern: str = "*.log") -> Iterator[Path]:
    yield from logdir.rglob(pattern)


def parse_one_log(path: Path) -> Tuple[Optional[str], List[Dict]]:
    text = path.read_text(errors="replace")
    identifier = extract_identifier(text)
    rows = parse_gnina_table_rows(text)
    return identifier, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--logdir", type=str, required=True, help="Directory containing gnina .log files")
    ap.add_argument("--out", type=str, required=True, help="Output file path (.parquet or .csv)")
    ap.add_argument("--format", type=str, default="parquet", choices=["parquet", "csv"], help="Output format")
    ap.add_argument("--pattern", type=str, default="*.log", help="Glob pattern to match log files (default: *.log)")
    ap.add_argument("--strict", action="store_true", help="If set, error on logs with missing id or zero rows")
    args = ap.parse_args()

    logdir = Path(args.logdir)
    if not logdir.exists():
        raise FileNotFoundError(f"logdir not found: {logdir}")

    all_rows: List[Dict] = []
    missing_id = 0
    missing_table = 0

    for p in iter_log_files(logdir, args.pattern):
        try:
            identifier, rows = parse_one_log(p)
        except Exception as e:
            # if a log is corrupted, either skip or die based on strict
            if args.strict:
                raise RuntimeError(f"Failed parsing {p}: {e}") from e
            continue

        if identifier is None:
            missing_id += 1
            if args.strict:
                raise ValueError(f"Could not extract identifier from: {p}")
            # still store rows with a placeholder id if they exist
            identifier = p.stem

        if not rows:
            missing_table += 1
            if args.strict:
                raise ValueError(f"No gnina table rows found in: {p}")
            continue

        for r in rows:
            r["identifier"] = identifier
            r["log_path"] = str(p)
            all_rows.append(r)

    df = pd.DataFrame(all_rows, columns=[
        "identifier",
        "pose",
        "affinity_kcal_mol",
        "intramol_kcal_mol",
        "cnn_pose_score",
        "cnn_affinity",
        "log_path",
    ])

    # nice deterministic ordering
    if not df.empty:
        df.sort_values(["identifier", "pose"], inplace=True, ignore_index=True)

    out_path = Path(args.out)

    if args.format == "parquet":
        # uses pyarrow if installed; fast + compact
        df.to_parquet(out_path, index=False)
    else:
        df.to_csv(out_path, index=False)

    print(f"parsed rows: {len(df):,}")
    print(f"unique identifiers: {df['identifier'].nunique() if not df.empty else 0:,}")
    print(f"logs missing identifier: {missing_id:,}")
    print(f"logs missing table rows: {missing_table:,}")
    print(f"wrote: {out_path}")


if __name__ == "__main__":
    main()
