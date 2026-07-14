#!/usr/bin/env python3
"""
Audit (and optionally fix) GRCh37 positions of the curated SNP panel.

Implements plan task 3.3 / 1.2.A: every rsID in the SNP database is resolved
against Ensembl GRCh37 (batch REST) and compared to the chrom/pos stored in the
CSV. Wrong positions make a SNP silently fail to match in PLINK --score, so the
corresponding trait drops out of the report — a correctness bug, not a nicety.

Usage:
    python audit_snp_positions.py                       # dry-run report
    python audit_snp_positions.py --json report.json    # + machine-readable report
    python audit_snp_positions.py --fix                 # rewrite CSV (keeps .bak)

Default input: prs_research_pipeline/data/snp_database_annotated.csv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ENSEMBL_GRCH37 = "https://grch37.rest.ensembl.org/variation/human"
BATCH_SIZE = 190          # Ensembl POST limit is 200 ids
CANONICAL = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}

DEFAULT_CSV = (
    Path(__file__).resolve().parents[2] / "data" / "snp_database_annotated.csv"
)


def norm_chrom(c: str) -> str:
    """'chr15' -> '15', 'chrMT' -> 'MT'."""
    c = (c or "").strip()
    if c.lower().startswith("chr"):
        c = c[3:]
    return "MT" if c.upper() in {"M", "MT"} else c.upper()


def ensembl_batch(ids: list[str]) -> dict:
    """POST a batch of rsIDs to Ensembl GRCh37; returns {rsid: variant_json}."""
    body = json.dumps({"ids": ids}).encode()
    req = urllib.request.Request(
        ENSEMBL_GRCH37,
        data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=45) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 429:  # rate limited
                wait = int(e.headers.get("Retry-After", 2 * (attempt + 1)))
                time.sleep(wait)
                continue
            raise
        except (urllib.error.URLError, TimeoutError):
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("Ensembl batch request failed after retries")


def grch37_mapping(variant: dict) -> tuple[str, int] | None:
    """Extract the canonical-chromosome GRCh37 (chrom, start) from a variant."""
    best = None
    for m in variant.get("mappings", []):
        if m.get("assembly_name") != "GRCh37":
            continue
        chrom = norm_chrom(m.get("seq_region_name", ""))
        if chrom in CANONICAL:  # prefer real chromosomes over patches/haplotypes
            return chrom, int(m["start"])
        best = best or (chrom, int(m["start"]))
    return best


def load_rows(path: Path) -> tuple[list[str], list[dict]]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        return reader.fieldnames, list(reader)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--json", type=Path, help="write machine-readable report here")
    ap.add_argument("--fix", action="store_true",
                    help="rewrite CSV with corrected positions (backup to .bak)")
    args = ap.parse_args()

    if not args.csv.exists():
        print(f"ERROR: CSV not found: {args.csv}", file=sys.stderr)
        return 1

    fieldnames, rows = load_rows(args.csv)
    rsids = sorted({r["rsid"].strip() for r in rows if r.get("rsid", "").strip()})
    print(f"Auditing {len(rsids)} unique rsIDs from {args.csv.name} "
          f"({len(rows)} rows) against Ensembl GRCh37...\n")

    resolved: dict[str, tuple[str, int] | None] = {}
    for i in range(0, len(rsids), BATCH_SIZE):
        chunk = rsids[i:i + BATCH_SIZE]
        data = ensembl_batch(chunk)
        for rsid in chunk:
            v = data.get(rsid)
            resolved[rsid] = grch37_mapping(v) if v else None
        print(f"  resolved {min(i + BATCH_SIZE, len(rsids))}/{len(rsids)}")
        time.sleep(0.3)  # be polite

    # Two classes of problem:
    #  - pos_error:  same chromosome, different position -> safe to auto-fix.
    #  - chrom_error: Ensembl maps the rsID to a DIFFERENT chromosome than the CSV
    #                 gene sits on. The rsID (or gene/trait label) is almost
    #                 certainly wrong; fixing only the position would relocate the
    #                 SNP into the wrong gene. NEVER auto-fix — flag for review.
    matches, pos_errors, chrom_errors, not_found = [], [], [], []
    for r in rows:
        rsid = r["rsid"].strip()
        csv_chrom, csv_pos = norm_chrom(r.get("chrom", "")), r.get("pos", "").strip()
        ref = resolved.get(rsid)
        if ref is None:
            not_found.append(rsid)
            continue
        ens_chrom, ens_pos = ref
        rec = {
            "rsid": rsid, "gene": r.get("gene", ""), "trait": r.get("trait_category", ""),
            "csv": f"chr{csv_chrom}:{csv_pos}", "ensembl_grch37": f"chr{ens_chrom}:{ens_pos}",
            "correct_chrom": f"chr{ens_chrom}", "correct_pos": str(ens_pos),
        }
        if csv_chrom != ens_chrom:
            chrom_errors.append(rec)
        else:
            try:
                same = int(csv_pos) == ens_pos
            except ValueError:
                same = False
            (matches if same else pos_errors).append(rec if not same else rsid)

    print("\n" + "=" * 78)
    print(f"  MATCH: {len(matches)}   "
          f"POS ERROR (same chrom, fixable): {len(pos_errors)}   "
          f"CHROM ERROR (needs review): {len(chrom_errors)}   "
          f"NOT FOUND: {len(not_found)}")
    print("=" * 78)

    if pos_errors:
        print("\nPOSITION ERRORS — same chromosome, wrong coordinate (auto-fixable):\n")
        print(f"  {'rsID':<14}{'gene':<10}{'CSV':<20}{'CORRECT (GRCh37)':<20}trait")
        for m in pos_errors:
            print(f"  {m['rsid']:<14}{m['gene']:<10}{m['csv']:<20}"
                  f"{m['ensembl_grch37']:<20}{m['trait']}")
    if chrom_errors:
        print("\n⚠  CHROMOSOME MISMATCHES — rsID maps to a DIFFERENT chromosome than "
              "the CSV gene.\n   The rsID or the gene/trait label is wrong. MANUAL "
              "review required (not auto-fixed):\n")
        print(f"  {'rsID':<14}{'gene':<10}{'CSV':<20}{'ENSEMBL rsID maps to':<22}trait")
        for m in chrom_errors:
            print(f"  {m['rsid']:<14}{m['gene']:<10}{m['csv']:<20}"
                  f"{m['ensembl_grch37']:<22}{m['trait']}")
    if not_found:
        print(f"\nNOT FOUND in Ensembl GRCh37 (merged/withdrawn rsIDs or indels): "
              f"{', '.join(not_found)}")

    if args.json:
        report = {"csv": str(args.csv), "n_rows": len(rows), "n_rsids": len(rsids),
                  "matches": len(matches), "pos_errors": pos_errors,
                  "chrom_errors": chrom_errors, "not_found": not_found}
        args.json.write_text(json.dumps(report, indent=2))
        print(f"\nReport written to {args.json}")

    if args.fix and pos_errors:
        backup = args.csv.with_suffix(args.csv.suffix + ".bak")
        backup.write_bytes(args.csv.read_bytes())
        fix = {m["rsid"]: (m["correct_chrom"], m["correct_pos"]) for m in pos_errors}
        for r in rows:
            if r["rsid"].strip() in fix:
                r["chrom"], r["pos"] = fix[r["rsid"].strip()]
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)
        print(f"\nFIXED {len(fix)} position error(s). Backup: {backup.name}")
        if chrom_errors:
            print(f"⚠  {len(chrom_errors)} chromosome mismatches were NOT fixed — "
                  "review manually.")

    return 1 if (pos_errors or chrom_errors) else 0


if __name__ == "__main__":
    sys.exit(main())
