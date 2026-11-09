#!/usr/bin/env python3
"""Sweep the magnitude of a uniform GTF coordinate shift and record what BuildGuard reports."""
from __future__ import annotations

import csv
import shutil
import sys
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import synthesize
from synthesize import write_fasta, write_gtf, write_vcf
from parse import parse_any
from detector import call_file, check_bundle

SHIFTS = [0, 10_000, 40_000, 60_000, 100_000, 500_000, 1_000_000, 2_000_000, 5_000_000, 10_000_000]
BUILDS = ["hg19", "hg38", "t2t"]
STYLES = ["ucsc", "ensembl", "ncbi"]
WORK = Path(tempfile.gettempdir()) / "buildguard_shift_sweep"


def main() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)
    rows = []
    for shift in SHIFTS:
        for build in BUILDS:
            for style in STYLES:
                for rep in range(3):
                    synthesize._REP = rep
                    d = WORK / f"{shift}_{build}_{style}_{rep}"
                    d.mkdir(parents=True)
                    write_fasta(d / "ref.fa", build, style)
                    write_gtf(d / "genes.gtf", build, style, shift_coords=shift)
                    write_vcf(d / "vars.vcf", build, style)
                    calls = [call_file(parse_any(p)) for p in sorted(d.glob("*"))]
                    rep_ = check_bundle(calls)
                    gtf = next(c for c in calls if c.kind == "gtf")
                    gtf_call = gtf.build or "unknown"
                    rows.append({
                        "shift_bp": shift, "build": build, "style": style, "rep": rep,
                        "gtf_call": "correct" if gtf_call == build else ("unknown" if gtf_call == "unknown" else "other_build"),
                        "severity": rep_.severity,
                    })
    out = ROOT / "results" / "shift_sweep.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{'shift':>10}  n  GTF correct/unknown/other   bundle ok/warning/error   flagged")
    for shift in SHIFTS:
        sub = [r for r in rows if r["shift_bp"] == shift]
        g = Counter(r["gtf_call"] for r in sub)
        s = Counter(r["severity"] for r in sub)
        flagged = (s["warning"] + s["error"]) / len(sub)
        print(f"{shift:>10} {len(sub):2d}  {g['correct']:2d}/{g['unknown']:2d}/{g['other_build']:2d}"
              f"                  {s['ok']:2d}/{s['warning']:2d}/{s['error']:2d}            {flagged:.2f}")
    print("by build at 2 Mb:", Counter((r["build"], r["gtf_call"], r["severity"]) for r in rows if r["shift_bp"] == 2_000_000))


if __name__ == "__main__":
    main()
