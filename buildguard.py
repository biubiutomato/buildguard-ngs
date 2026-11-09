#!/usr/bin/env python3
"""CLI: inspect FASTA / GTF / VCF files for genome-build and contig-name mismatch."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from detector import call_file, check_bundle
from parse import parse_any


def main() -> int:
    ap = argparse.ArgumentParser(description="Detect genome-build / contig naming mismatches")
    ap.add_argument("files", nargs="+", help="FASTA, GTF/GFF, or VCF")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args()
    calls = []
    for f in args.files:
        p = Path(f)
        if not p.exists():
            print(f"missing: {p}", file=sys.stderr)
            return 2
        calls.append(call_file(parse_any(p)))
    rep = check_bundle(calls)
    if args.json:
        print(json.dumps({
            "severity": rep.severity,
            "predicted_build": rep.predicted_build,
            "consistent": rep.consistent,
            "files": [
                {"path": c.path, "kind": c.kind, "build": c.build,
                 "naming": c.naming, "votes": c.votes, "warnings": c.warnings,
                 "evidence": c.evidence}
                for c in calls
            ],
            "issues": rep.issues,
        }, indent=2))
        return 0 if rep.severity != "error" else 1
    print(f"bundle: {rep.severity}  predicted_build={rep.predicted_build}")
    for c in calls:
        print(f"  {c.kind:5}  build={c.build or '?':5}  naming={c.naming:8}  {Path(c.path).name}")
        for e in c.evidence[:4]:
            print(f"         evidence: {e}")
        for w in c.warnings:
            print(f"         warn: {w}")
    for i in rep.issues:
        print(f"  [{i['severity']}] {i['code']}: {i['detail']}")
    return 0 if rep.severity != "error" else 1


if __name__ == "__main__":
    raise SystemExit(main())
