#!/usr/bin/env python3
"""Run BuildGuard on every combination of the public reference files in real_cases/."""
from __future__ import annotations

import csv
import itertools
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
from parse import parse_any
from detector import call_file, check_bundle

R = ROOT / "real_cases"
FASTA = {"hg19": R / "hg19.ucsc.fa", "hg38": R / "hg38.ucsc.fa", "t2t": R / "t2t.hs1.fa"}
GTF = {"hg19": R / "gencode.v19.hg19.gtf", "hg38": R / "gencode.v44.hg38.gtf"}
VCF = {"hg38": R / "1kGP.hg38.chr21.header.vcf", None: None}


def strip_header(src: Path, dst: Path) -> Path:
    """Copy a GTF without its '#' header lines, removing any free-text build label."""
    dst.write_text("".join(l for l in src.read_text().splitlines(True) if not l.startswith("#")))
    return dst


def main() -> None:
    tmp = Path(tempfile.mkdtemp(prefix="buildguard_real_"))
    rows = []
    for (fb, fa), (gb, gtf), (vb, vcf), stripped in itertools.product(FASTA.items(), GTF.items(), VCF.items(), [False, True]):
        g = strip_header(gtf, tmp / f"{gb}.noheader.gtf") if stripped else gtf
        files = [fa, g] + ([vcf] if vcf else [])
        calls = [call_file(parse_any(p)) for p in files]
        rep = check_bundle(calls)
        builds = {fb, gb} | ({vb} if vb else set())
        expected = "error" if len(builds) > 1 else "no error"
        observed = "error" if rep.severity == "error" else "no error"
        rows.append({
            "fasta": fb, "gtf": gb, "vcf": vb or "-", "gtf_header": "stripped" if stripped else "kept",
            "fasta_call": calls[0].build, "gtf_call": calls[1].build or "unknown",
            "vcf_call": (calls[2].build or "unknown") if vcf else "-",
            "severity": rep.severity, "expected": expected, "correct": observed == expected,
            "codes": "|".join(sorted({i["code"] for i in rep.issues})),
        })
    out = ROOT / "results" / "real_bundles.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['fasta']:5} {r['gtf']:5} {r['vcf']:5} {r['gtf_header']:8} -> fasta={r['fasta_call']} gtf={r['gtf_call']} vcf={r['vcf_call']} "
              f"{r['severity']:7} expected={r['expected']:8} {'OK' if r['correct'] else 'WRONG'} {r['codes']}")
    print(f"correct: {sum(r['correct'] for r in rows)}/{len(rows)}")


if __name__ == "__main__":
    main()
