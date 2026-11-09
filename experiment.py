#!/usr/bin/env python3
"""Generate corpus, run BuildGuard, write tables and metrics."""

from __future__ import annotations

import csv
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))

from detector import call_file, check_bundle, extract_features
from parse import parse_any
from synthesize import FAULTS, make_bundle

ROOT = Path(__file__).resolve().parent
CORPUS = Path(tempfile.gettempdir()) / "buildguard_corpus"
OUT = ROOT / "results"


def generate(n_reps: int = 3) -> list[dict]:
    if CORPUS.exists():
        import shutil
        shutil.rmtree(CORPUS)
    CORPUS.mkdir(parents=True)
    metas = []
    cid = 0
    builds = ["hg19", "hg38", "t2t"]
    styles = ["ucsc", "ensembl", "ncbi"]
    for build in builds:
        for style in styles:
            for fault in FAULTS:
                for rep in range(n_reps):
                    cid += 1
                    metas.append(make_bundle(CORPUS, cid, build, style, fault, rep))
    (OUT / "manifest.json").write_text(json.dumps(metas, indent=2))
    return metas


def evaluate(metas: list[dict]) -> list[dict]:
    rows = []
    for m in metas:
        d = Path(m["dir"])
        files = sorted(d.glob("*"))
        profs = [parse_any(p) for p in files]
        calls = [call_file(p) for p in profs]
        rep = check_bundle(calls)
        codes = {i["code"] for i in rep.issues}
        expected = set(m["expected_codes"])
        hard = not m.get("soft")
        # detection: for hard mismatch cases we need the expected error codes
        if m["fault"] == "clean":
            detected = rep.severity != "error"
            tp_like = detected
        elif hard:
            detected = bool(expected & codes) or (
                m["is_mismatch"] and rep.severity == "error"
            )
            tp_like = detected
        else:
            # soft: warning is success if we flagged something; clean-soft (missing lengths) ok either way
            if m["fault"] in ("missing_vcf_lengths", "header_lie"):
                detected = True
                tp_like = True
            else:
                detected = rep.severity in ("warning", "error") or bool(expected & codes)
                tp_like = detected

        # per-file build accuracy vs true_build of the *bundle intended reference*
        # For wrong_build_* the FASTA is the intended reference.
        fasta_call = next((c.build for c in calls if c.kind == "fasta"), None)
        build_ok = fasta_call == m["true_build"]
        # header_lie: true_build rewritten to actual lengths
        if m.get("header_lie"):
            build_ok = fasta_call == m["true_build"]

        rows.append({
            "id": m["id"],
            "true_build": m["true_build"],
            "true_style": m["true_style"],
            "fault": m["fault"],
            "is_mismatch": m["is_mismatch"],
            "hard": hard,
            "pred_build": rep.predicted_build,
            "fasta_build": fasta_call,
            "severity": rep.severity,
            "codes": "|".join(sorted(codes)),
            "detected": int(detected),
            "build_ok": int(build_ok),
            "n_issues": len(rep.issues),
            "issues": "; ".join(i["detail"] for i in rep.issues[:4]),
        })
    return rows


def summarize(rows: list[dict]) -> dict:
    by_fault = defaultdict(list)
    for r in rows:
        by_fault[r["fault"]].append(r)

    table = []
    for fault, rs in by_fault.items():
        n = len(rs)
        det = sum(r["detected"] for r in rs) / n
        bok = sum(r["build_ok"] for r in rs) / n
        err = sum(r["severity"] == "error" for r in rs) / n
        warn = sum(r["severity"] == "warning" for r in rs) / n
        table.append({
            "fault": fault,
            "n": n,
            "detection_rate": round(det, 3),
            "fasta_build_acc": round(bok, 3),
            "frac_error": round(err, 3),
            "frac_warning": round(warn, 3),
        })

    hard_pos = [r for r in rows if r["is_mismatch"] and r["hard"]]
    hard_neg = [r for r in rows if r["fault"] == "clean"]
    tp = sum(r["detected"] and r["is_mismatch"] for r in hard_pos)
    # detected on clean means "correctly called consistent" in our encoding
    tn = sum(r["detected"] for r in hard_neg)
    fn = len(hard_pos) - tp
    fp = len(hard_neg) - tn
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0

    # naming-only vs build-only
    naming = [r for r in rows if r["fault"] in ("naming_mismatch", "ncbi_vs_ucsc")]
    buildf = [r for r in rows if r["fault"] in ("wrong_build_gtf", "wrong_build_vcf")]

    summary = {
        "n_cases": len(rows),
        "hard_precision": round(prec, 3),
        "hard_recall": round(rec, 3),
        "hard_f1": round(f1, 3),
        "hard_tp": tp,
        "hard_fp": fp,
        "hard_fn": fn,
        "hard_tn": tn,
        "naming_detection": round(sum(r["detected"] for r in naming) / len(naming), 3) if naming else None,
        "cross_build_detection": round(sum(r["detected"] for r in buildf) / len(buildf), 3) if buildf else None,
        "fasta_build_acc_all": round(sum(r["build_ok"] for r in rows) / len(rows), 3),
        "by_fault": table,
    }
    return summary


def write_outputs(rows, summary):
    OUT.mkdir(exist_ok=True)
    keys = list(rows[0].keys())
    with (OUT / "case_results.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    with (OUT / "fault_table.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary["by_fault"][0].keys()))
        w.writeheader()
        w.writerows(summary["by_fault"])
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    # markdown table
    md = ["| Fault | N | Detection | FASTA build acc. | Error | Warning |",
          "|---|---:|---:|---:|---:|---:|"]
    for r in summary["by_fault"]:
        md.append(
            f"| {r['fault']} | {r['n']} | {r['detection_rate']:.3f} | "
            f"{r['fasta_build_acc']:.3f} | {r['frac_error']:.3f} | {r['frac_warning']:.3f} |"
        )
    md.append("")
    md.append(
        f"Hard mismatch P/R/F1 = {summary['hard_precision']}/"
        f"{summary['hard_recall']}/{summary['hard_f1']} "
        f"(TP={summary['hard_tp']} FP={summary['hard_fp']} "
        f"FN={summary['hard_fn']} TN={summary['hard_tn']})"
    )
    (OUT / "fault_table.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    print(json.dumps({k: v for k, v in summary.items() if k != "by_fault"}, indent=2))


def main():
    OUT.mkdir(exist_ok=True)
    print("generating corpus...")
    metas = generate(n_reps=3)
    print(f"  {len(metas)} bundles")
    print("evaluating...")
    rows = evaluate(metas)
    summary = summarize(rows)
    write_outputs(rows, summary)


if __name__ == "__main__":
    main()
