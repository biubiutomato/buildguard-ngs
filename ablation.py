#!/usr/bin/env python3
"""Reproduce the evidence-source ablation reported in Table II."""
from __future__ import annotations

import csv
import shutil
from pathlib import Path
import tempfile
from collections import Counter

from detector import (
    FileCall, _naming_of, _length_votes, _coord_oob, _landmark_votes,
    _mt_flavor, check_bundle,
)
from fingerprints import BUILDS, NCBI_TO_CHR
from parse import parse_any
from synthesize import FAULTS, make_bundle

ROOT = Path(tempfile.gettempdir()) / 'buildguard_ablation_corpus'
OUT = Path(__file__).resolve().parent / 'results' / 'ablation.csv'


def call_file_with_weights(prof, *, w_len=3.0, w_acc=4.0, w_lm=1.0,
                           w_tok=0.5, w_oob=2.0):
    naming = _naming_of(prof)
    lv, n_hits, lev = _length_votes(prof)
    lm, lmev = _landmark_votes(prof)
    scores = {b: 0.0 for b in BUILDS}
    warnings = []

    for b, v in lv.items():
        scores[b] += w_len * v
    for b, v in lm.items():
        scores[b] += w_lm * v
    if prof.header_mentions_build:
        scores[prof.header_mentions_build] += w_tok
    for c in prof.contigs:
        if c.raw_name in NCBI_TO_CHR:
            _, b = NCBI_TO_CHR[c.raw_name]
            scores[b] += w_acc
    for b in BUILDS:
        n_oob, _ = _coord_oob(prof, b)
        scores[b] -= w_oob * n_oob

    if prof.alt_like:
        warnings.append('alt/decoy/unplaced contig name(s)')
    styles = {c.style for c in prof.contigs if c.style in ('ucsc','ensembl','ncbi')}
    if naming == 'mixed' or len(styles) > 1:
        warnings.append('mixed contig naming styles in one file')
        naming = 'mixed' if len(styles) > 1 else naming
    if prof.missing_contig_lengths and prof.kind == 'vcf':
        warnings.append('VCF has no ##contig length lines; build call is weaker')

    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    best, best_s = ranked[0]
    second_s = ranked[1][1]
    build = best if best_s > 0 and best_s > second_s else (best if best_s >= 2 else None)
    if best_s > 0 and best_s == second_s:
        build = None
        warnings.append('ambiguous length/landmark votes')

    mt = _mt_flavor(prof)
    if mt == 'yoruba' and build in ('hg38','t2t'):
        warnings.append('chrM length 16571 with a GRCh38/T2T nuclear build')
    if mt == 'rcrs' and build == 'hg19':
        warnings.append('chrM length 16569 with an hg19 nuclear build')
    if build:
        n_oob_best, _ = _coord_oob(prof, build)
        if n_oob_best:
            warnings.append(f'{n_oob_best} coordinate(s) exceed {build} contig length')

    return FileCall(
        path=prof.path, kind=prof.kind, build=build, naming=naming,
        votes=scores, evidence=(lev + lmev)[:12], warnings=warnings,
        n_length_hits=n_hits, n_oob=0, mt_flavor=mt,
    )


def generate_cases():
    if ROOT.exists():
        shutil.rmtree(ROOT)
    ROOT.mkdir(parents=True)
    metas = []
    cid = 0
    for build in ('hg19','hg38','t2t'):
        for style in ('ucsc','ensembl','ncbi'):
            for fault in FAULTS:
                for rep in range(3):
                    cid += 1
                    metas.append(make_bundle(ROOT, cid, build, style, fault, rep))
    return metas


def evaluate(metas, weights):
    rows = []
    for m in metas:
        profs = [parse_any(p) for p in sorted(Path(m['dir']).glob('*'))]
        calls = [call_file_with_weights(p, **weights) for p in profs]
        report = check_bundle(calls)
        codes = {i['code'] for i in report.issues}
        expected = set(m['expected_codes'])
        hard = not m.get('soft')
        if m['fault'] == 'clean':
            detected = report.severity != 'error'
        elif hard:
            detected = bool(expected & codes) or (m['is_mismatch'] and report.severity == 'error')
        else:
            detected = True

        fasta = next((c.build for c in calls if c.kind == 'fasta'), None)
        rows.append((m, detected, fasta == m['true_build']))

    pos = [r for r in rows if r[0]['is_mismatch'] and not r[0].get('soft')]
    neg = [r for r in rows if r[0]['fault'] == 'clean']
    tp = sum(r[1] for r in pos)
    tn = sum(r[1] for r in neg)
    fn = len(pos) - tp
    fp = len(neg) - tn
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    fasta_acc = sum(r[2] for r in rows) / len(rows)
    return precision, recall, fasta_acc


def main():
    metas = generate_cases()
    base = dict(w_len=3.0, w_acc=4.0, w_lm=1.0, w_tok=0.5, w_oob=2.0)
    configs = [
        ('none (full score)', {}),
        ('header token', {'w_tok': 0}),
        ('accessions', {'w_acc': 0}),
        ('landmark genes', {'w_lm': 0}),
        ('contig lengths', {'w_len': 0}),
        ('lengths + accessions', {'w_len': 0, 'w_acc': 0}),
        ('lengths + accessions + landmarks', {'w_len': 0, 'w_acc': 0, 'w_lm': 0}),
    ]
    OUT.parent.mkdir(exist_ok=True)
    with OUT.open('w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['evidence_removed','precision','recall','fasta_build_acc'])
        for name, changes in configs:
            weights = base.copy(); weights.update(changes)
            p, r, a = evaluate(metas, weights)
            w.writerow([name, f'{p:.2f}', f'{r:.2f}', f'{a:.2f}'])
            print(f'{name:34s} P={p:.2f} R={r:.2f} FASTA={a:.2f}')
    print(f'Wrote {OUT}')


if __name__ == '__main__':
    main()
