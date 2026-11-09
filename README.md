# BuildGuard

Code and data for:

> Ye Ye. **BuildGuard: Detecting Genome-Build and Contig-Naming Mismatches in NGS Input Files.**
> IEEE BIBM 2026 Workshops (AI-BioNGS, paper S15206).

BuildGuard is a lightweight pre-flight checker that reads a FASTA, a GTF/GFF, and a VCF and reports whether they
agree on the human genome build (hg19 / GRCh37, hg38 / GRCh38, T2T-CHM13v2.0) and on the contig-naming style
(UCSC, Ensembl, NCBI/RefSeq). It uses primary-contig lengths, RefSeq accessions, naming patterns, and four
landmark genes, and every call is reported with the evidence behind it. It does not convert coordinates.

## Usage

```bash
python buildguard.py ref.fa genes.gtf vars.vcf
python buildguard.py --json ref.fa genes.gtf vars.vcf
```

The exit code is 1 when the bundle severity is `error` (a build or naming mismatch), so the check can stop a workflow.
Outcomes: `error` = cross-file mismatch; `warning` = a hazard that does not contradict the build
(missing VCF lengths, alt/decoy contigs, out-of-range coordinates, mitochondrial-length disagreement) or a file
whose build is `unknown`; `ok` = no issue.

## Reproduce the paper

No dependencies beyond the Python standard library (Python ≥ 3.9).

```bash
python experiment.py    # synthetic corpus (324 bundles) -> results/summary.json, fault_table.*, case_results.csv  (Table II)
python ablation.py      # evidence ablation -> results/ablation.csv                                                  (Table III)
python shift_sweep.py   # shifted-GTF sweep -> results/shift_sweep.csv                                               (Table I)
python real_bundles.py  # 24 bundles from public files -> results/real_bundles.csv
```

Synthetic corpora are written to the system temporary directory. Generation is deterministic.

## Files

| Path | Contents |
|---|---|
| `buildguard.py` | Command-line entry point |
| `parse.py` | Streaming FASTA / GTF / VCF parsers |
| `fingerprints.py` | Contig lengths, accessions, and landmark windows per build |
| `detector.py` | Per-file build scoring and bundle-level checks |
| `synthesize.py` | Fault-injection corpus generator |
| `real_cases/` | Public metadata extracts: UCSC `chrom.sizes` for hg19/hg38/hs1 written as FASTA headers, GENCODE v19 and v44 landmark-gene rows, and the 1000 Genomes 30× GRCh38 `chr21` VCF header |
| `results/` | Outputs of the scripts above |

## Known limitations

- The synthetic generator and the detector share the same length and landmark tables; perfect scores on the
  synthetic hard classes show the rules work as designed, not real-world accuracy.
- The T2T-CHM13 landmark windows are approximate (derived from hg38 positions).
- Uniform annotation shifts smaller than the ±50 kb landmark window, or in annotations named by RefSeq accession,
  are not detected (see `shift_sweep.py`).

## Citation

```bibtex
@inproceedings{ye2026buildguard,
  author    = {Ye, Ye},
  title     = {BuildGuard: Detecting Genome-Build and Contig-Naming Mismatches in NGS Input Files},
  booktitle = {Proc. IEEE International Conference on Bioinformatics and Biomedicine (BIBM) Workshops},
  year      = {2026}
}
```

## License

Code: MIT (see `LICENSE`). Files in `real_cases/` are excerpts of UCSC, GENCODE, and 1000 Genomes resources and
remain subject to their providers' terms.
