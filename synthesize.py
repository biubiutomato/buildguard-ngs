"""Synthesize miniature FASTA / GTF / VCF files and inject controlled faults."""

from __future__ import annotations

import random
import zlib

# Replicate index; set by make_bundle so that replicates use different seeds.
_REP = 0


def _seed(*parts) -> int:
    """Deterministic seed (independent of Python's per-process hash randomization)."""
    return zlib.crc32("|".join(map(str, parts + (_REP,))).encode())
from pathlib import Path

from fingerprints import BUILD_LEN, NCBI_TO_CHR

# reverse NCBI map: (build, chrom) -> accession
CHR_TO_NCBI: dict[tuple[str, str], str] = {}
for acc, (chrom, build) in NCBI_TO_CHR.items():
    CHR_TO_NCBI[(build, chrom)] = acc

PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# short dummy sequence; length is declared in header, body is truncated for speed
BODY = "ACGTN" * 40  # 200 bp


def _name(chrom: str, style: str, build: str) -> str:
    bare = chrom.replace("chr", "")
    if style == "ucsc":
        return chrom
    if style == "ensembl":
        return "MT" if chrom == "chrM" else bare
    if style == "ncbi":
        acc = CHR_TO_NCBI.get((build, chrom))
        if acc:
            return acc
        if build == "t2t" and chrom == "chrM":
            return "CP068254.1"
        return chrom
    return chrom


def write_fasta(path: Path, build: str, style: str, chroms: list[str] | None = None,
                mention_build: bool = False, include_alt: bool = False,
                corrupt_length: dict[str, int] | None = None,
                drop_mt: bool = False) -> None:
    table = BUILD_LEN[build]
    chroms = list(chroms or PRIMARY)
    if drop_mt:
        chroms = [c for c in chroms if c != "chrM"]
    lines = []
    for c in chroms:
        raw = _name(c, style, build)
        ln = (corrupt_length or {}).get(c, table[c])
        extra = f" LN:{ln}"
        if mention_build:
            extra += f" build:{build}"
        lines.append(f">{raw}{extra}")
        lines.append(BODY)
    if include_alt:
        lines.append(">chr1_KI270706v1_random LN:43750")
        lines.append(BODY)
        lines.append(">chrUn_GL000220v1 LN:161802")
        lines.append(BODY)
    path.write_text("\n".join(lines) + "\n")


def write_gtf(path: Path, build: str, style: str, n_genes: int = 40,
              mention_build: bool = False, shift_coords: int = 0,
              mix_style: bool = False, oob: bool = False) -> None:
    rng = random.Random(_seed(build, style, n_genes))
    table = BUILD_LEN[build]
    genes = ["TP53", "BRCA1", "EGFR", "MYC", "BRCA2", "KRAS", "PIK3CA", "PTEN"]
    # approximate gene windows per build (reuse landmarks + jitter)
    from fingerprints import LANDMARKS
    rows = []
    if mention_build:
        rows.append(f"#!genome-build {build}")
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    for i in range(n_genes):
        if i < len(genes) and genes[i] in LANDMARKS[build]:
            chrom, lo, hi = LANDMARKS[build][genes[i]]
            gname = genes[i]
        else:
            chrom = rng.choice(chroms)
            lo = rng.randint(10_000, max(20_000, table[chrom] // 4))
            hi = lo + rng.randint(800, 4000)
            gname = f"GENE{i:03d}"
        start = max(1, lo + shift_coords)
        end = hi + shift_coords
        if oob and i == 0:
            start = table[chrom] + 500
            end = start + 200
        st = style
        if mix_style and i % 5 == 0:
            st = "ensembl" if style == "ucsc" else "ucsc"
        raw = _name(chrom, st, build)
        attr = f'gene_id "{gname}"; gene_name "{gname}"; transcript_id "{gname}.1";'
        rows.append(f"{raw}\thavana\tgene\t{start}\t{end}\t.\t+\t.\t{attr}")
        rows.append(f"{raw}\thavana\texon\t{start}\t{min(end, start+300)}\t.\t+\t.\t{attr}")
    path.write_text("\n".join(rows) + "\n")


def write_vcf(path: Path, build: str, style: str, n_var: int = 30,
              with_contig: bool = True, mention_build: bool = False,
              drop_length: bool = False, include_alt: bool = False,
              shift_coords: int = 0, oob: bool = False,
              mt_swap: bool = False) -> None:
    rng = random.Random(_seed(build, style, n_var))
    table = dict(BUILD_LEN[build])
    if mt_swap:
        table["chrM"] = 16571 if table["chrM"] == 16569 else 16569
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrM"]
    rows = ["##fileformat=VCFv4.2"]
    if mention_build:
        rows.append(f"##reference={build}")
    if with_contig:
        for c in chroms:
            raw = _name(c, style, build)
            if drop_length:
                rows.append(f"##contig=<ID={raw}>")
            else:
                rows.append(f"##contig=<ID={raw},length={table[c]}>")
        if include_alt:
            rows.append("##contig=<ID=chr1_KI270706v1_random,length=43750>")
    rows.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for i in range(n_var):
        chrom = chroms[i % len(chroms)]
        lim = table[chrom]
        pos = 100 + (i * 17) % max(100, lim // 20)
        pos = max(1, pos + shift_coords)
        if oob and i == 0:
            pos = lim + 50
        raw = _name(chrom, style, build)
        rows.append(f"{raw}\t{pos}\trs{i}\tA\tG\t60\tPASS\t.")
    path.write_text("\n".join(rows) + "\n")


FAULTS = (
    "clean",
    "wrong_build_gtf",       # GTF from another build
    "wrong_build_vcf",
    "naming_mismatch",       # FASTA ucsc, GTF ensembl
    "ncbi_vs_ucsc",
    "missing_vcf_lengths",
    "mt_flavor_swap",
    "coord_oob",
    "alt_in_fasta_only",
    "header_lie",            # header says hg38 but lengths are hg19
    "mixed_names_in_gtf",
    "shifted_gtf",           # large coord shift mimicking wrong build without length
)


def make_bundle(root: Path, case_id: int, build: str, style: str, fault: str, rep: int = 0) -> dict:
    global _REP
    _REP = rep
    d = root / f"case_{case_id:04d}_{build}_{style}_{fault}"
    d.mkdir(parents=True, exist_ok=True)
    other = {"hg19": "hg38", "hg38": "t2t", "t2t": "hg19"}[build]
    other2 = {"hg19": "t2t", "hg38": "hg19", "t2t": "hg38"}[build]
    alt_style = {"ucsc": "ensembl", "ensembl": "ucsc", "ncbi": "ucsc"}[style]

    meta = {
        "id": case_id,
        "true_build": build,
        "true_style": style,
        "fault": fault,
        "dir": str(d),
        "is_mismatch": fault != "clean",
        "expected_codes": [],
    }

    if fault == "clean":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", build, style, mention_build=True)
        write_vcf(d / "vars.vcf", build, style, mention_build=True)
        meta["expected_codes"] = []
    elif fault == "wrong_build_gtf":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", other, style)
        write_vcf(d / "vars.vcf", build, style)
        meta["expected_codes"] = ["BUILD_MISMATCH"]
    elif fault == "wrong_build_vcf":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", build, style)
        write_vcf(d / "vars.vcf", other2, style)
        meta["expected_codes"] = ["BUILD_MISMATCH"]
    elif fault == "naming_mismatch":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", build, alt_style)
        write_vcf(d / "vars.vcf", build, style)
        meta["expected_codes"] = ["NAMING_MISMATCH"]
    elif fault == "ncbi_vs_ucsc":
        write_fasta(d / "ref.fa", build, "ncbi")
        write_gtf(d / "genes.gtf", build, "ucsc")
        write_vcf(d / "vars.vcf", build, "ensembl")
        meta["expected_codes"] = ["NAMING_MISMATCH"]
    elif fault == "missing_vcf_lengths":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", build, style)
        write_vcf(d / "vars.vcf", build, style, drop_length=True)
        meta["expected_codes"] = []  # warning only; still consistent
        meta["is_mismatch"] = False
        meta["soft"] = True
    elif fault == "mt_flavor_swap":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", build, style)
        write_vcf(d / "vars.vcf", build, style, mt_swap=True)
        meta["expected_codes"] = ["MT_FLAVOR_MISMATCH"]
        meta["soft"] = True
    elif fault == "coord_oob":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", build, style, oob=True)
        write_vcf(d / "vars.vcf", build, style, oob=True)
        meta["expected_codes"] = ["FILE_WARNING"]
        meta["soft"] = True
    elif fault == "alt_in_fasta_only":
        write_fasta(d / "ref.fa", build, style, include_alt=True)
        write_gtf(d / "genes.gtf", build, style)
        write_vcf(d / "vars.vcf", build, style)
        meta["expected_codes"] = ["FILE_WARNING"]
        meta["soft"] = True
    elif fault == "header_lie":
        # lengths tell the truth (other), header token lies (build)
        write_fasta(d / "ref.fa", other, style, mention_build=True)
        # rewrite header token
        txt = (d / "ref.fa").read_text().replace(f"build:{other}", f"build:{build}")
        (d / "ref.fa").write_text(txt)
        write_gtf(d / "genes.gtf", other, style)
        write_vcf(d / "vars.vcf", other, style, mention_build=True)
        meta["true_build"] = other  # files are actually `other`
        meta["expected_codes"] = []
        meta["is_mismatch"] = False
        meta["soft"] = True
        meta["header_lie"] = True
    elif fault == "mixed_names_in_gtf":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", build, style, mix_style=True)
        write_vcf(d / "vars.vcf", build, style)
        meta["expected_codes"] = ["FILE_WARNING"]
        meta["soft"] = True
    elif fault == "shifted_gtf":
        write_fasta(d / "ref.fa", build, style)
        write_gtf(d / "genes.gtf", build, style, shift_coords=2_000_000)
        write_vcf(d / "vars.vcf", build, style)
        meta["expected_codes"] = ["BUILD_UNKNOWN", "FILE_WARNING", "BUILD_MISMATCH"]
        meta["soft"] = True
    else:
        raise ValueError(fault)
    return meta
