"""Lightweight parsers for FASTA headers, GTF seqnames, and VCF contig metadata."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from fingerprints import GENBANK_TO_CHR, NCBI_TO_CHR


@dataclass
class Contig:
    raw_name: str
    canon: str  # chr1 / chrX / chrM / other
    length: int | None = None
    style: str = "unknown"  # ucsc | ensembl | ncbi | genbank | mixed | other
    extra: str = ""


@dataclass
class FileProfile:
    path: str
    kind: str  # fasta | gtf | gff | vcf
    contigs: list[Contig] = field(default_factory=list)
    n_records: int = 0
    has_chr_prefix: int = 0
    no_chr_prefix: int = 0
    ncbi_names: int = 0
    alt_like: int = 0
    header_mentions_build: str | None = None
    max_coord: dict[str, int] = field(default_factory=dict)
    gene_hits: dict[str, list[tuple[str, int]]] = field(default_factory=dict)
    missing_contig_lengths: bool = False
    notes: list[str] = field(default_factory=list)


_CHR_RE = re.compile(r"^(?:chr)?([0-9]{1,2}|[XxYy]|[Mm][Tt]?)$")
_BUILD_RE = re.compile(
    r"(hg19|hg38|hg37|grch37|grch38|b37|b38|hs37d5|t2t|chm13|hs1)",
    re.IGNORECASE,
)
_GENE_RE = re.compile(r'gene_name "([^"]+)"')


def canon_name(raw: str) -> tuple[str, str]:
    raw = raw.strip()
    if raw in NCBI_TO_CHR:
        return NCBI_TO_CHR[raw][0], "ncbi"
    if raw in GENBANK_TO_CHR:
        return GENBANK_TO_CHR[raw], "genbank"
    # strip versioned accession leftovers
    key = raw.split()[0]
    if key in NCBI_TO_CHR:
        return NCBI_TO_CHR[key][0], "ncbi"
    if key in GENBANK_TO_CHR:
        return GENBANK_TO_CHR[key], "genbank"
    m = _CHR_RE.match(key)
    if m:
        tok = m.group(1).upper()
        if tok in ("M", "MT"):
            chrom = "chrM"
        elif tok in ("X", "Y"):
            chrom = "chr" + tok
        else:
            chrom = "chr" + str(int(tok))
        style = "ucsc" if key.lower().startswith("chr") else "ensembl"
        return chrom, style
    return key, "other"


def _scan_build_token(text: str) -> str | None:
    m = _BUILD_RE.search(text)
    if not m:
        return None
    t = m.group(1).lower()
    if t in ("hg19", "hg37", "grch37", "b37", "hs37d5"):
        return "hg19"
    if t in ("hg38", "grch38", "b38"):
        return "hg38"
    if t in ("t2t", "chm13", "hs1"):
        return "t2t"
    return None


def parse_fasta(path: Path, max_headers: int = 500) -> FileProfile:
    prof = FileProfile(path=str(path), kind="fasta")
    headers: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                headers.append(line[1:].rstrip("\n"))
                if len(headers) >= max_headers:
                    break
    # second pass for lengths if file is tiny synthetic
    lengths: dict[str, int] = {}
    cur = None
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    lengths[cur] = n
                cur = line[1:].split()[0]
                n = 0
            else:
                n += len(line.strip())
        if cur is not None:
            lengths[cur] = n

    for h in headers:
        raw = h.split()[0]
        canon, style = canon_name(raw)
        body_ln = lengths.get(raw)
        ln = None
        m = re.search(r"LN:(\d+)", h)
        if m:
            ln = int(m.group(1))
        m2 = re.search(r"length[=:](\d+)", h, re.I)
        if m2 and ln is None:
            ln = int(m2.group(1))
        if ln is None:
            ln = body_ln
        c = Contig(raw_name=raw, canon=canon, length=ln, style=style, extra=h)
        prof.contigs.append(c)
        if style == "ucsc":
            prof.has_chr_prefix += 1
        elif style == "ensembl":
            prof.no_chr_prefix += 1
        elif style == "ncbi":
            prof.ncbi_names += 1
        if any(p.lower() in raw.lower() for p in ("_alt", "_random", "decoy", "chrUn", "KI27")):
            prof.alt_like += 1
        b = _scan_build_token(h)
        if b:
            prof.header_mentions_build = b
    prof.n_records = len(prof.contigs)
    return prof


def parse_gtf(path: Path, max_lines: int = 20000) -> FileProfile:
    prof = FileProfile(path=str(path), kind="gtf")
    seen: dict[str, Contig] = {}
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                b = _scan_build_token(line)
                if b:
                    prof.header_mentions_build = b
                continue
            n += 1
            if n > max_lines:
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            raw, start, end = parts[0], int(parts[3]), int(parts[4])
            canon, style = canon_name(raw)
            if raw not in seen:
                seen[raw] = Contig(raw_name=raw, canon=canon, style=style)
            prof.max_coord[canon] = max(prof.max_coord.get(canon, 0), end)
            gm = _GENE_RE.search(parts[8] if len(parts) > 8 else "")
            if gm:
                g = gm.group(1)
                prof.gene_hits.setdefault(g, []).append((canon, start))
            if style == "ucsc":
                prof.has_chr_prefix += 1
            elif style == "ensembl":
                prof.no_chr_prefix += 1
            elif style == "ncbi":
                prof.ncbi_names += 1
    prof.contigs = list(seen.values())
    prof.n_records = n
    prof.missing_contig_lengths = True
    return prof


def parse_vcf(path: Path, max_body: int = 5000) -> FileProfile:
    prof = FileProfile(path=str(path), kind="vcf")
    seen: dict[str, Contig] = {}
    body = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith("##contig="):
                raw = ""
                ln = None
                m = re.search(r"ID=([^,>]+)", line)
                if m:
                    raw = m.group(1)
                m2 = re.search(r"length=(\d+)", line)
                if m2:
                    ln = int(m2.group(1))
                canon, style = canon_name(raw)
                seen[raw] = Contig(raw_name=raw, canon=canon, length=ln, style=style)
                if style == "ucsc":
                    prof.has_chr_prefix += 1
                elif style == "ensembl":
                    prof.no_chr_prefix += 1
                elif style == "ncbi":
                    prof.ncbi_names += 1
                continue
            if line.startswith("#"):
                b = _scan_build_token(line)
                if b:
                    prof.header_mentions_build = b
                continue
            body += 1
            if body > max_body:
                break
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            raw = parts[0]
            pos = int(parts[1])
            canon, style = canon_name(raw)
            if raw not in seen:
                seen[raw] = Contig(raw_name=raw, canon=canon, style=style)
            prof.max_coord[canon] = max(prof.max_coord.get(canon, 0), pos)
    prof.contigs = list(seen.values())
    prof.n_records = body
    if not any(c.length for c in prof.contigs):
        prof.missing_contig_lengths = True
    return prof


def parse_any(path: Path) -> FileProfile:
    suf = path.suffix.lower()
    name = path.name.lower()
    if suf in (".fa", ".fasta", ".fna") or name.endswith(".fa") or "fasta" in name:
        return parse_fasta(path)
    if suf in (".gtf", ".gff", ".gff3"):
        return parse_gtf(path)
    if suf in (".vcf",) or name.endswith(".vcf"):
        return parse_vcf(path)
    # sniff
    with path.open() as fh:
        first = fh.readline()
    if first.startswith(">"):
        return parse_fasta(path)
    if first.startswith("##fileformat=VCF") or first.startswith("#CHROM"):
        return parse_vcf(path)
    return parse_gtf(path)
