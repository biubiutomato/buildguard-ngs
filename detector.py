"""Rule engine + lightweight scored classifier for genome-build / file mismatch."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from fingerprints import (
    ALT_PATTERNS,
    BUILD_LEN,
    BUILDS,
    LANDMARKS,
    LEN_INDEX,
    NCBI_TO_CHR,
    PRIMARY,
)
from parse import FileProfile


@dataclass
class FileCall:
    path: str
    kind: str
    build: str | None
    naming: str  # ucsc | ensembl | ncbi | mixed | unknown
    votes: dict[str, float] = field(default_factory=dict)
    evidence: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    n_length_hits: int = 0
    n_oob: int = 0
    mt_flavor: str | None = None  # rcrs | yoruba | unknown


@dataclass
class BundleReport:
    files: list[FileCall]
    consistent: bool
    predicted_build: str | None
    issues: list[dict]
    severity: str  # ok | warning | error


def _naming_of(prof: FileProfile) -> str:
    styles = Counter(c.style for c in prof.contigs if c.style != "other")
    if not styles:
        return "unknown"
    if len(styles) > 1:
        # ignore tiny minorities
        top, n = styles.most_common(1)[0]
        if n >= 0.8 * sum(styles.values()):
            return top
        return "mixed"
    return styles.most_common(1)[0][0]


def _length_votes(prof: FileProfile) -> tuple[Counter, int, list[str]]:
    votes: Counter = Counter()
    hits = 0
    ev = []
    for c in prof.contigs:
        if c.length is None:
            continue
        matches = LEN_INDEX.get(c.length, [])
        if not matches:
            # try nearby ±1 for MT / off-by-one
            for d in (1, -1):
                matches = LEN_INDEX.get(c.length + d, [])
                if matches:
                    ev.append(f"near-length {c.raw_name}={c.length} (~{matches})")
                    break
        if not matches:
            continue
        hits += 1
        for build, chrom in matches:
            votes[build] += 1
            if c.canon != "other" and chrom != c.canon:
                ev.append(f"length of {c.raw_name} matches {build}:{chrom} but name is {c.canon}")
    return votes, hits, ev


def _coord_oob(prof: FileProfile, assumed: str) -> tuple[int, list[str]]:
    table = BUILD_LEN[assumed]
    n = 0
    ev = []
    for chrom, mx in prof.max_coord.items():
        lim = table.get(chrom)
        if lim is None:
            continue
        if mx > lim:
            n += 1
            ev.append(f"{chrom} coordinate {mx} exceeds {assumed} length {lim}")
    return n, ev


def _landmark_votes(prof: FileProfile) -> tuple[Counter, list[str]]:
    votes: Counter = Counter()
    ev = []
    for gene, hits in prof.gene_hits.items():
        for build, spec in LANDMARKS.items():
            if gene not in spec:
                continue
            exp_chr, lo, hi = spec[gene]
            for chrom, start in hits:
                if chrom == exp_chr and lo - 50_000 <= start <= hi + 50_000:
                    votes[build] += 1
                    ev.append(f"landmark {gene}@{chrom}:{start} ~ {build}")
    return votes, ev


def _mt_flavor(prof: FileProfile) -> str | None:
    for c in prof.contigs:
        if c.canon == "chrM" and c.length:
            if c.length == 16571:
                return "yoruba"
            if c.length == 16569:
                return "rcrs"
    return None


def call_file(prof: FileProfile) -> FileCall:
    naming = _naming_of(prof)
    lv, n_hits, lev = _length_votes(prof)
    lm, lmev = _landmark_votes(prof)
    scores = {b: 0.0 for b in BUILDS}
    ev = list(lev) + list(lmev)
    warn = []

    # length votes dominate
    for b, v in lv.items():
        scores[b] += 3.0 * v
    for b, v in lm.items():
        scores[b] += 1.0 * v

    if prof.header_mentions_build:
        scores[prof.header_mentions_build] += 0.5
        ev.append(f"header token {prof.header_mentions_build}")

    # NCBI accession is a hard vote
    for c in prof.contigs:
        if c.raw_name in NCBI_TO_CHR:
            _, b = NCBI_TO_CHR[c.raw_name]
            scores[b] += 4.0
            ev.append(f"NCBI accession {c.raw_name} => {b}")

    # out-of-bound coords as negative evidence
    oob_total = 0
    for b in BUILDS:
        n_oob, oev = _coord_oob(prof, b)
        if n_oob:
            scores[b] -= 2.0 * n_oob
            if b == max(scores, key=scores.get):
                ev.extend(oev[:3])
        oob_total += n_oob

    # alt contig presence is not a build vote but a warning
    if prof.alt_like:
        warn.append(f"{prof.alt_like} alt/decoy/unplaced contig name(s)")

    styles = {c.style for c in prof.contigs if c.style in ("ucsc", "ensembl", "ncbi")}
    if naming == "mixed" or len(styles) > 1:
        warn.append("mixed contig naming styles in one file")
        naming = "mixed" if len(styles) > 1 else naming

    if prof.missing_contig_lengths and prof.kind == "vcf":
        warn.append("VCF has no ##contig length lines; build call is weaker")

    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    best, best_s = ranked[0]
    second_s = ranked[1][1]
    build = best if best_s > 0 and best_s > second_s else (best if best_s >= 2 else None)
    if best_s > 0 and best_s == second_s:
        build = None
        warn.append("ambiguous length/landmark votes")

    mt = _mt_flavor(prof)
    if mt == "yoruba" and build in ("hg38", "t2t"):
        warn.append("chrM length 16571 (NC_001807) with a GRCh38/T2T nuclear build")
    if mt == "rcrs" and build == "hg19":
        warn.append("chrM length 16569 (rCRS) with an hg19 nuclear build — common but inconsistent")

    if build:
        n_oob_best, oev_best = _coord_oob(prof, build)
        if n_oob_best:
            warn.append(f"{n_oob_best} coordinate(s) exceed {build} contig length")
            ev.extend(oev_best[:3])

    return FileCall(
        path=prof.path,
        kind=prof.kind,
        build=build,
        naming=naming,
        votes=scores,
        evidence=ev[:12],
        warnings=warn,
        n_length_hits=n_hits,
        n_oob=oob_total,
        mt_flavor=mt,
    )


def check_bundle(calls: list[FileCall]) -> BundleReport:
    issues = []
    builds = [c.build for c in calls if c.build]
    names = [c.naming for c in calls if c.naming not in ("unknown",)]

    # cross-file build
    if len(set(builds)) > 1:
        issues.append({
            "code": "BUILD_MISMATCH",
            "severity": "error",
            "detail": "files vote for different assemblies: "
            + ", ".join(f"{c.kind}:{c.build}" for c in calls),
        })
    if any(c.build is None for c in calls) and builds:
        issues.append({
            "code": "BUILD_UNKNOWN",
            "severity": "warning",
            "detail": "at least one file could not be assigned a build",
        })

    # naming
    core_names = [n for n in names if n in ("ucsc", "ensembl", "ncbi")]
    if len(set(core_names)) > 1:
        issues.append({
            "code": "NAMING_MISMATCH",
            "severity": "error",
            "detail": "contig naming styles differ across files: "
            + ", ".join(f"{c.kind}:{c.naming}" for c in calls),
        })

    # MT
    mts = {c.mt_flavor for c in calls if c.mt_flavor}
    if len(mts) > 1:
        issues.append({
            "code": "MT_FLAVOR_MISMATCH",
            "severity": "warning",
            "detail": f"mitochondrial length flavors mixed: {mts}",
        })

    for c in calls:
        for w in c.warnings:
            issues.append({
                "code": "FILE_WARNING",
                "severity": "warning",
                "detail": f"{c.kind}: {w}",
            })

    # majority build
    pred = None
    if builds:
        pred = Counter(builds).most_common(1)[0][0]

    sev = "ok"
    if any(i["severity"] == "error" for i in issues):
        sev = "error"
    elif issues:
        sev = "warning"
    return BundleReport(
        files=calls,
        consistent=(sev == "ok"),
        predicted_build=pred,
        issues=issues,
        severity=sev,
    )


def extract_features(prof: FileProfile) -> dict[str, float]:
    """Numeric features for the scored group classifier / ablation."""
    lv, n_hits, _ = _length_votes(prof)
    lm, _ = _landmark_votes(prof)
    naming = _naming_of(prof)
    return {
        "n_contigs": float(len(prof.contigs)),
        "n_length_hits": float(n_hits),
        "vote_hg19": float(lv.get("hg19", 0)),
        "vote_hg38": float(lv.get("hg38", 0)),
        "vote_t2t": float(lv.get("t2t", 0)),
        "lm_hg19": float(lm.get("hg19", 0)),
        "lm_hg38": float(lm.get("hg38", 0)),
        "lm_t2t": float(lm.get("t2t", 0)),
        "has_chr": float(naming == "ucsc"),
        "ensembl": float(naming == "ensembl"),
        "ncbi": float(naming == "ncbi"),
        "mixed_name": float(naming == "mixed"),
        "alt_like": float(prof.alt_like),
        "header_build": float(1.0 if prof.header_mentions_build else 0.0),
        "missing_ln": float(1.0 if prof.missing_contig_lengths else 0.0),
    }
