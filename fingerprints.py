"""Canonical contig fingerprints for human genome builds.

Lengths are primary assembled chromosomes only (no alt/decoy/unplaced).
Sources: GRCh37.p13 / hg19, GRCh38.p14 / hg38, T2T-CHM13v2.0 / hs1.
"""

from __future__ import annotations

BUILDS = ("hg19", "hg38", "t2t")

# UCSC-style names
HG19_LEN = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
    "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
    "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
    "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
    "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566,
    "chrM": 16571,  # historical Yoruba NC_001807; rCRS is 16569
}

HG38_LEN = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468, "chrX": 156040895, "chrY": 57227415,
    "chrM": 16569,  # rCRS
}

T2T_LEN = {
    "chr1": 248387328, "chr2": 242696752, "chr3": 201105948, "chr4": 193574945,
    "chr5": 182045439, "chr6": 172126628, "chr7": 160567428, "chr8": 146259331,
    "chr9": 150617247, "chr10": 134758134, "chr11": 135127769, "chr12": 133324548,
    "chr13": 113566686, "chr14": 101161492, "chr15": 99753195, "chr16": 96330374,
    "chr17": 84276897, "chr18": 80542538, "chr19": 61707364, "chr20": 66210255,
    "chr21": 45090682, "chr22": 51324926, "chrX": 154259566, "chrY": 62460029,
    "chrM": 16569,
}

BUILD_LEN = {"hg19": HG19_LEN, "hg38": HG38_LEN, "t2t": T2T_LEN}

# NCBI RefSeq accessions (primary)
NCBI_TO_CHR = {
    # GRCh37
    "NC_000001.10": ("chr1", "hg19"), "NC_000002.11": ("chr2", "hg19"),
    "NC_000003.11": ("chr3", "hg19"), "NC_000004.11": ("chr4", "hg19"),
    "NC_000005.9": ("chr5", "hg19"), "NC_000006.11": ("chr6", "hg19"),
    "NC_000007.13": ("chr7", "hg19"), "NC_000008.10": ("chr8", "hg19"),
    "NC_000009.11": ("chr9", "hg19"), "NC_000010.10": ("chr10", "hg19"),
    "NC_000011.9": ("chr11", "hg19"), "NC_000012.11": ("chr12", "hg19"),
    "NC_000013.10": ("chr13", "hg19"), "NC_000014.8": ("chr14", "hg19"),
    "NC_000015.9": ("chr15", "hg19"), "NC_000016.9": ("chr16", "hg19"),
    "NC_000017.10": ("chr17", "hg19"), "NC_000018.9": ("chr18", "hg19"),
    "NC_000019.9": ("chr19", "hg19"), "NC_000020.10": ("chr20", "hg19"),
    "NC_000021.8": ("chr21", "hg19"), "NC_000022.10": ("chr22", "hg19"),
    "NC_000023.10": ("chrX", "hg19"), "NC_000024.9": ("chrY", "hg19"),
    "NC_001807.4": ("chrM", "hg19"),
    # GRCh38
    "NC_000001.11": ("chr1", "hg38"), "NC_000002.12": ("chr2", "hg38"),
    "NC_000003.12": ("chr3", "hg38"), "NC_000004.12": ("chr4", "hg38"),
    "NC_000005.10": ("chr5", "hg38"), "NC_000006.12": ("chr6", "hg38"),
    "NC_000007.14": ("chr7", "hg38"), "NC_000008.11": ("chr8", "hg38"),
    "NC_000009.12": ("chr9", "hg38"), "NC_000010.11": ("chr10", "hg38"),
    "NC_000011.10": ("chr11", "hg38"), "NC_000012.12": ("chr12", "hg38"),
    "NC_000013.11": ("chr13", "hg38"), "NC_000014.9": ("chr14", "hg38"),
    "NC_000015.10": ("chr15", "hg38"), "NC_000016.10": ("chr16", "hg38"),
    "NC_000017.11": ("chr17", "hg38"), "NC_000018.10": ("chr18", "hg38"),
    "NC_000019.10": ("chr19", "hg38"), "NC_000020.11": ("chr20", "hg38"),
    "NC_000021.9": ("chr21", "hg38"), "NC_000022.11": ("chr22", "hg38"),
    "NC_000023.11": ("chrX", "hg38"), "NC_000024.10": ("chrY", "hg38"),
    "NC_012920.1": ("chrM", "hg38"),
    # T2T-CHM13v2.0
    "NC_060925.1": ("chr1", "t2t"), "NC_060926.1": ("chr2", "t2t"),
    "NC_060927.1": ("chr3", "t2t"), "NC_060928.1": ("chr4", "t2t"),
    "NC_060929.1": ("chr5", "t2t"), "NC_060930.1": ("chr6", "t2t"),
    "NC_060931.1": ("chr7", "t2t"), "NC_060932.1": ("chr8", "t2t"),
    "NC_060933.1": ("chr9", "t2t"), "NC_060934.1": ("chr10", "t2t"),
    "NC_060935.1": ("chr11", "t2t"), "NC_060936.1": ("chr12", "t2t"),
    "NC_060937.1": ("chr13", "t2t"), "NC_060938.1": ("chr14", "t2t"),
    "NC_060939.1": ("chr15", "t2t"), "NC_060940.1": ("chr16", "t2t"),
    "NC_060941.1": ("chr17", "t2t"), "NC_060942.1": ("chr18", "t2t"),
    "NC_060943.1": ("chr19", "t2t"), "NC_060944.1": ("chr20", "t2t"),
    "NC_060945.1": ("chr21", "t2t"), "NC_060946.1": ("chr22", "t2t"),
    "NC_060947.1": ("chrX", "t2t"), "NC_060948.1": ("chrY", "t2t"),
}

# GenBank CP accessions used in some T2T FASTA headers
GENBANK_TO_CHR = {
    "CP068277.2": "chr1", "CP068276.2": "chr2", "CP068275.2": "chr3",
    "CP068274.2": "chr4", "CP068273.2": "chr5", "CP068272.2": "chr6",
    "CP068271.2": "chr7", "CP068270.2": "chr8", "CP068269.2": "chr9",
    "CP068268.2": "chr10", "CP068267.2": "chr11", "CP068266.2": "chr12",
    "CP068265.2": "chr13", "CP068264.2": "chr14", "CP068263.2": "chr15",
    "CP068262.2": "chr16", "CP068261.2": "chr17", "CP068260.2": "chr18",
    "CP068259.2": "chr19", "CP068258.2": "chr20", "CP068257.2": "chr21",
    "CP068256.2": "chr22", "CP068255.2": "chrX", "CP086569.2": "chrY",
    "CP068254.1": "chrM",
}

# Landmark genes used when contig lengths are absent (GTF-only / headerless VCF)
# (name, chrom, start_approx, end_approx) — mid-exon anchors, not exact CDS
LANDMARKS = {
    "hg19": {
        "TP53": ("chr17", 7571720, 7590868),
        "BRCA1": ("chr17", 41196312, 41277500),
        "EGFR": ("chr7", 55086725, 55275031),
        "MYC": ("chr8", 128748315, 128753680),
    },
    "hg38": {
        "TP53": ("chr17", 7661779, 7687550),
        "BRCA1": ("chr17", 43044295, 43125483),
        "EGFR": ("chr7", 55019017, 55211628),
        "MYC": ("chr8", 127735434, 127741434),
    },
    "t2t": {
        # approximate T2T coords (offset from hg38 by known synteny; used only as weak prior)
        "TP53": ("chr17", 7680000, 7710000),
        "BRCA1": ("chr17", 43100000, 43200000),
        "EGFR": ("chr7", 55050000, 55250000),
        "MYC": ("chr8", 127900000, 128000000),
    },
}

ALT_PATTERNS = (
    "_random", "_alt", "_fix", "Un_", "chrUn", "HLA-", "decoy",
    "KI27", "GL00", "JTFH", "KQ98", "KN53",
)


def length_index() -> dict[int, list[tuple[str, str]]]:
    """Map exact length -> list of (build, chrom)."""
    idx: dict[int, list[tuple[str, str]]] = {}
    for build, table in BUILD_LEN.items():
        for chrom, ln in table.items():
            idx.setdefault(ln, []).append((build, chrom))
    return idx


LEN_INDEX = length_index()
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
