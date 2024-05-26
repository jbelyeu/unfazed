
HOM_ALT = 3
GT_UNKNOWN = 2 #this isn't used, just stated for clarity
HET = 1
HOM_REF = 0
SEX_KEY = {"male": 1, "female": 2}
VCF_TYPES = ["vcf", "vcf.gz", "bcf"]
SV_TYPES = ["DEL", "DUP", "INV", "CNV", "DUP:TANDEM", "DEL:ME", "CPX", "CTX"]
SNV_TYPES = ["POINT","SNV","INDEL"]
LABELS = ["chrom", "start", "end", "kid", "vartype"]
QUIET_MODE = False
# pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
CIGAR_MAP = {
    0: "M",
    1: "I",
    2: "D",
    3: "N",
    4: "S",
    5: "H",
    6: "P",
    7: "=",
    8: "X",
    9: "B",
}
# https://www.ncbi.nlm.nih.gov/grc/human
grch37_par1 = {
    "x": [10001, 2781479],
    "y": [10001, 2781479],
}
grch37_par2 = {
    "x": [155701383, 156030895],
    "y": [56887903, 57217415],
}

grch38_par1 = {
    "x": [60001, 2699520],
    "y": [10001, 2649520],
}

grch38_par2 = {
    "x": [154931044, 155260560],
    "y": [59034050, 59363566],
}


def get_prefix(vcf):
    chrom_prefix = ""
    for var in vcf:
        if "chr" in var.CHROM.lower():
            chrom_prefix = var.CHROM[:3]
        return chrom_prefix
    return ""

