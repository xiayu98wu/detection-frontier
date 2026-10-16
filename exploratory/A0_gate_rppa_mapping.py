"""Track A gate check.

Can the 458 RPPA peptide targets be mapped to gene symbols, and how many genes
have all three layers (CNV + mRNA + protein)?

RPPA peptide_target naming is messy: antibody-clone suffixes (BAP1C4), phospho
sites (ACC_pS79), aliases (1433BETA). This measures how much survives cleaning.
Also checks which layers are actually quantitative in the CNV matrix.
"""
import gzip, os, collections, statistics
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")

# ---------------- RPPA ----------------
rppa = {}
with open(os.path.join(DM, "MESO_rppa.tsv"), encoding="utf-8") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    samps = hdr[1:]
    for l in f:
        c = l.rstrip("\n").split("\t")
        vals = []
        for v in c[1:]:
            try:
                vals.append(float(v))
            except Exception:
                vals.append(None)
        rppa[c[0]] = vals
print("RPPA peptides: %d  samples: %d" % (len(rppa), len(samps)))

# ---------------- expression ----------------
with gzip.open(os.path.join(DM, "MESO_expr_full.tsv.gz"), "rt", encoding="utf-8", errors="replace") as g:
    h = g.readline().rstrip("\n").split("\t")
    ecols = [c.strip('"')[:12] for c in h[1:]]
    expr_genes = set()
    for line in g:
        expr_genes.add(line.split("\t", 1)[0].strip('"'))
print("expression genes: %d" % len(expr_genes))

# ---------------- CNV ----------------
cnv_genes = set()
with gzip.open(os.path.join(DM, "MESO_cnv_matrix.tsv.gz"), "rt", encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    ccols = [c.strip('"') for c in h[1:]]
    for l in f:
        cnv_genes.add(l.split("\t", 1)[0])
print("CNV genes: %d   CNV samples: %d" % (len(cnv_genes), len(ccols)))

# ---------------- clean RPPA names ----------------
PHOS = ("_pS", "_pT", "_pY", "_p", "P_p", "_p70", "_pS65")
ALIAS = {
    "1433BETA": "YWHAB", "1433EPSILON": "YWHAE", "1433ZETA": "YWHAZ",
    "4EBP1": "EIF4EBP1", "53BP1": "TP53BP1", "ACC1": "ACACA", "ACC": "ACACA",
    "AR": "AR", "BCL2": "BCL2", "BCLXL": "BCL2L1", "CD31": "PECAM1",
    "CD45": "PTPRC", "CD49B": "ITGA2", "CDK1": "CDK1", "CHK1": "CHEK1",
    "CHK2": "CHEK2", "CLAUDIN7": "CLDN7", "COLLAGENVI": "COL6A1",
    "CYCLINB1": "CCNB1", "CYCLIND1": "CCND1", "CYCLINE1": "CCNE1",
    "DDR1": "DDR1", "E_CADHERIN": "CDH1", "ECADHERIN": "CDH1",
    "EGFR": "EGFR", "ERALPHA": "ESR1", "ERALPHA_P": "ESR1", "ERBB2": "ERBB2",
    "ERBB3": "ERBB3", "FIBRONECTIN": "FN1", "GAPDH": "GAPDH",
    "GATA3": "GATA3", "GSK3ALPHABETA": "GSK3A", "HER2": "ERBB2",
    "HER3": "ERBB3", "HSP70": "HSPA1A", "IGFBP2": "IGFBP2",
    "INSULIN_R": "INSR", "JAK2": "JAK2", "KERATIN18": "KRT18",
    "KERATIN19": "KRT19", "KIT": "KIT", "KU80": "XRCC5", "LCK": "LCK",
    "MAPK": "MAPK1", "MEK1": "MAP2K1", "MET": "MET", "MIG6": "ERRFI1",
    "MSH2": "MSH2", "MSH6": "MSH6", "MTOR": "MTOR", "MYH9": "MYH9",
    "NCADHERIN": "CDH2", "NF2": "NF2", "NFKB": "NFKB1", "NOTCH1": "NOTCH1",
    "P21": "CDKN1A", "P27": "CDKN1B", "P38": "MAPK14", "P53": "TP53",
    "P62": "SQSTM1", "P70S6K": "RPS6KB1", "P90RSK": "RPS6KA1",
    "PARP": "PARP1", "PDL1": "CD274", "PD1": "PDCD1", "PI3KP110ALPHA": "PIK3CA",
    "PR": "PGR", "PTEN": "PTEN", "RB": "RB1", "S6": "RPS6",
    "SMAD1": "SMAD1", "SMAD3": "SMAD3", "SMAD4": "SMAD4", "SNAIL": "SNAI1",
    "SRC": "SRC", "STAT1": "STAT1", "STAT3": "STAT3", "STAT5ALPHA": "STAT5A",
    "STING": "TMEM173", "SYNAPTOPHYSIN": "SYP", "TAZ": "TAFAZZIN",
    "TFRC": "TFRC", "TSC1": "TSC1", "TSC2": "TSC2", "TWIST": "TWIST1",
    "VIMENTIN": "VIM", "WT1": "WT1", "XRCC1": "XRCC1", "YAP": "YAP1",
    "CD3": "CD3E", "CD4": "CD4", "CD8": "CD8A", "CD20": "MS4A1",
    "CD68": "CD68", "CD163": "CD163", "CD11C": "ITGAX", "CD14": "CD14",
    "CD44": "CD44", "CDH3": "CDH3", "CDH2": "CDH2",
    "GRANZYMEB": "GZMB", "GRANZYME-B": "GZMB", "HLA-DQA1": "HLA-DQA1",
    "IDO": "IDO1", "ZAP70": "ZAP70", "ZAP-70": "ZAP70", "SYK": "SYK",
    "VAV1": "VAV1", "HMHA1": "ARHGAP45", "PREX1": "PREX1",
    "BAP1C4": "BAP1", "BAP1": "BAP1",
}


def clean(name):
    n = name.strip().upper()
    if n in ALIAS:
        return ALIAS[n]
    # strip common clone suffixes (letter+digits at end)
    base = n
    for suf in ["C4", "C2", "C1", "EP", "AB", "-B", "_B"]:
        if base.endswith(suf) and len(base) > len(suf) + 2:
            base = base[: -len(suf)]
            break
    if base in ALIAS:
        return ALIAS[base]
    return base


total, phospho, mapped = [], [], {}
for name in rppa:
    is_phos = any(p in name for p in ("_pS", "_pT", "_pY")) or name.endswith("_P")
    if is_phos:
        phospho.append(name)
        continue
    total.append(name)
    mapped[name] = clean(name)
print("\ntotal-protein peptides: %d   phospho peptides: %d" % (len(total), len(phospho)))

genes = set(mapped.values())
print("distinct gene symbols after cleaning: %d" % len(genes))

in_e = {g for g in genes if g in expr_genes}
in_c = {g for g in genes if g in cnv_genes}
three = in_e & in_c
print("\n  mappable to expression : %d (%.0f%%)" % (len(in_e), 100 * len(in_e) / len(genes)))
print("  mappable to CNV        : %d (%.0f%%)" % (len(in_c), 100 * len(in_c) / len(genes)))
print("  ALL THREE layers       : %d" % len(three))

print("\n  unmapped to expression (%d):" % len(genes - in_e))
print("   ", sorted(genes - in_e)[:40])

print("\n  three-layer gene list:")
print("   ", ", ".join(sorted(three)))

# sanity: are the key genes there?
print("\n  key genes:")
for g in ["BAP1", "NF2", "CDKN2A", "MTAP", "CD274", "VIM", "TP53", "PTEN", "RB1"]:
    print("    %-8s expr=%-5s cnv=%-5s prot=%-5s" % (g, g in expr_genes, g in cnv_genes, g in genes))

# how many RPPA samples overlap CNV/expression samples
print("\n  sample overlap:")
print("    RPPA samples        : %d" % len(samps))
print("    in expression       : %d" % len([s for s in samps if s in set(ecols)]))
