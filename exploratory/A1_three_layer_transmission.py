"""A1: three-layer transmission (CNV -> mRNA -> protein) in TCGA-MESO.

For the 214 genes present in all three layers, across the 62 shared samples:
  a = rho(CNV, mRNA)      genomic -> transcript transmission
  b = rho(mRNA, protein)  transcript -> protein transmission
  c = rho(CNV, protein)   overall genomic -> protein transmission

Transmission ratio  T = c / a   (fraction of the cis-CNV effect surviving to protein)

Genes are only informative where CNV actually varies; a separate analysis handles
the mRNA-proteome correlation for all 214.

Controls built in:
  - CNV-varying genes vs CNV-stable genes (does the signal require CNV variance?)
  - the same statistic computed on a random gene subset for calibration
"""
import gzip, os, collections, statistics, math
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")
OUT = os.path.join(_ROOT, "results")
os.makedirs(OUT, exist_ok=True)

# ---------- load RPPA ----------
rppa = {}
with open(os.path.join(DM, "MESO_rppa.tsv"), encoding="utf-8") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    rsamps = hdr[1:]
    for l in f:
        c = l.rstrip("\n").split("\t")
        rppa[c[0]] = [float(v) if v != "NA" else None for v in c[1:]]

# ---------- RPPA name cleaning (same map as the gate script, condensed) ----------
ALIAS = {
    "1433BETA": "YWHAB", "1433EPSILON": "YWHAE", "1433ZETA": "YWHAZ",
    "4EBP1": "EIF4EBP1", "53BP1": "TP53BP1", "ACC1": "ACACA", "ACC": "ACACA",
    "BCLXL": "BCL2L1", "CD31": "PECAM1", "CD45": "PTPRC", "CD49B": "ITGA2",
    "CHK1": "CHEK1", "CHK2": "CHEK2", "CLAUDIN7": "CLDN7", "COLLAGENVI": "COL6A1",
    "CYCLINB1": "CCNB1", "CYCLIND1": "CCND1", "CYCLINE1": "CCNE1",
    "E_CADHERIN": "CDH1", "ECADHERIN": "CDH1", "ERALPHA": "ESR1",
    "FIBRONECTIN": "FN1", "GSK3ALPHABETA": "GSK3A", "HER2": "ERBB2",
    "HER3": "ERBB3", "HSP70": "HSPA1A", "INSULIN_R": "INSR",
    "KERATIN18": "KRT18", "KERATIN19": "KRT19", "KU80": "XRCC5",
    "MAPK": "MAPK1", "MEK1": "MAP2K1", "MIG6": "ERRFI1", "NCADHERIN": "CDH2",
    "NFKB": "NFKB1", "P21": "CDKN1A", "P27": "CDKN1B", "P38": "MAPK14",
    "P53": "TP53", "P62": "SQSTM1", "P70S6K": "RPS6KB1", "P90RSK": "RPS6KA1",
    "PARP": "PARP1", "PDL1": "CD274", "PD1": "PDCD1", "PI3KP110ALPHA": "PIK3CA",
    "RB": "RB1", "S6": "RPS6", "SNAIL": "SNAI1", "STAT5ALPHA": "STAT5A",
    "STING": "TMEM173", "SYNAPTOPHYSIN": "SYP", "TAZ": "TAFAZZIN",
    "TWIST": "TWIST1", "VIMENTIN": "VIM", "CD3": "CD3E", "CD8": "CD8A",
    "CD20": "MS4A1", "CD11C": "ITGAX", "GRANZYMEB": "GZMB",
    "GRANZYME-B": "GZMB", "IDO": "IDO1", "ZAP-70": "ZAP70", "ZAP70": "ZAP70",
    "HMHA1": "ARHGAP45", "BAP1C4": "BAP1", "CD134": "TNFRSF4",
    "CONNEXIN-43": "GJA1", "CMYC": "MYC", "CMET": "MET", "BETACATENIN": "CTNNB1",
    "B-ACTIN": "ACTB", "AKT": "AKT1", "AKTP_S473": "AKT1", "AKTP_T308": "AKT1",
    "B7-H3": "CD276", "B7-H4": "VTCN1", "BECLIN": "BECN1", "BIM": "BCL2L11",
    "BIP-GRP78": "HSPA5", "CABL": "ABL1", "CAVEOLIN1": "CAV1", "CASPASE3": "CASP3",
    "CASPASE7": "CASP7", "CASPASE8": "CASP8", "CGAS": "MB21D1", "CIAP": "BIRC2",
    "C-IAP2": "BIRC3", "CKIT": "KIT", "COX-IV": "COX4I1", "COX2": "PTGS2",
    "ANNEXIN1": "ANXA1", "ANNEXINVII": "ANXA7", "ATP5A": "ATP5F1A",
    "AURORA": "AURKA", "AURORA-A": "AURKA", "ADAR1": "ADAR", "ACECS1": "ACACA",
}


def clean(name):
    n = name.strip().upper()
    if n in ALIAS:
        return ALIAS[n]
    for suf in ["C4", "C2", "C1"]:
        if n.endswith(suf) and len(n) > len(suf) + 2 and n[: -len(suf)] in ALIAS:
            return ALIAS[n[: -len(suf)]]
    return n


prot = {}   # gene -> {sample: value}
for name, vals in rppa.items():
    if any(p in name for p in ("_pS", "_pT", "_pY")):
        continue
    g = clean(name)
    d = {}
    for s, v in zip(rsamps, vals):
        if v is not None:
            d[s] = v
    if len(d) >= 50:
        prot[g] = d

# ---------- load expression ----------
with gzip.open(os.path.join(DM, "MESO_expr_full.tsv.gz"), "rt", encoding="utf-8", errors="replace") as g:
    h = g.readline().rstrip("\n").split("\t")
    ecols = [c.strip('"')[:12] for c in h[1:]]
    expr = {}
    for line in g:
        c = line.rstrip("\n").split("\t")
        expr[c[0].strip('"')] = [float(x) for x in c[1:]]

# ---------- load CNV ----------
cnv = {}
with gzip.open(os.path.join(DM, "MESO_cnv_matrix.tsv.gz"), "rt", encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    ccols = [c.strip('"') for c in h[1:]]
    for l in f:
        c = l.rstrip("\n").split("\t")
        d = {}
        for i, b in enumerate(ccols):
            if c[i + 1] != "NA":
                d[b] = int(float(c[i + 1]))
        cnv[c[0]] = d

# ---------- shared samples: all three layers ----------
SHARED = [s for s in ecols if s in rsamps and s in set(ccols)]
print("shared samples (expression + RPPA + CNV): %d" % len(SHARED))
eidx = {b: i for i, b in enumerate(ecols)}
GENES = sorted(set(prot) & set(expr) & set(cnv))
print("three-layer genes: %d" % len(GENES))


def spearman(x, y):
    n = len(x)
    if n < 8:
        return None
    def rk(v):
        o = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[o[j + 1]] == v[o[i]]:
                j += 1
            a = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[o[k]] = a
            i = j + 1
        return r
    rx, ry = rk(x), rk(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


rows = []
for g in GENES:
    cn = [cnv[g].get(s) for s in SHARED]
    ex = [expr[g][eidx[s]] for s in SHARED]
    pr = [prot[g].get(s) for s in SHARED]
    ok = [i for i in range(len(SHARED)) if cn[i] is not None and pr[i] is not None]
    if len(ok) < 15:
        continue
    a = spearman([cn[i] for i in ok], [ex[i] for i in ok])
    b = spearman([ex[i] for i in ok], [pr[i] for i in ok])
    c = spearman([cn[i] for i in ok], [pr[i] for i in ok])
    vsd = statistics.pstdev([cn[i] for i in ok])
    rows.append({
        "gene": g, "n": len(ok), "cnv_sd": vsd,
        "r_cnv_mrna": a, "r_mrna_prot": b, "r_cnv_prot": c,
    })

print("genes with >=15 usable samples: %d" % len(rows))

# ============ global summary ============
def summ(vals):
    v = [x for x in vals if x is not None]
    if not v:
        return "n/a"
    return "n=%3d  mean=%+.3f  median=%+.3f  min=%+.3f  max=%+.3f" % (
        len(v), statistics.mean(v), statistics.median(v), min(v), max(v))

print("\n" + "=" * 90)
print("A1. LAYER-TO-LAYER CORRELATION (Spearman), TCGA-MESO, %d genes" % len(rows))
print("=" * 90)
print("  rho(CNV, mRNA)     : %s" % summ([r["r_cnv_mrna"] for r in rows]))
print("  rho(mRNA, protein) : %s" % summ([r["r_mrna_prot"] for r in rows]))
print("  rho(CNV,  protein) : %s" % summ([r["r_cnv_prot"] for r in rows]))

# CNV-varying vs stable
VAR = [r for r in rows if r["cnv_sd"] >= 0.3]
STABLE = [r for r in rows if r["cnv_sd"] < 0.1]
print("\n  --- split by CNV variance ---")
print("  CNV-varying (sd>=0.3) n=%d : rho(CNV,mRNA)=%s" % (len(VAR), summ([r["r_cnv_mrna"] for r in VAR])))
print("                        rho(mRNA,prot)=%s" % summ([r["r_mrna_prot"] for r in VAR]))
print("                        rho(CNV,prot) =%s" % summ([r["r_cnv_prot"] for r in VAR]))
print("  CNV-stable  (sd<0.1)  n=%d : rho(mRNA,prot)=%s" % (len(STABLE), summ([r["r_mrna_prot"] for r in STABLE])))

# ============ transmission ratio ============
print("\n" + "=" * 90)
print("A1b. TRANSMISSION RATIO  T = rho(CNV,prot) / rho(CNV,mRNA)")
print("     (fraction of the cis copy-number effect that survives to protein)")
print("=" * 90)
TR = []
for r in VAR:
    a, c = r["r_cnv_mrna"], r["r_cnv_prot"]
    if a is not None and c is not None and a > 0.15:
        r2 = dict(r)
        r2["T"] = c / a
        TR.append(r2)
print("  genes with rho(CNV,mRNA) > 0.15: %d" % len(TR))
if TR:
    tt = [r["T"] for r in TR]
    print("  T distribution: mean=%+.3f median=%+.3f  <0.5: %d  >0.8: %d" % (
        statistics.mean(tt), statistics.median(tt),
        sum(1 for t in tt if t < 0.5), sum(1 for t in tt if t > 0.8)))
    TR.sort(key=lambda r: -r["T"])
    print("\n  %-10s %8s %10s %10s %8s" % ("gene", "cnv_sd", "r(cnv,mrna)", "r(cnv,prot)", "T"))
    print("  " + "-" * 50)
    print("  --- most CONDUCTED (highest T) ---")
    for r in TR[:12]:
        print("  %-10s %8.2f %10.3f %10.3f %8.2f" % (r["gene"], r["cnv_sd"], r["r_cnv_mrna"], r["r_cnv_prot"], r["T"]))
    print("  --- most BUFFERED (lowest T) ---")
    for r in TR[-12:]:
        print("  %-10s %8.2f %10.3f %10.3f %8.2f" % (r["gene"], r["cnv_sd"], r["r_cnv_mrna"], r["r_cnv_prot"], r["T"]))

TR.sort(key=lambda r: -r["r_mrna_prot"] if r["r_mrna_prot"] is not None else 0)
print("\n  --- strongest mRNA-protein coupling (top 12) ---")
print("  %-10s %12s %12s %10s %8s" % ("gene", "r(mrna,prot)", "r(cnv,mrna)", "r(cnv,prot)", "T"))
for r in TR[:12]:
    print("  %-10s %12.3f %12.3f %10.3f %8.2f" % (r["gene"], r["r_mrna_prot"], r["r_cnv_mrna"], r["r_cnv_prot"], r["T"]))
print("  --- weakest mRNA-protein coupling (bottom 12) ---")
for r in TR[-12:]:
    print("  %-10s %12.3f %12.3f %10.3f %8.2f" % (r["gene"], r["r_mrna_prot"], r["r_cnv_mrna"], r["r_cnv_prot"], r["T"]))

# ============ save ============
with open(os.path.join(OUT, "A1_three_layer.tsv"), "w", encoding="utf-8") as f:
    f.write("gene\tn_samples\tcnv_sd\trho_cnv_mrna\trho_mrna_prot\trho_cnv_prot\tT\n")
    for r in sorted(rows, key=lambda x: x["gene"]):
        T = ""
        a, c = r["r_cnv_mrna"], r["r_cnv_prot"]
        if a and c is not None and a > 0.15:
            T = "%.3f" % (c / a)
        f.write("%s\t%d\t%.3f\t%s\t%s\t%s\t%s\n" % (
            r["gene"], r["n"], r["cnv_sd"],
            "" if r["r_cnv_mrna"] is None else "%.3f" % r["r_cnv_mrna"],
            "" if r["r_mrna_prot"] is None else "%.3f" % r["r_mrna_prot"],
            "" if r["r_cnv_prot"] is None else "%.3f" % r["r_cnv_prot"], T))
print("\n[saved] A1_three_layer.tsv")
