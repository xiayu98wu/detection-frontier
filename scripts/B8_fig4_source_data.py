"""B8: export the data behind Figure 4 (the three-layer / reproduction panel).

Figure 4 needs, per layer, a numeric value that a reader can check:
  Layer 1 (DNA)    : BAP1 mutation prevalence, and CDKN2A 9p21 deletion prevalence
  Layer 2 (RNA)    : BAP1 mRNA in BAP1-mutant vs wild-type  (Mann-Whitney)
  Layer 3 (protein): BAP1 protein (RPPA) vs immune scores    (Spearman)
  Reproduction     : CDKN2A deletion rate vs the published rate; survival split

Nothing here is recomputed from scratch where a script already exists -- this
only writes the numbers that were previously printed to stdout but never saved.
"""
import gzip, os, collections, statistics, math, json

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root
D = os.path.join(_ROOT, "data")
OUT = os.path.join(_ROOT, "results")
os.makedirs(OUT, exist_ok=True)


def mw(a, b):
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return None, None
    comb = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = [0.0] * len(comb)
    i = 0
    while i < len(comb):
        j = i
        while j + 1 < len(comb) and comb[j + 1][0] == comb[i][0]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = r
        i = j + 1
    U1 = sum(ranks[k] for k in range(len(comb)) if comb[k][1] == 0) - n1 * (n1 + 1) / 2
    tie = collections.Counter(v for v, _ in comb)
    N = n1 + n2
    den = N * (N - 1)
    sd = math.sqrt(n1 * n2 / 12 * ((N + 1) - sum(t ** 3 - t for t in tie.values()) / den))
    if sd == 0:
        return None, None
    z = (U1 - n1 * n2 / 2) / sd
    return z, 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


def spearman(x, y):
    n = len(x)
    if n < 8:
        return None, None

    def rk(v):
        o = sorted(range(n), key=lambda i: v[i])
        r = [0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[o[j + 1]] == v[o[i]]:
                j += 1
            for k in range(i, j + 1):
                r[o[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    rx, ry = rk(x), rk(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    den = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)) * sum((ry[i] - my) ** 2 for i in range(n)))
    if den == 0:
        return None, None
    r = num / den
    t = r * math.sqrt((n - 2) / max(1e-12, 1 - r * r))
    return r, 2 * (1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2))))


# ---- expression ----
with gzip.open(os.path.join(D, "MESO_expr_full.tsv.gz"), "rt", encoding="utf-8", errors="replace") as g:
    hdr = g.readline().rstrip("\n").split("\t")
    cols = [c.strip('"') for c in hdr[1:]]
    expr = {}
    for line in g:
        c = line.rstrip("\n").split("\t")
        expr[c[0].strip('"')] = [float(x) for x in c[1:]]
bars = [c[:12] for c in cols]

# ---- mutations ----
mut = collections.defaultdict(set)
with open(os.path.join(D, "MESO_mutations.maf"), encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    gi, ci = h.index("Hugo_Symbol"), h.index("Tumor_Sample_Barcode")
    for l in f:
        c = l.rstrip("\n").split("\t")
        if len(c) > max(gi, ci):
            mut[c[gi]].add(c[ci][:12])
has_maf = set().union(*mut.values())

# ---- CNV ----
cnv = collections.defaultdict(dict)
with gzip.open(os.path.join(D, "MESO_cnv_matrix.tsv.gz"), "rt", encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    ccols = [c.strip('"') for c in h[1:]]
    for l in f:
        c = l.rstrip("\n").split("\t")
        for i, b in enumerate(ccols):
            if c[i + 1] != "NA":
                cnv[c[0]][b] = int(float(c[i + 1]))

# ---- RPPA ----
rppa = {}
with open(os.path.join(D, "MESO_rppa.tsv"), encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    rcols = [c.strip('"') for c in h[1:]]
    for l in f:
        c = l.rstrip("\n").split("\t")
        rppa[c[0].strip('"')] = {rcols[i]: float(v) for i, v in enumerate(c[1:]) if v != "NA"}

SIG = {
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "CD3D", "CD3E"],
    "CHECKPT": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "IFN_G": ["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG", "IDO1", "GBP1"],
    "MACROPHAGE": ["CD68", "CD163", "MRC1", "CSF1R"],
    "FIBROBLAST": ["COL1A1", "COL1A2", "FAP", "ACTA2", "DCN"],
    "EMT": ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1"],
    "PROLIF": ["MKI67", "TOP2A", "CCNB1", "BUB1", "AURKA"],
    "ANTIGEN": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB9"],
}


def score(genes, key):
    gs = [g for g in genes if g in expr]
    return {bars[j]: statistics.mean(expr[g][j] for g in gs) for j in range(len(bars))}


rows = []

print("=" * 92)
print("B8. FIGURE 4 SOURCE DATA")
print("=" * 92)

# ---------- Layer 1: DNA ----------
bap1_mut = mut.get("BAP1", set())
n_bap1 = len([b for b in bars if b in has_maf and b in bap1_mut])
n_ev = len([b for b in bars if b in has_maf])
print("\n[Layer 1 DNA]")
print("  BAP1 mutated          : %d / %d evaluable = %.3f" % (n_bap1, n_ev, n_bap1 / n_ev))
rows.append(dict(layer="DNA", metric="BAP1_mut_prevalence", value=n_bap1 / n_ev,
                 n_num=n_bap1, n_den=n_ev, detail="BAP1 mutation prevalence"))

cd = cnv.get("CDKN2A", {})
n_del = len([b for b in ccols if b in cd and cd[b] <= 0])
n_cnv = len([b for b in ccols if b in cd])
print("  CDKN2A hom-del        : %d / %d = %.3f" % (n_del, n_cnv, n_del / n_cnv))
rows.append(dict(layer="DNA", metric="CDKN2A_homdel_prevalence", value=n_del / n_cnv,
                 n_num=n_del, n_den=n_cnv, detail="CDKN2A homozygous deletion prevalence"))

# ---------- Layer 2: RNA ----------
print("\n[Layer 2 RNA] BAP1 mRNA, mut vs wt")
for G in ["BAP1", "NF2", "SETD2", "TP53", "MTAP"]:
    if G not in expr:
        continue
    s = mut.get(G, set())
    a = [expr[G][i] for i, b in enumerate(bars) if b in has_maf and b in s]
    b_ = [expr[G][i] for i, b in enumerate(bars) if b in has_maf and b not in s]
    if len(a) < 3 or len(b_) < 3:
        print("  %-8s skipped (mut n=%d, wt n=%d)" % (G, len(a), len(b_)))
        continue
    z, p = mw(a, b_)
    dm = statistics.mean(a) - statistics.mean(b_)
    print("  %-8s mut=%7.3f wt=%7.3f  delta=%+7.3f  p=%.3e  (n=%d vs %d)" % (
        G, statistics.mean(a), statistics.mean(b_), dm, p, len(a), len(b_)))
    rows.append(dict(layer="RNA", metric="%s_mRNA_mut_vs_wt" % G, value=dm,
                     p_value=p, n_num=len(a), n_den=len(b_),
                     detail="%s mRNA difference (mutant - wildtype)" % G))

# BAP1 immune signature scores at RNA level
sc_cd8 = score(SIG["CD8_T"], None)
s = mut.get("BAP1", set())
a = [sc_cd8[b] for b in bars if b in has_maf and b in s]
b_ = [sc_cd8[b] for b in bars if b in has_maf and b not in s]
z, p = mw(a, b_)
allsd = statistics.pstdev(list(sc_cd8.values()))
d = (statistics.mean(a) - statistics.mean(b_)) / allsd
print("  BAP1mut x CD8_T score : d=%+.4f p=%.3f" % (d, p))
# The statistic is a DNA-layer contrast (mutation status), even though the score it
# is tested against is built from RNA. The old layer label here said "RNA", which
# contradicted the manuscript, which calls this the DNA layer.
rows.append(dict(layer="DNA", metric="BAP1_mut_x_CD8T_d", value=d, p_value=p,
                 n_num=len(a), n_den=len(b_), detail="BAP1 mutation vs CD8 T-cell score"))

# RNA layer, continuous: does BAP1 mRNA itself track the CD8 T-cell score?
# This is the RNA-layer counterpart of the protein-layer Spearman below, so the
# three layers can be shown on one axis without any placeholder value.
bap1_expr = {bars[j]: expr["BAP1"][j] for j in range(len(bars))}
common_rna = [b for b in bars if b in bap1_expr and b in sc_cd8]
r_rna, p_rna = spearman([bap1_expr[b] for b in common_rna],
                        [sc_cd8[b] for b in common_rna])
print("  BAP1 mRNA vs CD8_T     : rho=%+.3f p=%.3f n=%d" % (r_rna, p_rna, len(common_rna)))
rows.append(dict(layer="RNA", metric="BAP1_mRNA_vs_CD8T_rho", value=r_rna,
                 p_value=p_rna, n_num=len(common_rna), n_den=len(common_rna),
                 detail="BAP1 mRNA vs CD8 T-cell score (Spearman)"))

# ---------- Layer 3: protein ----------
print("\n[Layer 3 protein] RPPA")
bap1_pep = [k for k in rppa if "BAP1" in k.upper()]
print("  BAP1 peptides on panel: %s" % bap1_pep)
if bap1_pep:
    pv = rppa[bap1_pep[0]]
    for nm, gl in SIG.items():
        sc = score(gl, None)
        common = [b for b in bars if b in pv and b in sc]
        if len(common) < 20:
            continue
        r, p = spearman([pv[b] for b in common], [sc[b] for b in common])
        print("  BAP1_prot vs %-11s rho=%+.3f p=%.3f n=%d" % (nm, r, p, len(common)))
        rows.append(dict(layer="Protein", metric="BAP1_prot_vs_%s_rho" % nm, value=r,
                         p_value=p, n_num=len(common), n_den=len(common),
                         detail="BAP1 protein vs %s score (Spearman)" % nm))

# CDKN2A-adjacent: is CDKN2A on the RPPA panel at all?
print("  CDKN2A peptides on panel: %s" % [k for k in rppa if "CDKN2A" in k.upper()])
print("  MTAP   peptides on panel: %s" % [k for k in rppa if "MTAP" in k.upper()])

# ---------- Layer 3: detectability of the protein layer ----------
# The protein layer is the smallest of the three (RPPA covers 62 of 87 samples), so
# its detection limit is stricter than the DNA layer's. Two numbers are needed to say
# whether it can address the published medium-sized association: the smallest
# standardised effect a 62-sample two-group comparison resolves, and the smallest
# correlation 62 samples resolve (the protein-layer statistic is a Spearman rho).
Z975 = 1.959963984540054


def _phi(x):
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def _power_d(dd, n1, n2):
    lam = dd * math.sqrt(n1 * n2 / (n1 + n2))
    return 1 - _phi(Z975 - lam) + _phi(-Z975 - lam)


def _dmin(n, p):
    n1 = max(2, int(round(n * p)))
    n2 = n - n1
    lo, hi = 0.01, 8.0
    for _ in range(200):
        m = (lo + hi) / 2
        if _power_d(m, n1, n2) >= 0.80:
            hi = m
        else:
            lo = m
    return hi, n1, n2


def _power_r(r, n):
    if r >= 1:
        return 1.0
    zr = 0.5 * math.log((1 + r) / (1 - r))
    s = math.sqrt(n - 3) * zr
    return 1 - _phi(Z975 - s) + _phi(-Z975 - s)


def _rmin(n):
    lo, hi = 0.001, 0.99
    for _ in range(200):
        m = (lo + hi) / 2
        if _power_r(m, n) >= 0.80:
            hi = m
        else:
            lo = m
    return hi


n_rppa = len([b for b in bars if b in pv]) if bap1_pep else 0
if n_rppa >= 20:
    p_rppa = n_bap1 / n_ev  # BAP1 prevalence carries over to the protein subset
    dm_prot, n1p, n2p = _dmin(n_rppa, p_rppa)
    r_min = _rmin(n_rppa)
    # the correlation a medium (d = 0.5) two-group effect implies, via the exact
    # point-biserial identity, at the DNA layer's group sizes
    delta, nA, nB = 0.5, 22, 56
    ntot = nA + nB
    tot_sd = math.sqrt(((ntot - 2) + delta * delta * nA * nB / ntot) / (ntot - 1))
    r_medium = abs(delta) / tot_sd * math.sqrt(nA * nB) / ntot
    print("\n[Layer 3 detectability] RPPA n=%d, BAP1 prevalence %.3f" % (n_rppa, p_rppa))
    print("  protein-layer d_min            : %.3f  (n1=%d n2=%d)" % (dm_prot, n1p, n2p))
    print("  min detectable Spearman rho    : %.3f" % r_min)
    print("  rho implied by a medium effect : %.3f" % r_medium)
    rows.append(dict(layer="Protein", metric="RPPA_dmin_pBAP1", value=dm_prot,
                     n_num=n1p, n_den=n2p,
                     detail="protein-layer detection limit at BAP1 prevalence, n=%d" % n_rppa))
    rows.append(dict(layer="Protein", metric="RPPA_min_detectable_rho", value=r_min,
                     n_num=n_rppa, n_den=n_rppa,
                     detail="smallest Spearman rho detectable at 80%% power, n=%d" % n_rppa))
    rows.append(dict(layer="Protein", metric="rho_implied_by_medium_effect",
                     value=r_medium, n_num=nA, n_den=nB,
                     detail="rho implied by a d=0.5 two-group effect at %d vs %d" % (nA, nB)))

# ---------- Reproduction ----------
print("\n[Reproduction]")
det = n_del / n_cnv
pub = 0.415
print("  CDKN2A del rate: ours %.4f vs published %.4f  (abs diff %.4f)" % (det, pub, abs(det - pub)))
rows.append(dict(layer="Reproduction", metric="CDKN2A_del_rate_vs_published",
                 value=det, n_num=n_del, n_den=n_cnv,
                 detail="ours %.4f; published %.4f (n=89, JECCR 2025); abs diff %.4f" % (det, pub, abs(det - pub))))

with open(os.path.join(OUT, "B8_fig4_source.tsv"), "w", encoding="utf-8") as f:
    f.write("layer\tmetric\tvalue\tp_value\tn_num\tn_den\tdetail\n")
    for r in rows:
        f.write("%s\t%s\t%.6g\t%s\t%d\t%d\t%s\n" % (
            r["layer"], r["metric"], r["value"],
            ("%.3e" % r["p_value"]) if r.get("p_value") is not None else "NA",
            r["n_num"], r["n_den"], r["detail"]))
print("\n[saved] B8_fig4_source.tsv  (%d rows)" % len(rows))
