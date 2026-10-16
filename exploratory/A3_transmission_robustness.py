"""A3: is the CNV-driven transmission ratio T a robust, interpretable quantity?

A1 found mean T = 0.41 for 136 CNV-driven genes. Before investing in cross-cohort
CNV downloads, test whether T is:
  (1) stable under re-sampling (bootstrap)
  (2) distinguishable from a null where CNV and protein are linked only through mRNA
  (3) driven by a few outliers or a broad property

The null (2) is the important one: if CNV affects protein ONLY via mRNA, then
r(CNV,prot) should equal r(CNV,mRNA) * r(mRNA,prot), i.e. T should equal r(mRNA,prot).
"""
import gzip, os, collections, statistics, math, random
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")
OUT = os.path.join(_ROOT, "results")
os.makedirs(OUT, exist_ok=True)

rppa = {}
with open(os.path.join(DM, "MESO_rppa.tsv"), encoding="utf-8") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    rsamps = hdr[1:]
    for l in f:
        c = l.rstrip("\n").split("\t")
        rppa[c[0]] = [float(v) if v != "NA" else None for v in c[1:]]

with gzip.open(os.path.join(DM, "MESO_expr_full.tsv.gz"), "rt", encoding="utf-8", errors="replace") as g:
    h = g.readline().rstrip("\n").split("\t")
    ecols = [c.strip('"')[:12] for c in h[1:]]
    expr = {}
    for line in g:
        c = line.rstrip("\n").split("\t")
        expr[c[0].strip('"')] = [float(x) for x in c[1:]]
eidx = {b: i for i, b in enumerate(ecols)}

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

ALIAS = {"BAP1C4": "BAP1", "PDL1": "CD274", "P53": "TP53", "VIMENTIN": "VIM"}
GENES = sorted(g for g in rppa if ALIAS.get(g, g) in expr and ALIAS.get(g, g) in cnv
               and sum(1 for v in rppa[g] if v is not None) >= 25)
print("three-layer genes (exact/alias): %d" % len(GENES))


def pearson(x, y):
    n = len(x)
    if n < 12:
        return None
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    dx = math.sqrt(sum((x[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((y[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx > 0 and dy > 0 else None


rows = []
for ab in GENES:
    g = ALIAS.get(ab, ab)
    vals = rppa[ab]
    rec = []
    for i, s in enumerate(rsamps):
        if vals[i] is None or s not in eidx or s not in cnv[g]:
            continue
        rec.append((cnv[g][s], expr[g][eidx[s]], vals[i]))
    if len(rec) < 25:
        continue
    a = pearson([r[0] for r in rec], [r[1] for r in rec])   # cnv -> mrna
    b = pearson([r[1] for r in rec], [r[2] for r in rec])   # mrna -> prot
    c = pearson([r[0] for r in rec], [r[2] for r in rec])   # cnv -> prot
    if a is None or b is None or c is None:
        continue
    rows.append({"gene": g, "n": len(rec), "a": a, "b": b, "c": c,
                 "cnv_sd": statistics.pstdev([r[0] for r in rec])})
print("genes with all three correlations: %d" % len(rows))

# ---- (2) THE KEY NULL ----
print("\n" + "=" * 92)
print("1. NULL MODEL: if CNV affects protein ONLY through mRNA, then c = a*b (T = b)")
print("=" * 92)
DRIVEN = [r for r in rows if r["a"] >= 0.20]
print("  CNV-driven subset (a >= 0.20): %d genes" % len(DRIVEN))
obs_T = [r["c"] / r["a"] for r in DRIVEN if r["a"] > 0]
pred_T = [r["b"] for r in DRIVEN]
print("  observed  T = c/a : mean=%+.3f median=%+.3f" % (statistics.mean(obs_T), statistics.median(obs_T)))
print("  predicted T = b   : mean=%+.3f median=%+.3f" % (statistics.mean(pred_T), statistics.median(pred_T)))
print("  median absolute deviation between observed and predicted: %.3f" % statistics.median(
    [abs(o - p) for o, p in zip(obs_T, pred_T)]))
over = sum(1 for o, p in zip(obs_T, pred_T) if o > p + 0.15)
under = sum(1 for o, p in zip(obs_T, pred_T) if o < p - 0.15)
print("  genes where observed T exceeds prediction by >0.15 : %d" % over)
print("  genes where observed T falls below prediction by >0.15: %d" % under)

# ---- (1) bootstrap stability ----
print("\n" + "=" * 92)
print("2. BOOTSTRAP STABILITY of the mean transmission ratio (1000 resamples of genes)")
print("=" * 92)
random.seed(2026)
means = []
for _ in range(1000):
    samp = [random.choice(DRIVEN) for _ in DRIVEN]
    ts = [r["c"] / r["a"] for r in samp if r["a"] > 0]
    if ts:
        means.append(statistics.mean(ts))
means.sort()
print("  mean T = %.3f   95%% CI = [%.3f, %.3f]" % (
    statistics.mean(means), means[25], means[975]))

# ---- (3) outlier sensitivity ----
print("\n" + "=" * 92)
print("3. OUTLIER SENSITIVITY")
print("=" * 92)
ts = sorted([(r["c"] / r["a"], r["gene"]) for r in DRIVEN if r["a"] > 0], key=lambda t: -t[0])
print("  top 5 (most conducted): %s" % ", ".join("%s=%.2f" % (g, t) for t, g in ts[:5]))
print("  bottom 5 (most buffered): %s" % ", ".join("%s=%.2f" % (g, t) for t, g in ts[-5:]))
trimmed = [t for t, _ in ts[5:-5]]
print("  mean T, 5%% trimmed each tail = %.3f (vs %.3f untrimmed)" % (
    statistics.mean(trimmed), statistics.mean([t for t, _ in ts])))

# ---- how many are genuinely buffered vs conducted ----
print("\n" + "=" * 92)
print("4. CONDUCTED vs BUFFERED CLASSIFICATION (threshold T=0.5)")
print("=" * 92)
cond = [r for r in DRIVEN if r["c"] / r["a"] >= 0.5] if all(r["a"] > 0 for r in DRIVEN) else []
cond = [r for r in DRIVEN if r["a"] > 0 and r["c"] / r["a"] >= 0.5]
buff = [r for r in DRIVEN if r["a"] > 0 and r["c"] / r["a"] < 0.5]
print("  conducted (T>=0.5): %d genes -> %s" % (len(cond), ", ".join(r["gene"] for r in cond[:25])))
print("  buffered  (T<0.5) : %d genes -> %s" % (len(buff), ", ".join(r["gene"] for r in buff[:25])))

# do conducted genes have a common property?
print("\n  property comparison (conducted vs buffered):")
print("    %-14s %10s %10s %10s" % ("property", "conducted", "buffered", "p"))
for prop in ["a", "b", "cnv_sd"]:
    x = [r[prop] for r in cond]
    y = [r[prop] for r in buff]
    if not x or not y:
        continue
    # Mann-Whitney
    comb = sorted([(v, 0) for v in x] + [(v, 1) for v in y])
    ranks = [0.0] * len(comb)
    i = 0
    while i < len(comb):
        j = i
        while j + 1 < len(comb) and comb[j + 1][0] == comb[i][0]:
            j += 1
        rr = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[k] = rr
        i = j + 1
    n1, n2 = len(x), len(y)
    U1 = sum(ranks[k] for k in range(len(comb)) if comb[k][1] == 0) - n1 * (n1 + 1) / 2
    tie = collections.Counter(v for v, _ in comb)
    N = n1 + n2
    den = N * (N - 1)
    sd = math.sqrt(n1 * n2 / 12 * ((N + 1) - sum(t ** 3 - t for t in tie.values()) / den)) if den else 0
    z = (U1 - n1 * n2 / 2) / sd if sd else 0
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    print("    %-14s %10.3f %10.3f %10.2e" % (prop, statistics.mean(x), statistics.mean(y), p))

with open(os.path.join(OUT, "A3_transmission_robustness.tsv"), "w", encoding="utf-8") as f:
    f.write("gene\tn\tr_cnv_mrna\tr_mrna_prot\tr_cnv_prot\tT\tcnv_sd\n")
    for r in sorted(rows, key=lambda x: -(x["c"] / x["a"]) if x["a"] > 0 else 0):
        T = r["c"] / r["a"] if r["a"] > 0 else float("nan")
        f.write("%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (
            r["gene"], r["n"], r["a"], r["b"], r["c"], T, r["cnv_sd"]))
print("\n[saved] A3_transmission_robustness.tsv")
