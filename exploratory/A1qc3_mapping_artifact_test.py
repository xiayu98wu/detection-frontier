"""A1-QC3: is the low mRNA-protein correlation real, or an artefact of loose
antibody->gene mapping?

Two checks:
  (1) exact-name subset: antibodies whose name IS the gene symbol, no alias mapping
  (2) permutation control: shuffle sample labels to establish the null
  (3) literature calibration: within-cohort mRNA-protein correlations are typically
     0.2-0.4; the 0.51 RPPA500 figure is a PAN-CANCER value (pools biological
     variance across 32 tumour types), so it is not the right comparator.
"""
import gzip, os, math, statistics, random
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")

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


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    dx = math.sqrt(sum((x[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((y[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx > 0 and dy > 0 else None


# ---------- (1) exact-name subset ----------
gene_set = set(expr)
exact = [ab for ab in rppa if ab in gene_set]
print("=" * 84)
print("1. EXACT-NAME antibodies (no alias mapping at all): %d of %d" % (len(exact), len(rppa)))
print("=" * 84)
cors = []
detail = []
for ab in exact:
    vals = rppa[ab]
    pairs = [(vals[i], expr[ab][eidx[s]]) for i, s in enumerate(rsamps) if vals[i] is not None and s in eidx]
    if len(pairs) < 25:
        continue
    r = pearson([a for a, _ in pairs], [b for _, b in pairs])
    if r is not None:
        cors.append(r)
        detail.append((ab, r, len(pairs)))
print("  n=%d antibodies with >=25 samples" % len(cors))
print("  mean r = %+.3f   median = %+.3f" % (statistics.mean(cors), statistics.median(cors)))
detail.sort(key=lambda t: -t[1])
print("\n  top 10:")
for ab, r, n in detail[:10]:
    print("    %-14s r=%+.3f  n=%d" % (ab, r, n))
print("  bottom 10:")
for ab, r, n in detail[-10:]:
    print("    %-14s r=%+.3f  n=%d" % (ab, r, n))

# ---------- (2) permutation null ----------
print("\n" + "=" * 84)
print("2. PERMUTATION NULL (shuffle sample labels within each antibody)")
print("=" * 84)
random.seed(2026)
null_means = []
for _ in range(200):
    nc = []
    for ab in exact:
        vals = rppa[ab]
        pr = [(vals[i], expr[ab][eidx[s]]) for i, s in enumerate(rsamps) if vals[i] is not None and s in eidx]
        if len(pr) < 25:
            continue
        a = [x for x, _ in pr]
        b = [y for _, y in pr]
        random.shuffle(b)
        r = pearson(a, b)
        if r is not None:
            nc.append(r)
    if nc:
        null_means.append(statistics.mean(nc))
print("  null mean r: %.4f  sd=%.4f  max |mean| = %.4f" % (
    statistics.mean(null_means), statistics.pstdev(null_means),
    max(abs(x) for x in null_means)))
print("  observed mean r = %+.3f" % statistics.mean(cors))
print("  -> %s" % ("observed is FAR outside the null: real biological coupling"
                   if abs(statistics.mean(cors)) > 5 * statistics.pstdev(null_means)
                   else "observation is within null range - no coupling detected"))

# ---------- (3) what does the literature actually say ----------
print("\n" + "=" * 84)
print("3. CALIBRATION: what is the right comparator?")
print("=" * 84)
print("  RPPA500 pan-cancer figure (0.51): computed by pooling samples ACROSS")
print("    32 tumour types. Between-cancer variance inflates the correlation.")
print("  Within-cohort mRNA-protein correlations reported in the literature")
print("    (CPTAC, TCGA per-cancer analyses) typically fall in 0.2-0.4.")
print("  So the correct test is: does MESO differ from OTHER COHORTS computed")
print("    with THIS SAME pipeline? (cross-cohort job running separately)")
print("  raw MESO mean r here = %+.3f" % statistics.mean(cors))

# ---------- (4) how does this compare to the CNV layer ----------
print("\n" + "=" * 84)
print("4. CONTEXT: within the same samples, how strong is CNV->mRNA?")
print("=" * 84)
cnv = {}
with gzip.open(os.path.join(DM, "MESO_cnv_matrix.tsv.gz"), "rt", encoding="utf-8") as f:
    hh = f.readline().rstrip("\n").split("\t")
    ccols = [c.strip('"') for c in hh[1:]]
    for l in f:
        c = l.rstrip("\n").split("\t")
        d = {}
        for i, b in enumerate(ccols):
            if c[i + 1] != "NA":
                d[b] = int(float(c[i + 1]))
        cnv[c[0]] = d
cm = []
for g in exact:
    if g not in cnv:
        continue
    pr = [(cnv[g][s], expr[g][eidx[s]]) for s in ecols if s in cnv[g]]
    if len(pr) < 25:
        continue
    r = pearson([a for a, _ in pr], [b for _, b in pr])
    if r is not None:
        cm.append(r)
print("  mean r(CNV, mRNA) for the same %d genes = %+.3f" % (len(cm), statistics.mean(cm)))
print("  mean r(mRNA, protein)                     = %+.3f" % statistics.mean(cors))
print("  -> the CNV->mRNA step is the stronger link by %.1fx" % (
    statistics.mean(cm) / statistics.mean(cors) if statistics.mean(cors) else float("nan")))
