"""Fixed specificity null + histology-stratified survival."""
import gzip, os, collections, statistics, math, random
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")

with gzip.open(os.path.join(DM, "MESO_expr_full.tsv.gz"), "rt", encoding="utf-8", errors="replace") as g:
    hdr = g.readline().rstrip("\n").split("\t")
    cols = [c.strip('"') for c in hdr[1:]]
    expr = {}
    for line in g:
        c = line.rstrip("\n").split("\t")
        expr[c[0].strip('"')] = [float(x) for x in c[1:]]
bars = [c[:12] for c in cols]

cnv = collections.defaultdict(dict)
with gzip.open(os.path.join(DM, "MESO_cnv_matrix.tsv.gz"), "rt", encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    cb = [c.strip('"') for c in h[1:]]
    for l in f:
        c = l.rstrip("\n").split("\t")
        for i, b in enumerate(cb):
            if c[i + 1] != "NA":
                cnv[c[0]][b] = int(float(c[i + 1]))

with open(os.path.join(DM, "MESO_clinicalMatrix.tsv"), encoding="utf-8", errors="replace") as f:
    ch = f.readline().rstrip("\n").split("\t")
    crows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
hi = ch.index("histological_type")
si = ch.index("sampleID")
hist = {r[si][:12]: r[hi] for r in crows if len(r) > hi}

DEL = {b: (1 if cnv["CDKN2A"].get(b) == 0 else 0) for b in bars}
EPI = [b for b in bars if "Epithelioid" in hist.get(b, "")]
NONEPI = [b for b in bars if "Epithelioid" not in hist.get(b, "")]

SIGS = {
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "CD3D", "CD3E"],
    "FIBROBLAST": ["COL1A1", "COL1A2", "FAP", "ACTA2", "DCN"],
    "EMT": ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1"],
    "PROLIF": ["MKI67", "TOP2A", "CCNB1", "BUB1", "AURKA"],
    "ANTIGEN": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB9"],
}


def score_vec(genes):
    gs = [g for g in genes if g in expr]
    return [statistics.mean(expr[g][j] for g in gs) for j in range(len(bars))]


def mw(a, b):
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return None, None
    comb = sorted([(v, 0) for v in a] + [(v, 1) for v in b])
    ranks = [0] * len(comb)
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
    den = (N * (N - 1))
    sd = math.sqrt(n1 * n2 / 12 * ((N + 1) - sum(t ** 3 - t for t in tie.values()) / den)) if den else 0
    if sd == 0:
        return None, None
    z = (U1 - n1 * n2 / 2) / sd
    return z, 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


# ---------- expression-matched null, fixed ----------
random.seed(2026)
gene_means = {g: statistics.mean(expr[g]) for g in expr}
ordered = sorted(gene_means, key=lambda g: gene_means[g])
N = len(ordered)
NBINS = 20
bin_of = {}
for idx, g in enumerate(ordered):
    bin_of[g] = min(NBINS - 1, int(NBINS * idx / N))


def null_set(real_genes, n=300):
    """Draw random genes with the same expression-bin composition as real_genes."""
    want = collections.Counter(bin_of[g] for g in real_genes if g in bin_of)
    pool = collections.defaultdict(list)
    for g in ordered:
        pool[bin_of[g]].append(g)
    out = []
    for _ in range(n):
        sel = []
        for bi, k in want.items():
            cand = [g for g in pool[bi] if g not in (real_genes or [])]
            if cand:
                sel += random.sample(cand, min(k, len(cand)))
        out.append(sel)
    return out


print("=" * 76)
print("SPECIFICITY NULL (expression-matched), within EPITHELIOID")
print("  observed z vs null distribution of |z| (n=300 random sets)")
print("=" * 76)
print("%-11s %8s %9s %8s %9s %8s" % ("signature", "obs_z", "null_mean", "null_max", "percentile", "null_p"))
A = [b for b in EPI if DEL[b] == 1]
B = [b for b in EPI if DEL[b] == 0]
for nm, gl in SIGS.items():
    gs = [g for g in gl if g in expr]
    sc = score_vec(gs)
    obs_z, _ = mw([sc[j] for j, b in enumerate(bars) if b in A], [sc[j] for j, b in enumerate(bars) if b in B])
    zs = []
    for gset in null_set(gs, 300):
        if len(gset) < 3:
            continue
        s2 = score_vec(gset)
        z2, _ = mw([s2[j] for j, b in enumerate(bars) if b in A], [s2[j] for j, b in enumerate(bars) if b in B])
        if z2 is not None:
            zs.append(z2)
    if not zs or obs_z is None:
        print("%-11s (null failed)" % nm)
        continue
    absz = [abs(z) for z in zs]
    pct = 100 * sum(1 for z in absz if z < abs(obs_z)) / len(absz)
    np_ = sum(1 for z in absz if z >= abs(obs_z)) / len(absz)
    print("%-11s %+8.2f %9.2f %8.2f %8.0f%% %8.3f" % (nm, obs_z, statistics.mean(absz), max(absz), pct, np_))

# ---------- survival stratified by histology ----------
print("\n" + "=" * 76)
print("SURVIVAL: 9p21 deletion, stratified by histology")
print("=" * 76)
sv = os.path.join(DM, "MESO_survival.txt")
with open(sv, encoding="utf-8", errors="replace") as f:
    sh = f.readline().rstrip("\n").split("\t")
    srows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
idx = {k: i for i, k in enumerate(sh)}
oi, ti = idx.get("OS"), idx.get("OS.time")
surv = {}
for r in srows:
    bc = r[0][:12]
    try:
        surv[bc] = (int(float(r[oi])), float(r[ti]))
    except Exception:
        pass
print("survival records: %d" % len(surv))


def logrank(ev):
    ev = sorted(ev, key=lambda t: t[1])
    n1 = sum(1 for _, _, g in ev if g == 1)
    n2 = len(ev) - n1
    if n1 < 3 or n2 < 3:
        return None
    O1 = E1 = V = 0.0
    at1, at2 = n1, n2
    for e, t, g in ev:
        tot = at1 + at2
        if tot <= 0:
            break
        if e == 1:
            p1 = at1 / tot
            O1 += (1 if g == 1 else 0)
            E1 += p1
            V += p1 * (1 - p1) * (tot - 1) / max(1e-9, tot)
        if g == 1:
            at1 -= 1
        else:
            at2 -= 1
    chi = (O1 - E1) ** 2 / V if V > 0 else 0
    p = math.exp(-chi / 2) if chi > 0 else 1.0
    return chi, p, O1, E1, n1, n2


for label, cohort in [("ALL (n=%d)" % len(bars), bars),
                      ("EPITHELIOID ONLY (n=%d)" % len(EPI), EPI)]:
    ev = [(surv[b][0], surv[b][1], DEL[b]) for b in cohort if b in surv and surv[b][1] > 0]
    r = logrank(ev)
    if r is None:
        print("\n%s: too few events" % label)
        continue
    chi, p, O1, E1, n1, n2 = r
    print("\n%s" % label)
    print("  del n=%d  intact n=%d" % (n1, n2))
    print("  log-rank chi2=%.3f  p=%.5f  (obs deaths del=%.1f, expected=%.1f)" % (chi, p, O1, E1))
    md = statistics.median([t for e, t, g in ev if g == 1])
    mi = statistics.median([t for e, t, g in ev if g == 0])
    print("  median OS: del=%.0f d   intact=%.0f d   ratio=%.2f" % (md, mi, md / mi))

# multivariable-ish: check histology distribution within epithelioid-only del/intact
print("\n  histology composition of deletion groups:")
for lab, cohort in [("ALL", bars), ("EPI", EPI)]:
    dd = collections.Counter(hist.get(b, "?") for b in cohort if DEL[b] == 1)
    ii = collections.Counter(hist.get(b, "?") for b in cohort if DEL[b] == 0)
    print("   %s del: %s" % (lab, dict(dd)))
    print("   %s int: %s" % (lab, dict(ii)))
