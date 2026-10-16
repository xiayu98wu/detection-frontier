import gzip, os, collections, statistics, math
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

print("=" * 76)
print("CONFOUND CHECK: is CDKN2A deletion independent of histology?")
print("=" * 76)
tt = collections.Counter()
td = collections.Counter()
for b in bars:
    h = hist.get(b, "?")
    tt[h] += 1
    if cnv["CDKN2A"].get(b) == 0:
        td[h] += 1
print("\n%-40s %4s %11s %7s" % ("histology", "n", "CDKN2A del", "rate"))
for h, n in tt.most_common():
    print("%-40s %4d %11d %6.1f%%" % (h[:40], n, td[h], 100 * td[h] / n))

SIGS = {
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "CD3D", "CD3E"],
    "CHECKPT": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "IFN_G": ["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG", "IDO1", "GBP1"],
    "FIBROBLAST": ["COL1A1", "COL1A2", "FAP", "ACTA2", "DCN"],
    "EMT": ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1"],
    "PROLIF": ["MKI67", "TOP2A", "CCNB1", "BUB1", "AURKA"],
    "ANTIGEN": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB9"],
}


def score(gl):
    gs = [g for g in gl if g in expr]
    return {bars[j]: statistics.mean(expr[g][j] for g in gs) for j in range(len(bars))}


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
    sd = math.sqrt(n1 * n2 / 12 * ((N + 1) - sum(t ** 3 - t for t in tie.values()) / (N * (N - 1))))
    if sd == 0:
        return None, None
    z = (U1 - n1 * n2 / 2) / sd
    return z, 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


EPI = [b for b in bars if "Epithelioid" in hist.get(b, "")]
print("\n" + "=" * 76)
print("CDKN2A HOMDEL WITHIN EPITHELIOID ONLY (n=%d)  <- confound removed" % len(EPI))
print("=" * 76)
A = [b for b in EPI if cnv["CDKN2A"].get(b) == 0]
B = [b for b in EPI if cnv["CDKN2A"].get(b) == 2]
print("  epithelioid: del n=%d   diploid n=%d" % (len(A), len(B)))
print("%-12s %9s %9s %8s %7s %10s" % ("signature", "del", "diploid", "delta", "z", "p"))
for nm, gl in SIGS.items():
    sc = score(gl)
    z, p = mw([sc[b] for b in A], [sc[b] for b in B])
    if z is None:
        continue
    flag = "  ***" if p < 0.001 else ("  **" if p < 0.01 else ("  *" if p < 0.05 else ""))
    ma, mb = statistics.mean([sc[b] for b in A]), statistics.mean([sc[b] for b in B])
    print("%-12s %9.3f %9.3f %+8.3f %+7.2f %10.2e%s" % (nm, ma, mb, ma - mb, z, p, flag))
