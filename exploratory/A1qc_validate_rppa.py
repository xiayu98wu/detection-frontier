"""A1-QC: validate the reconstructed RPPA matrix before trusting any conclusion.

1. Pearson vs Spearman (the published pan-cancer benchmark 0.51 is Pearson)
2. positive control: CDKN2A CNV -> mRNA (must be strongly negative; not on RPPA)
3. cross-check against Xena's official MESO RPPA matrix
"""
import gzip, os, math, statistics, urllib.request, io
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")

# ---------- my reconstructed RPPA ----------
rppa = {}
with open(os.path.join(DM, "MESO_rppa.tsv"), encoding="utf-8") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    rsamps = hdr[1:]
    for l in f:
        c = l.rstrip("\n").split("\t")
        rppa[c[0]] = [float(v) if v != "NA" else None for v in c[1:]]

# ---------- expression / CNV ----------
with gzip.open(os.path.join(DM, "MESO_expr_full.tsv.gz"), "rt", encoding="utf-8", errors="replace") as g:
    h = g.readline().rstrip("\n").split("\t")
    ecols = [c.strip('"')[:12] for c in h[1:]]
    expr = {}
    for line in g:
        c = line.rstrip("\n").split("\t")
        expr[c[0].strip('"')] = [float(x) for x in c[1:]]

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


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    dx = math.sqrt(sum((x[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((y[i] - my) ** 2 for i in range(n)))
    if dx == 0 or dy == 0:
        return None
    return num / (dx * dy)


def spearman(x, y):
    n = len(x)
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
    return pearson(rx, ry)


print("=" * 88)
print("QC-1. CDKN2A positive control: CNV -> mRNA (not on RPPA, tests the CNV matrix itself)")
print("=" * 88)
for gene in ["CDKN2A", "CDKN2B", "MTAP", "BAP1", "NF2"]:
    if gene not in cnv or gene not in expr:
        print("  %-8s missing" % gene)
        continue
    idx = {b: i for i, b in enumerate(ecols)}
    pairs = [(cnv[gene][b], expr[gene][idx[b]]) for b in ecols if b in cnv[gene]]
    pr = pearson([x for x, _ in pairs], [y for _, y in pairs])
    sp = spearman([x for x, _ in pairs], [y for _, y in pairs])
    print("  %-8s n=%2d  pearson=%+.3f  spearman=%+.3f" % (gene, len(pairs), pr, sp))

# ---------- Xena official RPPA ----------
print("\n" + "=" * 88)
print("QC-2. cross-check against Xena official MESO RPPA")
print("=" * 88)
B = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"
xena = None
for u in ["TCGA.MESO.sampleMap%2FRPPA",
          "TCGA.MESO.sampleMap%2FMESO_RPPA",
          "TCGA.MESO.sampleMap%2Fmeso_rppa"]:
    try:
        req = urllib.request.Request(B + u, headers={"User-Agent": "Mozilla/5.0"})
        raw = urllib.request.urlopen(req, timeout=60).read()
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        txt = raw.decode("utf-8", "replace")
        first = txt.split("\n", 1)[0]
        if "sample" in first.lower() or "\t" in first:
            xena = txt
            print("  [ok] %s  (%.1f KB)" % (u, len(raw) / 1024))
            break
    except Exception as e:
        print("  [--] %-40s %s" % (u, type(e).__name__))

if xena:
    lines = [l for l in xena.split("\n") if l.strip()]
    xh = lines[0].split("\t")
    print("  Xena RPPA header cols: %d, first 3: %s" % (len(xh), xh[:3]))
    xdata = {}
    for l in lines[1:]:
        c = l.split("\t")
        xdata[c[0]] = c[1:]
    print("  Xena RPPA rows (antibodies): %d" % len(xdata))
    # compare overlapping antibodies
    mine = set(rppa)
    theirs = set(xdata)
    ov = sorted(mine & theirs)
    print("  antibody overlap: %d" % len(ov))
    # correlate shared samples
    xs = xh[1:]
    xb = [s[:12] for s in xs]
    ncmp = 0
    cors = []
    for ab in ov[:400]:
        mv = rppa[ab]
        tv = xdata[ab]
        pairs = []
        for i, s in enumerate(xs):
            b = s[:12]
            if b in rsamps and tv[i] not in ("", "NA"):
                try:
                    pairs.append((mv[rsamps.index(b)], float(tv[i])))
                except Exception:
                    pass
        pairs = [(a, c) for a, c in pairs if a is not None]
        if len(pairs) >= 20:
            r = pearson([a for a, _ in pairs], [c for _, c in pairs])
            if r is not None:
                cors.append(r)
                ncmp += 1
    if cors:
        print("  antibodies compared: %d   mean Pearson(mine, Xena) = %.4f" % (ncmp, statistics.mean(cors)))
        print("  min=%.3f  max=%.3f  frac>0.95=%.2f" % (min(cors), max(cors),
              sum(1 for r in cors if r > 0.95) / len(cors)))
else:
    print("  Xena RPPA not reachable; skipping cross-check")

# ---------- QC-3: Pearson vs Spearman for mRNA-protein ----------
print("\n" + "=" * 88)
print("QC-3. mRNA-protein coupling: Pearson vs Spearman (benchmark 0.51 is Pearson)")
print("=" * 88)
ALIAS = {"BAP1C4": "BAP1", "PDL1": "CD274", "P53": "TP53", "VIMENTIN": "VIM",
         "ECADHERIN": "CDH1", "E_CADHERIN": "CDH1", "NCADHERIN": "CDH2",
         "CYCLIND1": "CCND1", "CYCLINB1": "CCNB1", "CYCLINE1": "CCNE1",
         "CD45": "PTPRC", "CD31": "PECAM1", "HER2": "ERBB2", "HER3": "ERBB3",
         "4EBP1": "EIF4EBP1", "RB": "RB1", "P21": "CDKN1A", "P27": "CDKN1B",
         "SNAIL": "SNAI1", "TWIST": "TWIST1", "IDO": "IDO1", "ZAP70": "ZAP70",
         "ZAP-70": "ZAP70", "GRANZYMEB": "GZMB", "CD8": "CD8A", "CD4": "CD4"}
SHARED = [s for s in ecols if s in rsamps]
eidx = {b: i for i, b in enumerate(ecols)}
ps, ss = [], []
for ab, vals in rppa.items():
    g = ALIAS.get(ab, ab)
    if g not in expr:
        continue
    pairs = []
    for s, v in zip(rsamps, vals):
        if v is not None and s in eidx:
            pairs.append((v, expr[g][eidx[s]]))
    if len(pairs) < 25:
        continue
    p = pearson([a for a, _ in pairs], [b for _, b in pairs])
    sp = spearman([a for a, _ in pairs], [b for _, b in pairs])
    if p is not None:
        ps.append(p)
        ss.append(sp)
print("  matched antibodies: %d" % len(ps))
print("  mean Pearson  r(mRNA, protein) = %+.3f  (median %+.3f)" % (statistics.mean(ps), statistics.median(ps)))
print("  mean Spearman r(mRNA, protein) = %+.3f  (median %+.3f)" % (statistics.mean(ss), statistics.median(ss)))
print("  Pan-cancer RPPA500 benchmark (Pearson) = +0.51")
