"""Track B module B4: consolidate every driver x phenotype test run in the
TCGA-MESO case study into one empirical table.

Purpose: quantify the multiplicity / winner's-curse picture. If a study scans
D drivers x S signatures in an underpowered cohort, how many "hits" appear by
chance, and what effect size would a reported hit actually imply?
"""
import gzip, os, collections, statistics, math, random
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")
OUT = os.path.join(_ROOT, "results")
os.makedirs(OUT, exist_ok=True)

# ---------------- load ----------------
with gzip.open(os.path.join(DM, "MESO_expr_full.tsv.gz"), "rt", encoding="utf-8", errors="replace") as g:
    h = g.readline().rstrip("\n").split("\t")
    ecols = [c.strip('"')[:12] for c in h[1:]]
    expr = {}
    for line in g:
        c = line.rstrip("\n").split("\t")
        expr[c[0].strip('"')] = [float(x) for x in c[1:]]

cnv = collections.defaultdict(dict)
with gzip.open(os.path.join(DM, "MESO_cnv_matrix.tsv.gz"), "rt", encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    ccols = [c.strip('"') for c in h[1:]]
    for l in f:
        c = l.rstrip("\n").split("\t")
        for i, b in enumerate(ccols):
            if c[i + 1] != "NA":
                cnv[c[0]][b] = int(float(c[i + 1]))

mut = collections.defaultdict(set)
with open(os.path.join(DM, "MESO_mutations.maf"), encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    gi = h.index("Hugo_Symbol")
    ci = h.index("Tumor_Sample_Barcode")
    for l in f:
        c = l.rstrip("\n").split("\t")
        if len(c) > max(gi, ci):
            mut[c[gi]].add(c[ci][:12])
has_maf = set().union(*mut.values())

with open(os.path.join(DM, "MESO_clinicalMatrix.tsv"), encoding="utf-8", errors="replace") as f:
    ch = f.readline().rstrip("\n").split("\t")
    crows = [l.rstrip("\n").split("\t") for l in f if l.strip()]
hi = ch.index("histological_type")
si = ch.index("sampleID")
hist = {r[si][:12]: r[hi] for r in crows if len(r) > hi}

SIGS = {
    "CD8_T": ["CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "CD3D", "CD3E"],
    "CHECKPT": ["PDCD1", "CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT"],
    "IFN_G": ["STAT1", "IRF1", "CXCL9", "CXCL10", "IFNG", "IDO1", "GBP1"],
    "MACROPHAGE": ["CD68", "CD163", "MRC1", "CSF1R"],
    "FIBROBLAST": ["COL1A1", "COL1A2", "FAP", "ACTA2", "DCN"],
    "EMT": ["VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "TWIST1", "ZEB1"],
    "PROLIF": ["MKI67", "TOP2A", "CCNB1", "BUB1", "AURKA"],
    "ANTIGEN": ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "PSMB9"],
}


def score(gl):
    gs = [g for g in gl if g in expr]
    return [statistics.mean(expr[g][j] for g in gs) for j in range(len(ecols))]


def mw_p(a, b):
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return None
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
    sd = math.sqrt(n1 * n2 / 12 * ((N + 1) - sum(t ** 3 - t for t in tie.values()) / den)) if den else 0
    if sd == 0:
        return None
    z = (U1 - n1 * n2 / 2) / sd
    return 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))


PHI = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))


def power_d(d, n1, n2, alpha=0.05):
    if n1 < 2 or n2 < 2:
        return 0.0
    ncp = d * math.sqrt(n1 * n2 / (n1 + n2))
    return (1 - PHI(1.959963985 - ncp)) + PHI(-1.959963985 - ncp)


def d_min(n, p, power=0.80):
    n1 = max(2, int(round(n * p)))
    n2 = max(2, n - n1)
    lo, hi = 0.01, 8.0
    if power_d(hi, n1, n2) < power:
        return float("nan")
    for _ in range(120):
        mid = (lo + hi) / 2
        if power_d(mid, n1, n2) < power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ---------------- driver definitions ----------------
DRIVERS = {
    "BAP1_mut":      {"type": "mut", "gene": "BAP1"},
    "NF2_mut":       {"type": "mut", "gene": "NF2"},
    "TP53_mut":      {"type": "mut", "gene": "TP53"},
    "SETD2_mut":     {"type": "mut", "gene": "SETD2"},
    "CDKN2A_del":    {"type": "cnv0", "gene": "CDKN2A"},
    "MTAP_del":      {"type": "cnv0", "gene": "MTAP"},
    "CDKN2A_loss":   {"type": "cnvle1", "gene": "CDKN2A"},
    "BAP1_cnvloss":  {"type": "cnvle1", "gene": "BAP1"},
    "NF2_cnvloss":   {"type": "cnvle1", "gene": "NF2"},
}

rows = []
for dname, spec in DRIVERS.items():
    if spec["type"] == "mut":
        s = mut.get(spec["gene"], set())
        grp = {b: (1 if (b in has_maf and b in s) else 0) for b in ecols}
        unknown = [b for b in ecols if b not in has_maf]
    else:
        thr = 0 if spec["type"] == "cnv0" else 1
        cd = cnv.get(spec["gene"], {})
        grp = {b: (1 if (b in cd and cd[b] <= thr) else 0) for b in ecols}
        unknown = [b for b in ecols if b not in cd]
    A = [b for b in ecols if grp[b] == 1]
    B = [b for b in ecols if grp[b] == 0 and b not in unknown]
    n1, n2 = len(A), len(B)
    if n1 < 5 or n2 < 5:
        continue
    p_prev = n1 / (n1 + n2)
    dmin = d_min(n1 + n2, p_prev)

    for sn, gl in SIGS.items():
        sc = score(gl)
        a = [sc[i] for i, b in enumerate(ecols) if b in A]
        b_ = [sc[i] for i, b in enumerate(ecols) if b in B]
        allsd = statistics.pstdev(sc)
        d = (statistics.mean(a) - statistics.mean(b_)) / allsd if allsd > 0 else 0
        pv = mw_p(a, b_)
        pw = power_d(abs(d), n1, n2)
        rows.append({
            "driver": dname, "signature": sn, "n_altered": n1, "n_wt": n2,
            "prevalence": p_prev, "d_min_80": dmin, "observed_d": d,
            "power_at_observed": pw, "p_value": pv,
        })

# ---------------- report ----------------
print("=" * 100)
print("B4. ALL DRIVER x PHENOTYPE TESTS IN THE TCGA-MESO CASE STUDY")
print("=" * 100)
print("%-14s %-11s %5s %5s %8s %9s %9s %10s" % (
    "driver", "signature", "n_alt", "n_wt", "d_min80", "obs_d", "power", "p"))
print("-" * 100)
for r in sorted(rows, key=lambda x: x["p_value"] if x["p_value"] is not None else 1):
    flag = ""
    if r["p_value"] is not None:
        if r["p_value"] < 0.05 / len(rows):
            flag = "  **bonf"
        elif r["p_value"] < 0.05:
            flag = "  *"
    print("%-14s %-11s %5d %5d %8.2f %+9.3f %8.1f%% %10.2e%s" % (
        r["driver"], r["signature"], r["n_altered"], r["n_wt"],
        r["d_min_80"], r["observed_d"], 100 * r["power_at_observed"],
        r["p_value"] if r["p_value"] is not None else float("nan"), flag))

ntest = len(rows)
bonf = 0.05 / ntest
nhit = sum(1 for r in rows if r["p_value"] is not None and r["p_value"] < 0.05)
nbonf = sum(1 for r in rows if r["p_value"] is not None and r["p_value"] < bonf)
print("\n  total tests            : %d" % ntest)
print("  Bonferroni threshold   : %.2e" % bonf)
print("  raw p<0.05             : %d (expected by chance if all null: %.1f)" % (nhit, 0.05 * ntest))
print("  surviving Bonferroni   : %d" % nbonf)

# how many tests are even capable of detecting a moderate effect?
MODERATE = 0.40   # typical effect size for driver-immune associations in the literature
cap = sum(1 for r in rows if r["d_min_80"] <= MODERATE)
print("\n  tests with d_min(80%%) <= %.2f (i.e. able to detect a moderate effect): %d / %d" % (
    MODERATE, cap, ntest))
print("  -> %.0f%% of this study's comparisons are structurally incapable of" % (100 * (1 - cap / ntest)))
print("     detecting an effect of the size the literature typically reports.")

# winner's curse: if you report the single best hit
best = min([r for r in rows if r["p_value"] is not None], key=lambda r: r["p_value"])
print("\n  best-scoring comparison: %s x %s  d=%.3f  p=%.2e  (d_min80=%.2f)" % (
    best["driver"], best["signature"], best["observed_d"], best["p_value"], best["d_min_80"]))
print("  -> a reported effect of d=%.2f from this n, taken at face value," % abs(best["observed_d"]))
print("     is %.1fx above what an unselected test of this size can resolve." % (abs(best["observed_d"]) / best["d_min_80"]))

with open(os.path.join(OUT, "B4_all_tests.tsv"), "w", encoding="utf-8") as f:
    f.write("driver\tsignature\tn_altered\tn_wt\tprevalence\td_min_80\tobserved_d\tpower\tp_value\n")
    for r in rows:
        f.write("%s\t%s\t%d\t%d\t%.3f\t%.3f\t%.3f\t%.4f\t%.3e\n" % (
            r["driver"], r["signature"], r["n_altered"], r["n_wt"], r["prevalence"],
            r["d_min_80"], r["observed_d"], r["power_at_observed"],
            r["p_value"] if r["p_value"] is not None else float("nan")))
print("\n[saved] B4_all_tests.tsv")
