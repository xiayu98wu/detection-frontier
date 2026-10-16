"""Track B modules B1+B2.

B1: power surface. For a two-group comparison (driver-altered vs wild-type) with
    total n and alteration prevalence p, what standardized effect size d is
    detectable at 80% power / alpha=0.05?

B2: instantiate the surface on the rare-tumour TCGA cohorts present locally.

This is the quantitative backbone of the methods paper: it converts "n is small"
into a specific, checkable number.
"""
import gzip, os, math, glob, statistics
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")
FULL = os.environ.get("TCGA_FULL_DIR",
                      os.path.join(_ROOT, "data", "TCGA_full"))
OUT = os.path.join(_ROOT, "results")
os.makedirs(OUT, exist_ok=True)

PHI = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))


def power_two_group(d, n1, n2, alpha=0.05):
    """Normal-approximation power of a two-sample test (equal variance)."""
    if n1 < 2 or n2 < 2:
        return 0.0
    zc = 1.959963985 if alpha == 0.05 else 1.644853627
    ncp = d * math.sqrt(n1 * n2 / (n1 + n2))
    return (1 - PHI(zc - ncp)) + PHI(-zc - ncp)


def d_min(n, p, power=0.80, alpha=0.05):
    """Smallest detectable d for total n with a group of proportion p."""
    n1 = max(2, int(round(n * p)))
    n2 = max(2, n - n1)
    lo, hi = 0.01, 8.0
    if power_two_group(hi, n1, n2, alpha) < power:
        return float("nan")
    for _ in range(200):
        mid = (lo + hi) / 2
        if power_two_group(mid, n1, n2, alpha) < power:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ============ B1: the surface ============
NGRID = [40, 50, 60, 70, 80, 90, 100, 120, 150, 200, 300, 400, 600, 800, 1200]
PGRID = [0.50, 0.40, 0.30, 0.25, 0.20, 0.15, 0.10]

print("=" * 96)
print("B1. MINIMUM DETECTABLE EFFECT SIZE d  (80% power, alpha=0.05, two-group)")
print("=" * 96)
hdr = "%-8s" % "n\\p" + "".join("%10s" % ("p=%.2f" % p) for p in PGRID)
print(hdr)
print("-" * len(hdr))
lines = ["n\t" + "\t".join("p=%.2f" % p for p in PGRID)]
for n in NGRID:
    row = "%-8d" % n + "".join("%10.2f" % d_min(n, p) for p in PGRID)
    print(row)
    lines.append("%d\t" % n + "\t".join("%.3f" % d_min(n, p) for p in PGRID))

# ============ B2: instantiate on rare TCGA cohorts ============
print("\n" + "=" * 96)
print("B2. INSTANTIATED ON TCGA COHORTS (sample counts from local full matrices)")
print("=" * 96)

RARE = ["ACC", "CHOL", "DLBC", "MESO", "UCS", "UVM", "PCPG",
        "TGCT", "THYM", "KICH", "CESC", "READ", "ESCA", "BLCA", "KIRP"]
rows = []
for code in sorted(os.listdir(FULL)) if os.path.isdir(FULL) else []:
    pass

found = {}
for code in RARE:
    f = os.path.join(FULL, code + "_full.tsv.gz")
    if not os.path.exists(f):
        f2 = os.path.join(FULL, code + "_HiSeqV2.gz")
        f = f2 if os.path.exists(f2) else None
    if not f:
        continue
    with gzip.open(f, "rt", encoding="utf-8", errors="replace") as g:
        h = g.readline().rstrip("\n").split("\t")
        cols = [c.strip('"') for c in h[1:]]
    tumour = [c for c in cols if len(c) >= 15 and c[13:15] == "01"]
    found[code] = len(tumour) if tumour else len(cols)

print("%-8s %8s %12s %12s %12s" % ("cohort", "n_tumour", "d_min(p=.30)", "d_min(p=.20)", "d_min(p=.10)"))
print("-" * 58)
inst = ["cohort\tn_tumour\tdmin_p30\tdmin_p20\tdmin_p10"]
for code, n in sorted(found.items(), key=lambda x: x[1]):
    d30, d20, d10 = d_min(n, 0.30), d_min(n, 0.20), d_min(n, 0.10)
    print("%-8s %8d %12.2f %12.2f %12.2f" % (code, n, d30, d20, d10))
    inst.append("%s\t%d\t%.3f\t%.3f\t%.3f" % (code, n, d30, d20, d10))

if not found:
    # The pan-cancer matrices live outside this repository. Without them B2 is
    # empty, and writing an empty file would silently destroy the cohort table
    # shipped in results/ -- which B5 depends on. Refuse to overwrite instead.
    print("\n  NOTE: no expression matrices found under\n    %s\n"
          "  Set TCGA_FULL_DIR to the directory holding <COHORT>_full.tsv.gz\n"
          "  (or <COHORT>_HiSeqV2.gz) to regenerate B2. The table shipped in\n"
          "  results/B2_cohort_instantiation.tsv is left untouched." % FULL)
    inst = None

# ============ where do real observations fall? ============
print("\n" + "=" * 96)
print("B3-prep. OBSERVED effect sizes from the MESO case study, vs what n=87 can detect")
print("=" * 96)
N_MESO = found.get("MESO", 87)
OBSERVED = [
    ("BAP1 mut vs CD8_T",            0.011, 22, 56),
    ("BAP1 mut vs IFNg",             0.037, 22, 56),
    ("BAP1 mut vs checkpoint",       0.084, 22, 56),
    ("NF2 mut vs CD8_T",             0.124, 19, 59),
    ("NF2 mut vs IFNg",              0.341, 19, 59),
    ("NF2 mut vs macrophage",        0.383, 19, 59),
    ("CDKN2A del vs EMT",            0.882, 36, 51),
    ("CDKN2A del vs PROLIF",         0.995, 36, 51),
    ("CDKN2A del vs FIBROBLAST",     0.804, 36, 51),
]
print("%-26s %8s %8s %10s %10s" % ("comparison", "n1", "n2", "observed d", "power @n"))
print("-" * 68)
obs_rows = ["comparison\tn_altered\tn_wt\tobserved_d\tpower"]
for label, d, n1, n2 in OBSERVED:
    pw = power_two_group(d, n1, n2)
    print("%-26s %8d %8d %10.3f %9.1f%%" % (label, n1, n2, d, 100 * pw))
    obs_rows.append("%s\t%d\t%d\t%.3f\t%.4f" % (label, n1, n2, d, pw))

print("\n  MESO cohort n=%d" % N_MESO)
print("  d_min at p=0.25 (BAP1-like prevalence 22/78): %.2f" % d_min(78, 0.25))
print("  d_min at p=0.46 (CDKN2A-like prevalence 36/87): %.2f" % d_min(87, 0.46))

with open(os.path.join(OUT, "B1_power_surface.tsv"), "w", encoding="utf-8") as f:
    f.write("\n".join(lines) + "\n")
if inst is not None:
    with open(os.path.join(OUT, "B2_cohort_instantiation.tsv"), "w", encoding="utf-8") as f:
        f.write("\n".join(inst) + "\n")
    print("\n[saved] B1_power_surface.tsv / B2_cohort_instantiation.tsv / B3_observed_effects.tsv")
else:
    print("\n[saved] B1_power_surface.tsv / B3_observed_effects.tsv"
          "  (B2 not overwritten)")
