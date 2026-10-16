"""Track B module B5: the applicability frontier.

Converts the power surface into a decision rule a researcher can use BEFORE
starting a study: given a cohort of size n and a driver of prevalence p, which
question classes are answerable?

Question classes, by the effect size the literature typically reports:
  large   d ~ 0.8-1.0  (e.g. a deletion that shifts a proliferation programme)
  medium  d ~ 0.5      (e.g. a driver that shifts an immune signature)
  small   d ~ 0.3      (e.g. most single-driver / single-signature associations)

Also reports the "required n" — how many samples a study would need to answer
each class at 80% power.
"""
import os, math, statistics, json
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

OUT = os.path.join(_ROOT, "results")
os.makedirs(OUT, exist_ok=True)
PHI = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))


def power(d, n1, n2, alpha=0.05):
    if n1 < 2 or n2 < 2:
        return 0.0
    ncp = d * math.sqrt(n1 * n2 / (n1 + n2))
    return (1 - PHI(1.959963985 - ncp)) + PHI(-1.959963985 - ncp)


def d_min(n, p, target=0.80):
    n1 = max(2, int(round(n * p)))
    n2 = max(2, n - n1)
    lo, hi = 0.01, 10.0
    if power(hi, n1, n2) < target:
        return float("nan")
    for _ in range(150):
        mid = (lo + hi) / 2
        if power(mid, n1, n2) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def n_required(d, p, target=0.80):
    """Smallest total n achieving target power for effect d at prevalence p."""
    for n in range(20, 20001, 2):
        n1 = max(2, int(round(n * p)))
        n2 = max(2, n - n1)
        if power(d, n1, n2) >= target:
            return n
    return None


CLASSES = [("large", 0.90), ("medium", 0.50), ("small", 0.30)]
PREV = [0.50, 0.30, 0.20, 0.10]

print("=" * 96)
print("B5. APPLICABILITY FRONTIER: minimum cohort n needed to detect each effect class")
print("     (80% power, alpha=0.05, two-group comparison)")
print("=" * 96)
print("%-10s %10s %12s %12s %12s %12s" % ("effect class", "d", *["p=%.2f" % p for p in PREV]))
print("-" * 96)
req_rows = []
for label, d in CLASSES:
    vals = [n_required(d, p) for p in PREV]
    print("%-10s %10.2f %12s %12s %12s %12s" % (
        label, d, *[("%d" % v if v else ">20000") for v in vals]))
    req_rows.append((label, d) + tuple(vals))

print("\n" + "=" * 96)
print("B5b. WHICH ONLY-EXISTING TCGA COHORTS CLEAR EACH BAR?")
print("=" * 96)
# Cohort sizes are read from the B2 output rather than repeated here, so that
# every number in the release traces back to one place. Run B1_power_surface.py
# (which writes B2_cohort_instantiation.tsv) first; the shipped release/results/
# already contains it.
_B2 = os.path.join(OUT, "B2_cohort_instantiation.tsv")
if not os.path.exists(_B2):
    raise SystemExit("B2_cohort_instantiation.tsv not found in %s -- "
                     "run B1_power_surface.py first." % OUT)
COHORTS = []
with open(_B2, encoding="utf-8-sig") as _f:
    next(_f)
    for _line in _f:
        _c = _line.rstrip("\n").split("\t")
        if len(_c) >= 2 and _c[0]:
            COHORTS.append((_c[0], int(_c[1])))
COHORTS.sort(key=lambda t: t[1])
print("%-8s %7s %14s %14s %14s" % ("cohort", "n", "large(d=0.9)", "medium(d=0.5)", "small(d=0.3)"))
print("-" * 62)
for code, n in COHORTS:
    marks = []
    for label, d in CLASSES:
        dm = d_min(n, 0.30)   # at a representative 30% prevalence
        marks.append("YES" if dm <= d else "no")
    print("%-8s %7d %14s %14s %14s" % (code, n, marks[0], marks[1], marks[2]))

print("\n  reading: entries are whether that cohort can detect the effect class")
print("  at p=0.30 with 80% power.")

print("\n" + "=" * 96)
print("B5c. THE MESO CASE, CONCRETELY")
print("=" * 96)
meso_n = 87
print("  cohort n = 87")
for label, d in CLASSES:
    n1 = int(round(meso_n * 0.30))
    n2 = meso_n - n1
    pw = power(d, n1, n2)
    dm = d_min(meso_n, 0.30)
    print("    %-7s d=%.2f  power=%.1f%%   %s" % (
        label, d, 100 * pw, "detectable" if dm <= d else "NOT detectable"))
print("\n  observed in the MESO case study:")
OBS = [("CDKN2A del x PROLIF", 0.898, 36, 51), ("CDKN2A del x EMT", 0.839, 36, 51),
       ("NF2 mut x IFN-g", 0.341, 19, 59), ("BAP1 mut x CD8_T", 0.011, 22, 56)]
for label, d, n1, n2 in OBS:
    print("    %-22s d=%.3f  power=%.1f%%" % (label, d, 100 * power(d, n1, n2)))

# ---- save ----
with open(os.path.join(OUT, "B5_required_n.tsv"), "w", encoding="utf-8") as f:
    f.write("effect_class\td\t" + "\t".join("n_required_p%.2f" % p for p in PREV) + "\n")
    for row in req_rows:
        f.write("%s\t%.2f\t" % (row[0], row[1]) +
                "\t".join(str(v) if v else ">20000" for v in row[2:]) + "\n")
with open(os.path.join(OUT, "B5_cohort_frontier.tsv"), "w", encoding="utf-8") as f:
    f.write("cohort\tn\tlarge_d0.9\tmedium_d0.5\tsmall_d0.3\n")
    for code, n in COHORTS:
        marks = []
        for label, d in CLASSES:
            marks.append("Y" if d_min(n, 0.30) <= d else "N")
        f.write("%s\t%d\t%s\n" % (code, n, "\t".join(marks)))
print("\n[saved] B5_required_n.tsv / B5_cohort_frontier.tsv")
