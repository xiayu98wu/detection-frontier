"""B6: figures-ready export for the methods paper.

Fig 3 is the pivotal figure: for every driver x phenotype test in the MESO case
study, plot |observed d| against d_min(n,p). The diagonal splits tests that COULD
have detected their own effect from those that could not.

Also computes the summary statistics quoted in the text.
"""
import os, math, statistics, collections
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

RES = os.path.join(_ROOT, "results")

# read the B4 table
rows = []
with open(os.path.join(RES, "B4_all_tests.tsv"), encoding="utf-8") as f:
    h = f.readline().rstrip("\n").split("\t")
    for l in f:
        c = l.rstrip("\n").split("\t")
        rows.append({
            "driver": c[0], "signature": c[1],
            "n_altered": int(c[2]), "n_wt": int(c[3]),
            "prevalence": float(c[4]), "d_min": float(c[5]),
            "d": float(c[6]), "power": float(c[7]), "p": float(c[8]),
        })
print("tests loaded: %d" % len(rows))

# ---------------- Fig 3 export ----------------
print("\n" + "=" * 88)
print("FIG 3 DATA: |d| vs d_min, per test")
print("=" * 88)
det = [r for r in rows if abs(r["d"]) >= r["d_min"]]
und = [r for r in rows if abs(r["d"]) < r["d_min"]]
print("  tests whose own effect EXCEEDS their detection floor (would be 'detectable'): %d" % len(det))
for r in sorted(det, key=lambda x: -abs(x["d"])):
    print("    %-14s x %-11s d=%+.3f  d_min=%.2f  p=%.1e" % (
        r["driver"], r["signature"], r["d"], r["d_min"], r["p"]))
print("  tests whose own effect is BELOW their detection floor: %d (%.0f%%)" % (
    len(und), 100 * len(und) / len(rows)))

# ---------------- which drivers are ever informative ----------------
print("\n" + "=" * 88)
print("INFORMATIVENESS BY DRIVER")
print("=" * 88)
by_d = collections.defaultdict(list)
for r in rows:
    by_d[r["driver"]].append(r)
print("%-14s %5s %9s %9s %11s" % ("driver", "n_tests", "max|d|", "d_min", "detectable"))
for d, rs in sorted(by_d.items(), key=lambda kv: -max(abs(r["d"]) for r in kv[1])):
    mx = max(abs(r["d"]) for r in rs)
    dm = rs[0]["d_min"]
    nd = sum(1 for r in rs if abs(r["d"]) >= r["d_min"])
    print("%-14s %5d %9.3f %9.2f %8d/%d" % (d, len(rs), mx, dm, nd, len(rs)))

# ---------------- the headline claim ----------------
print("\n" + "=" * 88)
print("HEADLINE NUMBERS FOR THE TEXT")
print("=" * 88)
cap = sum(1 for r in rows if r["d_min"] <= 0.40)
print("  tests capable of detecting a moderate effect (d_min <= 0.40): %d / %d" % (cap, len(rows)))
mx_all = max(abs(r["d"]) for r in rows)
print("  largest effect observed anywhere in the case study: d = %.3f" % mx_all)
print("  smallest detection floor anywhere in the case study: d_min = %.3f" % min(r["d_min"] for r in rows))
print("  ratio: the study's best signal is %.1fx its best detection floor" % (
    mx_all / min(r["d_min"] for r in rows)))

# how many tests would have been needed to detect their own effect?
print("\n  Of the %d tests that came out non-significant (p>=0.05):" % sum(1 for r in rows if r["p"] >= 0.05))
ns = [r for r in rows if r["p"] >= 0.05]
print("    %d had |d| below their own detection floor -> uninformative negative" % sum(
    1 for r in ns if abs(r["d"]) < r["d_min"]))
print("    %d had |d| above their floor but still p>=0.05 (small-sample noise)" % sum(
    1 for r in ns if abs(r["d"]) >= r["d_min"]))

# ---------------- export for plotting ----------------
out = os.path.join(RES, "B6_fig3_points.tsv")
with open(out, "w", encoding="utf-8") as f:
    f.write("driver\tsignature\tn_altered\tn_wt\tabs_d\td_min\tpower\tp_value\tdetectable\n")
    for r in sorted(rows, key=lambda x: abs(x["d"])):
        f.write("%s\t%s\t%d\t%d\t%.4f\t%.3f\t%.4f\t%.3e\t%s\n" % (
            r["driver"], r["signature"], r["n_altered"], r["n_wt"],
            abs(r["d"]), r["d_min"], r["power"], r["p"],
            "yes" if abs(r["d"]) >= r["d_min"] else "no"))
print("\n[saved] B6_fig3_points.tsv")

# ---------------- Fig 2 export: cohort frontier ----------------
print("\n" + "=" * 88)
print("FIG 2 DATA: cohort reachability")
print("=" * 88)
print("  (see B5_cohort_frontier.tsv; reproduced here with effect classes)")
print("  %-8s %6s %10s %10s %10s" % ("cohort", "n", "d_min", "can_d0.5", "can_d0.3"))
import glob
cohorts = []
with open(os.path.join(RES, "B2_cohort_instantiation.tsv"), encoding="utf-8") as f:
    f.readline()
    for l in f:
        c = l.rstrip("\n").split("\t")
        cohorts.append((c[0], int(c[1]), float(c[2]), float(c[3]), float(c[4])))
for code, n, d30, d20, d10 in sorted(cohorts, key=lambda x: x[1]):
    print("  %-8s %6d %10.2f %10s %10s" % (
        code, n, d30,
        "yes" if d30 <= 0.50 else "no",
        "yes" if d30 <= 0.30 else "no"))
