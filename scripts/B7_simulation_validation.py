"""B7: SIMULATION VALIDATION of the normal-approximation power surface.

Fully vectorised with numpy (rank computation, tie correction, welch t) so the
whole grid runs in seconds instead of hours.

Why this is needed: BMC Med Res Methodol reviewers will ask whether the closed-form
d_min(n,p) matches empirical power under conditions that are NOT ideal:
  (a) scores are gene-signature means, not normal
  (b) group sizes are unbalanced and set by a driver's prevalence
  (c) d standardised by the pooled SD over the WHOLE cohort (manuscript definition)
  (d) the study's tests are Mann-Whitney (rank-based), with tie correction

Protocol: generate from a raw shift delta -> compute the REALIZED d under the
manuscript definition -> compare empirical power at that realized d against the
closed form. Bisect on the realized d to get the empirical 80% threshold.

REPRODUCIBILITY NOTE
--------------------
Seeds are derived from zlib.crc32 of the scenario label, NOT from Python's
built-in hash(). hash() is salted per process for strings (PEP 456), so a
seed built from hash() differs on every run and the calibration would not be
reproducible. This was a real defect in an earlier version of this script.

GROUP-SIZE CONFIGURATIONS
-------------------------
Three of the four configurations are the ACTUAL group sizes of driver strata in
the TCGA-MESO case study, taken from results/B4_all_tests.tsv:
    CDKN2A_del   36 : 51   (p = 0.414)  <- the association that WAS resolvable
    CDKN2A_loss  48 : 39   (p = 0.552)  <- any-copy loss, high prevalence
    BAP1_mut     22 : 56   (p = 0.282)  <- the association that was NOT resolvable
The fourth is a balanced 180:180 design. It has NO real-cohort analogue and is
included only as a large-n scaling reference; it is labelled as synthetic
everywhere it appears.

An earlier version of this script used the pairs 29:23 and 36:33 under the
labels n=87/120/300. Those pairs had no derivation and did not sum to the
stated n. They have been removed.
"""
import os, math, sys, zlib
import numpy as np
from scipy.stats import rankdata
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

OUT = os.path.join(_ROOT, "results")
os.makedirs(OUT, exist_ok=True)
PHI = lambda x: 0.5 * (1 + math.erf(x / math.sqrt(2)))
ZC = 1.959963985


def stable_seed(*parts):
    """Deterministic across processes (unlike built-in hash())."""
    return zlib.crc32("|".join(str(p) for p in parts).encode("utf-8")) % (2 ** 31)


def formula_power(d, n1, n2):
    if n1 < 2 or n2 < 2:
        return 0.0
    ncp = d * math.sqrt(n1 * n2 / (n1 + n2))
    return (1 - PHI(ZC - ncp)) + PHI(-ZC - ncp)


def formula_dmin(n1, n2, target=0.80):
    lo, hi = 0.01, 8.0
    for _ in range(120):
        mid = (lo + hi) / 2
        if formula_power(mid, n1, n2) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def mw_pvalue_vec(A, B):
    """Vectorised two-sided Mann-Whitney with tie correction.
    A: (nrep, n1), B: (nrep, n2) -> returns p-values (nrep,)."""
    nrep, n1 = A.shape
    n2 = B.shape[1]
    comb = np.concatenate([A, B], axis=1)
    ranks = rankdata(comb, axis=1)          # average ranks, ties handled
    U1 = ranks[:, :n1].sum(axis=1) - n1 * (n1 + 1) / 2.0

    # tie correction: sum(t^3 - t) per row
    s = np.sort(comb, axis=1)
    n = s.shape[1]
    newrun = np.ones_like(s, dtype=bool)
    newrun[:, 1:] = s[:, 1:] != s[:, :-1]
    runidx = np.cumsum(newrun, axis=1) - 1
    maxrun = int(runidx.max()) + 1
    counts = np.zeros((nrep, maxrun), dtype=np.int64)
    rows = np.repeat(np.arange(nrep), n)
    np.add.at(counts, (rows, runidx.ravel()), 1)
    tie = (counts.astype(np.float64) ** 3 - counts).sum(axis=1)

    N = n1 + n2
    den = N * (N - 1)
    var = n1 * n2 / 12.0 * ((N + 1) - tie / den)
    var = np.where(var <= 0, np.nan, var)
    z = (U1 - n1 * n2 / 2.0) / np.sqrt(var)
    p = 2 * (1 - 0.5 * (1 + np.vectorize(math.erf)(np.abs(z) / math.sqrt(2))))
    return p, z


def welch_pvalue_vec(A, B):
    n1, n2 = A.shape[1], B.shape[1]
    v1 = A.var(axis=1, ddof=1)
    v2 = B.var(axis=1, ddof=1)
    se = np.sqrt(v1 / n1 + v2 / n2)
    se = np.where(se <= 0, np.nan, se)
    t = (A.mean(axis=1) - B.mean(axis=1)) / se
    return 2 * (1 - 0.5 * (1 + np.vectorize(math.erf)(np.abs(t) / math.sqrt(2))))


def realized_d_vec(A, B):
    """Manuscript definition: mean difference / SD pooled over the whole cohort."""
    allv = np.concatenate([A, B], axis=1)
    sd = allv.std(axis=1, ddof=0)
    sd = np.where(sd <= 0, np.nan, sd)
    return (A.mean(axis=1) - B.mean(axis=1)) / sd


# ---------------- generators: (nrep, n) arrays ----------------
def gen(kind, n1, n2, delta, nrep, rng):
    if kind == "ideal normal":
        A = rng.normal(delta, 1.0, (nrep, n1))
        B = rng.normal(0.0, 1.0, (nrep, n2))
    elif kind == "gene-set score":
        # A gene-set score is a mean over G marker genes. Averaging 40 normals is
        # still normal, so to test the thing that actually matters (a SKEWED score
        # distribution) the per-gene noise is log-normal and the resulting score is
        # rescaled to SD = 1, which makes delta equal d directly.
        G = 40
        M, S = 1.6487212707, 0.3417365957      # mean and SD of mean-of-40 lognormal(0,1)
        A = (rng.lognormal(0.0, 1.0, (nrep, n1, G)).mean(axis=2) - M) / S
        B = (rng.lognormal(0.0, 1.0, (nrep, n2, G)).mean(axis=2) - M) / S
        A = A + delta                          # effect carried by group A (positive d)
    elif kind == "dropout-contaminated":
        # Equal within-group variance; 6% of samples per group carry a left-shifted
        # outlier, which creates the heavy left tail and tied ranks seen in assay
        # data. NOT heteroscedastic: both groups share the same variance by design.
        A = rng.normal(delta, 0.9, (nrep, n1))
        B = rng.normal(0.0, 0.9, (nrep, n2))
        mA = rng.random((nrep, n1)) < 0.06
        mB = rng.random((nrep, n2)) < 0.06
        A = np.where(mA, A - rng.exponential(2.0, (nrep, n1)), A)
        B = np.where(mB, B - rng.exponential(2.0, (nrep, n2)), B)
    else:
        raise ValueError(kind)
    return A, B


def calibrate(kind, n1, n2, delta, nrep, seed):
    rng = np.random.default_rng(seed)
    A, B = gen(kind, n1, n2, delta, nrep, rng)
    p_mw, _ = mw_pvalue_vec(A, B)
    p_w = welch_pvalue_vec(A, B)
    dbar = float(np.nanmean(realized_d_vec(A, B)))
    return {
        "d_realized": dbar,
        "power_mw": float(np.nanmean(p_mw < 0.05)),
        "power_welch": float(np.nanmean(p_w < 0.05)),
        "power_formula": formula_power(dbar, n1, n2),
    }


def empirical_dmin(kind, n1, n2, target=0.80, tol=0.03, nrep=1200, nrep_final=6000):
    """Empirical 80%-power threshold, returned IN d UNITS.

    The bisection runs over the generating shift delta (monotone in d), but the
    closed-form d_min is expressed in units of the whole-cohort standardised
    effect. Returning the bisected delta and subtracting it from the formula
    would compare two different scales; at delta ~0.6 they differ by ~6%.
    The final threshold is therefore converted by evaluating the realized d at
    the converged delta with a larger replicate count.
    """
    lo, hi = 0.0, 4.0
    while hi - lo > tol:
        mid = (lo + hi) / 2
        r = calibrate(kind, n1, n2, mid, nrep,
                      stable_seed("dmin", kind, n1, n2, round(mid, 4)))
        if r["power_mw"] < target:
            lo = mid
        else:
            hi = mid
    delta_star = (lo + hi) / 2
    r = calibrate(kind, n1, n2, delta_star, nrep_final,
                  stable_seed("dminfinal", kind, n1, n2, round(delta_star, 4)))
    return r["d_realized"], delta_star


KINDS = ["ideal normal", "gene-set score", "dropout-contaminated"]

# (label, n1, n2, provenance) -- see module docstring
CONFIGS = [
    ("CDKN2A_del",  36,  51, "TCGA-MESO real stratum (p=0.414)"),
    ("CDKN2A_loss", 48,  39, "TCGA-MESO real stratum (p=0.552)"),
    ("BAP1_mut",    22,  56, "TCGA-MESO real stratum (p=0.282)"),
    ("balanced_ref", 180, 180, "SYNTHETIC balanced reference, no cohort analogue"),
]
DELTAS = [0.0, 0.3, 0.6, 0.9, 1.2]

print("=" * 112, flush=True)
print("B7. SIMULATION VALIDATION OF THE POWER SURFACE", flush=True)
print("=" * 112, flush=True)
print("group-size configurations:", flush=True)
for lab, a, b, prov in CONFIGS:
    print("  %-14s %4d : %-4d  n=%-4d  %s" % (lab, a, b, a + b, prov), flush=True)

rows = []
for kind in KINDS:
    print("\n" + "-" * 112, flush=True)
    print("scenario: %s" % kind, flush=True)
    print("-" * 112, flush=True)
    print("%-14s %5s %5s %7s %8s %9s %9s %9s %9s" % (
        "config", "n1", "n2", "delta", "d_real", "pow(MW)", "pow(t)", "pow(form)", "err"), flush=True)
    for (lab, n1, n2, _prov) in CONFIGS:
        for delta in DELTAS:
            seed = stable_seed(kind, lab, n1, n2, delta)
            r = calibrate(kind, n1, n2, delta, 2000, seed)
            err = r["power_mw"] - r["power_formula"]
            print("%-14s %5d %5d %7.2f %8.3f %9.3f %9.3f %9.3f %+9.3f" % (
                lab, n1, n2, delta, r["d_realized"], r["power_mw"],
                r["power_welch"], r["power_formula"], err), flush=True)
            rows.append(dict(kind=kind, config=lab, n1=n1, n2=n2, n=n1 + n2,
                             delta=delta, **r))

errs = np.array([r["power_mw"] - r["power_formula"] for r in rows])
print("\n" + "=" * 112, flush=True)
print("CALIBRATION SUMMARY: empirical MW power minus closed-form power", flush=True)
print("=" * 112, flush=True)
print("  comparisons            : %d" % len(errs), flush=True)
print("  mean error             : %+.4f" % errs.mean(), flush=True)
print("  median error           : %+.4f" % float(np.median(errs)), flush=True)
print("  max |error|            : %.4f" % np.abs(errs).max(), flush=True)
print("  95th pct |error|       : %.4f" % np.percentile(np.abs(errs), 95), flush=True)
print("  errors within +/-0.05  : %d/%d" % (int((np.abs(errs) <= 0.05).sum()), len(errs)), flush=True)
print("  errors within +/-0.10  : %d/%d" % (int((np.abs(errs) <= 0.10).sum()), len(errs)), flush=True)

print("\n" + "=" * 112, flush=True)
print("d_min CALIBRATION (empirical 80%%-power threshold vs closed form, BOTH in d units)", flush=True)
print("=" * 112, flush=True)
print("%-14s %-6s %5s %5s %5s %11s %13s %13s %9s" % (
    "config", "mech", "n1", "n2", "n", "delta*", "d_min(emp)", "d_min(form)", "error"), flush=True)
dmin_rows = []
for kind in KINDS:
    for (lab, n1, n2, _prov) in CONFIGS:
        emp_d, delta_star = empirical_dmin(kind, n1, n2)
        fm = formula_dmin(n1, n2)
        print("%-14s %-6s %5d %5d %5d %11.3f %13.3f %13.3f %+9.3f" % (
            lab, kind[:6], n1, n2, n1 + n2, delta_star, emp_d, fm, emp_d - fm), flush=True)
        dmin_rows.append((kind, lab, n1, n2, delta_star, emp_d, fm))

with open(os.path.join(OUT, "B7_simulation.tsv"), "w", encoding="utf-8") as f:
    f.write("scenario\tconfig\tn1\tn2\tn_total\tdelta\td_realized\tpower_mw\tpower_welch\tpower_formula\terror\n")
    for r in rows:
        f.write("%s\t%s\t%d\t%d\t%d\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%+.4f\n" % (
            r["kind"], r["config"], r["n1"], r["n2"], r["n"], r["delta"], r["d_realized"],
            r["power_mw"], r["power_welch"], r["power_formula"],
            r["power_mw"] - r["power_formula"]))
with open(os.path.join(OUT, "B7_dmin_calibration.tsv"), "w", encoding="utf-8") as f:
    f.write("scenario\tconfig\tn1\tn2\tn_total\tdelta_star\td_min_empirical\td_min_formula\terror\n")
    for k, lab, n1, n2, ds, emp, fm in dmin_rows:
        f.write("%s\t%s\t%d\t%d\t%d\t%.3f\t%.3f\t%.3f\t%+.3f\n" % (
            k, lab, n1, n2, n1 + n2, ds, emp, fm, emp - fm))
errs_d = np.array([emp - fm for (_k, _l, _a, _b, _ds, emp, fm) in dmin_rows])
print("\n  d_min error: max |err| = %.3f, mean %+.4f, within +/-0.03: %d/%d, within +/-0.05: %d/%d" % (
    np.abs(errs_d).max(), errs_d.mean(),
    int((np.abs(errs_d) <= 0.03).sum()), len(errs_d),
    int((np.abs(errs_d) <= 0.05).sum()), len(errs_d)), flush=True)
print("\n[saved] B7_simulation.tsv / B7_dmin_calibration.tsv", flush=True)
