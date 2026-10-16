"""A2-minimal: the decisive cheap test.

Is TCGA-MESO's mRNA-protein correlation (0.31) different from other cohorts
computed with the identical pipeline?

Only needs RPPA + expression for the comparison cohorts (no CNV download).

If MESO ~= other cohorts -> 0.31 is normal, no MESO-specific "decoupling" story.
If MESO << other cohorts  -> a real cohort-specific observation worth pursuing.
"""
import json, urllib.request, os, io, gzip, math, statistics, sys
from concurrent.futures import ThreadPoolExecutor
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")
FULL = os.environ.get("TCGA_FULL_DIR",
                      os.path.join(_ROOT, "data", "TCGA_full"))
CACHE = os.path.join(DM, "rppa_cache")
os.makedirs(CACHE, exist_ok=True)

COHORTS = sys.argv[1:] if len(sys.argv) > 1 else ["LUAD", "KIRC", "HNSC", "UCEC"]


def gdc(endpoint, payload):
    req = urllib.request.Request("https://api.gdc.cancer.gov/" + endpoint,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    return json.load(urllib.request.urlopen(req, timeout=120))


def fetch(fid, dest, tries=4):
    if os.path.exists(dest) and os.path.getsize(dest) > 3000:
        return True
    for t in range(tries):
        try:
            b = urllib.request.urlopen("https://api.gdc.cancer.gov/data/" + fid, timeout=240).read()
            if len(b) < 3000:
                raise ValueError("small")
            open(dest, "wb").write(b)
            return True
        except Exception:
            pass
    return False


def load_expr(code):
    for nm in [code + "_full.tsv.gz", code + "_HiSeqV2.gz"]:
        p = os.path.join(FULL, nm)
        if os.path.exists(p):
            with gzip.open(p, "rt", encoding="utf-8", errors="replace") as g:
                h = g.readline().rstrip("\n").split("\t")
                cols = [c.strip('"')[:12] for c in h[1:]]
                d = {}
                for line in g:
                    c = line.rstrip("\n").split("\t")
                    d[c[0].strip('"')] = [float(x) for x in c[1:]]
            return d, cols
    return None, None


def pearson(x, y):
    n = len(x)
    if n < 20:
        return None
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    dx = math.sqrt(sum((x[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((y[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx > 0 and dy > 0 else None


def cohort_rppa(code, ecols):
    """Return {antibody: [values aligned to ecols]}, downloading what is needed."""
    cache = os.path.join(CACHE, code + "_rppa.json")
    if os.path.exists(cache):
        with open(cache) as f:
            d = json.load(f)
        return d
    flt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-" + code]}},
        {"op": "in", "content": {"field": "data_category", "value": ["Proteome Profiling"]}}]}
    hits = gdc("files", {"filters": flt, "size": 3000, "format": "JSON",
                         "fields": "file_id,file_name"}).get("data", {}).get("hits", [])
    want = set(ecols)
    todo = [(h["file_id"], h["file_name"].split("_RPPA")[0][:12]) for h in hits
            if h["file_name"].split("_RPPA")[0][:12] in want]
    print("  %-5s %4d RPPA files, %4d match expression" % (code, len(hits), len(todo)), flush=True)
    if not todo:
        return {}
    ddir = os.path.join(CACHE, code)
    os.makedirs(ddir, exist_ok=True)
    with ThreadPoolExecutor(max_workers=8) as ex:
        list(ex.map(lambda t: fetch(t[0], os.path.join(ddir, t[0] + ".tsv")), todo))

    # per-barcode value maps
    per = {ab: {} for ab in []}
    rows = {}
    for fid, bc in todo:
        p = os.path.join(ddir, fid + ".tsv")
        if not os.path.exists(p):
            continue
        try:
            with open(p, encoding="utf-8", errors="replace") as f:
                lines = [l for l in f if l.strip()]
            h = lines[0].split("\t")
            pi = h.index("peptide_target")
            ei = h.index("protein_expression")
            for l in lines[1:]:
                c = l.split("\t")
                if len(c) > max(pi, ei):
                    try:
                        rows.setdefault(c[pi], {})[bc] = float(c[ei])
                    except Exception:
                        pass
        except Exception:
            continue
    with open(cache, "w") as f:
        json.dump(rows, f)
    return rows


print("=" * 92)
print("A2-minimal. r(mRNA, protein), identical pipeline, per cohort")
print("=" * 92)
print("%-8s %8s %9s %10s %10s %9s" % ("cohort", "n_samples", "n_antib", "mean_r", "median_r", "n_exact"))
print("-" * 92)
summary = []
for code in ["MESO"] + COHORTS:
    if code == "MESO":
        expr, ecols = load_expr("MESO")
        if expr is None:
            d = os.path.join(DM, "MESO_expr_full.tsv.gz")
            with gzip.open(d, "rt", encoding="utf-8", errors="replace") as g:
                h = g.readline().rstrip("\n").split("\t")
                ecols = [c.strip('"')[:12] for c in h[1:]]
                expr = {}
                for line in g:
                    c = line.rstrip("\n").split("\t")
                    expr[c[0].strip('"')] = [float(x) for x in c[1:]]
        rp = {}
        with open(os.path.join(DM, "MESO_rppa.tsv"), encoding="utf-8") as f:
            hh = f.readline().rstrip("\n").split("\t")
            rs = hh[1:]
            for l in f:
                c = l.rstrip("\n").split("\t")
                d = {}
                for s, v in zip(rs, c[1:]):
                    if v != "NA":
                        d[s] = float(v)
                rp[c[0]] = d
    else:
        expr, ecols = load_expr(code)
        if expr is None:
            print("%-8s (no local expression matrix)" % code)
            continue
        rp = cohort_rppa(code, ecols)
        if not rp:
            print("%-8s (no RPPA)" % code)
            continue

    eidx = {b: i for i, b in enumerate(ecols)}
    exact = [ab for ab in rp if ab in expr]
    cors = []
    for ab in exact:
        vals = rp[ab]
        pr = [(vals[s], expr[ab][eidx[s]]) for s in vals if s in eidx]
        if len(pr) < 20:
            continue
        r = pearson([a for a, _ in pr], [b for _, b in pr])
        if r is not None:
            cors.append(r)
    if not cors:
        print("%-8s (no usable antibodies)" % code)
        continue
    ns = len(set().union(*[set(v) for v in list(rp.values())[:50]])) if rp else 0
    print("%-8s %8d %9d %+10.3f %+10.3f %9d" % (
        code, ns, len(cors), statistics.mean(cors), statistics.median(cors), len(exact)))
    summary.append((code, len(cors), statistics.mean(cors)))

print("\n" + "-" * 92)
if len(summary) > 1:
    others = [s[2] for s in summary if s[0] != "MESO"]
    meso = [s[2] for s in summary if s[0] == "MESO"]
    if others and meso:
        print("  MESO            mean r = %+.3f" % meso[0])
        print("  other cohorts   mean r = %+.3f  (sd %.3f across %d cohorts)" % (
            statistics.mean(others), statistics.pstdev(others) if len(others) > 1 else 0, len(others)))
        d = meso[0] - statistics.mean(others)
        print("  difference      = %+.3f" % d)
        print("\n  -> %s" % (
            "MESO mRNA-protein coupling is NOT unusual; 0.31 is the norm for single cohorts"
            if abs(d) < 0.06 else
            "MESO differs from other cohorts - potentially MESO-specific"))
