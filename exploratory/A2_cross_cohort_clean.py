"""A2-clean: per-cohort r(mRNA, protein) using the already-downloaded RPPA cache.

Barcode recovery: the file name on disk is the GDC file_id, so the barcode must
come from the GDC manifest, not from the path. We rebuild the manifest here.
"""
import json, urllib.request, os, gzip, math, statistics, glob
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")
FULL = os.environ.get("TCGA_FULL_DIR",
                      os.path.join(_ROOT, "data", "TCGA_full"))
CACHE = os.path.join(DM, "rppa_cache")


def gdc(endpoint, payload):
    req = urllib.request.Request("https://api.gdc.cancer.gov/" + endpoint,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    return json.load(urllib.request.urlopen(req, timeout=120))


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


def parse_cohort(code, ecols):
    """Rebuild {antibody: {barcode: value}} from the on-disk cache + GDC manifest."""
    ddir = os.path.join(CACHE, code)
    if not os.path.isdir(ddir):
        return None
    files = glob.glob(os.path.join(ddir, "*.tsv"))
    if not files:
        return None
    # manifest: file_id -> barcode
    flt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-" + code]}},
        {"op": "in", "content": {"field": "data_category", "value": ["Proteome Profiling"]}}]}
    hits = gdc("files", {"filters": flt, "size": 3000, "format": "JSON",
                         "fields": "file_id,file_name"}).get("data", {}).get("hits", [])
    fmap = {h["file_id"]: h["file_name"].split("_RPPA")[0][:12] for h in hits}
    eset = set(ecols)
    out = {}
    used = 0
    for p in files:
        fid = os.path.basename(p)[:-4]
        bc = fmap.get(fid)
        if bc is None or bc not in eset:
            continue
        used += 1
        with open(p, encoding="utf-8", errors="replace") as f:
            lines = [l for l in f if l.strip()]
        h = lines[0].rstrip("\n").split("\t")
        pi = h.index("peptide_target")
        ei = h.index("protein_expression")
        for l in lines[1:]:
            c = l.rstrip("\n").split("\t")
            if len(c) > max(pi, ei):
                try:
                    out.setdefault(c[pi], {})[bc] = float(c[ei])
                except Exception:
                    pass
    print("      parsed %d/%d cache files into %d antibodies" % (used, len(files), len(out)), flush=True)
    return out


print("=" * 92)
print("A2-clean. r(mRNA, protein) across cohorts, identical pipeline")
print("=" * 92)
print("%-8s %9s %8s %10s %10s %9s" % ("cohort", "n_samples", "n_antib", "mean_r", "median_r", "n_exact"))
print("-" * 92)
summary = []

# ---- MESO from its own file ----
meso_expr, meso_cols = None, None
p = os.path.join(DM, "MESO_expr_full.tsv.gz")
with gzip.open(p, "rt", encoding="utf-8", errors="replace") as g:
    h = g.readline().rstrip("\n").split("\t")
    meso_cols = [c.strip('"')[:12] for c in h[1:]]
    meso_expr = {}
    for line in g:
        c = line.rstrip("\n").split("\t")
        meso_expr[c[0].strip('"')] = [float(x) for x in c[1:]]
meso_rp = {}
with open(os.path.join(DM, "MESO_rppa.tsv"), encoding="utf-8") as f:
    hh = f.readline().rstrip("\n").split("\t")
    rs = hh[1:]
    for l in f:
        c = l.rstrip("\n").split("\t")
        d = {}
        for s, v in zip(rs, c[1:]):
            if v != "NA":
                d[s] = float(v)
        meso_rp[c[0]] = d

for code in ["MESO", "LUAD", "KIRC", "HNSC", "UCEC"]:
    if code == "MESO":
        expr, ecols, rp = meso_expr, meso_cols, meso_rp
    else:
        expr, ecols = load_expr(code)
        if expr is None:
            print("%-8s  (no expression matrix)" % code)
            continue
        rp = parse_cohort(code, ecols)
        if not rp:
            print("%-8s  (cache incomplete - still downloading)" % code)
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
        print("%-8s  (no usable antibodies)" % code)
        continue
    ns = len(set().union(*[set(v) for v in rp.values()])) if rp else 0
    print("%-8s %9d %8d %+10.3f %+10.3f %9d" % (
        code, ns, len(cors), statistics.mean(cors), statistics.median(cors), len(exact)))
    summary.append((code, statistics.mean(cors), len(cors)))

print("\n" + "-" * 92)
if len(summary) >= 2:
    meso = [s[1] for s in summary if s[0] == "MESO"]
    other = [s[1] for s in summary if s[0] != "MESO"]
    if meso and other:
        print("  MESO          r = %+.3f" % meso[0])
        print("  other cohorts r = %+.3f  (range %+.3f to %+.3f, %d cohorts)" % (
            statistics.mean(other), min(other), max(other), len(other)))
        d = meso[0] - statistics.mean(other)
        sd = statistics.pstdev(other) if len(other) > 1 else 0
        print("  difference    = %+.3f" % d)
        if sd > 0:
            print("  MESO is %.1f SD from the other-cohort mean" % (d / sd))
        print("\n  -> %s" % (
            "MESO is within the normal range: no cohort-specific decoupling"
            if abs(d) < max(0.06, 1.5 * sd) else
            "MESO differs from other cohorts - cohort-specific observation"))
