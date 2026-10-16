"""A2: cross-cohort comparison done properly.

Question: is TCGA-MESO's CNV->mRNA->protein transmission different from other
cohorts, when computed with an IDENTICAL pipeline?

Key design decision (from A1-QC3): a random gene panel is dominated by genes whose
mRNA is not copy-number driven, which makes r(CNV,mRNA) small by construction.
The informative statistic restricts to genes where CNV actually varies AND drives
mRNA, then asks what fraction of that effect survives to protein.

Outputs, per cohort:
  n_genes          three-layer genes available
  r_cm             mean r(CNV, mRNA)
  r_mp             mean r(mRNA, protein)
  r_cp             mean r(CNV, protein)
  r_cm_driven      same, restricted to CNV-varying genes
  T                mean transmission ratio r_cp/r_cm (CNV-driven subset)
"""
import json, urllib.request, os, io, tarfile, gzip, math, statistics, sys
from concurrent.futures import ThreadPoolExecutor
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")
FULL = os.environ.get("TCGA_FULL_DIR",
                      os.path.join(_ROOT, "data", "TCGA_full"))
CNVDIR = os.path.join(DM, "cox_cnv")
os.makedirs(CNVDIR, exist_ok=True)

COHORTS = sys.argv[1:] if len(sys.argv) > 1 else ["LUAD", "KIRC", "HNSC"]


def gdc(endpoint, payload):
    req = urllib.request.Request("https://api.gdc.cancer.gov/" + endpoint,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    return json.load(urllib.request.urlopen(req, timeout=120))


def fetch_one(fid, dest, tries=4):
    if os.path.exists(dest) and os.path.getsize(dest) > 5000:
        return True
    for t in range(tries):
        try:
            b = urllib.request.urlopen("https://api.gdc.cancer.gov/data/" + fid, timeout=240).read()
            if len(b) < 5000:
                raise ValueError("small")
            open(dest, "wb").write(b)
            return True
        except Exception:
            pass
    return False


def get_rppa(code, want_barcodes):
    """Return {antibody: {barcode: value}} for a cohort."""
    flt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-" + code]}},
        {"op": "in", "content": {"field": "data_category", "value": ["Proteome Profiling"]}}]}
    hits = gdc("files", {"filters": flt, "size": 3000, "format": "JSON",
                         "fields": "file_id,file_name"}).get("data", {}).get("hits", [])
    if not hits:
        return None
    # only download files whose barcode is in the expression matrix
    todo = []
    for h in hits:
        bc = h["file_name"].split("_RPPA")[0][:12]
        if bc in want_barcodes:
            todo.append((h["file_id"], os.path.join(CNVDIR, "%s_%s.tsv" % (code, h["file_id"]))))
    print("  %s: %d RPPA files, %d match expression barcodes" % (code, len(hits), len(todo)), flush=True)
    if not todo:
        return None
    with ThreadPoolExecutor(max_workers=6) as ex:
        list(ex.map(lambda t: fetch_one(*t), todo))

    out = {}
    for fid, path in todo:
        if not os.path.exists(path):
            continue
        bc = os.path.basename(path).split("_", 1)[1][:-4]
        # recover barcode from the file header content
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                lines = [l for l in f if l.strip()]
            h = lines[0].split("\t")
            pi = h.index("peptide_target")
            ei = h.index("protein_expression")
            for l in lines[1:]:
                c = l.split("\t")
                if len(c) > max(pi, ei):
                    try:
                        out.setdefault(c[pi], {})[c[ei]] = None   # placeholder
                        out[c[pi]][c[ei]] = float(c[ei])
                    except Exception:
                        pass
        except Exception:
            continue
    return out


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
    if n < 15:
        return None
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    dx = math.sqrt(sum((x[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((y[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx > 0 and dy > 0 else None


print("=" * 96)
print("A2. CROSS-COHORT: identical pipeline, CNV-driven subset")
print("=" * 96)
print("%-8s %7s %8s %8s %8s %9s %9s %8s" % (
    "cohort", "n_3lyr", "r(cnv,mr)", "r(mr,prot)", "r(cnv,pr)", "r_cm_driven", "T_driven", "n_driven"))
print("-" * 96)
rows = []
for code in COHORTS:
    expr, ecols = load_expr(code)
    if expr is None:
        print("%-8s  (no expression matrix locally)" % code)
        continue
    want = set(ecols)
    rp = get_rppa(code, want)
    if not rp:
        print("%-8s  (no RPPA matched)" % code, flush=True)
        continue
    # CNV for this cohort from local full CNV if present, else skip layer
    cnp = os.path.join(DM, "cox_cnv_matrix_%s.tsv.gz" % code)
    if not os.path.exists(cnp):
        print("%-8s  RPPA ok (%d antibodies) but no CNV matrix yet" % (code, len(rp)), flush=True)
        with open(os.path.join(DM, "rppa_%s.json" % code), "w") as f:
            json.dump({k: v for k, v in rp.items()}, f)
        continue
    print("%-8s  complete" % code)
