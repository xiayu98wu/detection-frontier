"""A1-QC2: is r(mRNA,protein)=0.31 MESO-specific, or an artefact of my pipeline?

Run the identical computation on other TCGA cohorts that have RPPA. If they return
~0.5 (the published pan-cancer benchmark), then MESO's low value is real.
"""
import json, urllib.request, os, io, tarfile, gzip, math, statistics, collections
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

DM = os.path.join(_ROOT, "data")
FULL = os.environ.get("TCGA_FULL_DIR",
                      os.path.join(_ROOT, "data", "TCGA_full"))

COHORTS = ["BRCA", "LUAD", "KIRC", "HNSC", "UCEC", "COAD"]


def gdc_files(project, category="Proteome Profiling"):
    flt = {"op": "and", "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": [project]}},
        {"op": "in", "content": {"field": "data_category", "value": [category]}}]}
    req = urllib.request.Request("https://api.gdc.cancer.gov/files",
        data=json.dumps({"filters": flt, "size": 2000, "format": "JSON",
                         "fields": "file_id,file_name"}).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"})
    return json.load(urllib.request.urlopen(req, timeout=120))["data"]["hits"]


def fetch_rppa(project):
    hits = gdc_files(project)
    if not hits:
        return None
    ids = [h["file_id"] for h in hits]
    # chunked download
    out = {}
    CH = 60
    for i in range(0, len(ids), CH):
        chunk = ids[i:i + CH]
        for t in range(4):
            try:
                raw = urllib.request.urlopen(urllib.request.Request(
                    "https://api.gdc.cancer.gov/data",
                    data=json.dumps({"ids": chunk}).encode(),
                    headers={"Content-Type": "application/json",
                             "User-Agent": "Mozilla/5.0"}), timeout=600).read()
                tf = tarfile.open(fileobj=io.BytesIO(raw))
                for m in tf.getmembers():
                    if not (m.isfile() and m.name.endswith(".tsv")):
                        continue
                    bc = os.path.basename(m.name).split("_RPPA")[0][:12]
                    txt = tf.extractfile(m).read().decode("utf-8", "replace")
                    lines = [l for l in txt.splitlines() if l.strip()]
                    h = lines[0].split("\t")
                    pi = h.index("peptide_target")
                    ei = h.index("protein_expression")
                    for l in lines[1:]:
                        c = l.split("\t")
                        if len(c) > max(pi, ei):
                            try:
                                out.setdefault(c[pi], {})[bc] = float(c[ei])
                            except Exception:
                                pass
                break
            except Exception:
                if t == 3:
                    return None
    return out


ALIAS = {"BAP1C4": "BAP1", "PDL1": "CD274", "P53": "TP53", "VIMENTIN": "VIM",
         "ECADHERIN": "CDH1", "E_CADHERIN": "CDH1", "NCADHERIN": "CDH2",
         "CYCLIND1": "CCND1", "CYCLINB1": "CCNB1", "CYCLINE1": "CCNE1",
         "CD45": "PTPRC", "CD31": "PECAM1", "HER2": "ERBB2", "HER3": "ERBB3",
         "4EBP1": "EIF4EBP1", "RB": "RB1", "P21": "CDKN1A", "P27": "CDKN1B",
         "SNAIL": "SNAI1", "TWIST": "TWIST1", "IDO": "IDO1",
         "GRANZYMEB": "GZMB", "CD8": "CD8A", "CD3": "CD3E", "CD20": "MS4A1",
         "ACC1": "ACACA", "ACC": "ACACA", "BCLXL": "BCL2L1", "CHK1": "CHEK1",
         "CHK2": "CHEK2", "HSP70": "HSPA1A", "MAPK": "MAPK1", "MEK1": "MAP2K1",
         "MIG6": "ERRFI1", "NFKB": "NFKB1", "P38": "MAPK14", "P62": "SQSTM1",
         "PARP": "PARP1", "RB1": "RB1", "S6": "RPS6", "STAT5ALPHA": "STAT5A"}


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    num = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    dx = math.sqrt(sum((x[i] - mx) ** 2 for i in range(n)))
    dy = math.sqrt(sum((y[i] - my) ** 2 for i in range(n)))
    return num / (dx * dy) if dx > 0 and dy > 0 else None


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


print("=" * 88)
print("r(mRNA, protein) ACROSS TCGA COHORTS - identical pipeline")
print("=" * 88)
print("%-8s %8s %10s %10s %10s" % ("cohort", "n_rppa", "n_genes", "mean_r", "median_r"))
print("-" * 52)
results = []
for code in COHORTS:
    expr, ecols = load_expr(code)
    if expr is None:
        print("%-8s  (no local expression matrix)" % code)
        continue
    rp = fetch_rppa("TCGA-" + code)
    if not rp:
        print("%-8s  (no RPPA via GDC)" % code)
        continue
    eidx = {b: i for i, b in enumerate(ecols)}
    rs = list(next(iter(rp.values())).keys())
    rs = [s for s in rs if s in eidx]
    cors = []
    for ab, vals in rp.items():
        g = ALIAS.get(ab, ab)
        if g not in expr:
            continue
        pairs = [(vals[s], expr[g][eidx[s]]) for s in rs if s in vals]
        if len(pairs) < 25:
            continue
        r = pearson([a for a, _ in pairs], [b for _, b in pairs])
        if r is not None:
            cors.append(r)
    if cors:
        print("%-8s %8d %10d %+10.3f %+10.3f" % (code, len(rs), len(cors),
              statistics.mean(cors), statistics.median(cors)))
        results.append((code, len(rs), statistics.mean(cors)))

print("\n  published pan-cancer RPPA500 benchmark (Pearson) = +0.510")
if results:
    print("  my pipeline, non-MESO cohorts mean = %+.3f" % statistics.mean([r[2] for r in results]))
