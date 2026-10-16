import json, urllib.request, os, io, collections, statistics, math, gzip, tarfile
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file
D = os.path.join(_ROOT, "data")
if False:
    pass
flt={"op":"and","content":[
 {"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-MESO"]}},
 {"op":"in","content":{"field":"data_category","value":["Proteome Profiling"]}}]}
req=urllib.request.Request("https://api.gdc.cancer.gov/files", data=json.dumps({"filters":flt,"size":300,"format":"JSON"}).encode(),
    headers={"Content-Type":"application/json","User-Agent":"Mozilla/5.0"})
hits=json.load(urllib.request.urlopen(req,timeout=120))["data"]["hits"]
raw=urllib.request.urlopen(urllib.request.Request("https://api.gdc.cancer.gov/data",
    data=json.dumps({"ids":[h["file_id"] for h in hits]}).encode(),
    headers={"Content-Type":"application/json","User-Agent":"Mozilla/5.0"}), timeout=600).read()
tf=tarfile.open(fileobj=io.BytesIO(raw))
prot={}; samp=[]
def f2(s):
    try: return float(s)
    except: return None
for m in tf.getmembers():
    if not (m.isfile() and m.name.endswith(".tsv")): continue
    barcode=os.path.basename(m.name).split("_RPPA")[0][:12]
    txt=tf.extractfile(m).read().decode("utf-8","replace")
    lines=[l for l in txt.splitlines() if l.strip()]
    h=lines[0].split("\t"); pi=h.index("peptide_target"); ei=h.index("protein_expression")
    samp.append(barcode)
    for l in lines[1:]:
        c=l.split("\t")
        if len(c)>max(pi,ei):
            v=f2(c[ei])
            if v is not None: prot.setdefault(c[pi],{})[barcode]=v
print(f"RPPA: {len(prot)} proteins x {len(samp)} samples")
with open(os.path.join(D,"MESO_rppa.tsv"),"w",encoding="utf-8",newline="") as f:
    f.write("peptide_target\t"+"\t".join(samp)+"\n")
    for k in prot:
        f.write(k+"\t"+"\t".join(str(prot[k].get(s,"NA")) for s in samp)+"\n")

with gzip.open(os.path.join(D,"MESO_expr_full.tsv.gz"),"rt",encoding="utf-8",errors="replace") as g:
    hdr=g.readline().rstrip("\n").split("\t"); cols=[c.strip('"') for c in hdr[1:]]
    expr={}
    for line in g:
        c=line.rstrip("\n").split("\t"); expr[c[0].strip('"')]=[float(x) for x in c[1:]]
bars=[c[:12] for c in cols]
common=[b for b in bars if b in set(samp)]
print(f"samples with BOTH RPPA and expression: {len(common)}")
cands=[k for k in prot if "BAP1" in k.upper()]
print("BAP1 RPPA peptides:", cands)

def score(genes):
    gs=[g for g in genes if g in expr]
    return {bars[j]: statistics.mean(expr[g][j] for g in gs) for j in range(len(bars))}
def spearman(x,y):
    n=len(x)
    if n<8: return None,None
    def rk(v):
        o=sorted(range(n), key=lambda i:v[i]); r=[0]*n; i=0
        while i<n:
            j=i
            while j+1<n and v[o[j+1]]==v[o[i]]: j+=1
            for k in range(i,j+1): r[o[k]]=(i+j)/2+1
            i=j+1
        return r
    rx,ry=rk(x),rk(y); mx,my=sum(rx)/n,sum(ry)/n
    num=sum((rx[i]-mx)*(ry[i]-my) for i in range(n))
    den=math.sqrt(sum((rx[i]-mx)**2 for i in range(n))*sum((ry[i]-my)**2 for i in range(n)))
    if den==0: return None,None
    r=num/den; t=r*math.sqrt((n-2)/max(1e-12,1-r*r))
    return r, 2*(1-0.5*(1+math.erf(abs(t)/math.sqrt(2))))

SIG={"CD8_T":["CD8A","CD8B","GZMA","GZMB","PRF1","CD3D","CD3E"],
     "CHECKPT":["PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT"],
     "IFN_G":["STAT1","IRF1","CXCL9","CXCL10","IFNG","IDO1","GBP1"],
     "MACROPHAGE":["CD68","CD163","MRC1","CSF1R"],"EMT":["VIM","CDH2","FN1","SNAI1","ZEB1"]}
print("\n"+"="*70)
print("BAP1 PROTEIN vs immune signatures (RPPA n=%d)" % len(common))
print("="*70)
if cands:
    bv=prot[cands[0]]
    for nm,gl in SIG.items():
        sc=score(gl); pairs=[(bv[b],sc[b]) for b in common if b in bv]
        r,p=spearman([x for x,_ in pairs],[y for _,y in pairs])
        print(f"  BAP1_prot vs {nm:12s} rho={r:+.3f}  p={p:.3e}  n={len(pairs)}")
else: print("  BAP1 not on panel")

print("\n"+"="*70)
print("UNBIASED: RPPA proteins tracking CD8_T (top 25 by |rho|)")
print("="*70)
sc=score(SIG["CD8_T"]); res=[]
for pname,pv in prot.items():
    pairs=[(pv[b],sc[b]) for b in common if b in pv and b in sc]
    if len(pairs)<20: continue
    r,p=spearman([x for x,_ in pairs],[y for _,y in pairs])
    if r is not None: res.append((r,p,pname,len(pairs)))
res.sort(key=lambda t:-abs(t[0]))
print(f"  {'protein':30s} {'rho':>7s} {'p':>10s} {'n':>4s}")
for r,p,pname,n in res[:25]:
    flag="  ***" if p<0.001 else ("  **" if p<0.01 else ("  *" if p<0.05 else ""))
    print(f"  {pname:30s} {r:+7.3f} {p:10.2e} {n:4d}{flag}")