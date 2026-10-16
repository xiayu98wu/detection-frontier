import gzip, os, collections, statistics, math
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file
DM = os.path.join(_ROOT, "data")
with gzip.open(os.path.join(DM,"MESO_expr_full.tsv.gz"),"rt",encoding="utf-8",errors="replace") as g:
    hdr=g.readline().rstrip("\n").split("\t"); cols=[c.strip('"') for c in hdr[1:]]
    expr={}
    for line in g:
        c=line.rstrip("\n").split("\t"); expr[c[0].strip('"')]=[float(x) for x in c[1:]]
bars=[c[:12] for c in cols]
cnv=collections.defaultdict(dict)
with gzip.open(os.path.join(DM,"MESO_cnv_matrix.tsv.gz"),"rt",encoding="utf-8") as f:
    h=f.readline().rstrip("\n").split("\t"); cb=[c.strip('"') for c in h[1:]]
    for l in f:
        c=l.rstrip("\n").split("\t"); gname=c[0]
        for i,b in enumerate(cb):
            v=c[i+1]
            if v!="NA": cnv[gname][b]=int(float(v))
mut=collections.defaultdict(set)
with open(os.path.join(DM,"MESO_mutations.maf"),encoding="utf-8") as f:
    h=f.readline().rstrip("\n").split("\t"); gi=h.index("Hugo_Symbol"); ci=h.index("Tumor_Sample_Barcode")
    for l in f:
        c=l.rstrip("\n").split("\t")
        if len(c)>max(gi,ci): mut[c[gi]].add(c[ci][:12])
with open(os.path.join(DM,"MESO_clinicalMatrix.tsv"),encoding="utf-8",errors="replace") as f:
    ch=f.readline().rstrip("\n").split("\t")
    crows=[l.rstrip("\n").split("\t") for l in f if l.strip()]
hi=ch.index("histological_type"); si=ch.index("sampleID")
hist={r[si][:12]:r[hi] for r in crows if len(r)>hi}
print(f"expr {len(expr)}x{len(bars)}  cnv genes={len(cnv)}  cnv bars={len(cb)}")
print(f"cnv bars in expr: {sum(1 for b in cb if b in set(bars))}/{len(cb)}")

SIGS={
 "CD8_T":   ["CD8A","CD8B","GZMA","GZMB","PRF1","CD3D","CD3E"],
 "CHECKPT": ["PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT"],
 "IFN_G":   ["STAT1","IRF1","CXCL9","CXCL10","IFNG","IDO1","GBP1"],
 "MACROPHAGE":["CD68","CD163","MRC1","CSF1R"],
 "FIBROBLAST":["COL1A1","COL1A2","FAP","ACTA2","DCN"],
 "EMT":     ["VIM","CDH2","FN1","SNAI1","SNAI2","TWIST1","ZEB1"],
 "PROLIF":  ["MKI67","TOP2A","CCNB1","BUB1","AURKA"],
 "ANTIGEN": ["HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","PSMB9"],
}
def score(gl):
    gs=[g for g in gl if g in expr]
    return {bars[j]: statistics.mean(expr[g][j] for g in gs) for j in range(len(bars))}
def mw(a,b):
    n1,n2=len(a),len(b)
    if n1<3 or n2<3: return None,None
    comb=sorted([(v,0) for v in a]+[(v,1) for v in b]); ranks=[0]*len(comb); i=0
    while i<len(comb):
        j=i
        while j+1<len(comb) and comb[j+1][0]==comb[i][0]: j+=1
        r=(i+j)/2+1
        for k in range(i,j+1): ranks[k]=r
        i=j+1
    U1=sum(ranks[k] for k in range(len(comb)) if comb[k][1]==0)-n1*(n1+1)/2
    tie=collections.Counter(v for v,_ in comb); N=n1+n2
    sd=math.sqrt(n1*n2/12*((N+1)-sum(t**3-t for t in tie.values())/(N*(N-1))))
    if sd==0: return None,None
    z=(U1-n1*n2/2)/sd
    return z, 2*(1-0.5*(1+math.erf(abs(z)/math.sqrt(2))))

def test(name, grp):
    A=[b for b in bars if grp.get(b)==1]
    B=[b for b in bars if grp.get(b)==0]
    print(f"\n### {name}:  altered n={len(A)}   intact n={len(B)}")
    print(f"{'signature':12s} {'alt':>9s} {'intact':>9s} {'delta':>8s} {'z':>7s} {'p':>10s}")
    out=[]
    for nm,gl in SIGS.items():
        sc=score(gl)
        a=[sc[b] for b in A]; b_=[sc[b] for b in B]
        z,p=mw(a,b_)
        if z is None: continue
        flag="  ***" if p<0.001 else ("  **" if p<0.01 else ("  *" if p<0.05 else ""))
        print(f"{nm:12s} {statistics.mean(a):9.3f} {statistics.mean(b_):9.3f} "
              f"{statistics.mean(a)-statistics.mean(b_):+8.3f} {z:+7.2f} {p:10.2e}{flag}")
        out.append((nm,p))
    return out

# 9p21 homozygous deletion
for gene in ["CDKN2A","MTAP"]:
    grp={b:(1 if cnv[gene].get(b)==0 else 0) for b in bars if b in cnv[gene]}
    test(f"{gene} HOMDEL (CN=0)", grp)

# any loss of CDKN2A (CN<=1) 
grp={b:(1 if cnv["CDKN2A"].get(b,99)<=1 else 0) for b in bars if b in cnv["CDKN2A"]}
test("CDKN2A any loss (CN<=1)", grp)

# BAP1 CNV loss
grp={b:(1 if cnv["BAP1"].get(b,99)<=1 else 0) for b in bars if b in cnv["BAP1"]}
test("BAP1 CNV loss (CN<=1)", grp)

# combined 9p21+BAP1
grp={b:(1 if (cnv["CDKN2A"].get(b,99)==0 and cnv["BAP1"].get(b,99)<=1) else 0) for b in bars}
test("9p21del AND BAP1loss", grp)