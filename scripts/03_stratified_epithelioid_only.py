import gzip, os, collections, statistics, math
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file
D = os.path.join(_ROOT, "data")
with gzip.open(os.path.join(D,"MESO_expr_full.tsv.gz"),"rt",encoding="utf-8",errors="replace") as g:
    hdr=g.readline().rstrip("\n").split("\t"); cols=[c.strip('"') for c in hdr[1:]]
    expr={}
    for line in g:
        c=line.rstrip("\n").split("\t"); expr[c[0].strip('"')]=[float(x) for x in c[1:]]
bars=[c[:12] for c in cols]
mut=collections.defaultdict(set)
with open(os.path.join(D,"MESO_mutations.maf"),encoding="utf-8") as f:
    h=f.readline().rstrip("\n").split("\t")
    gi=h.index("Hugo_Symbol"); ci=h.index("Tumor_Sample_Barcode")
    for l in f:
        c=l.rstrip("\n").split("\t")
        if len(c)>max(gi,ci): mut[c[gi]].add(c[ci][:12])
has_maf=set().union(*mut.values())
with open(os.path.join(D,"MESO_clinicalMatrix.tsv"),encoding="utf-8",errors="replace") as f:
    ch=f.readline().rstrip("\n").split("\t")
    crows=[l.rstrip("\n").split("\t") for l in f if l.strip()]
hi=ch.index("histological_type"); si=ch.index("sampleID")
hist={r[si][:12]: r[hi] for r in crows if len(r)>hi}
print("histology map:", len(hist))

SIGS={
 "CD8_T":   ["CD8A","CD8B","GZMA","GZMB","PRF1","CD3D","CD3E"],
 "CHECKPT": ["PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT"],
 "IFN_G":   ["STAT1","IRF1","CXCL9","CXCL10","IFNG","IDO1","GBP1"],
 "MACROPHAGE":["CD68","CD163","MRC1","CSF1R"],
 "FIBROBLAST":["COL1A1","COL1A2","FAP","ACTA2","DCN"],
 "EMT":     ["VIM","CDH2","FN1","SNAI1","SNAI2","TWIST1","ZEB1"],
 "PROLIF":  ["MKI67","TOP2A","CCNB1","BUB1","AURKA"],
}
def score(genes):
    gs=[g for g in genes if g in expr]
    return [statistics.mean(expr[g][j] for g in gs) for j in range(len(bars))]
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

EPI = [i for i,b in enumerate(bars) if "Epithelioid" in hist.get(b,"")]
print(f"\nEpithelioid-only samples: {len(EPI)}")

print("\n" + "="*72)
print("HISTOLOGY EFFECT (all 87) - epithelioid vs non-epithelioid")
print("="*72)
NE = [i for i,b in enumerate(bars) if "Epithelioid" not in hist.get(b,"")]
print(f"{'signature':12s} {'epi_mean':>9s} {'non_mean':>9s} {'delta':>8s} {'z':>7s} {'p':>9s}")
for nm,genes in SIGS.items():
    sc=score(genes)
    z,p=mw([sc[i] for i in EPI],[sc[i] for i in NE])
    flag="  ***" if p<0.001 else ("  **" if p<0.01 else ("  *" if p<0.05 else ""))
    print(f"{nm:12s} {statistics.mean([sc[i] for i in EPI]):9.3f} {statistics.mean([sc[i] for i in NE]):9.3f} "
          f"{statistics.mean([sc[i] for i in EPI])-statistics.mean([sc[i] for i in NE]):+8.3f} {z:+7.2f} {p:9.2e}{flag}")

print("\n" + "="*72)
print("BAP1 / NF2 WITHIN EPITHELIOID ONLY")
print("="*72)
for GENE in ["BAP1","NF2"]:
    s=set(mut.get(GENE,set()))
    a=[i for i in EPI if bars[i] in has_maf and bars[i] in s]
    b_=[i for i in EPI if bars[i] in has_maf and bars[i] not in s]
    print(f"\n### {GENE}: mut n={len(a)}  wt n={len(b_)}")
    print(f"{'signature':12s} {'mut':>9s} {'wt':>9s} {'delta':>8s} {'z':>7s} {'p':>9s}")
    for nm,genes in SIGS.items():
        sc=score(genes)
        z,p=mw([sc[i] for i in a],[sc[i] for i in b_])
        flag="  ***" if p<0.001 else ("  **" if p<0.01 else ("  *" if p<0.05 else ""))
        print(f"{nm:12s} {statistics.mean([sc[i] for i in a]):9.3f} {statistics.mean([sc[i] for i in b_]):9.3f} "
              f"{statistics.mean([sc[i] for i in a])-statistics.mean([sc[i] for i in b_]):+8.3f} {z:+7.2f} {p:9.2e}{flag}")