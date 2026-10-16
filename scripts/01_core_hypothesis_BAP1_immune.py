import gzip, os, collections, random, math, statistics
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file
D = os.path.join(_ROOT, "data")

# ---- expression ----
with gzip.open(os.path.join(D,"MESO_expr_full.tsv.gz"),"rt",encoding="utf-8",errors="replace") as g:
    hdr=g.readline().rstrip("\n").split("\t"); cols=[c.strip('"') for c in hdr[1:]]
    expr={}
    for line in g:
        c=line.rstrip("\n").split("\t"); expr[c[0].strip('"')]=[float(x) for x in c[1:]]
bars=[c[:12] for c in cols]
print(f"expr: {len(expr)} genes x {len(cols)} samples")

# ---- mutations ----
mut=collections.defaultdict(set)
with open(os.path.join(D,"MESO_mutations.maf"),encoding="utf-8") as f:
    h=f.readline().rstrip("\n").split("\t"); gi=h.index("Hugo_Symbol"); ci=h.index("Tumor_Sample_Barcode")
    for l in f:
        c=l.rstrip("\n").split("\t")
        if len(c)>max(gi,ci): mut[c[gi]].add(c[ci][:12])
print(f"maf: samples={len(set().union(*mut.values()))}")

NSAMP=len(cols)
def status(gene):
    s=set(mut.get(gene,set()))
    return [1 if b in s else 0 for b in bars]
# Only test samples with a MAF entry (else mut status unknown)
has_maf=set().union(*mut.values())
known=[b in has_maf for b in bars]
print(f"samples with MAF: {sum(known)}/{NSAMP}")

# ---- signatures ----
SIGS={
 "CD8_T":   ["CD8A","CD8B","GZMA","GZMB","PRF1","CD3D","CD3E"],
 "CHECKPT": ["PDCD1","CD274","CTLA4","LAG3","HAVCR2","TIGIT"],
 "IFN_G":   ["STAT1","IRF1","CXCL9","CXCL10","IFNG","IDO1","GBP1"],
 "CYTOLYTIC":["GZMA","PRF1"],
 "MACROPHAGE":["CD68","CD163","MRC1","CSF1R"],
 "FIBROBLAST":["COL1A1","COL1A2","FAP","ACTA2","DCN"],
 "EMT":     ["VIM","CDH2","FN1","SNAI1","SNAI2","TWIST1","ZEB1"],
 "PROLIF":  ["MKI67","TOP2A","CCNB1","BUB1","AURKA"],
}
def zscore(genes):
    gs=[g for g in genes if g in expr]
    if len(gs)<2: return None
    vecs=[expr[g] for g in gs]
    z=[]
    for j in range(NSAMP):
        vals=[v[j] for v in vecs]
        m=statistics.mean(vals); sd=statistics.pstdev(vals) or 1
        z.append((vals[j]-m)/sd)   # per-gene z is degenerate; use mean instead below
    return None
def score(genes):
    gs=[g for g in genes if g in expr]
    if len(gs)<2: return None
    return [statistics.mean(expr[g][j] for g in gs) for j in range(NSAMP)]

def mannwhitney(a,b):
    # rank-based, exact-ish normal approx
    n1,n2=len(a),len(b)
    if n1<3 or n2<3: return None,None
    comb=[(v,0) for v in a]+[(v,1) for v in b]
    comb.sort()
    # average ranks
    ranks={}
    i=0
    while i<len(comb):
        j=i
        while j+1<len(comb) and comb[j+1][0]==comb[i][0]: j+=1
        r=(i+j)/2+1
        for k in range(i,j+1): ranks.setdefault(k,r)
        i=j+1
    R1=sum(ranks[k] for k in range(len(comb)) if comb[k][1]==0)
    U1=R1-n1*(n1+1)/2
    mu=n1*n2/2
    tie=collections.Counter(v for v,_ in comb)
    N=n1+n2
    tiecorr=sum(t**3-t for t in tie.values())
    sd=math.sqrt(n1*n2/12*((N+1)-tiecorr/(N*(N-1))))
    if sd==0: return None,None
    z=(U1-mu)/sd
    return z, 2*(1-0.5*(1+math.erf(abs(z)/math.sqrt(2))))

print("\n" + "="*76)
print("BAP1 mutation status vs phenotype (tumours with MAF, n=%d mut / %d wt)" %
      (sum(s and k for s,k in zip(status("BAP1"),known)), 0))
print("="*76)
for GENE in ["BAP1","NF2","TP53","SETD2"]:
    st=status(GENE)
    grp=[(s,k) for s,k in zip(st,known) if k]
    a_idx=[i for i,(s,k) in enumerate(zip(st,known)) if k and s==1]
    b_idx=[i for i,(s,k) in enumerate(zip(st,known)) if k and s==0]
    print(f"\n### {GENE}:  mutant n={len(a_idx)}   wildtype n={len(b_idx)}")
    print(f"{'signature':12s} {'mut_mean':>9s} {'wt_mean':>9s} {'delta':>8s} {'z':>7s} {'p':>9s}")
    for nm,genes in SIGS.items():
        sc=score(genes)
        if sc is None: continue
        a=[sc[i] for i in a_idx]; b=[sc[i] for i in b_idx]
        z,p=mannwhitney(a,b)
        if z is None: continue
        flag="  ***" if p<0.001 else ("  **" if p<0.01 else ("  *" if p<0.05 else ""))
        print(f"{nm:12s} {statistics.mean(a):9.3f} {statistics.mean(b):9.3f} {statistics.mean(a)-statistics.mean(b):+8.3f} {z:+7.2f} {p:9.2e}{flag}")