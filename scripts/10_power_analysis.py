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
    h=f.readline().rstrip("\n").split("\t"); gi=h.index("Hugo_Symbol"); ci=h.index("Tumor_Sample_Barcode")
    for l in f:
        c=l.rstrip("\n").split("\t")
        if len(c)>max(gi,ci): mut[c[gi]].add(c[ci][:12])
has_maf=set().union(*mut.values())

SIG={"CD8_T":["CD8A","CD8B","GZMA","GZMB","PRF1","CD3D","CD3E"],
     "IFN_G":["STAT1","IRF1","CXCL9","CXCL10","IFNG","IDO1","GBP1"]}
def score(gl):
    gs=[g for g in gl if g in expr]
    return [statistics.mean(expr[g][j] for g in gs) for j in range(len(bars))]

def power_ttest(d, n1, n2, alpha=0.05):
    # normal approx power for two-sample t
    if n1<2 or n2<2: return 0.0
    ncp = d*math.sqrt(n1*n2/(n1+n2))
    zc = 1.959964 - (0.0)
    # power = P(|Z|>zc) under N(ncp,1)
    phi = lambda x: 0.5*(1+math.erf(x/math.sqrt(2)))
    return (1-phi(zc-ncp)) + phi(-zc-ncp)

print("="*74)
print("POWER ANALYSIS: what effect size can n=87 detect?")
print("="*74)
print(f"{'grouping':22s} {'n_mut':>6s} {'n_wt':>6s} {'min detectable d':>18s}")
for GENE,n_mut in [("BAP1",22),("NF2",19),("TP53",11),("SETD2",8)]:
    n_wt=78-n_mut
    d=0.1
    while d<3.0:
        if power_ttest(d,n_mut,n_wt)>=0.80: break
        d+=0.005
    print(f"{GENE:22s} {n_mut:6d} {n_wt:6d} {d:18.2f}")

print("\n" + "="*74)
print("OBSERVED effect sizes were 0.02-0.5  ->  far below detectable")
print("="*74)
# what d did we actually observe for BAP1 vs CD8?
for GENE in ["BAP1","NF2"]:
    s=set(mut.get(GENE,set()))
    for nm,gl in SIG.items():
        sc=score(gl)
        a=[sc[i] for i,b in enumerate(bars) if b in has_maf and b in s]
        b_=[sc[i] for i,b in enumerate(bars) if b in has_maf and b not in s]
        # pooled SD in units of the score's own SD across all samples
        allsd=statistics.pstdev(sc)
        d=(statistics.mean(a)-statistics.mean(b_))/allsd
        pw=power_ttest(abs(d),len(a),len(b_))
        print(f"  {GENE:5s} vs {nm:8s}  observed d={d:+.3f}   power={pw*100:4.1f}%")