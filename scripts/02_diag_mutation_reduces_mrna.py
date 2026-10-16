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
mut=collections.defaultdict(set); cls=collections.defaultdict(set)
with open(os.path.join(D,"MESO_mutations.maf"),encoding="utf-8") as f:
    h=f.readline().rstrip("\n").split("\t")
    gi=h.index("Hugo_Symbol"); ci=h.index("Tumor_Sample_Barcode")
    vi=h.index("Variant_Classification") if "Variant_Classification" in h else None
    for l in f:
        c=l.rstrip("\n").split("\t")
        if len(c)>max(gi,ci):
            mut[c[gi]].add(c[ci][:12])
            if vi is not None: cls[c[gi]].add(c[vi])
has_maf=set().union(*mut.values())

print("=== BAP1 mutation classes in MESO ===")
for k,v in collections.Counter(cls.get("BAP1",[])).items(): print(f"  {k:28s} {v}")
print("=== NF2 mutation classes ===")
for k,v in collections.Counter(cls.get("NF2",[])).items(): print(f"  {k:28s} {v}")

def mw(a,b):
    n1,n2=len(a),len(b)
    if n1<3 or n2<3: return None,None
    comb=sorted([(v,0) for v in a]+[(v,1) for v in b])
    ranks=[0]*len(comb); i=0
    while i<len(comb):
        j=i
        while j+1<len(comb) and comb[j+1][0]==comb[i][0]: j+=1
        r=(i+j)/2+1
        for k in range(i,j+1): ranks[k]=r
        i=j+1
    R1=sum(ranks[k] for k in range(len(comb)) if comb[k][1]==0)
    U1=R1-n1*(n1+1)/2; mu=n1*n2/2
    tie=collections.Counter(v for v,_ in comb); N=n1+n2
    sd=math.sqrt(n1*n2/12*((N+1)-sum(t**3-t for t in tie.values())/(N*(N-1))))
    if sd==0: return None,None
    z=(U1-mu)/sd
    return z, 2*(1-0.5*(1+math.erf(abs(z)/math.sqrt(2))))

print("\n=== does the mutation reduce mRNA? ===")
print(f"{'gene':10s} {'mut_mean':>9s} {'wt_mean':>9s} {'delta':>8s} {'z':>7s} {'p':>9s}")
for GENE in ["BAP1","NF2","SETD2","TP53"]:
    if GENE not in expr: continue
    s=set(mut.get(GENE,set()))
    a=[expr[GENE][i] for i,b in enumerate(bars) if b in has_maf and b in s]
    b_=[expr[GENE][i] for i,b in enumerate(bars) if b in has_maf and b not in s]
    z,p=mw(a,b_)
    print(f"{GENE:10s} {statistics.mean(a):9.3f} {statistics.mean(b_):9.3f} {statistics.mean(a)-statistics.mean(b_):+8.3f} {z:+7.2f} {p:9.2e}")

print("\n=== BAP1 co-occurrence with other drivers ===")
bp=mut.get("BAP1",set())
for g2 in ["NF2","TP53","SETD2","LATS2"]:
    s2=mut.get(g2,set())
    both=len(bp&s2); a=len(bp-s2); b=len(s2-bp); d=len(has_maf-bp-s2)
    print(f"  BAP1 x {g2:6s}  both={both:2d}  BAP1only={a:2d}  {g2}only={b:2d}  neither={d:2d}")