import tarfile, io, os, gzip, collections
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file
D = os.path.join(_ROOT, "data")
raw = open(os.path.join(D,"MESO_maf.tar.gz"),"rb").read()
tf = tarfile.open(fileobj=io.BytesIO(raw))
hdr=None; rows=[]
for m in tf.getmembers():
    if not m.isfile() or not m.name.endswith(".maf.gz"): continue
    b = tf.extractfile(m).read()
    if b[:2]==b"\x1f\x8b": b = gzip.decompress(b)
    lines=[l for l in b.decode("utf-8","replace").splitlines() if l and not l.startswith("#")]
    if not lines: continue
    h=lines[0].split("\t")
    if "Hugo_Symbol" not in h: continue
    if hdr is None: hdr=h
    rows += lines[1:]
print(f"MAF files merged, rows={len(rows)}")
with open(os.path.join(D,"MESO_mutations.maf"),"w",encoding="utf-8",newline="") as f:
    f.write("\t".join(hdr)+"\n"+"\n".join(rows)+"\n")

gi=hdr.index("Hugo_Symbol"); ci=hdr.index("Tumor_Sample_Barcode")
vi=hdr.index("Variant_Classification") if "Variant_Classification" in hdr else None
cnt=collections.defaultdict(set)
for l in rows:
    c=l.split("\t")
    if len(c)>max(gi,ci): cnt[c[gi]].add(c[ci][:12])
ns=len(set().union(*cnt.values()))
print(f"unique samples: {ns}")
print("\n=== TOP 25 MUTATED ===")
for g,s in sorted(cnt.items(),key=lambda x:-len(x[1]))[:25]:
    print(f"  {g:12s} {len(s):3d}  {100*len(s)/ns:5.1f}%")
print("\n=== TARGETED (design premise) ===")
for k in ["BAP1","NF2","CDKN2A","MTAP","TP53","SETD2","PBRM1","LATS2","TRAF7","TTN","SF3B1","KMT2D"]:
    n=len(cnt.get(k,[]))
    print(f"  {k:9s} {n:3d}  {100*n/ns:5.1f}%")