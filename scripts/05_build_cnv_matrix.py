import json, urllib.request, os, glob, gzip, collections
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file
DM = os.path.join(_ROOT, "data")
D = os.path.join(DM,"MESO_cnv_ascat3")
flt={"op":"and","content":[
 {"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-MESO"]}},
 {"op":"in","content":{"field":"data_type","value":["Gene Level Copy Number"]}},
 {"op":"in","content":{"field":"analysis.workflow_type","value":["ASCAT3"]}}]}
req=urllib.request.Request("https://api.gdc.cancer.gov/files",
    data=json.dumps({"filters":flt,"size":200,"format":"JSON",
      "fields":"file_id,associated_entities.entity_submitter_id"}).encode(),
    headers={"Content-Type":"application/json","User-Agent":"Mozilla/5.0"})
hits=json.load(urllib.request.urlopen(req,timeout=120))["data"]["hits"]
fid2bar={}
for h in hits:
    bars=set()
    for e in (h.get("associated_entities") or []):
        sid=e.get("entity_submitter_id","")
        if len(sid)>=15 and sid[13:15]=="01": bars.add(sid[:12])
    if len(bars)==1: fid2bar[h["file_id"]]=bars.pop()
    else: print("  ambiguous:", h["file_id"], bars)
print(f"clean tumour mapping: {len(fid2bar)}/87")

# expression barcodes
with gzip.open(os.path.join(DM,"MESO_expr_full.tsv.gz"),"rt",encoding="utf-8",errors="replace") as g:
    hdr=g.readline().rstrip("\n").split("\t"); ecols=[c.strip('"')[:12] for c in hdr[1:]]
eset=set(ecols)
print(f"expression barcodes: {len(eset)} (unique {len(set(ecols))})")
inter=set(fid2bar.values()) & eset
print(f"CNV tumour barcodes also in expression: {len(inter)}")

# rebuild CNV matrix keyed by tumour barcode
data=collections.defaultdict(dict); got=[]
for p in sorted(glob.glob(os.path.join(D,"*.tsv"))):
    fid=os.path.basename(p)[:-4]
    bc=fid2bar.get(fid)
    if bc is None: continue
    got.append(bc)
    with open(p,encoding="utf-8",errors="replace") as f:
        h=f.readline().rstrip("\n").split("\t")
        ci=h.index("copy_number"); gi=h.index("gene_name")
        for l in f:
            c=l.rstrip("\n").split("\t")
            if len(c)>max(ci,gi) and c[ci].strip():
                data[c[gi]][bc]=int(float(c[ci]))
print(f"CNV matrix: {len(data)} genes x {len(got)} tumours")
bcs=sorted(got)
with gzip.open(os.path.join(DM,"MESO_cnv_matrix.tsv.gz"),"wt",encoding="utf-8",newline="") as f:
    f.write("gene\t"+"\t".join(bcs)+"\n")
    for g in sorted(data):
        f.write(g+"\t"+"\t".join(str(data[g].get(b,"NA")) for b in bcs)+"\n")
print("[saved] MESO_cnv_matrix.tsv.gz (tumour-only, barcode-keyed)")

print("\n=== 9p21 structure among 87 TUMOURS ===")
pairs=collections.Counter()
for b in bcs:
    pairs[(data["CDKN2A"].get(b), data["MTAP"].get(b))]+=1
for k,v in sorted(pairs.items(),key=lambda x:-x[1])[:8]:
    print(f"  CDKN2A={k[0]} MTAP={k[1]}  n={v}")
hom=sum(1 for b in bcs if data["CDKN2A"].get(b)==0)
print(f"\n  CDKN2A homozygous deletion (CN=0): {hom}/87 = {100*hom/87:.1f}%")
print(f"  MTAP   homozygous deletion (CN=0): {sum(1 for b in bcs if data['MTAP'].get(b)==0)}/87")