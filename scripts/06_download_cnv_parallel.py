import json, urllib.request, os, time
from concurrent.futures import ThreadPoolExecutor, as_completed
D = os.path.join(_ROOT, "data")
OUT = os.path.join(D, "MESO_cnv_ascat3")
os.makedirs(OUT, exist_ok=True)
req=urllib.request.Request("https://api.gdc.cancer.gov/files",
    data=json.dumps({"filters":{"op":"and","content":[
     {"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-MESO"]}},
     {"op":"in","content":{"field":"data_type","value":["Gene Level Copy Number"]}},
     {"op":"in","content":{"field":"analysis.workflow_type","value":["ASCAT3"]}}]},
     "size":200,"format":"JSON"}).encode(),
    headers={"Content-Type":"application/json","User-Agent":"Mozilla/5.0"})
hits=json.load(urllib.request.urlopen(req,timeout=120))["data"]["hits"]
print(f"files to fetch: {len(hits)}", flush=True)

def get(h):
    fid=h["file_id"]; dest=os.path.join(OUT, fid+".tsv")
    if os.path.exists(dest) and os.path.getsize(dest)>10000: return fid,"cached"
    for t in range(4):
        try:
            b=urllib.request.urlopen(f"https://api.gdc.cancer.gov/data/{fid}", timeout=180).read()
            if len(b)<10000: raise ValueError("too small")
            open(dest,"wb").write(b)
            return fid, f"ok {len(b)//1024}KB"
        except Exception as e:
            if t==3: return fid, f"FAIL {type(e).__name__}"
            time.sleep(3)
import threading
import os
_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file
lock=threading.Lock(); cnt=[0]
with ThreadPoolExecutor(max_workers=6) as ex:
    futs={ex.submit(get,h):h for h in hits}
    for fu in as_completed(futs):
        fid,st=fu.result()
        with lock:
            cnt[0]+=1
            if cnt[0]%10==0 or "FAIL" in st: print(f"  {cnt[0]}/{len(hits)}  {st}", flush=True)
ok=len([f for f in os.listdir(OUT) if f.endswith(".tsv")])
print(f"DONE: {ok}/{len(hits)} files on disk", flush=True)