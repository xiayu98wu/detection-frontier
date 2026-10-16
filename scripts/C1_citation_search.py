"""Targeted citation search for the Reproducibility-framed manuscript.
Uses PubMed E-utilities (no API key needed for low volume).
"""
import json, urllib.request, urllib.parse, time, os, sys

# PubMed titles contain non-breaking spaces and other characters that the default
# Windows console codepage (cp936/cp1252) cannot encode, which aborted a query
# mid-run during development. Force UTF-8 on stdout so the search is portable.
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repository root, resolved from this file

OUT = os.path.join(_ROOT, "results")
os.makedirs(OUT, exist_ok=True)
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def esearch(term, retmax=6, sort="relevance"):
    q = urllib.parse.urlencode({"db": "pubmed", "term": term, "retmax": retmax,
                                "retmode": "json", "sort": sort})
    url = BASE + "esearch.fcgi?" + q
    return json.load(urllib.request.urlopen(url, timeout=45))["esearchresult"]


def esummary(ids):
    if not ids:
        return {}
    q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "json"})
    return json.load(urllib.request.urlopen(BASE + "esummary.fcgi?" + q, timeout=45)).get("result", {})


QUERIES = [
    ("power / underpowered", '(underpowered OR "statistical power") AND (cancer OR oncology) AND (reproducib* OR replicat*)'),
    ("TCGA rare cancer", '(TCGA OR "Cancer Genome Atlas") AND ("rare cancer" OR "rare tumour" OR "rare tumor")'),
    ("absence of evidence", '"absence of evidence" AND (association OR biomarker) AND (power OR sample size)'),
    ("winner's curse", "\"winner's curse\" AND (genetic association OR biomarker)"),
    ("study size effect estimates", '"study size" AND ("effect size" OR "effect estimates") AND (bias OR overestimat*)'),
    ("sample size justification", '"sample size" AND (justification OR reporting) AND (guideline OR practice) AND (biomedical OR clinical research)'),
    ("mesothelioma BAP1 immune", 'BAP1 AND mesothelioma AND (immune OR immunotherapy)'),
    ("mesothelioma CDKN2A prognosis", 'CDKN2A AND mesothelioma AND (prognosis OR survival OR deletion)'),
    ("multiplicity testing cancer genomics", '(multiplicity OR "multiple testing") AND (cancer genomics OR TCGA OR bioinformatics) AND (false discovery OR reporting)'),
    ("preregistration observational", '(preregistration OR "pre-registration") AND (observational OR cohort) AND (bias OR transparency)'),
    ("p-value problem epidemiology", '"p-value" AND (epidemiology OR "observational research") AND (misinterpret* OR limit*)'),
    ("underpowered replication failure", '(replication OR reproducibility) AND (failure OR crisis) AND (underpowered OR "low power")'),
]

results = {}
print("=" * 96)
print("CITATION SEARCH (PubMed E-utilities)")
print("=" * 96)
for label, term in QUERIES:
    try:
        r = esearch(term, retmax=8)
        ids = r.get("idlist", [])
        tot = r.get("count", "?")
        print("\n### %s   (total hits: %s)" % (label, tot))
        print("    query: %s" % term[:110])
        if not ids:
            print("    (no hits)")
            continue
        s = esummary(ids)
        entries = []
        for i in ids:
            d = s.get(i, {})
            if not d:
                continue
            title = d.get("title", "")
            journal = d.get("source", "")
            year = (d.get("pubdate", "") or "")[:4]
            auth = d.get("sortfirstauthor", "")
            pmid = d.get("uid", "")
            entries.append({"pmid": pmid, "title": title, "journal": journal,
                            "year": year, "author": auth})
            print("    - [%s] %s %s (%s) PMID:%s" % (year, title[:88], journal, auth, pmid))
        results[label] = entries
        time.sleep(0.4)
    except Exception as e:
        print("  ERROR %s: %s" % (label, e))

with open(os.path.join(OUT, "citations.json"), "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False, indent=2)
print("\n[saved] citations.json")
print("total queries: %d, total entries: %d" % (len(results), sum(len(v) for v in results.values())))
