#!/usr/bin/env python3
"""Re-run the prespecified PubMed candidate screen.

Requires the environment variable NCBI_EMAIL (and optionally NCBI_API_KEY).
Without network access or credentials, the cached record list
(data/screen/records.csv, search date recorded in data/screen/search_log.txt)
is left unchanged. The screen identifies candidate worked examples only; it is
not a systematic review and is not used to estimate prevalence.
"""
import csv
import json
import os
import sys
import urllib.parse
import urllib.request
from datetime import date

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCREEN = os.path.join(ROOT, "data", "screen")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def call(endpoint, **params):
    params.update(tool="pairedbounds", email=os.environ["NCBI_EMAIL"])
    if os.environ.get("NCBI_API_KEY"):
        params["api_key"] = os.environ["NCBI_API_KEY"]
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)


def main():
    if not os.environ.get("NCBI_EMAIL"):
        print("NCBI_EMAIL not set; keeping cached data/screen/records.csv")
        return 0
    query = open(os.path.join(SCREEN, "query.txt")).read().strip()
    res = call("esearch.fcgi", db="pubmed", term=query, retmax=500, retmode="json")
    ids = res["esearchresult"]["idlist"]
    summ = call("esummary.fcgi", db="pubmed", id=",".join(ids), retmode="json")["result"] if ids else {}
    rows = []
    for i in ids:
        s = summ[i]
        doi = next((a["value"] for a in s.get("articleids", []) if a["idtype"] == "doi"), "")
        rows.append({"pmid": i, "pubdate": s.get("pubdate", ""), "journal": s.get("source", ""),
                     "first_author": (s.get("authors") or [{"name": ""}])[0]["name"],
                     "title": s.get("title", ""), "doi": doi})
    with open(os.path.join(SCREEN, "records.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["pmid"])
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(SCREEN, "search_log.txt"), "w") as f:
        f.write(f"database: PubMed (E-utilities)\nsearch_date: {date.today().isoformat()}\n"
                f"records: {len(rows)}\n")
    print(f"{len(rows)} records written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
