#!/usr/bin/env python3
"""Refresh the bibliographic metadata of the census papers from the arXiv API.

This is the one script here that needs the network, and it is not part of `make all`: its output,
data/census_papers.csv, is committed, so every other step (and every consumer of the bibliography)
runs offline. Run it only when data/published_spectra.csv gains an arXiv ID.

Resumes from the committed CSV and fetches only IDs missing from it, so a re-run after adding one
paper costs one request. The arXiv API asks for no more than one request every three seconds; the
batch size and delays here stay well inside that.

Reads  data/published_spectra.csv (for the arXiv IDs) and data/census_papers.csv (as cache).
Writes data/census_papers.csv.
Pure standard library, but requires outbound HTTPS.
"""
import csv, os, re, sys, time, urllib.request, urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def _p(*a): return os.path.join(ROOT, *a)

API = "https://export.arxiv.org/api/query"
UA = {"User-Agent": "search-budget/1.0 (https://github.com/rodem-hep/search-budget)"}
BATCH, GAP = 25, 6
COLS = ["arxiv", "title", "journal", "doi", "year"]

ids = []
for r in csv.DictReader(open(_p("data", "published_spectra.csv"))):
    for a in r["arxiv"].split():
        if a not in ids:
            ids.append(a)

cache = {}
if os.path.exists(_p("data", "census_papers.csv")):
    cache = {r["arxiv"]: r for r in csv.DictReader(open(_p("data", "census_papers.csv")))}

todo = [i for i in ids if i not in cache]
print(f"{len(ids)} arXiv IDs in the census; {len(todo)} to fetch")


def _field(block, tag):
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", block, re.S)
    return re.sub(r"\s+", " ", m.group(1)).strip() if m else ""


for start in range(0, len(todo), BATCH):
    chunk = todo[start:start + BATCH]
    q = urllib.parse.urlencode({"id_list": ",".join(chunk), "max_results": len(chunk)})
    for attempt in range(6):
        try:
            xml = urllib.request.urlopen(
                urllib.request.Request(f"{API}?{q}", headers=UA), timeout=120).read().decode()
            break
        except Exception as e:                                   # 429/503 under load are routine
            wait = 20 * (attempt + 1)
            print(f"  retry {attempt + 1} in {wait}s: {e}")
            time.sleep(wait)
    else:
        sys.exit(f"arXiv API unreachable for batch starting at {start}; nothing written")

    for e in re.findall(r"<entry>(.*?)</entry>", xml, re.S):
        m = re.search(r"abs/(.+?)v?\d*$", _field(e, "id"))
        if not m or m.group(1) not in chunk:
            continue
        cache[m.group(1)] = {
            "arxiv": m.group(1),
            "title": _field(e, "title"),
            "journal": _field(e, "arxiv:journal_ref"),
            "doi": _field(e, "arxiv:doi"),
            "year": _field(e, "published")[:4],
        }
    print(f"  batch {start // BATCH + 1}: {len(cache)}/{len(ids)} resolved")
    time.sleep(GAP)

missing = [i for i in ids if i not in cache]
if missing:
    sys.exit(f"unresolved after fetching: {missing}")

with open(_p("data", "census_papers.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    for i in ids:
        w.writerow({k: cache[i][k] for k in COLS})

print(f"wrote data/census_papers.csv ({len(ids)} papers, "
      f"{sum(1 for i in ids if cache[i]['journal'])} with a journal reference)")
