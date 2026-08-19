import csv
import os
import re
import sys
import time
import urllib.parse
import urllib.request

API = "https://export.arxiv.org/api/query"
UA = {"User-Agent": "search-budget/1.0 (https://github.com/rodem-hep/search-budget)"}
BATCH, GAP, RETRIES = 25, 6, 6


def census_ids(path):
    ids = []
    for row in csv.DictReader(open(path)):
        for a in row["arxiv"].split():
            if a not in ids:
                ids.append(a)
    return ids


def cached(path):
    if not os.path.exists(path):
        return {}
    return {r["arxiv"]: r for r in csv.DictReader(open(path))}


def field(block, tag):
    m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", block, re.S)
    return re.sub(r"\s+", " ", m.group(1)).strip() if m else ""


def _query(chunk):
    q = urllib.parse.urlencode({"id_list": ",".join(chunk), "max_results": len(chunk)})
    for attempt in range(RETRIES):
        try:
            request = urllib.request.Request(f"{API}?{q}", headers=UA)
            return urllib.request.urlopen(request, timeout=120).read().decode()
        except Exception as exc:
            wait = 20 * (attempt + 1)
            print(f"  retry {attempt + 1} in {wait}s: {exc}")
            time.sleep(wait)
    return None


def harvest(ids, cache, extract, batch=BATCH, gap=GAP):
    todo = [i for i in ids if i not in cache]
    print(f"{len(ids)} arXiv IDs in the census; {len(todo)} to fetch")
    for start in range(0, len(todo), batch):
        chunk = todo[start:start + batch]
        xml = _query(chunk)
        if xml is None:
            sys.exit(f"arXiv API unreachable for batch starting at {start}; nothing written")
        for entry in re.findall(r"<entry>(.*?)</entry>", xml, re.S):
            m = re.search(r"abs/(.+?)v?\d*$", field(entry, "id"))
            if not m or m.group(1) not in chunk:
                continue
            cache[m.group(1)] = extract(m.group(1), entry)
        print(f"  batch {start // batch + 1}: {len(cache)}/{len(ids)} resolved")
        time.sleep(gap)
    missing = [i for i in ids if i not in cache]
    if missing:
        sys.exit(f"unresolved after fetching: {missing}")
    return cache
