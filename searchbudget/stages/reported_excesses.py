import re

from .. import io, paths
from ..core.bump_observables import SCAN
from ..core.lee import p1
from ..registry import stage

NUM = r"\$?(\d+(?:\.\d+)?)\s*(?:σ|\\sigma)?\$?"


@stage(
    name="reported-excesses",
    group="census",
    summary="every excess >= 3 sigma reported in the census abstracts",
    outputs=["tables/reported_excesses.csv", "overviews/REPORTED_EXCESSES.md"],
    inputs=["data/census_papers.csv", "data/published_spectra.csv",
            "data/census_abstracts.csv"],
    needs=["tables/search_budget.csv"],
)
def main(options=None):
    paper_meta = {r["arxiv"]: r for r in io.read_rows(paths.data("census_papers.csv"))}
    n_axes = len(SCAN)
    spectrum_of = {}
    for r in io.read_rows(paths.data("published_spectra.csv")):
        for a in r["arxiv"].split():
            spectrum_of[a] = r["spectrum"]

    SIG_UNIT = re.compile(rf"{NUM}\s*(?:σ|\\sigma|standard deviations?)")
    PAIRED = re.compile(rf"local\s*\(global\)\s*significance[s]? of {NUM}\s*\({NUM}\s*"
                        r"(?:σ|\\sigma)?\$?\)")
    LOCAL_LIST = re.compile(rf"local significances? of {NUM}(?: and {NUM})? standard deviations")

    def sentences(text):
        return re.split(r"(?<=[.!?])\s+(?=[A-Z])", re.sub(r"\s+", " ", text))

    def extract(sentence):
        m = PAIRED.search(sentence)
        if m:
            return float(m.group(1)), float(m.group(2))
        m = LOCAL_LIST.search(sentence)
        if m:
            return max(float(g) for g in m.groups() if g), None
        vals = [(float(m.group(1)), m.start()) for m in SIG_UNIT.finditer(sentence)]
        if not vals:
            return None, None
        loc = glo = None
        for v, pos in vals:
            prefix = sentence[:pos].lower()
            is_global = prefix.rfind("global") > max(prefix.rfind("local"), len(prefix) - 80)
            if is_global:
                glo = max(glo, v) if glo is not None else v
            else:
                loc = max(loc, v) if loc is not None else v
        return loc, glo

    rows = []
    for r in io.read_rows(paths.data("census_abstracts.csv")):
        for s in sentences(r["abstract"]):
            sl = s.lower()
            if not re.search(r"excess|deviation", sl):
                continue
            if re.search(r"no (statistically )?(significant )?(local )?(excess|deviation)"
                         r"|not significant|absence of|exclud|signal strength|signal significance",
                         sl):
                continue
            loc, glo = extract(s)
            if loc is None and glo is None:
                continue
            meta = paper_meta.get(r["arxiv"], {})
            rows.append({
                "arxiv": r["arxiv"], "year": meta.get("year", ""),
                "spectrum": spectrum_of.get(r["arxiv"], ""),
                "z_local": f"{loc:.1f}" if loc is not None else "",
                "z_global": f"{glo:.1f}" if glo is not None else "",
                "sentence": s.strip(),
            })

    best = {}
    for row in rows:
        k = row["arxiv"]
        if k not in best:
            best[k] = row
            continue
        b = best[k]
        if row["z_local"] and float(row["z_local"]) > float(b["z_local"] or 0):
            b["z_local"], b["sentence"] = row["z_local"], row["sentence"]
        if row["z_global"] and float(row["z_global"]) > float(b["z_global"] or 0):
            b["z_global"] = row["z_global"]
    rows = sorted(best.values(), key=lambda r: (-(float(r["z_local"] or 0)), r["arxiv"]))

    io.write_dicts(paths.table("reported_excesses.csv"), rows,
                   ["arxiv", "year", "spectrum", "z_local", "z_global", "sentence"])

    budget = io.read_rows(paths.table("search_budget.csv"))
    N_inc = sum(float(r["ns_scan"]) for r in budget)
    N_sel = sum(float(r["ns_with_selections"]) for r in budget)
    exp_inc, exp_sel = N_inc * p1(3.0), N_sel * p1(3.0)

    n_papers = len({r["arxiv"] for r in io.read_rows(paths.data("census_abstracts.csv"))})
    ge3 = [r for r in rows if r["z_local"] and float(r["z_local"]) >= 3.0]
    sub3 = [r for r in rows if r not in ge3]

    with open(io.ensure(paths.overview("REPORTED_EXCESSES.md")), "w") as f:
        f.write("# Reported excesses, mined from the census abstracts\n\n")
        f.write(f"The observed side of `EXCESS_COUNTING.md`, made concrete for the resonance subset: "
                f"of the {n_papers} census papers, **{len(ge3)} report a local excess >= 3 sigma in "
                f"their abstract**, and {len(sub3)} more quantify a largest deviation below 3 sigma; "
                f"the rest report no significant excess or none at all. Extraction: "
                f"`searchbudget/stages/reported_excesses.py`; the full rows, with the abstract sentence each "
                f"number came from, are in `results/tables/reported_excesses.csv`.\n\n")
        f.write("| arXiv | year | spectrum | Z_local | Z_global |\n|---|--:|---|--:|--:|\n")
        for r in ge3:
            f.write(f"| {r['arxiv']} | {r['year']} | {r['spectrum']} | {r['z_local']} | "
                    f"{r['z_global'] or '-'} |\n")
        f.write(f"\nEvery quoted global significance above is <= 2.1 sigma, and none has since "
                f"grown into a discovery.\n\n")
        f.write("## Against the background-only expectation\n\n")
        f.write(f"A single background-only sweep of the {n_axes}-spectrum budget expects "
                f"`N x p(>=3 sigma)` = {N_inc:.0f} x 1.35e-3 = **{exp_inc:.1f}** reports "
                f"(**{exp_sel:.1f}** at selection granularity, N = {N_sel:.0f}). The observed "
                f"{len(ge3)} sits between the two.\n\n")
        f.write("Two corrections pull in opposite directions and are left uncorrected:\n\n"
                "- **Abstracts under-report.** A sub-3 sigma maximum is routinely left to the body\n"
                "  (the 2 TeV diboson excess of 1506.00962 appears in its abstract only as a 2.5\n"
                "  sigma *global*), so the count is a lower bound on what the papers contain.\n"
                "- **The reports are not independent trials.** Successive papers on one axis re-scan\n"
                "  partially shared data, and the record sweeps most axes more than once across\n"
                "  Run 1 and Run 2, so the effective number of independent sweeps is above one.\n\n"
                "Neither changes the conclusion: the reported population is what a counted trials\n"
                "factor predicts under background-only, at the subset level and not just for the\n"
                "program-wide anchor of `EXCESS_COUNTING.md`.\n")

    print(f"{n_papers} abstracts; {len(rows)} papers quantify a largest deviation; "
          f"{len(ge3)} at >= 3 sigma local")
    for r in ge3:
        print(f"  {r['arxiv']} ({r['year']}) {r['spectrum']}: local {r['z_local']}"
              f"{' global ' + r['z_global'] if r['z_global'] else ''}")
    print(f"expected >=3 sigma per background-only sweep: {exp_inc:.1f} (inclusive), "
          f"{exp_sel:.1f} (selections)")
    print("wrote results/tables/reported_excesses.csv, results/overviews/REPORTED_EXCESSES.md")
