import re

from .. import io, paths
from ..registry import stage

# A charged (search, axis) pair counts as model independent when any of the
# search's papers states a model-independent result in its own abstract: a
# generic Gaussian-shape limit, a model-agnostic or anomaly-detection scan, or
# a cross-section limit the search itself declares model independent.
EVIDENCE = re.compile(
    r"model[- ]independent|model[- ]agnostic|anomaly detection|"
    r"generic gaussian|gaussian[- ]shaped?|"
    r"gaussian[- ]?(?:signal|resonance|contribution|model|line)", re.I)


@stage(
    name="model-independence",
    group="census",
    summary="which charged (search, axis) pairs carry a model-independent result",
    outputs=["tables/model_independence.csv"],
    inputs=["data/published_spectra.csv", "data/census_abstracts.csv"],
    needs=["tables/census_budget.csv"],
)
def main(options=None):
    abstracts = {r["arxiv"]: r["abstract"]
                 for r in io.read_rows(paths.data("census_abstracts.csv"))}
    spectra = {r["spectrum"]: r for r in io.read_rows(paths.data("published_spectra.csv"))}
    charged = [r for r in io.read_rows(paths.table("census_budget.csv"))
               if float(r["n_s"]) > 0]

    evidence = {}
    for name, row in spectra.items():
        hits = []
        for a in row["arxiv"].split():
            m = EVIDENCE.search(abstracts.get(a, ""))
            if m:
                hits.append((a, m.group(0)))
        evidence[name] = hits

    io.write_rows(
        paths.table("model_independence.csv"),
        ["spectrum", "budget_axis", "model_independent", "evidence_arxiv", "evidence_phrase"],
        [[r["spectrum"], r["budget_axis"], "yes" if evidence[r["spectrum"]] else "no",
          " ".join(a for a, _ in evidence[r["spectrum"]]),
          "; ".join(p for _, p in evidence[r["spectrum"]])]
         for r in charged])

    mi = [r for r in charged if evidence[r["spectrum"]]]
    io.note(f"of the {len(charged)} charged (search, axis) pairs, {len(mi)} from "
            f"{len({r['spectrum'] for r in mi})} of the {len({r['spectrum'] for r in charged})} "
            f"searches carry a model-independent result; the other {len(charged) - len(mi)} "
            f"are interpreted only in benchmark models")
