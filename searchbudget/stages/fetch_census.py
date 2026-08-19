import re

from .. import arxiv, io, paths
from ..core.public_obs_map import LITERATURE
from ..registry import stage

META_COLS = ["arxiv", "title", "journal", "doi", "year"]
ABSTRACT_COLS = ["arxiv", "abstract"]
MODEL_COLS = ["arxiv", "title", "authors", "journal", "doi", "year"]


def _ids():
    return arxiv.census_ids(paths.data("published_spectra.csv"))


@stage(
    name="fetch-census-meta",
    group="fetch",
    summary="refresh the bibliographic details of the census papers from the arXiv API",
    outputs=["data/census_papers.csv"],
    inputs=["data/published_spectra.csv"],
    network=True,
)
def fetch_meta(options=None):
    out = paths.data("census_papers.csv")
    ids = _ids()

    def extract(aid, entry):
        return {"arxiv": aid,
                "title": arxiv.field(entry, "title"),
                "journal": arxiv.field(entry, "arxiv:journal_ref"),
                "doi": arxiv.field(entry, "arxiv:doi"),
                "year": arxiv.field(entry, "published")[:4]}

    cache = arxiv.harvest(ids, arxiv.cached(out), extract)
    io.write_dicts(out, [{k: cache[i][k] for k in META_COLS} for i in ids], META_COLS)
    print(f"wrote data/census_papers.csv ({len(ids)} papers, "
          f"{sum(1 for i in ids if cache[i]['journal'])} with a journal reference)")


@stage(
    name="fetch-census-abstracts",
    group="fetch",
    summary="refresh the abstracts of the census papers from the arXiv API",
    outputs=["data/census_abstracts.csv"],
    inputs=["data/published_spectra.csv"],
    network=True,
)
def fetch_abstracts(options=None):
    out = paths.data("census_abstracts.csv")
    ids = _ids()

    def extract(aid, entry):
        return {"arxiv": aid, "abstract": arxiv.field(entry, "summary")}

    cache = arxiv.harvest(ids, arxiv.cached(out), extract)
    io.write_dicts(out, [{k: cache[i][k] for k in ABSTRACT_COLS} for i in ids], ABSTRACT_COLS)
    print(f"wrote data/census_abstracts.csv ({len(ids)} abstracts)")


@stage(
    name="fetch-model-meta",
    group="fetch",
    summary="refresh the bibliographic details of the literature-sourced model papers",
    outputs=["data/model_papers.csv"],
    network=True,
)
def fetch_model_meta(options=None):
    out = paths.data("model_papers.csv")
    ids = []
    for refs in LITERATURE.values():
        for a in refs:
            if a not in ids:
                ids.append(a)

    def extract(aid, entry):
        names = re.findall(r"<name>(.*?)</name>", entry, re.S)
        return {"arxiv": aid,
                "title": arxiv.field(entry, "title"),
                "authors": "; ".join(re.sub(r"\s+", " ", n).strip() for n in names),
                "journal": arxiv.field(entry, "arxiv:journal_ref"),
                "doi": arxiv.field(entry, "arxiv:doi"),
                "year": arxiv.field(entry, "published")[:4]}

    cache = arxiv.harvest(ids, arxiv.cached(out), extract)
    io.write_dicts(out, [{k: cache[i][k] for k in MODEL_COLS} for i in ids], MODEL_COLS)
    print(f"wrote data/model_papers.csv ({len(ids)} papers behind {len(LITERATURE)} "
          f"literature-sourced model classes)")
