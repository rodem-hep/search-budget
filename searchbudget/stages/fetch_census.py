from .. import arxiv, io, paths
from ..registry import stage

META_COLS = ["arxiv", "title", "journal", "doi", "year"]
ABSTRACT_COLS = ["arxiv", "abstract"]


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
