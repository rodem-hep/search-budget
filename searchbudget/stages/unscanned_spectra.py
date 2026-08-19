from .. import io, paths
from ..core import combinatorial_budget as CB
from ..core.catalogue import canonical_order, models_by_spectrum
from ..core.scan_alphabet import WIDE_ARGS
from ..registry import stage
from .scaled_scan import MOTIVATED, canon

PAIR_AXIS = "m(tt)/m(jj)"

# The pair-produced axis resolved into the scan's own units, one leg spectrum per
# (mass group, category): label, leg composition, category content in the wide alphabet.
LEGS = [
    ("m(jj) pair-leg in 4j",   "jj", {"j": 4}),
    ("m(jj) pair-leg in 2t2j", "jj", {"j": 2, "t": 2}),
    ("m(tt) pair-leg in 2t2j", "tt", {"j": 2, "t": 2}),
    ("m(tt) pair-leg in 4t",   "tt", {"t": 4}),
]

# The one leg a published search covers: the paired-dijet resonance search scans the jj leg
# in the four-jet category (the census charges it to the inclusive m(jj) axis).
LEG_COVERED = {"m(jj) pair-leg in 4j": "Pair-produced dijet resonances (4 jets)"}

EXPECTED = 14


def _cat_content(cat):
    out = {}
    core = cat.split("_")[0]
    for i in range(0, len(core), 2):
        if int(core[i + 1]):
            out[core[i]] = int(core[i + 1])
    return out


@stage(
    name="unscanned-spectra",
    group="budget",
    summary="the model-motivated spectra no published ATLAS search scans, in the scan's units",
    outputs=["tables/unscanned_spectra.csv"],
    needs=["tables/census_budget.csv"],
)
def main(options=None):
    by_obs = models_by_spectrum()
    ranked = canonical_order(by_obs)
    published = {r["budget_axis"] for r in io.read_rows(paths.table("census_budget.csv"))
                 if r["budget_axis"] != "-" and float(r["n_s"]) > 0}

    s = CB.enumerate_scan(**WIDE_ARGS)
    fittable_comps = set(s.by_type)
    leg_fit = {label: any(r.n_s > 0 and canon(r.group) == grp and _cat_content(r.cat) == content
                          for r in s.rows)
               for label, grp, content in LEGS}

    rows_out = []
    for o in ranked:
        if o in published:
            continue
        if o == PAIR_AXIS:
            for label, _grp, _content in LEGS:
                covered = LEG_COVERED.get(label, "")
                rows_out.append([label, o, "scanned" if covered else "unscanned", covered,
                                 "yes" if leg_fit[label] else "no", len(by_obs[o])])
            continue
        comps = MOTIVATED.get(o)
        fit = bool(comps) and any(canon(c) in fittable_comps for c in comps)
        rows_out.append([o, o, "unscanned", "", "yes" if fit else "no", len(by_obs[o])])

    unscanned = [r for r in rows_out if r[2] == "unscanned"]
    fittable = [r for r in unscanned if r[4] == "yes"]
    if len(unscanned) != EXPECTED:
        raise SystemExit(f"{len(unscanned)} unscanned spectra, expected {EXPECTED}; "
                         "if the change is intended, update EXPECTED and the paper")

    io.write_rows(paths.table("unscanned_spectra.csv"),
                  ["spectrum", "budget_axis", "status", "covered_by", "fittable_in_scan",
                   "n_model_classes"], rows_out)

    print(f"model-motivated spectra with no published ATLAS scan, in the scan's "
          f"(category, mass-group) units: {len(unscanned)}, of which {len(fittable)} fittable")
    for r in rows_out:
        note = f"  [scanned: {r[3]}]" if r[2] == "scanned" else \
               ("" if r[4] == "yes" else "  (not fittable)")
        if r[2] == "unscanned" or note:
            print(f"  {r[0]:28s} {r[1]:14s}{note}")
    print("wrote results/tables/unscanned_spectra.csv")
