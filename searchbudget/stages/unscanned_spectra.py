import collections

from .. import io, paths
from ..core import combinatorial_budget as CB
from ..core.catalogue import canonical_order, models_by_spectrum
from ..core.public_obs_map import LITERATURE
from ..viz.labels import mathify, textsafe
from ..core import yield_model as YM
from ..core.scan_alphabet import DATASETS, WIDE_ARGS
from ..registry import stage
from .scaled_scan import MOTIVATED, canon, cat_content, predicted_axes

PAIR_AXIS = "m(tt)/m(jj)"


def _model_tex(m):
    """A model class as the appendix tables name it: escaped, and cited where the class comes
    from the literature sweep rather than the FeynRules/UFO database."""
    cite = "~\\cite{" + ",".join(f"lit:{a}" for a in LITERATURE[m]) + "}" if m in LITERATURE else ""
    for a, b in (("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_")):
        m = m.replace(a, b)
    return textsafe(m) + cite

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

# The headline is the count that can actually be fitted in the combinatorial scan: a spectrum the
# scan cannot hold a fit in is not one of its 4,438. The total is kept alongside it.
EXPECTED_FITTABLE = 19
EXPECTED_TOTAL = 21


@stage(
    name="unscanned-spectra",
    group="budget",
    summary="the model-motivated spectra no published ATLAS search scans, in the scan's units",
    outputs=["tables/unscanned_spectra.csv", "tables/unscanned_scan_units.csv",
             "tables/new_spectra.csv", "tex/new_spectra.tex"],
    needs=["tables/census_budget.csv"],
)
def main(options=None):
    by_obs = models_by_spectrum()
    ranked = canonical_order(by_obs)
    published = {r["budget_axis"] for r in io.read_rows(paths.table("census_budget.csv"))
                 if r["budget_axis"] != "-" and float(r["n_s"]) > 0}

    s_scan = CB.enumerate_scan(**WIDE_ARGS)
    fittable_comps = set(s_scan.by_type)
    leg_fit = {label: any(r.n_s > 0 and canon(r.group) == grp and cat_content(r.cat) == content
                          for r in s_scan.rows)
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
    if (len(fittable), len(unscanned)) != (EXPECTED_FITTABLE, EXPECTED_TOTAL):
        raise SystemExit(f"{len(fittable)} fittable of {len(unscanned)} unscanned spectra, "
                         f"expected {EXPECTED_FITTABLE} of {EXPECTED_TOTAL}; if the change is "
                         "intended, update the two constants and the paper")

    io.write_rows(paths.table("unscanned_spectra.csv"),
                  ["spectrum", "budget_axis", "status", "covered_by", "fittable_in_scan",
                   "n_model_classes"], rows_out)

    # The same question in the scan's own units, where a spectrum is one axis in one category and
    # counts only when some model predicts that final state with the mass on its resonant
    # sub-system (`scaled_scan.predicted_axes`). Priced on the same dataset as the headline scan,
    # Run 2+3, and counted the way that scan counts, `split` spectra per row.
    _n_ref = YM.N_REF
    try:
        YM.N_REF = _n_ref * DATASETS[-1][1]
        rows_scan = [r for r in CB.enumerate_scan(**WIDE_ARGS).rows if r.n_s > 0]
    finally:
        YM.N_REF = _n_ref

    n_fit = n_mot = n_mot_unscanned = 0
    hit_axes, hit_axes_unscanned, hit_states = set(), set(), set()
    new_rows, new_looks = [], 0.0
    for r in rows_scan:
        n_fit += r.split
        axes = predicted_axes(r)
        if not axes:
            continue
        n_mot += r.split
        hit_axes |= axes
        hit_states.add((canon(r.group), r.cat))
        fresh = {a for a in axes if a not in published}
        if len(fresh) == len(axes):          # no published search on any axis that predicts it
            n_mot_unscanned += r.split
            hit_axes_unscanned |= fresh
            new_looks += r.n_s
            new_rows.append([canon(r.group), r.cat.split("_")[0], r.split, f"{r.n_s:.1f}",
                             "; ".join(sorted(fresh)),
                             "; ".join(sorted({m for a in fresh for m in by_obs[a]}))])

    io.write_rows(paths.table("new_spectra.csv"),
                  ["mass_group", "category", "histograms", "n_s", "budget_axes", "model_classes"],
                  sorted(new_rows, key=lambda x: (-float(x[3]), x[0], x[1])))

    # the same list for the paper, one row per axis
    TEX = {"e": "e", "m": r"{\mu}", "T": r"{\tau}", "g": r"{\gamma}", "j": "j", "b": "b",
           "t": "t", "V": "V", "H": "H", "Z": "Z"}

    def _tex(counts):
        return "$" + "".join(("" if n == 1 else str(n)) + TEX[k] for k, n in
                             sorted(counts.items(), key=lambda kv: "emTgjbtVHZ".index(kv[0]))) + "$"

    def _fs(cat):
        return _tex(cat_content(cat))

    def _grp(g):
        return "$m(" + _tex(collections.Counter(g)).strip("$") + ")$"

    by_axis = collections.defaultdict(lambda: [0, 0.0, set(), set(), set()])
    for grp, cat, split, ns, axes, models in new_rows:
        e = by_axis[axes]
        e[0] += split
        e[1] += float(ns)
        e[2].add(_fs(cat))
        e[3].add(_grp(grp))
        e[4] |= set(models.split("; "))
    RAG = ">{\\raggedright\\arraybackslash}"
    COL = (f"@{{}}l{RAG}p{{0.14\\textwidth}}{RAG}p{{0.10\\textwidth}}rr"
           f"{RAG}p{{0.30\\textwidth}}@{{}}")
    with open(io.ensure(paths.tex("new_spectra.tex")), "w") as f:
        f.write("% Generated by searchbudget/stages/unscanned_spectra.py. "
                "Do not edit: regenerate instead.\n")
        f.write(f"\\begin{{tabular}}{{{COL}}}\n\\toprule\n")
        f.write("axis & mass & final state(s) & spectra & looks & model classes that predict it"
                " \\\\\n\\midrule\n")
        for axis, (nh, looks, cats, masses, models) in sorted(by_axis.items(),
                                                              key=lambda kv: -kv[1][1]):
            f.write(f"{mathify(axis)} & {', '.join(sorted(masses))} & {', '.join(sorted(cats))} "
                    f"& {nh} & {looks:.0f} & {', '.join(_model_tex(m) for m in sorted(models))} \\\\\n")
        f.write("\\midrule\n"
                f"total & & & {n_mot_unscanned} & {new_looks:.0f} & \\\\\n"
                "\\bottomrule\n\\end{tabular}\n")

    io.write_rows(paths.table("unscanned_scan_units.csv"),
                  ["quantity", "spectra", "share_of_fittable", "axes"],
                  [["new spectra: predicted, and on an axis ATLAS has not scanned",
                    n_mot_unscanned, f"{n_mot_unscanned / n_fit:.4f}", len(hit_axes_unscanned)],
                   ["... of the predicted spectra, the category being one a model produces and "
                    "the mass its resonant sub-system", n_mot, f"{n_mot / n_fit:.4f}",
                    len(hit_axes)],
                   ["... of the fittable spectra of the combinatorial scan", n_fit, "1.000", ""],
                   ["no model predicts that final state", n_fit - n_mot,
                    f"{1 - n_mot / n_fit:.4f}", ""]])
    print(f"NEW SPECTRA the scan would add: {n_mot_unscanned} of its {n_fit} fittable spectra "
          f"({100 * n_mot_unscanned / n_fit:.2f}%) are a final state a public model predicts, with "
          f"the mass on its resonant sub-system, and sit on an axis no published ATLAS search "
          f"scans; they carry {new_looks:,.0f} looks over {len(hit_axes_unscanned)} axes: "
          f"{', '.join(sorted(hit_axes_unscanned))}")
    print(f"  for scale, {n_mot} of {n_fit} ({100 * n_mot / n_fit:.1f}%) are predicted at all, "
          f"over {len(hit_axes)} axes and {len(hit_states)} (final state, mass) pairs; the rest of "
          f"the scan is territory no model points at")
    print("wrote results/tables/unscanned_scan_units.csv, results/tables/new_spectra.csv")

    print(f"model-motivated spectra with no published ATLAS scan that the combinatorial scan can "
          f"fit, in its (category, mass-group) units: {len(fittable)}; "
          f"{len(unscanned)} counting the {len(unscanned) - len(fittable)} it cannot fit")
    for r in rows_out:
        note = f"  [scanned: {r[3]}]" if r[2] == "scanned" else \
               ("" if r[4] == "yes" else "  (not fittable)")
        if r[2] == "unscanned" or note:
            print(f"  {r[0]:28s} {r[1]:14s}{note}")
    print("wrote results/tables/unscanned_spectra.csv")
