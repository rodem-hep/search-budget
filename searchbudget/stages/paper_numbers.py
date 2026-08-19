import math

from .. import io, paths
from ..core import yield_model as YM
from ..core.bump_observables import res
from ..core.lee import merge_segments, p1, z_global as zg, z_local as z5
from ..registry import stage


@stage(
    name="paper-numbers",
    group="paper",
    summary="every number the paper quotes, generated into the README block",
    outputs=["README.md"],
    needs=["tables/search_budget.csv", "tables/search_budget_selections.csv",
           "tables/published_census.csv", "tables/census_budget.csv",
           "tables/scan_summary.csv", "tables/lens_scan.csv", "tables/priority_scan.csv",
           "tables/model_independence.csv",
           "tables/budget_uncertainty.csv", "tables/two_body_matrix.csv",
           "tables/ab_split_scan.csv", "tables/ab_split_toys.csv",
           "tables/ab_guard_toys.csv", "tables/estimator_defects.csv",
           "tables/selection_rules.csv"],
)
def main(options=None):
    START = "<!-- paper-numbers:start -->"
    END = "<!-- paper-numbers:end -->"


    def load(name):
        return io.read_rows(paths.table(name))


    budget = load("search_budget.csv")
    selections = load("search_budget_selections.csv")
    census_families = [r for r in load("published_census.csv") if r["family"] != "TOTAL"]
    census = load("census_budget.csv")
    scan = {r["dataset"]: r for r in load("scan_summary.csv")}
    lenses = [r for r in load("lens_scan.csv") if r["dataset"].startswith("Run 2+3")]
    unc = {(r["source"], r["direction"]): r for r in load("budget_uncertainty.csv")}
    twobody = load("two_body_matrix.csv")
    ab = {r["basis"]: r for r in load("ab_split_scan.csv")}
    defects = {r["quantity"]: r["value"] for r in load("estimator_defects.csv")}
    rules = load("selection_rules.csv")
    toys = load("ab_split_toys.csv")
    guard = load("ab_guard_toys.csv")
    mi_rows = load("model_independence.csv")

    ns = {r["observable"]: float(r["ns_scan"]) for r in budget}
    N_model = sum(ns.values())
    N_sel = sum(float(r["ns_with_selections"]) for r in selections)
    n_sel_channels = sum(int(r["n_event_selections"]) for r in selections)
    sharp = sorted(ns.values(), reverse=True)

    pairs = len(census)
    priced = [r for r in census if float(r["n_s"]) > 0]
    N_census = sum(float(r["n_s"]) for r in priced)
    entries = {r["spectrum"] for r in census}
    carry = {r["spectrum"] for r in priced}
    fixed = sum(1 for r in priced if abs(float(r["n_s"]) - 1) < 1e-9)
    axes_covered = {r["budget_axis"] for r in priced if r["budget_axis"] != "-"}
    off_axis = [r for r in priced if r["budget_axis"] == "-"]
    unpriced = [r for r in census if float(r["n_s"]) == 0]
    papers = sum(int(r["papers"]) for r in census_families)
    stale = sum(int(r["stale"]) for r in census_families)
    run3 = sum(int(r["run3"]) for r in census_families)


    def union_N():
        seg = {}
        for r in priced:
            a = r["budget_axis"]
            if a == "-" or not r["window_GeV"] or r["window_GeV"] == "fixed":
                continue
            for s in r["window_GeV"].split("+"):
                lo, hi = (float(x) for x in s.split("-"))
                seg.setdefault(a, []).append((lo, hi))
        tot = 0.0
        for a, segs in seg.items():
            tot += sum(math.log(hi / lo) / res(a) for lo, hi in merge_segments(segs))
        return tot + sum(float(r["n_s"]) for r in off_axis)


    N_union = union_N()

    r23 = scan["Run 2+3, ~400 fb-1"]
    N_scan = float(r23["N_trials"])
    N_lens = float(r23["lensed_N_trials"])
    reach = 9407
    views = sum(int(x["views_kept"]) for x in lenses)
    thin = sum(int(x["views_too_thin"]) for x in lenses)
    reach = views + thin
    lens_looks = sum(float(x["looks_kept"]) for x in lenses)

    gaps = [r for r in twobody if r["status"] != "scanned"]
    gap_looks = sum(float(r["ns"]) for r in gaps)
    gap_lo = min(gaps, key=lambda r: float(r["ns"]))
    gap_hi = max(gaps, key=lambda r: float(r["ns"]))
    ee = next(r for r in twobody if r["object_1"] == "e" and r["object_2"] == "e")

    lensed_ab = ab["lensed scan"]
    at5 = {r["estimator"]: r for r in rules if r["mu"] == "5"}


    def env(source, key):
        r = unc.get((source, "envelope")) or unc[(source, "total")]
        v = r[key]
        if not v:
            return "--"
        return "/".join(f"{float(x):+.2f}" for x in v.split("/"))


    def rows(*tuples):
        out = ["| quantity | value | reproduced by |", "|---|---|---|"]
        out += [f"| {q} | {v} | {s} |" for q, v, s in tuples]
        return "\n".join(out)


    blocks = []

    blocks.append("### The ladder (Table 1)\n\n" + rows(
        ("spectra the public model space populates", f"**{len(budget)}**", "`search_budget.csv`"),
        ("the model space, published windows",
         f"`N = {N_model:,.0f}` → `Z_local = {z5(N_model):.2f}` "
         f"{env('total', 'dZ_model_space')}", "`search_budget.csv`, `budget_uncertainty.csv`"),
        ("the same at published event-selection granularity",
         f"{n_sel_channels} channels, `N = {N_sel:,.0f}` → `{z5(N_sel):.2f}`",
         "`search_budget_selections.csv`"),
        ("the published ATLAS program",
         f"{len(priced)} charged (search, axis) pairs, "
         f"`N = {N_census:,.0f}` → `{z5(N_census):.2f}`", "`census_budget.csv`"),
        ("a fully combinatorial scan, Run 2+3",
         f"{int(r23['fittable_spectra']):,} fittable spectra, `N = {N_scan:,.0f}` → "
         f"`{z5(N_scan):.2f}` {env('total', 'dZ_combinatorial_scan')}", "`scan_summary.csv`"),
        ("... with one event-level selection at a time",
         f"{int(r23['lensed_histograms']):,} histograms, `N = {N_lens:,.0f}` → "
         f"`{z5(N_lens):.2f}` {env('total', 'dZ_scan_with_lenses')}", "`scan_summary.csv`"),
        ("published record → lensed scan",
         f"a factor {N_lens/N_census:.0f} in `N`, worth {z5(N_lens)-z5(N_census):+.2f}σ", "derived"),
        ("the whole ladder", f"a factor {N_lens/N_model:.0f} in `N`, worth "
         f"{z5(N_lens)-z5(N_model):+.2f}σ", "derived"),
        ("model space → combinatorial scan",
         f"{z5(N_scan)-z5(N_model):+.2f}σ, {env('total', 'dZ_difference')} as a difference",
         "`budget_uncertainty.csv`"),
        ("a second experiment scanning the same axes",
         f"doubles `N`, worth {z5(2*N_model)-z5(N_model):+.2f}σ", "derived"),
    ))

    blocks.append("### Counting the looks (Section 2)\n\n" + rows(
        ("a histogram is fittable if it holds",
         f"≥ {YM.MIN_EVENTS:.0f} events and ≥ {YM.MIN_BINS} elements of ≥ 1 event", "`yield_model.py`"),
        ("the declared background",
         f"`n(m) = 10^{math.log10(YM.N_REF):.0f}·W·(m/{YM.M_REF:.0f} GeV)^(1−{YM.P:.0f})"
         f"·(r/{YM.R_REF})` per element", "`yield_model.py`"),
        ("the dataset", "140 fb⁻¹ at 13 TeV plus ~2× that at 13.6 TeV, i.e. ×3 the anchor and "
         f"×{float(r23['mass_reach_scale']):.2f} in mass reach", "`scan_summary.csv`"),
        ("what one look *is*, the leading uncertainty",
         f"`N` from 0.5·N to N·Z/√(2π): {env('the definition of one look', 'dZ_model_space')} on the "
         f"model space, {env('the definition of one look', 'dZ_combinatorial_scan')} on the scan, "
         f"{env('the definition of one look', 'dZ_difference')} on the difference",
         "`budget_uncertainty.csv`"),
        ("the resolution scale", f"every `r` ×2 either way: "
         f"{env('mass resolution, scale', 'dZ_model_space')} / "
         f"{env('mass resolution, scale', 'dZ_combinatorial_scan')} / "
         f"{env('mass resolution, scale', 'dZ_difference')}", "`budget_uncertainty.csv`"),
        ("... drawn per channel instead of in common",
         f"{env('mass resolution, per channel', 'dZ_model_space')}", "`budget_uncertainty.csv`"),
        ("Rice up-crossings against the element count",
         f"a factor `Z/√(2π)` = {z5(N_model)/math.sqrt(2*math.pi):.1f} at this bar", "derived"),
    ))

    blocks.append("### The model-space budget (Section 3)\n\n" + rows(
        ("cheapest and dearest spectrum",
         f"{min(ns.values()):.0f} looks (`{min(ns, key=ns.get)}`) to {max(ns.values()):.0f} "
         f"(`{max(ns, key=ns.get)}`)", "`search_budget.csv`"),
        ("what the five sharpest axes carry", f"{sum(sharp[:5])/N_model:.0%} of `N`",
         "`search_budget.csv`"),
        ("dropping the largest contributor",
         f"`Z_local` {z5(N_model):.2f} → {z5(N_model - max(ns.values())):.2f}", "`search_budget.csv`"),
        ("event selections instead of inclusive spectra",
         f"{n_sel_channels} channels, `{z5(N_sel):.2f}`", "`search_budget_selections.csv`"),
        ("the spectra the most model classes point at",
         ", ".join(f"`{r['observable']}` ({r['n_models_public']})" for r in
                   sorted(budget, key=lambda r: -int(r["n_models_public"]))[:3]),
         "`model_spectrum_map.csv`"),
    ))

    exp3, exp5 = N_census * p1(3.0), N_census * p1(5.0)
    blocks.append("### The published ATLAS program (Section 4)\n\n" + rows(
        ("the record", f"**{papers}** papers (2010–2026) in **{len(entries)}** catalogued searches over "
         f"{len(axes_covered) + len(off_axis) + len(unpriced)} bump observables",
         "`published_census.csv`, `census_budget.csv`"),
        ("(search, axis) pairs", f"{pairs}, of which **{len(priced)}** carry a chargeable range, over "
         f"{len(carry)} of the {len(entries)} searches", "`census_budget.csv`"),
        ("pairs charged a single look", f"{fixed} (a fixed mass, not a scan)", "`census_budget.csv`"),
        ("pairs with nothing chargeable", f"{len(unpriced)}", "`census_budget.csv`"),
        ("priced entry by entry", f"`N = {N_census:,.0f}` → `Z_local = {z5(N_census):.2f}`",
         "`census_budget.csv`"),
        ("priced once per axis, over the union of its ranges",
         f"`N = {N_union:,.0f}` → `{z5(N_union):.2f}`, within "
         f"{abs(N_union-N_model)/N_model:.1%} of the model space", "`census_budget.csv`"),
        ("recency", f"{stale} of the {len(entries)} have no paper since 2019, {run3} carry a "
         "published Run-3 result", "`published_census.csv`"),
        ("model independence",
         f"**{sum(1 for r in mi_rows if r['model_independent'] == 'yes')}** of the "
         f"{len(mi_rows)} charged pairs carry a model-independent result, from "
         f"{len({r['spectrum'] for r in mi_rows if r['model_independent'] == 'yes'})} of the "
         f"{len(carry)} searches", "`model_independence.csv`"),
        ("excesses a background-only sweep expects",
         f"{exp3:.0f} local ≥3σ and {exp5:.4f} spurious ≥5σ over the census",
         "`REPORTED_EXCESSES.md`, `EXCESS_COUNTING.md`"),
        ("what a local 5σ is worth globally",
         f"{zg(N_census):.1f}σ over the published record, {zg(N_scan):.1f}σ over the combinatorial "
         f"scan, nothing at all over the lensed scan", "derived"),
    ))

    blocks.append("### The combinatorial scan (Section 5)\n\n" + rows(
        ("the alphabet", f"ten object types, masses of ≤4 objects: "
         f"{int(r23['combinations']):,} combinations in {int(r23['categories']):,} categories",
         "`scan_summary.csv`"),
        ("what statistics leaves",
         f"**{int(r23['fittable_spectra']):,}** fittable, `N = {N_scan:,.0f}` → `{z5(N_scan):.2f}`",
         "`scan_summary.csv`"),
        ("two-body share of the survivors",
         f"{int(r23['k2_spectra'])/int(r23['fittable_spectra']):.0%}", "`scan_summary.csv`"),
        ("tier (i), every motivated composition once",
         f"{int(r23['tier0_spectra'])} spectra, `N = {float(r23['tier0_looks']):,.0f}` → "
         f"`{r23['tier0_Z_local']}`", "`scan_summary.csv`"),
        ("tier (ii), those compositions elsewhere",
         f"{int(r23['tier1_spectra']):,} spectra, `N = {float(r23['tier1_looks']):,.0f}`",
         "`scan_summary.csv`"),
        ("tier (iii), nothing motivates it",
         f"{int(r23['tier2_spectra']):,} spectra, `N = {float(r23['tier2_looks']):,.0f}`",
         "`scan_summary.csv`"),
        ("the share theory motivates",
         f"{int(r23['motivated_spectra']):,} of {int(r23['fittable_spectra']):,} spectra "
         f"({float(r23['motivated_frac']):.0%}), over {int(r23['compositions_motivated'])} of the "
         f"{int(r23['compositions_fittable'])} fittable compositions",
         "`scan_summary.csv`, per composition in `priority_scan.csv`"),
        ("selection lenses", f"{reach:,} conceivable views, {thin:,} ruled out by statistics, "
         f"{views:,} survive carrying {lens_looks:,.0f} looks", "`lens_scan.csv`"),
        ("the lensed scan", f"{int(r23['lensed_histograms']):,} histograms "
         f"({views/int(r23['fittable_spectra']):.2f} extra per spectrum), `N = {N_lens:,.0f}` → "
         f"`{z5(N_lens):.2f}`, so a lens costs {z5(N_lens)-z5(N_scan):+.2f}σ", "`scan_summary.csv`"),
        ("where the relation loses meaning",
         f"`25 − 2 ln N < 0` at `N = {math.exp(12.5)/1e5:.1f}·10^5`; the lensed scan is past it",
         "derived"),
        ("the yield anchor ×0.01 and ×100",
         f"{int(r23['anchor_x0.01_spectra']):,} to {int(r23['anchor_x100_spectra']):,} spectra, "
         f"`Z_local` {r23['anchor_x0.01_Z']} to {r23['anchor_x100_Z']}", "`scan_summary.csv`"),
        ("the same requirement applied to published axes",
         f"would keep {r23['published_axes_gated']} of {r23['published_axes_total']} axes and "
         f"{float(r23['published_looks_gated']):,.0f} of {float(r23['published_looks_total']):,.0f} "
         "looks, so exempting them is conservative", "`scan_summary.csv`"),
        ("the two-body grid", f"{len(gaps)} unscanned pairs cost {gap_looks:,.0f} looks, from "
         f"{float(gap_lo['ns']):.0f} (`{gap_lo['object_1']}{gap_lo['object_2']}`) to "
         f"{float(gap_hi['ns']):.0f} (`{gap_hi['object_1']}{gap_hi['object_2']}`); closing all of them "
         f"takes `N` {N_model:,.0f} → {N_model+gap_looks:,.0f} and `Z_local` {z5(N_model):.2f} → "
         f"{z5(N_model+gap_looks):.2f}", "`two_body_matrix.csv`"),
        ("the costliest compositions", ", ".join(
            f"`{c.split(':')[0]}` ({float(c.split(':')[1]):,.0f})"
            for c in r23["costliest_compositions"].split("; ")[:4]),
         "`scan_summary.csv`, per composition in `priority_scan.csv`"),
        ("the most heavily scanned pair",
         f"`{ee['object_1']}{ee['object_2']}`, {float(ee['ns']):.0f} looks over {ee['n_axes']} axes",
         "`two_body_matrix.csv`"),
    ))

    lens_tab = ["| lens | applies when | efficiency | spectra it could reach | views kept | looks |",
                "|---|---|--:|--:|--:|--:|"]
    for x in lenses:
        lens_tab.append(f"| {x['lens']} | {x['requirement']}"
                        + (f", window capped at {x['window_cap_GeV']} GeV" if x['window_cap_GeV']
                           else "")
                        + f" | {x['efficiency']} | {int(x['spectra_reached']):,} | "
                          f"{int(x['views_kept']):,} | {float(x['looks_kept']):,.0f} |")
    blocks.append("### The four selection lenses, Run 2+3 (Section 5)\n\nOne at a time, never a "
                  "product of two (`lens_scan.csv`):\n\n" + "\n".join(lens_tab))

    blocks.append("### Two-stage unblinding (Section 6)\n\n" + rows(
        ("the basis it is priced on", f"the lensed scan, `N = {float(lensed_ab['N_trials']):,.0f}`, "
         f"single stage exactly corrected `{lensed_ab['Z_single_stage']}`", "`ab_split_scan.csv`"),
        ("the optimised split", f"reach `{lensed_ab['reach_w3']}` at `w = 3` "
         f"({float(lensed_ab['reach_w3'])-float(lensed_ab['Z_single_stage']):+.2f}σ) and "
         f"`{lensed_ab['reach_w1']}` with pinned windows "
         f"({float(lensed_ab['reach_w1'])-float(lensed_ab['Z_single_stage']):+.2f}σ)",
         "`ab_split_scan.csv`"),
        ("the naive 50/50 split", f"reach `{float(lensed_ab['reach_5050_zcut3']):.2f}` "
         f"({float(lensed_ab['cost_5050']):+.2f}σ)", "`ab_split_scan.csv`"),
        ("the pre-registered list at `Z_cut = 3`",
         f"{int(lensed_ab['k_zcut3'])} windows, claim bar `{lensed_ab['claim_bar_zcut3']}`",
         "`ab_split_scan.csv`"),
        ("background-only toys", f"in 2·10⁴ toys the best confirmation reaches "
         f"`{float(lensed_ab['toy_best_zcut3']):.2f}` against that bar: "
         f"{int(lensed_ab['toy_false_claims_zcut3'])} false claims", "`ab_split_scan.csv`"),
        ("break-even trials inflation", f"`R* = {lensed_ab['R_star_w3']}` at `w = 3`, "
         f"`{lensed_ab['R_star_w1']}` with pinned windows", "`ab_split_scan.csv`"),
        ("the measured defect rate of a published scanner",
         f"{defects['published rate, analytical-function histograms']} of analytical-function and "
         f"{defects['published rate, simulated Standard Model histograms']} of simulated histograms "
         f"flag above 5σ; "
         f"{defects['spurious 5 sigma candidates returned'].replace('+/-', '±')} spurious candidates on "
         f"{defects['application histograms, no signal injected']} histograms",
         "`estimator_defects.csv`"),
        ("... per look, against the Gaussian tail",
         f"{defects['measured defect rate per look']} against {p1(5.0):.2e}, an inflation of "
         f"**{defects['inflation over the Gaussian tail']}×**", "`estimator_defects.csv`"),
        ("... against the break-even",
         f"{defects['margin of the measured inflation over break-even']}× past it, so the split is the "
         "more sensitive procedure", "`estimator_defects.csv`"),
        ("the coherent fraction the split cannot kill",
         f"{defects['candidates surviving the correlation veto']} candidates, "
         f"**{defects['coherent fraction'].replace('+/-', '±')}**", "`estimator_defects.csv`"),
        ("toy validation of the reach", "; ".join(
            f"{t['procedure'].split(',')[0]} {t['reach_toys']} against {t['reach_analytic']} analytic"
            for t in toys) + f" (on the model space, `N = {float(toys[0]['N_trials']):,.0f}`)",
         "`ab_split_toys.csv`"),
        ("the procedure on one toy spectrum", "; ".join(
            f"{g['toy']}: `Z_A = {g['Z_A']}` → `Z_B = {g['Z_B']}`" for g in guard)
         + f", against a claim bar of {guard[0]['claim_bar']}", "`ab_guard_toys.csv`"),
        ("the split's own exposure to that inflation",
         f"the list grows {defects['pre-registered list at Z_cut = 3, nominal']} → "
         f"{float(defects['... if the inflation persists at Z_cut']):,.0f} and the bar "
         f"{defects['claim bar, nominal list']} → {defects['claim bar, inflated list']}",
         "`estimator_defects.csv`"),
    ))

    sel = ["| estimator | argmax | threshold | BH | thr − argmax | nominal `q` | P(argmax *is* the signal) |",
           "|---|--:|--:|--:|--:|--:|--:|"]
    for r in rules:
        if r["mu"] != "5":
            continue
        sel.append(f"| {r['estimator']} | {r['argmax'] or '*untunable*'} | {r['threshold']} | "
                   f"{r['bh'] or '*unreachable*'} | {r['thr_minus_argmax'] or '--'} | "
                   f"{r['bh_nominal_q'] or '--'} | {r['argmax_is_signal']} |")
    blocks.append("### Selecting candidates (Section 7)\n\nConfirmation probability for a 5σ signal, "
                  f"every rule held to the same {float(at5['PERFECT estimator']['argmax_budget']):.2e} "
                  "false-confirmation budget, in per cent (`selection_rules.csv`):\n\n"
                  + "\n".join(sel))

    u = ["| source | varied over | model space | scan | difference |", "|---|---|--:|--:|--:|"]
    for (src, direction), r in unc.items():
        if direction != "envelope" and src != "total":
            continue
        if direction not in ("envelope", "total"):
            continue
        u.append(f"| {src} | {r['varied_over']} | {r['dZ_model_space'] or '--'} | "
                 f"{r['dZ_combinatorial_scan'] or '--'} | {r['dZ_difference'] or '--'} |")
    blocks.append("### The uncertainty budget (Appendix)\n\nEvery declared input, moved over the range "
                  "it is known to (`budget_uncertainty.csv`):\n\n" + "\n".join(u))

    y = ["| object | `F` | symmetric channel | `r` | events/element at 1 TeV | one-event mass |",
         "|---|--:|---|--:|--:|--:|"]
    for k, ch in YM.SYM.items():
        if ch is None:
            y.append(f"| `{k}` | {YM.F[k]:g} | -- | -- | -- | -- |")
            continue
        r, w = res(ch), YM.weight(k + k)
        y.append(f"| `{k}` | {YM.F[k]:g} | `{ch}` | {r:.3f} | {YM.per_element(YM.M_REF, r, w):,.1f} | "
                 f"{YM.one_event_mass(r, w):,.0f} GeV |")
    blocks.append("### The yield model (Appendix)\n\nCalibrated object by object on the Run-2 anchor "
                  "(`yield_model.py`, run it for this table):\n\n" + "\n".join(y))

    body = ("\n\n".join(blocks) + "\n\nEvery row above is written by "
            "[`paper_numbers.py`](searchbudget/stages/paper_numbers.py) from the committed tables, so `make all` "
            "keeps it and the paper in step; nothing here is typed by hand.")

    text = open(paths.README).read()
    i, j = text.index(START), text.index(END)
    io.write_text(paths.README, text[:i] + START + "\n\n" + body + "\n\n" + text[j:])
    print(f"wrote the paper-numbers block into README.md ({len(body.splitlines())} lines, "
          f"{len(blocks)} sections)")
