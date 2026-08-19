import collections
import math

from .. import io, paths
from ..core import combinatorial_budget as CB
from ..core import yield_model as YM
from ..core.bump_observables import n_s, res, scan_segments, z_local_for_global5
from ..core.scan_alphabet import (BASE, DATASETS, LENSES, MET_KEY, ORDER_WIDE,
                                  SIGMA_BASE_DERIVED, SIGMA_WIDE, WIDE, WIDE_ARGS, lens_views)
from ..registry import stage

TRIALS_BUDGET = 5.0e5

MOTIVATED = {
    "m(gammagamma)": ("gg",),      "m(egamma)": ("eg",),     "m(mugamma)": ("gm",),
    "m(jgamma)": ("gj",),          "m(Vgamma)": ("Vg",),
    "m(ee)": ("ee",),              "m(ee) SS": ("ee",),      "m(ee) (Zd)": ("ee",),
    "m(mumu)": ("mm",),            "m(mumu) SS": ("mm",),    "m(mumu) (Zd)": ("mm",),
    "m(emu) LFV": ("em",),         "m(emu) SS": ("em",),
    "m(tautau)": ("TT",),          "m(etau) LFV": ("eT",),   "m(mutau) LFV": ("Tm",),
    "m(ej)": ("ej",),              "m(muj)": ("jm",),        "m(tauj)": ("Tj",),
    "m(eb)": ("be",),              "m(mub)": ("bm",),        "m(taub)": ("Tb",),
    "m(eZ)": ("eZ",),              "m(muZ)": ("mZ",),
    "m(jj)": ("jj",),              "m(bj)": ("bj",),         "m(bb)": ("bb",),
    "m(cb) dijet": ("bj",),        "m(3j)": ("jjj",),
    "m(eejj)": ("eejj",),          "m(mumujj)": ("jjmm",),
    "m(VV)": ("VV",),              "m(Vh)": ("HV",),         "m(HH)": ("HH",),
    "m(tt)": ("tt",),              "m(tt)/m(jj)": ("tt", "jj"),
    "m(tb)": ("bt",),              "m(tW)": ("Vt",),         "m(Wb)": ("Vb",),
    "m(Ht)": ("Ht",),              "m(ttZ)/m(Zt)": ("Vt", "Zt"),
    "multilepton": ("eee", "eem", "emm", "mmm"),
    "m(multi)": None,
    "mT(ev)": None, "mT(muv)": None, "mT(taunu)": None,
}


def canon(c):
    return "".join(sorted(c))


MOTIVATED_COMPS = {canon(c) for cs in MOTIVATED.values() if cs for c in cs}


def tiers(rows):
    sp = [r for r in rows if r.n_s > 0]
    best = {}
    for i, r in enumerate(sp):
        k = canon(r.group)
        if k in MOTIVATED_COMPS and (k not in best or r.w > sp[best[k]].w):
            best[k] = i
    once = set(best.values())
    tier = [0 if i in once else 1 if canon(r.group) in MOTIVATED_COMPS else 2
            for i, r in enumerate(sp)]
    return sp, tier, best


@stage(
    name="scaled-scan",
    group="scan",
    summary="the ten-object scan: tiers, selection lenses, and the per-dataset summary",
    outputs=["tables/scaled_scan.csv", "tables/priority_scan.csv", "tables/lens_scan.csv",
             "tables/scan_summary.csv", "tables/scaled_scan.txt"],
)
def main(options=None):
    with io.captured(paths.table("scaled_scan.txt")):
        print("object alphabet: resolution and yield factor from the published symmetric channel")
        for k in ORDER_WIDE:
            ch = YM.SYM[k]
            print(f"  {k}  {WIDE[k][0]:12s} sigma = {SIGMA_WIDE[k]:.3f}  yield {YM.F[k]:8.1e}   "
                  f"({ch if ch else 'baseline value, no symmetric channel'})")
        print(f"  {MET_KEY}  {WIDE[MET_KEY][0]:12s} {'':14s} yield {YM.F[MET_KEY]:8.1e}   "
              f"(a category split, never a mass)")
        print(f"\nmodel-motivated compositions: {len(MOTIVATED_COMPS)} from "
              f"{sum(1 for v in MOTIVATED.values() if v)} of the {len(MOTIVATED)} axes "
              f"(no composition for {', '.join(k for k, v in MOTIVATED.items() if v is None)})")

        pub_axes = {a: cs for a, cs in MOTIVATED.items() if cs}

        def gated_published():
            kept, N_gated, N_full = 0, 0.0, 0.0
            for a, cs in pub_axes.items():
                r, w = res(a), max(YM.weight(c) for c in cs)
                N_full += sum(n_s(lo, hi, r) for lo, hi in scan_segments(a))
                g = 0.0
                for lo, hi in scan_segments(a):
                    _hs, ns, _ev, fits = YM.gate(lo, hi, r, w)
                    g += ns if fits else 0.0
                N_gated += g
                kept += g > 0
            return kept, N_gated, N_full

        pub_ok, pub_N_gated, pub_N = gated_published()
        print(f"the requirement is never applied to a published window. On the {len(pub_axes)} "
              f"published axes this alphabet can form it would leave {pub_ok} of them and "
              f"N = {pub_N_gated:,.0f} of {pub_N:,.0f}, so exempting them is the conservative "
              f"choice")
        print()

        VARIANTS = [
            ("five objects, lepton trigger", BASE),
            ("five objects, derived sigma",  {**BASE, "sigma": SIGMA_BASE_DERIVED}),
            ("ten objects, any trigger",     WIDE_ARGS),
        ]

        runs = []
        for label, kw in VARIANTS:
            s = CB.enumerate_scan(**kw)
            runs.append((label, kw, s))
            print(f"=== {label}  ({len(kw['order'])} object types, K<={kw['kmax']} objects per "
                  f"mass, <={kw['nobj']} per category, "
                  f"{'lepton required' if kw.get('trig', CB.TRIG) else 'any trigger'})")
            CB.report(s)
            print()

        full = runs[-1][2]

        spectra, tier_list, best = tiers(full.rows)
        ranked = sorted(((tier_list[i], -r.w, r.n_s, canon(r.group), r.cat, r.split)
                         for i, r in enumerate(spectra)),
                        key=lambda x: (x[0], x[1], x[2]))

        kept_looks, kept_n = collections.Counter(), collections.Counter()
        drop_looks, drop_n = collections.Counter(), collections.Counter()
        tier_N = collections.Counter()
        tier_n = collections.Counter()
        N_sel, n_sel, cut_rate, stopped = 0.0, 0, None, False
        kept_cats, tier0_cats, tier_axes = set(), set(), collections.Counter()
        for tier, negrate, ns, key, cat, nh in ranked:
            tier_N[tier] += ns
            tier_n[tier] += nh
            tier_axes[tier] += 1
            if tier == 0:
                tier0_cats.add(cat)
            if not stopped and N_sel + ns <= TRIALS_BUDGET:
                N_sel += ns
                n_sel += nh
                kept_looks[key] += ns
                kept_n[key] += nh
                kept_cats.add(cat)
                cut_rate = -negrate
            else:
                stopped = True
                drop_looks[key] += ns
                drop_n[key] += nh
        tierA_N = tier_N[0] + tier_N[1]

        print("=== full ten-object scan")
        print(f"spectra {full.n_hist:,}   N = {full.N:,.0f}   "
              f"Z_local = {z_local_for_global5(full.N):.2f}")
        print(f"tier 0, every motivated composition once: {tier_n[0]:6,d} spectra, "
              f"N = {tier_N[0]:10,.0f} ({100*tier_N[0]/full.N:4.1f} % of the scan), "
              f"Z_local = {z_local_for_global5(tier_N[0]):.2f}")
        print(f"tier 1, those compositions elsewhere   : {tier_n[1]:6,d} spectra, "
              f"N = {tier_N[1]:10,.0f} ({100*tier_N[1]/full.N:4.1f} %)")
        print(f"tier 2, the rest                       : {tier_n[2]:6,d} spectra, "
              f"N = {tier_N[2]:10,.0f} ({100*tier_N[2]/full.N:4.1f} %)")
        print(f"tiers 0+1 together: N = {tierA_N:,.0f}, which is "
              f"{tierA_N/TRIALS_BUDGET:.2f} times the budget on its own")
        NEW_TYPES = "TgtVH"
        _nn = sum(v for c, v in full.by_type.items() if any(k in c for k in NEW_TYPES))
        _nl = sum(v for c, v in full.looks.items() if any(k in c for k in NEW_TYPES))
        print(f"spectra reaching a new object type: {_nn:,} of {full.n_hist:,} "
              f"({100*_nn/full.n_hist:.0f} %), carrying {100*_nl/full.N:.0f} % of N; "
              f"two-body groups are {100*full.by_size[2]/full.n_hist:.0f} % of the spectra")
        _lost = sorted(MOTIVATED_COMPS - set(full.by_type))
        print(f"motivated compositions with no fittable histogram in any category: "
              f"{', '.join(_lost) if _lost else 'none'}")
        print()
        print(f"=== priority prefix that fits N <= {TRIALS_BUDGET:,.0f}")
        print(f"selected {n_sel:,} of {full.n_hist:,} spectra ({100*n_sel/full.n_hist:.1f} %) over "
              f"{len(kept_n)} of {len(full.by_type)} compositions and {len(kept_cats):,} of "
              f"{full.n_cat:,} categories")
        print(f"tier 0 is {tier_axes[0]} axes in {len(tier0_cats)} categories, "
              f"{tier_n[0]} histograms")
        if stopped:
            print(f"of the {tier_n[1]:,} tier-1 spectra, {n_sel - tier_n[0]:,} fit "
                  f"({100*(n_sel - tier_n[0])/tier_n[1]:.0f} %)")
            print(f"the cut lands at a category yield of {cut_rate:.1e}: thinner categories are "
                  f"dropped")
        else:
            print(f"nothing is cut: every fittable spectrum fits inside the budget, so the "
                  f"statistics requirement binds first and the priority order never has to be "
                  f"applied")
        print(f"N = {N_sel:,.0f} ({100*N_sel/full.N:.1f} % of the full scan), "
              f"Z_local = {z_local_for_global5(N_sel):.2f} "
              f"(band {z_local_for_global5(N_sel*0.5):.2f}-{z_local_for_global5(N_sel*2):.2f})")
        print(f"a local 5 sigma is then worth Z_global = "
              f"{math.sqrt(max(25.0 - 2*math.log(N_sel), 0)):.2f} sigma")
        print()

        print("yield anchor scaled (the model's own uncertainty): spectra, N, Z_local")
        _n_ref = YM.N_REF
        for _s in (0.01, 1.0, 100.0):
            YM.N_REF = _n_ref * _s
            _r = CB.enumerate_scan(**WIDE_ARGS, collect=False)
            print(f"  x{_s:<6g} {_r.n_hist:7,d} of {_r.n_hist + _r.n_thin:7,d} histograms, "
                  f"N = {_r.N:11,.0f}, Z_local = {z_local_for_global5(_r.N):.2f}")
        YM.N_REF = _n_ref
        print()

        print("object type    spectra in scan  selected   kept looks   what survives")
        for k in ORDER_WIDE:
            tot = sum(v for c, v in full.by_type.items() if k in c)
            sel = sum(v for c, v in kept_n.items() if k in c)
            lk = sum(v for c, v in kept_looks.items() if k in c)
            comps = sorted((c for c in kept_n if k in c), key=lambda c: -kept_looks[c])
            extra = [c for c in comps if c not in MOTIVATED_COMPS]
            if not comps:
                what = "DROPPED ENTIRELY"
            else:
                what = f"{len(comps)} compositions: " + " ".join(comps[:6]) + \
                       (" ..." if len(comps) > 6 else "") + \
                       (f"  (+{len(extra)} unmotivated)" if extra else "  (all motivated)")
            print(f"  {k} {WIDE[k][0]:12s} {tot:9,d} {sel:9,d} {lk:12,.0f}   {what}")
        print()

        print("compositions dropped entirely, costliest first:")
        gone = [(c, drop_looks[c], drop_n[c]) for c in drop_looks if c not in kept_n]
        for c, lk, n in sorted(gone, key=lambda x: -x[1])[:15]:
            print(f"  {c:6s} {n:6,d} spectra, {lk:10,.0f} looks")
        print(f"  ... {len(gone)} of {len(full.by_type)} compositions dropped entirely")
        print()
        print("costliest compositions that are kept:")
        for c in sorted(kept_looks, key=lambda c: -kept_looks[c])[:15]:
            tag = "motivated" if c in MOTIVATED_COMPS else "by rate"
            print(f"  {c:6s} {kept_n[c]:6,d} of {full.by_type[c]:6,d} spectra, "
                  f"{kept_looks[c]:9,.0f} looks   ({tag})")
        print()

        items, tier_of = [], {}
        lens_n, lens_N = collections.Counter(), collections.Counter()
        lens_thin = collections.Counter()
        for i, r in enumerate(spectra):
            key = canon(r.group)
            tier = tier_list[i]
            tier_of[i] = (tier, key)
            items.append((tier, -r.w, 0, r.n_s, key, r.cat, r.split, ""))
        for i, lk, li, lns, fits in lens_views(spectra):
            r, (tier, key) = spectra[i], tier_of[i]
            if not fits:
                lens_thin[lk] += r.split
                continue
            items.append((max(tier, 1), -r.w, li, lns, key, r.cat, r.split, lk))
            lens_n[lk] += r.split
            lens_N[lk] += lns
        items.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        lensed_n = full.n_hist + sum(lens_n.values())
        lensed_N = full.N + sum(lens_N.values())

        L_N, L_views, L_spec, L_stop = 0.0, 0, 0, False
        L_comps, L_cats = set(), set()
        L_lens_n, L_lens_N = collections.Counter(), collections.Counter()
        for tier, negrate, li, ns, key, cat, nh, lk in items:
            if L_stop or L_N + ns > TRIALS_BUDGET:
                L_stop = True
                continue
            L_N += ns
            L_views += nh
            L_comps.add(key)
            L_cats.add(cat)
            if li == 0:
                L_spec += nh
            else:
                L_lens_n[lk] += nh
                L_lens_N[lk] += ns

        print("selection lenses on the same axes, one at a time:")
        for key, label, rule, _ok, cap in LENSES:
            print(f"  {label:20s} {rule:34s} eff {YM.LENS_EFF[key]:5.3f} {lens_n[key]:7,d} views "
                  f"{lens_N[key]:9,.0f} looks, {lens_thin[key]:6,d} too thin"
                  f"{'' if cap is None else f' (window capped at {cap:.0f} GeV)'}")
        print(f"  {'all four':20s} {'':34s} {'':9s} {sum(lens_n.values()):7,d} views "
              f"{sum(lens_N.values()):9,.0f} looks, {sum(lens_thin.values()):6,d} too thin")
        print(f"the ten-object scan with lenses: {lensed_n:,} histograms "
              f"({lensed_n/full.n_hist:.1f} per spectrum), N = {lensed_N:,.0f}, "
              f"Z_local = {z_local_for_global5(lensed_N):.2f}")
        print()
        print(f"=== priority prefix with lenses that fits N <= {TRIALS_BUDGET:,.0f}")
        print(f"selected {L_spec:,} of {full.n_hist:,} spectra ({100*L_spec/full.n_hist:.1f} %) "
              f"plus {L_views - L_spec:,} lens views of them: {L_views:,} histograms "
              f"({100*L_views/lensed_n:.1f} % of the lensed scan) over {len(L_comps)} compositions "
              f"and {len(L_cats):,} of {full.n_cat:,} categories")
        print(f"N = {L_N:,.0f}, Z_local = {z_local_for_global5(L_N):.2f}")
        for key, label, _rule, _ok, _cap in LENSES:
            frac = 100.0 * L_lens_n[key] / lens_n[key] if lens_n[key] else 0.0
            print(f"  {label:20s} {L_lens_n[key]:7,d} of {lens_n[key]:7,d} views kept "
                  f"({frac:4.1f} %), {L_lens_N[key]:9,.0f} looks")
        if L_spec < n_sel:
            print(f"the same budget with no lens at all reaches {n_sel:,} spectra, so the lenses "
                  f"are paid for in coverage: {n_sel - L_spec:,} fewer axes-in-categories for "
                  f"{L_views - L_spec:,} lens views of the ones that remain")
        else:
            print(f"the lenses cost nothing in coverage: axes and lens views together still fit "
                  f"the budget, and {sum(lens_thin.values()):,} further views are ruled out by "
                  f"statistics alone")
        print()

        dataset_rows, summary_rows = [], []
        print(f"dataset scaled (the anchor is {DATASETS[0][0]}): spectra, N, Z_local, mass reach")
        for _label, _s in DATASETS:
            YM.N_REF = _n_ref * _s
            _b = CB.enumerate_scan(**BASE)
            _r = CB.enumerate_scan(**WIDE_ARGS)
            _sp, _tier, _best_i = tiers(_r.rows)
            _best = {canon(_sp[i].group): _sp[i] for i in _best_i.values()}
            _tn, _tN = collections.Counter(), collections.Counter()
            for i, x in enumerate(_sp):
                _tn[_tier[i]] += x.split
                _tN[_tier[i]] += x.n_s
            _ln_n, _ln_N, _ln_thin = (collections.Counter(), collections.Counter(),
                                      collections.Counter())
            for i, lk, _li, lns, fits in lens_views(_sp):
                if fits:
                    _ln_n[lk] += _sp[i].split
                    _ln_N[lk] += lns
                else:
                    _ln_thin[lk] += _sp[i].split
            _ln, _lN = _r.n_hist + sum(_ln_n.values()), _r.N + sum(_ln_N.values())
            _t0_n, _t0_N = _tn[0], _tN[0]
            _mot_n, _mot_N = _tn[0] + _tn[1], _tN[0] + _tN[1]
            dataset_rows.append((_label, _s, _b, _r, _ln, _lN, _best,
                                 dict(_ln_n), dict(_ln_N), dict(_ln_thin)))
            print(f"  {_label} (x{_s:g} the anchor, "
                  f"x{_s ** (1.0 / (YM.P - 1.0)):.2f} in mass reach)")
            print(f"    five objects   {_b.n_hist:6,d} of {_b.n_hist + _b.n_thin:6,d} histograms, "
                  f"N = {_b.N:9,.0f}, Z_local = {z_local_for_global5(_b.N):.2f}")
            print(f"    ten objects    {_r.n_hist:6,d} of {_r.n_hist + _r.n_thin:6,d} histograms, "
                  f"N = {_r.N:9,.0f}, Z_local = {z_local_for_global5(_r.N):.2f}")
            print(f"    with lenses    {_ln:6,d} histograms, "
                  f"N = {_lN:9,.0f}, Z_local = {z_local_for_global5(_lN):.2f}")
            print(f"    motivated once {_t0_n:6,d} spectra in "
                  f"{len({x.cat for x in _best.values()})} categories over {len(_best)} "
                  f"compositions, N = {_t0_N:9,.0f}, "
                  f"Z_local = {z_local_for_global5(_t0_N):.2f}")
            print(f"    tier 0 / 1 / 2 spectra : {_tn[0]:6,d} {_tn[1]:6,d} {_tn[2]:6,d}")
            print(f"    tier 0 / 1 / 2 looks   : {_tN[0]:9,.0f} {_tN[1]:9,.0f} {_tN[2]:9,.0f}")
            print(f"    on a motivated axis    : {_mot_n:,} of {_r.n_hist:,} spectra "
                  f"({100*_mot_n/_r.n_hist:.0f} %), {100*_mot_N/_r.N:.0f} % of N, over "
                  f"{len(_best)} of {len(_r.by_type)} object multisets")
            print(f"    by group size          : " + ", ".join(
                f"k={_k}: {_v:,} ({100*_v/_r.n_hist:.0f} %)"
                for _k, _v in sorted(_r.by_size.items())))
            print(f"    costliest compositions : " + ", ".join(
                f"{_c} {_r.looks[_c]:,.0f}"
                for _c in sorted(_r.looks, key=lambda c: -_r.looks[c])[:6]))
            for _key, _lab, _rule, _ok, _cap in LENSES:
                print(f"    lens {_lab:20s} {_ln_n[_key]:7,d} views {_ln_N[_key]:9,.0f} looks, "
                      f"{_ln_thin[_key]:6,d} too thin")
            _pk, _pg, _pf = gated_published()
            print(f"    published windows, if they were gated too: {_pk} of {len(pub_axes)} axes, "
                  f"N = {_pg:,.0f} of {_pf:,.0f}")
            print(f"    yield anchor x0.01 / x100 (this dataset): ", end="")
            _anchor = {}
            for _f in (0.01, 100.0):
                YM.N_REF = _n_ref * _s * _f
                _y = CB.enumerate_scan(**WIDE_ARGS, collect=False)
                _anchor[_f] = (_y.n_hist, z_local_for_global5(_y.N))
                print(f"{_y.n_hist:,} spectra Z = {z_local_for_global5(_y.N):.2f}   ", end="")
            print()
            YM.N_REF = _n_ref * _s
            summary_rows.append({
                "dataset": _label, "luminosity_scale": f"{_s:g}",
                "mass_reach_scale": f"{_s ** (1.0 / (YM.P - 1.0)):.2f}",
                "combinations": _r.n_hist + _r.n_thin, "categories": _r.n_cat,
                "fittable_spectra": _r.n_hist, "N_trials": f"{_r.N:.0f}",
                "Z_local": f"{z_local_for_global5(_r.N):.2f}",
                "lensed_histograms": _ln, "lensed_N_trials": f"{_lN:.0f}",
                "lensed_Z_local": f"{z_local_for_global5(_lN):.2f}",
                "tier0_spectra": _tn[0], "tier1_spectra": _tn[1], "tier2_spectra": _tn[2],
                "tier0_looks": f"{_tN[0]:.0f}", "tier1_looks": f"{_tN[1]:.0f}",
                "tier2_looks": f"{_tN[2]:.0f}",
                "tier0_Z_local": f"{z_local_for_global5(_tN[0]):.2f}",
                "motivated_spectra": _mot_n, "motivated_frac": f"{_mot_n/_r.n_hist:.3f}",
                "compositions_fittable": len(_r.by_type), "compositions_motivated": len(_best),
                "k2_spectra": _r.by_size.get(2, 0), "k3_spectra": _r.by_size.get(3, 0),
                "k4_spectra": _r.by_size.get(4, 0),
                "published_axes_gated": _pk, "published_axes_total": len(pub_axes),
                "published_looks_gated": f"{_pg:.0f}", "published_looks_total": f"{_pf:.0f}",
                "anchor_x0.01_spectra": _anchor[0.01][0], "anchor_x0.01_Z": f"{_anchor[0.01][1]:.2f}",
                "anchor_x100_spectra": _anchor[100.0][0], "anchor_x100_Z": f"{_anchor[100.0][1]:.2f}",
                "costliest_compositions": "; ".join(
                    f"{_c}:{_r.looks[_c]:.0f}" for _c in
                    sorted(_r.looks, key=lambda c: -_r.looks[c])[:6])})
        YM.N_REF = _n_ref
        print()

    def row(label, kw, cats, nsp, ncomp, N):
        return [label, len(kw["order"]), kw["kmax"], kw["nobj"],
                "lepton" if kw.get("trig", CB.TRIG) else "any", cats, nsp, ncomp, f"{N:.0f}",
                f"{z_local_for_global5(N):.2f}", f"{z_local_for_global5(N*0.5):.2f}",
                f"{z_local_for_global5(N*2):.2f}"]

    scaled = [row(label, kw, s.n_cat, s.n_hist, len(s.by_type), s.N) for label, kw, s in runs]
    scaled.append(row("... model-motivated compositions once each", VARIANTS[-1][1],
                      len(tier0_cats), tier_n[0], len(MOTIVATED_COMPS & set(full.by_type)),
                      tier_N[0]))
    scaled.append(row(f"... prioritised to N <= {TRIALS_BUDGET:.0e}", VARIANTS[-1][1],
                      len(kept_cats), n_sel, len(kept_n), N_sel))
    scaled.append(row("ten objects, with selection lenses", VARIANTS[-1][1], full.n_cat, lensed_n,
                      len(full.by_type), lensed_N))
    scaled.append(row(f"... prioritised to N <= {TRIALS_BUDGET:.0e}, lenses included",
                      VARIANTS[-1][1], len(L_cats), L_views, len(L_comps), L_N))
    for _label, _s, _b, _r, _ln, _lN, _best, _lnn, _lnN, _lnt in dataset_rows:
        if _s == 1.0:
            continue
        scaled.append(row(f"five objects, lepton trigger ({_label})", VARIANTS[0][1], _b.n_cat,
                          _b.n_hist, len(_b.by_type), _b.N))
        scaled.append(row(f"ten objects, any trigger ({_label})", VARIANTS[-1][1], _r.n_cat,
                          _r.n_hist, len(_r.by_type), _r.N))
        scaled.append(row(f"ten objects, with selection lenses ({_label})", VARIANTS[-1][1],
                          _r.n_cat, _ln, len(_r.by_type), _lN))
        scaled.append(row(f"... model-motivated compositions once each ({_label})",
                          VARIANTS[-1][1], len({x.cat for x in _best.values()}),
                          sum(x.split for x in _best.values()), len(_best),
                          sum(x.n_s for x in _best.values())))
    io.write_rows(paths.table("scaled_scan.csv"),
                  ["scan", "types", "K_max", "objects_per_category", "trigger", "categories",
                   "spectra", "compositions", "N_trials", "Z_local", "Z_lo", "Z_hi"], scaled)

    io.write_rows(
        paths.table("priority_scan.csv"),
        ["composition", "K", "tier", "spectra_total", "spectra_kept", "looks_total", "looks_kept"],
        [[c, len(c), "motivated" if c in MOTIVATED_COMPS else "rate", full.by_type[c],
          kept_n.get(c, 0), f"{full.looks[c]:.0f}", f"{kept_looks.get(c, 0.0):.0f}"]
         for c in sorted(full.by_type, key=lambda c: (len(c), -full.looks[c]))])

    lens_rows = []
    for _label, _s, _b, _r, _ln, _lN, _best, _lnn, _lnN, _lnt in dataset_rows:
        for key, label, rule, _ok, cap in LENSES:
            kept, thin = _lnn.get(key, 0), _lnt.get(key, 0)
            lens_rows.append([_label, label, rule, "" if cap is None else f"{cap:.0f}",
                              f"{YM.LENS_EFF[key]:g}", kept + thin, kept, thin,
                              f"{_lnN.get(key, 0.0):.0f}"])
    io.write_rows(paths.table("lens_scan.csv"),
                  ["dataset", "lens", "requirement", "window_cap_GeV", "efficiency",
                   "spectra_reached", "views_kept", "views_too_thin", "looks_kept"], lens_rows)

    io.write_dicts(paths.table("scan_summary.csv"), summary_rows)

    io.note(f"\nwrote results/tables/scaled_scan.csv "
            f"({len(runs) + 4 + 4 * (len(DATASETS) - 1)} rows), priority_scan.csv "
            f"({len(full.by_type)} compositions) and lens_scan.csv "
            f"({len(LENSES) * len(DATASETS)} rows) and scan_summary.csv "
            f"({len(summary_rows)} datasets)")
