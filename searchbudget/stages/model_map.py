from .. import io, paths
from ..core.bump_observables import ns_scan
from ..core.catalogue import canonical_order, models_by_spectrum
from ..core.public_obs_map import NSEL, NSEL_DEFAULT, PUBLIC_OBS, nsel
from ..registry import stage
from ..viz.labels import mathify, textsafe


def _tex(s):
    for a, b in (("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_")):
        s = s.replace(a, b)
    return textsafe(s)


@stage(
    name="model-map",
    group="budget",
    summary="per spectrum: the model classes behind it and the event selections counted",
    outputs=["tex/model_map_appendix.tex", "tables/model_spectrum_map.csv"],
)
def main(options=None):
    SEL_TEX = {
        "m(ee)":         "single dielectron scan (barrel/endcap are fit categories, combined)",
        "m(mumu)":       r"single dimuon scan (combined charge/$\eta$ categories)",
        "m(ee) (Zd)":    r"prompt $e$ lepton-jet / $4e$ scan",
        "m(mumu) (Zd)":  r"prompt muon lepton-jet / $4\mu$ scan",
        "m(emu) LFV":    r"LFV $e\mu$",
        "m(etau) LFV":   r"LFV $e\tau$ (hadronic $\tau$)",
        "m(mutau) LFV":  r"LFV $\mu\tau$ (hadronic $\tau$)",
        "m(ee) SS":      r"$H^{++}$ same-sign $ee$",
        "m(emu) SS":     r"$H^{++}$ same-sign $e\mu$",
        "m(mumu) SS":    r"$H^{++}$ same-sign $\mu\mu$",
        "multilepton":   r"$3\ell$ / $4\ell$ signal regions by flavour and charge "
                         r"(axis deliberately unsplit)",
        "m(gammagamma)": "spin-0 (ggF) + spin-2/VBF (converted/unconverted categories)",
        "m(ej)":         r"LQ pair $\to eq$: all four lepton--jet pairings scanned "
                         r"(arXiv:2006.05872)",
        "m(muj)":        r"LQ pair $\to \mu q$: all four lepton--jet pairings scanned "
                         r"(arXiv:2006.05872)",
        "m(tauj)":       r"LQ $\to \tau q$: $\tau_{\mathrm{had}}$ + $\tau_{\mathrm{lep}}$ selections",
        "m(eb)":         r"RPV stop $\to be$ + $b$-tagged LQ leg (both pairings)",
        "m(mub)":        r"RPV stop $\to b\mu$ + $b$-tagged LQ leg (both pairings)",
        "m(taub)":       r"LQ$_3 \to \tau b$: $\tau_{\mathrm{had}}$ + $\tau_{\mathrm{lep}}$",
        "m(bj)":         r"$b^* \to bg$: $b$-tagged leading + subleading jet pairing",
        "m(eZ)":         r"trilepton $e+Z$ resonance",
        "m(muZ)":        r"trilepton $\mu+Z$ resonance",
        "mT(ev)":        r"$W' \to e\nu$",
        "mT(muv)":       r"$W' \to \mu\nu$",
        "mT(taunu)":     r"$W' \to \tau\nu$ (single hadronic-$\tau$ channel)",
        "m(jj)":         "low-mass (TLA/ISR) + high-mass inclusive dijet",
        "m(cb) dijet":   r"light $H^+ \to cb$ in $t\bar t$: hadronic + semileptonic tag",
        "m(3j)":         "RPV/sgluon three-jet resonance (single SR)",
        "m(bb)":         r"$h \to aa \to 4b$ low-mass + $A/Y \to bb$ resolved + boosted",
        "m(jgamma)":     r"$q^* \to q\gamma$ (single photon+jet SR)",
        "m(tt)":         r"$t\bar t$ resonance: $\ell$+jets resolved, $\ell$+jets boosted, "
                         r"all-hadronic, dilepton",
        "m(tt)/m(jj)":   r"coloron pair $\to tt$ or $jj$ decay legs",
        "m(tb)":         r"$W'/H^+ \to tb$: 0-lepton, 1-lepton $\times$ $b$-tag",
        "m(VV)":         r"$WW/WZ/ZZ$ $\times$ $qqqq$ / $\ell\nu qq$ / $\ell\ell qq$ / "
                         r"$\nu\nu qq$ / $\ell\nu\ell\nu$, ggF+VBF",
        "m(Vh)":         r"$Wh/Zh$ $\times$ 0/1/2-lepton $\times$ $h \to bb$ resolved/boosted",
        "m(HH)":         r"$bbbb$ resolved + boosted, $bb\tau\tau$, $bb\gamma\gamma$, $bbVV$, "
                         r"multilepton",
        "m(ttZ)/m(Zt)":  r"VLQ single/pair $T \to tZ$ selections",
        "m(multi)":      "QBH / multijet (single high-multiplicity SR)",
        "m(tautau)":     r"$H/Z' \to \tau\tau$: $\tau_e\tau_h$, $\tau_\mu\tau_h$, $\tau_h\tau_h$",
        "m(eejj)":       r"$W_R \to eejj$",
        "m(mumujj)":     r"$W_R \to \mu\mu jj$",
        "m(Vgamma)":     r"$X \to Z\gamma$: $Z \to \ell\ell$ and boosted $Z \to qq$ (+ $W\gamma$)",
        "m(egamma)":     r"excited electron $e^*$",
        "m(mugamma)":    r"excited muon $\mu^*$",
        "m(tW)":         r"$b^* \to tW$: leptonic + hadronic top",
        "m(Wb)":         r"single VLQ $T/Y \to Wb$ (1-lepton SR)",
        "m(Ht)":         r"single VLQ $T \to Ht$ ($h \to bb$ tagged)",
    }

    by_obs = models_by_spectrum()
    order = canonical_order(by_obs)
    ranked = sorted(order, key=lambda o: -ns_scan(o))

    if len(ranked) != 46:
        raise SystemExit(f"{len(ranked)} spectra carry a model, expected 46")

    rows = [(o, ns_scan(o), sorted(by_obs[o], key=lambda m: (m.lower(), m))) for o in ranked]

    pairs = sum(len(ms) for _, _, ms in rows)
    widest, n_widest = max(((o, len(ms)) for o, _, ms in rows), key=lambda t: t[1])

    undoc = [o for o in ranked if o not in NSEL]
    if undoc:
        raise SystemExit(f"no event-selection note for {undoc}; NSEL_DEFAULT={NSEL_DEFAULT} would hide it")
    unset = [o for o in ranked if o not in SEL_TEX]
    if unset:
        raise SystemExit(f"no typeset event-selection note for {unset}; extend SEL_TEX")
    n_sel = sum(nsel(o) for o in ranked)
    big_n = sum(nsel(o) * ns_scan(o) for o in ranked)
    big_n_tex = f"{big_n:,.0f}".replace(",", r"\,")
    if n_sel != 94:
        raise SystemExit(f"event selections sum to {n_sel}, expected 94")

    with open(io.ensure(paths.tex("model_map_appendix.tex")), "w") as f:
        f.write("% Generated by searchbudget/stages/model_map.py. Do not edit: regenerate instead.\n")
        f.write("\\section{The Model-to-Spectrum Map}\n\\label{app:modelmap}\n\n")
        f.write(f"""Table~\\ref{{tab:modelmap}} is the input the count of Section~\\ref{{sec:budget}} is built
from: each canonical spectrum the budget scans, the public model classes predicting a resonance in
it, and the published event selections that scan it. Rows follow Figure~\\ref{{fig:budget}}, most
looks first. The model classes come from the public FeynRules
database~\\cite{{feynrules,feynrules_db}}, which distributes BSM Lagrangians in the UFO
format~\\cite{{ufo}} that the LHC generators read, together with the published search record. One
entry is a class of models sharing a decay topology rather than a single Lagrangian, so
\\emph{{2HDM}} covers the general, type-II and CP-violating implementations together: two
implementations peaking in the same spectrum are tested by the same search. Classes whose only
signature is non-resonant (mono-$X$, anomalous $\\mathrm{{d}}E/\\mathrm{{d}}x$, displaced-only decays)
are absent, populating no mass spectrum and buying no trials. The map is many-to-many, and its
off-diagonal weight is what makes breadth cheap: {len(PUBLIC_OBS)} classes make {pairs}
class--spectrum pairs over {len(rows)} spectra, {pairs / len(rows):.1f} per spectrum on average and
{n_widest} on {mathify(widest)}.

A channel $c_s$ is a distinct event selection producing its own bump spectrum: a $b$-tag category, a
boost regime, a sub-decay mode, or an ambiguous object pairing the search histograms both ways. A fit
category combined into a single limit does not count, and $c_s$ is a multiplicity \\emph{{within}} one
spectrum, never reaching across to a different final state. That column sums to {n_sel} and
$\\Nsig = \\sum_s c_s n_s = {big_n_tex}$, the selections-level entry of
Table~\\ref{{tab:granularity}}. If a search family does combine its channels, counting each as an
independent scan over-counts the trials, which is why that level is an upper bracket on the inclusive
count rather than a replacement for it.

\\begingroup\\small\\setlength{{\\tabcolsep}}{{4pt}}
\\begin{{longtable}}{{@{{}}>{{\\raggedright\\arraybackslash}}p{{0.145\\textwidth}}rr%
>{{\\raggedright\\arraybackslash}}p{{0.375\\textwidth}}%
>{{\\raggedright\\arraybackslash}}p{{0.305\\textwidth}}@{{}}}}
\\caption{{The {len(rows)} spectra of the budget, in the order of Figure~\\ref{{fig:budget}}: the public
model classes predicting a resonance in each, and the published event selections counted against
it.}}\\label{{tab:modelmap}}\\\\
\\toprule
spectrum & $n_s$ & $c_s$ & public model classes & event selections \\\\
\\midrule
\\endfirsthead
\\caption[]{{\\emph{{continued.}}}}\\\\
\\toprule
spectrum & $n_s$ & $c_s$ & public model classes & event selections \\\\
\\midrule
\\endhead
\\midrule
\\multicolumn{{2}}{{@{{}}l}}{{total}} & {n_sel} & & \\\\
\\bottomrule
\\endlastfoot
""")
        for o, ns, ms in rows:
            f.write(f"{mathify(o)} & {ns:.0f} & {nsel(o)} & {', '.join(_tex(m) for m in ms)} "
                    f"& {SEL_TEX[o]} \\\\\n")
        f.write("\\end{longtable}\n\\endgroup\n")

    io.write_rows(paths.table("model_spectrum_map.csv"),
                  ["observable", "ns_scan", "n_model_classes", "model_classes",
                   "n_selections", "ns_x_selections", "selections"],
                  [[o, f"{ns:.1f}", len(ms), "; ".join(ms),
                    nsel(o), f"{nsel(o) * ns:.1f}", NSEL[o][1]] for o, ns, ms in rows])

    print(f"wrote results/tex/model_map_appendix.tex ({len(rows)} spectra, "
          f"{len(PUBLIC_OBS)} model classes, {pairs} pairs; "
          f"{n_sel} event selections, N = {big_n:,.0f})")
    print(f"wrote results/tables/model_spectrum_map.csv")
