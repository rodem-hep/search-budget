#!/usr/bin/env python3
"""Two appendix tables that document the two inputs of the budget, per spectrum: which public
model classes point at each spectrum, and how the 46 inclusive spectra become 94 event selections.

The first is the inverse of the coverage matrix drawn by budget_plots.py, over the same 46 spectra:
both build their spectrum set from PUBLIC_OBS through canon(), and rows here follow the ranking of
the budget figure, so a row and a bar line up. Spectrum labels come from obs_labels.mathify, so the
table and the figure axis spell a spectrum the same way.

The second writes out NSEL, whose per-spectrum justification strings are the provenance of the 94.
Both the 94 and N = sum(c_s n_s) = 6597 are asserted here, so the table and the document's
granularity table cannot drift apart.

Reads  scripts/bump_observables.py, scripts/public_obs_map.py (modules; no data files).
Writes results/tex/model_map_appendix.tex, results/tables/model_spectrum_map.csv.
Pure standard library.
"""
import collections, csv, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
from bump_observables import canon, CANON_ORDER, ns_scan
from public_obs_map import PUBLIC_OBS, NSEL, NSEL_DEFAULT, nsel
from obs_labels import mathify, textsafe

def _p(*a): return os.path.join(ROOT, *a)

def _tex(s):
    for a, b in (("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_")):
        s = s.replace(a, b)
    return textsafe(s)

# The print form of the NSEL notes. Those are ASCII because they also go to the terminal reports;
# set in the same words here, typeset. Keyed on the observable, and every spectrum must appear.
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

# Same construction as budget_plots.py, so the two cannot disagree on the spectrum set.
by_obs = collections.defaultdict(set)
for m, obss in PUBLIC_OBS.items():
    for o in obss:
        by_obs[canon(o)].add(m)
order = [o for o in CANON_ORDER if o == canon(o) and o in by_obs]
order += [o for o in sorted(by_obs) if o not in order]
ranked = sorted(order, key=lambda o: -ns_scan(o))     # the order of the budget figure

if len(ranked) != 46:
    raise SystemExit(f"{len(ranked)} spectra carry a model, expected 46")

# the model sets are sets, so their order is hash-randomised per process: sort them
rows = [(o, ns_scan(o), sorted(by_obs[o], key=lambda m: (m.lower(), m))) for o in ranked]

pairs = sum(len(ms) for _, _, ms in rows)
widest, n_widest = max(((o, len(ms)) for o, _, ms in rows), key=lambda t: t[1])

# The selections level. Every spectrum must be named explicitly: a row falling through to
# NSEL_DEFAULT would make the printed column unable to reproduce its own total.
undoc = [o for o in ranked if o not in NSEL]
if undoc:
    raise SystemExit(f"no event-selection note for {undoc}; NSEL_DEFAULT={NSEL_DEFAULT} would hide it")
unset = [o for o in ranked if o not in SEL_TEX]
if unset:
    raise SystemExit(f"no typeset event-selection note for {unset}; extend SEL_TEX")
n_sel = sum(nsel(o) for o in ranked)
big_n = sum(nsel(o) * ns_scan(o) for o in ranked)
big_n_tex = f"{big_n:,.0f}".replace(",", r"\,")     # LaTeX thin space, as in the document's tables
if n_sel != 94:
    raise SystemExit(f"event selections sum to {n_sel}, expected 94")

# ------------------------------------------------------------------ the appendix
with open(_p("results", "tex", "model_map_appendix.tex"), "w") as f:
    f.write("% Generated by scripts/export_model_map_tex.py. Do not edit: regenerate instead.\n")
    f.write("\\section{The Model-to-Spectrum Map}\n\\label{app:modelmap}\n\n")
    f.write(f"""Table~\\ref{{tab:modelmap}} is the input behind the count of Section~\\ref{{sec:budget}}:
each canonical spectrum the budget scans, with the public model classes that predict a resonance in
it. Rows follow Figure~\\ref{{fig:budget}}, most looks first, and $n_s$ restates that figure's bar as
a number.

The model classes are assembled from the public FeynRules model
database~\\cite{{feynrules,feynrules_db}}, which distributes BSM Lagrangians in the UFO
format~\\cite{{ufo}} that the LHC generators read, together with the published search record. One
entry is a class of models sharing a decay topology rather than a single Lagrangian, so
\\emph{{2HDM}} covers the general, type-II and CP-violating implementations together; that is the
granularity at which the distinction matters here, since two implementations peaking in the same
spectrum are tested by the same search. Model classes whose only signature is non-resonant (mono-$X$
missing energy, anomalous $\\mathrm{{d}}E/\\mathrm{{d}}x$, displaced-only decays) are absent, because
they populate no mass spectrum and so buy no trials.

The map is deliberately many-to-many, and its off-diagonal weight is what makes breadth cheap:
{len(PUBLIC_OBS)} model classes make {pairs} class--spectrum pairs over {len(rows)} spectra, so a
spectrum carries {pairs / len(rows):.1f} model classes on average and {mathify(widest)} carries
{n_widest}. Adding a model class to a spectrum already listed costs no trials at all; only opening a
row costs its $n_s$.

\\begin{{longtable}}{{@{{}}>{{\\raggedright\\arraybackslash}}p{{0.17\\textwidth}}r>{{\\raggedright\\arraybackslash}}p{{0.61\\textwidth}}@{{}}}}
\\caption{{The public model classes predicting a resonance in each of the {len(rows)} spectra of the
budget, in the order of Figure~\\ref{{fig:budget}}.}}\\label{{tab:modelmap}}\\\\
\\toprule
spectrum & $n_s$ & public model classes \\\\
\\midrule
\\endfirsthead
\\caption[]{{\\emph{{continued.}}}}\\\\
\\toprule
spectrum & $n_s$ & public model classes \\\\
\\midrule
\\endhead
\\bottomrule
\\endlastfoot
""")
    for o, ns, ms in rows:
        f.write(f"{mathify(o)} & {ns:.0f} & {', '.join(_tex(m) for m in ms)} \\\\\n")
    f.write("\\end{longtable}\n")

    # ---------------------------------------------------------------- the 94
    f.write(f"""
\\section{{From {len(rows)} Spectra to {n_sel} Event Selections}}
\\label{{app:nsel}}

Table~\\ref{{tab:granularity}} quotes {n_sel} published event selections against the {len(rows)}
inclusive spectra, and Table~\\ref{{tab:nsel}} is where that number comes from. For each spectrum it
gives the multiplicity $c_s$ of distinct event selections that a published search family scans along
that same mass axis, and names them. The column sums to {n_sel}, and
$\\Nsig = \\sum_s c_s n_s = {big_n_tex}$ is the selections-level entry of
Table~\\ref{{tab:granularity}}.

What counts as a channel here is a distinct event selection that produces its own bump spectrum: a
$b$-tag category, a boost regime, a sub-decay mode, or an object pairing that is ambiguous and which
the search therefore histograms both ways. A fit category that is statistically combined into a
single limit does not count, and $c_s$ is a multiplicity \\emph{{within}} one spectrum, so it never
reaches across to a different final state: an electron and a muon are different objects with
different resolutions, so $m(ee)$ and $m(\\mu\\mu)$ are two rows of Table~\\ref{{tab:modelmap}} with
a window and an $n_s$ each. Where a search family does combine its channels, treating each as an
independent scan over-counts the trials, which is the conservative direction for a discovery
threshold and is why the selections level is quoted as an upper bracket on the inclusive count
rather than as a replacement for it.

\\begin{{longtable}}{{@{{}}>{{\\raggedright\\arraybackslash}}p{{0.17\\textwidth}}rr>{{\\raggedright\\arraybackslash}}p{{0.45\\textwidth}}@{{}}}}
\\caption{{The published event selections behind the {n_sel} channels of
Table~\\ref{{tab:granularity}}: $c_s$ per spectrum, and what they are. Same row order as
Table~\\ref{{tab:modelmap}}.}}\\label{{tab:nsel}}\\\\
\\toprule
spectrum & $n_s$ & $c_s$ & the event selections counted \\\\
\\midrule
\\endfirsthead
\\caption[]{{\\emph{{continued.}}}}\\\\
\\toprule
spectrum & $n_s$ & $c_s$ & the event selections counted \\\\
\\midrule
\\endhead
\\midrule
\\multicolumn{{2}}{{@{{}}l}}{{total}} & \\textbf{{{n_sel}}} & \\\\
\\bottomrule
\\endlastfoot
""")
    for o, ns, _ in rows:
        f.write(f"{mathify(o)} & {ns:.0f} & {nsel(o)} & {SEL_TEX[o]} \\\\\n")
    f.write("\\end{longtable}\n")

# ------------------------------------------------------------------ the machine-readable form
with open(_p("results", "tables", "model_spectrum_map.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["observable", "ns_scan", "n_model_classes", "model_classes",
                "n_selections", "ns_x_selections", "selections"])
    for o, ns, ms in rows:
        w.writerow([o, f"{ns:.1f}", len(ms), "; ".join(ms),
                    nsel(o), f"{nsel(o) * ns:.1f}", NSEL[o][1]])

print(f"wrote results/tex/model_map_appendix.tex ({len(rows)} spectra, "
      f"{len(PUBLIC_OBS)} model classes, {pairs} pairs; "
      f"{n_sel} event selections, N = {big_n:,.0f})")
print(f"wrote results/tables/model_spectrum_map.csv")
