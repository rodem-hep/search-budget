#!/usr/bin/env python3
"""Which reachable (<=4 object) invariant-mass compositions have a published ATLAS bump hunt?

Reads results/tables/combinatorial_budget.csv and prints the gap; the Makefile captures it as
results/tables/composition_gap.txt. The coverage mapping is deliberately generous, so the
uncovered list is a lower bound -- see docs/METHOD_NOTES.md.

Also writes results/tex/composition_appendix.tex: the same question laid out composition by
composition, so the summary counts a document quotes can be checked row by row rather than taken
on trust. Every composition appears there, covered or not, with what covers it and what it costs.
"""
import os, csv, collections, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from obs_labels import mathify

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(os.path.join(ROOT, "results/tables/combinatorial_budget.csv"))))
comps = collections.Counter(r["group"] for r in rows)

def norm(c):   # canonical multiset string, e.g. 'jbe' -> 'bej'
    return "".join(sorted(c))

# --- published ATLAS bump observables expressed as object compositions ------------------
# (e=electron, m=muon, j=light jet, b=b-jet, Z=leptonic Z)
COVERED = {}
def add(comp, obs):
    COVERED.setdefault(norm(comp), obs)

# m(ll) high-mass Drell-Yan; m(ll') LFV
for c in ["ee","mm"]:            add(c, "m(ll) high-mass dilepton")
add("em", "m(ll') LFV emu")
# m(lj) leptoquark  (lepton + jet, and the b-flavoured LQ variants)
for c in ["ej","mj"]:            add(c, "m(lj) leptoquark")
for c in ["eb","mb"]:            add(c, "m(lb) 3rd-gen leptoquark (LQ3 -> b tau/b mu)")
# dijet family
add("jj", "m(jj) dijet")
add("bb", "m(bb) b-tagged dijet")
add("jb", "m(jj) dijet (b-tag inclusive / b*)")
add("jjj", "m(3j) RPV gluino trijet")
# W_R -> lljj  and diboson semileptonic
for c in ["eejj","mmjj","emjj"]: add(c, "m(lljj) W_R Keung-Senjanovic")
add("jjjj", "m(4j) sgluon / paired dijet")
add("bbbb", "m(HH)->4b / X->HH")
add("bbjj", "m(HH) bbWW / bbjj")
# lepton + leptonic Z: genuine bump hunt, but 8 TeV ONLY (1506.01291, JHEP 09 (2015) 108).
# Never repeated at Run-2 luminosity -> counted as covered, but see NEVER_EXAMINED.md B.2.
# (The scanned observable is the Z-constrained trilepton mass; loose-3l comps eem/emm stay uncovered.)
for c in ["eZ","mZ"]:            add(c, "m(lZ) VLL/heavy-N -- 8 TeV only, arXiv:1506.01291")
# Zb: single+pair VLQ B->Zb, m(Zb) final discriminant -- 8 TeV ONLY (1409.5500). Abandoned class.
add("bZ", "m(Zb) VLQ B->Zb -- 8 TeV only, arXiv:1409.5500")
# heavy ZZ->4l: m(4l) is the fitted discriminant (2009.14791) -> 2 loose leptons + leptonic Z.
for c in ["eeZ","mmZ","emZ"]:    add(c, "m(4l) heavy ZZ->llll, arXiv:2009.14791")
# single/pair VLQ T->Zt (bump_observables.SCAN m(ttZ)/m(Zt), 1-4 TeV): leptonic Z + hadronic top.
add("jjbZ", "m(Zt) VLQ T->Zt (SCAN window 1-4 TeV)")
# Vh, VV with a leptonic Z
add("Zbb", "m(Vh), Z+h(bb)")
add("Zjj", "m(VV) semileptonic ZV -> lljj")
add("ZZ",  "m(VV) -> 4l / ZZ")
add("Zj",  "m(VV)/m(Zj) Z+jet resonance")
# tt / tb resonances -- top decays to b+jj, so the published m(tt), m(tb) cover these
add("bjj", "m(t...) hadronic top leg within m(tt)/m(tb)")
add("bbjj_tt", "m(tt) resolved")   # already have bbjj
add("bj", "m(tb) W'->tb (b + hadronic leg)")

covered_set = set(COVERED)

uncovered, covered_hits = [], []
for c, n in sorted(comps.items(), key=lambda kv: -kv[1]):
    (covered_hits if norm(c) in covered_set else uncovered).append((c, n, len(c)))

print("distinct reachable compositions : %d" % len(comps))
print("  with a published ATLAS bump hunt: %d" % len(covered_hits))
print("  with NONE                       : %d" % len(uncovered))
print()
byK = collections.Counter(k for _,_,k in uncovered)
byKall = collections.Counter(len(c) for c in comps)
print("K : uncovered / total")
for k in sorted(byKall):
    print("  K=%d : %3d / %3d" % (k, byK.get(k,0), byKall[k]))
print()

# Same coverage question weighted by what the scan actually histograms, and by what it costs:
# a composition is one final state, but it recurs across categories and object indices.
sp   = collections.Counter(); spc = collections.Counter()
lk   = collections.Counter(); lkc = collections.Counter()
for r in rows:
    k, w, ns = len(r["group"]), int(r["charge_split"]), float(r["n_s"])
    sp[k] += w; lk[k] += ns
    if norm(r["group"]) in covered_set:
        spc[k] += w; lkc[k] += ns
pc = lambda a, b: 100.0 * a / b if b else 0.0
print("covered fraction  compositions / spectra / looks")
for k in sorted(sp):
    print("  K=%d : %3d/%3d = %3.0f%%   %4d/%4d = %3.0f%%   %6.0f/%6.0f = %3.0f%%"
          % (k, byKall[k]-byK.get(k,0), byKall[k], pc(byKall[k]-byK.get(k,0), byKall[k]),
             spc[k], sp[k], pc(spc[k], sp[k]), lkc[k], lk[k], pc(lkc[k], lk[k])))
print("  all : %3d/%3d = %3.0f%%   %4d/%4d = %3.0f%%   %6.0f/%6.0f = %3.0f%%"
      % (len(covered_hits), len(comps), pc(len(covered_hits), len(comps)),
         sum(spc.values()), sum(sp.values()), pc(sum(spc.values()), sum(sp.values())),
         sum(lkc.values()), sum(lk.values()), pc(sum(lkc.values()), sum(lk.values()))))
print()
print("=== K=2 compositions with NO published ATLAS bump hunt ===")
for c,n,k in uncovered:
    if k==2: print("   %-6s (appears in %d exclusive categories)" % (c,n))
print()
print("=== K=3 uncovered (first 30) ===")
print("   " + ", ".join(c for c,n,k in uncovered if k==3))
print()
print("=== K=4 uncovered ===")
print("   " + ", ".join(c for c,n,k in uncovered if k==4))

# ---------------------------------------------------------------- the appendix
# Aggregate on the canonical multiset: two spellings of one composition would otherwise be two rows.
cats, looks = collections.Counter(), collections.defaultdict(float)
for r in rows:
    cats[norm(r["group"])] += 1
    looks[norm(r["group"])] += float(r["n_s"])
if len(cats) != len(comps):
    raise SystemExit(f"{len(comps)} raw groups collapse to {len(cats)} compositions; "
                     f"the summary counts double-count the difference")
if (len(cats), len(covered_hits), len(uncovered)) != (82, 25, 57):
    raise SystemExit(f"compositions {len(cats)}, covered {len(covered_hits)}, "
                     f"uncovered {len(uncovered)}; expected 82, 25, 57")

ORDER = "emjbZ"
GLYPH = {"e": "e", "m": r"\mu ", "j": "j", "b": "b", "Z": "Z"}

def comp_tex(c):
    return "$" + "".join(GLYPH[ch] for ch in sorted(c, key=ORDER.index)) + "$"

def obs_tex(s):
    s = mathify(s).replace("W_R", "$W_R$").replace("->", r"$\to$")
    if "_" in s.replace("$W_R$", ""):
        raise SystemExit(f"unescaped underscore in {s!r}; extend obs_tex")
    return s

out = os.path.join(ROOT, "results", "tex", "composition_appendix.tex")
with open(out, "w") as f:
    f.write("% Generated by scripts/composition_gap.py. Do not edit: regenerate instead.\n")
    f.write("\\section{The Combinatorial Scan, Composition by Composition}\n")
    f.write("\\label{app:compositions}\n\n")
    f.write(f"""Every one of the {len(cats)} object compositions the combinatorial scan of
Section~\\ref{{sec:combinatorial}} reaches, with whether a published search scans a mass built from
those object types. This is the row-by-row form of Table~\\ref{{tab:compgap}}: the
{len(covered_hits)} covered and {len(uncovered)} uncovered entries below are that table's columns,
so its counts can be checked rather than taken on trust.

\\emph{{cat.}} is the number of exclusive multiplicity categories the composition occurs in, and
\\emph{{looks}} the independent looks the scan spends on it, summed over those categories. Coverage
is assigned as generously as can be defended: any published ATLAS search scanning a mass of those
object types counts, and a published $V$ may stand for a $Z$, a $jj$ pair or an $\\ell\\ell$ pair. The
uncovered set is therefore a lower bound on the gap. Two entries are covered only by an $8\\TeV$
search and are marked. Nothing at $K=2$ is uncovered, which is a statement about this five-object
alphabet and not about two-body space, as Table~\\ref{{tab:twobody}} shows.

\\begin{{longtable}}{{@{{}}lrr>{{\\raggedright\\arraybackslash}}p{{0.52\\textwidth}}@{{}}}}
\\caption{{The {len(cats)} reachable object compositions of the combinatorial scan, with the
published search that covers each or a dash where none does.}}\\label{{tab:compositions}}\\\\
\\toprule
composition & cat. & looks & published scan covering it \\\\
\\midrule
\\endfirsthead
\\caption[]{{\\emph{{continued.}}}}\\\\
\\toprule
composition & cat. & looks & published scan covering it \\\\
\\midrule
\\endhead
\\bottomrule
\\endlastfoot
""")
    for k in sorted({len(c) for c in cats}):
        ck = [c for c in cats if len(c) == k]
        nc = sum(1 for c in ck if c in covered_set)
        f.write(f"\\multicolumn{{4}}{{@{{}}l}}{{\\textbf{{$K={k}$}}: {len(ck)} compositions, "
                f"{nc} covered, {len(ck) - nc} not}} \\\\[1pt]\n")
        for c in sorted(ck, key=lambda c: (c not in covered_set, -looks[c])):
            f.write(f"\\quad {comp_tex(c)} & {cats[c]} & {looks[c]:.0f} & "
                    f"{obs_tex(COVERED[c]) if c in covered_set else '---'} \\\\\n")
        f.write("\\addlinespace[3pt]\n")
    f.write("\\end{longtable}\n")

# stdout is the report the Makefile captures, so progress goes to stderr
print(f"wrote results/tex/composition_appendix.tex ({len(cats)} compositions)", file=sys.stderr)
