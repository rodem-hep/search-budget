import collections

from .. import io, paths
from ..registry import stage
from ..viz.labels import mathify

ORDER = "emjbZ"
GLYPH = {"e": "e", "m": r"\mu ", "j": "j", "b": "b", "Z": "Z"}
EXPECTED = (82, 25, 57)


def _covered():
    covered = {}

    def add(comp, obs):
        covered.setdefault("".join(sorted(comp)), obs)

    for c in ["ee", "mm"]:
        add(c, "m(ll) high-mass dilepton")
    add("em", "m(ll') LFV emu")
    for c in ["ej", "mj"]:
        add(c, "m(lj) leptoquark")
    for c in ["eb", "mb"]:
        add(c, "m(lb) 3rd-gen leptoquark (LQ3 -> b tau/b mu)")
    add("jj", "m(jj) dijet")
    add("bb", "m(bb) b-tagged dijet")
    add("jb", "m(jj) dijet (b-tag inclusive / b*)")
    add("jjj", "m(3j) RPV gluino trijet")
    for c in ["eejj", "mmjj", "emjj"]:
        add(c, "m(lljj) W_R Keung-Senjanovic")
    add("jjjj", "m(4j) sgluon / paired dijet")
    add("bbbb", "m(HH)->4b / X->HH")
    add("bbjj", "m(HH) bbWW / bbjj")
    for c in ["eZ", "mZ"]:
        add(c, "m(lZ) VLL/heavy-N -- 8 TeV only, arXiv:1506.01291")
    add("bZ", "m(Zb) VLQ B->Zb -- 8 TeV only, arXiv:1409.5500")
    for c in ["eeZ", "mmZ", "emZ"]:
        add(c, "m(4l) heavy ZZ->llll, arXiv:2009.14791")
    add("jjbZ", "m(Zt) VLQ T->Zt (SCAN window 1-4 TeV)")
    add("Zbb", "m(Vh), Z+h(bb)")
    add("Zjj", "m(VV) semileptonic ZV -> lljj")
    add("ZZ", "m(VV) -> 4l / ZZ")
    add("Zj", "m(VV)/m(Zj) Z+jet resonance")
    add("bjj", "m(t...) hadronic top leg within m(tt)/m(tb)")
    add("bbjj_tt", "m(tt) resolved")
    add("bj", "m(tb) W'->tb (b + hadronic leg)")
    return covered


@stage(
    name="composition-gap",
    group="scan",
    summary="which object compositions of the scan a published search already covers",
    outputs=["tables/composition_gap.txt", "tex/composition_appendix.tex"],
    needs=["tables/combinatorial_budget.csv"],
)
def main(options=None):
    rows = io.read_rows(paths.table("combinatorial_budget.csv"))
    comps = collections.Counter(r["group"] for r in rows)

    def norm(c):
        return "".join(sorted(c))

    COVERED = _covered()
    covered_set = set(COVERED)

    uncovered, covered_hits = [], []
    for c, n in sorted(comps.items(), key=lambda kv: -kv[1]):
        (covered_hits if norm(c) in covered_set else uncovered).append((c, n, len(c)))

    with io.captured(paths.table("composition_gap.txt")):
        print("distinct reachable compositions : %d" % len(comps))
        print("  with a published ATLAS bump hunt: %d" % len(covered_hits))
        print("  with NONE                       : %d" % len(uncovered))
        print()
        byK = collections.Counter(k for _, _, k in uncovered)
        byKall = collections.Counter(len(c) for c in comps)
        print("K : uncovered / total")
        for k in sorted(byKall):
            print("  K=%d : %3d / %3d" % (k, byK.get(k, 0), byKall[k]))
        print()

        sp = collections.Counter()
        spc = collections.Counter()
        lk = collections.Counter()
        lkc = collections.Counter()
        thin = collections.Counter()
        for r in rows:
            k, w, ns = len(r["group"]), int(r["charge_split"]), float(r["n_s"])
            if ns == 0.0:
                thin[k] += w
                continue
            sp[k] += w
            lk[k] += ns
            if norm(r["group"]) in covered_set:
                spc[k] += w
                lkc[k] += ns

        def pc(a, b):
            return 100.0 * a / b if b else 0.0

        print("too thin to fit (no window with 25 elements of >=1 event): %d histograms"
              % sum(thin.values()))
        print("covered fraction  compositions / spectra / looks")
        for k in sorted(sp):
            print("  K=%d : %3d/%3d = %3.0f%%   %4d/%4d = %3.0f%%   %6.0f/%6.0f = %3.0f%%"
                  % (k, byKall[k]-byK.get(k, 0), byKall[k], pc(byKall[k]-byK.get(k, 0), byKall[k]),
                     spc[k], sp[k], pc(spc[k], sp[k]), lkc[k], lk[k], pc(lkc[k], lk[k])))
        print("  all : %3d/%3d = %3.0f%%   %4d/%4d = %3.0f%%   %6.0f/%6.0f = %3.0f%%"
              % (len(covered_hits), len(comps), pc(len(covered_hits), len(comps)),
                 sum(spc.values()), sum(sp.values()), pc(sum(spc.values()), sum(sp.values())),
                 sum(lkc.values()), sum(lk.values()), pc(sum(lkc.values()), sum(lk.values()))))
        print()
        print("=== K=2 compositions with NO published ATLAS bump hunt ===")
        for c, n, k in uncovered:
            if k == 2:
                print("   %-6s (appears in %d exclusive categories)" % (c, n))
        print()
        print("=== K=3 uncovered (first 30) ===")
        print("   " + ", ".join(c for c, n, k in uncovered if k == 3))
        print()
        print("=== K=4 uncovered ===")
        print("   " + ", ".join(c for c, n, k in uncovered if k == 4))

    cats, looks = collections.Counter(), collections.defaultdict(float)
    for r in rows:
        cats[norm(r["group"])] += 1
        looks[norm(r["group"])] += float(r["n_s"])
    if len(cats) != len(comps):
        raise SystemExit(f"{len(comps)} raw groups collapse to {len(cats)} compositions; "
                         f"the summary counts double-count the difference")
    if (len(cats), len(covered_hits), len(uncovered)) != EXPECTED:
        raise SystemExit(f"compositions {len(cats)}, covered {len(covered_hits)}, "
                         f"uncovered {len(uncovered)}; expected {EXPECTED[0]}, {EXPECTED[1]}, "
                         f"{EXPECTED[2]}")

    def comp_tex(c):
        return "$" + "".join(GLYPH[ch] for ch in sorted(c, key=ORDER.index)) + "$"

    def obs_tex(s):
        s = mathify(s).replace("W_R", "$W_R$").replace("->", r"$\to$")
        if "_" in s.replace("$W_R$", ""):
            raise SystemExit(f"unescaped underscore in {s!r}; extend obs_tex")
        return s

    out = paths.tex("composition_appendix.tex")
    io.ensure(out)
    with open(out, "w") as f:
        f.write("% Generated by searchbudget/stages/composition_gap.py. "
                "Do not edit: regenerate instead.\n")
        f.write("\\section{The Combinatorial Scan, Composition by Composition}\n")
        f.write("\\label{app:compositions}\n\n")
        f.write(f"""Every one of the {len(cats)} object compositions the combinatorial scan of
Section~\\ref{{sec:combinatorial}} reaches, with whether a published search scans a mass built from
those object types. This is the row-by-row form of Table~\\ref{{tab:compgap}}: the
{len(covered_hits)} covered and {len(uncovered)} uncovered entries below are that table's columns,
so its counts can be checked rather than taken on trust.

\\emph{{cat.}} is the number of exclusive multiplicity categories the composition occurs in, and
\\emph{{looks}} the independent looks the scan spends on it, summed over those categories; zero looks
means no category of it holds a histogram with the events to support a fit. Coverage
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

    io.note(f"wrote results/tex/composition_appendix.tex ({len(cats)} compositions)")
