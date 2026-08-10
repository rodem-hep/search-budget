#!/usr/bin/env python3
"""Which reachable (<=4 object) invariant-mass compositions have a published ATLAS bump hunt?

Reads results/tables/combinatorial_budget.csv and prints the gap; the Makefile captures it as
results/tables/composition_gap.txt. The coverage mapping is deliberately generous, so the
uncovered list is a lower bound -- see docs/METHOD_NOTES.md.
"""
import os, csv, collections

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
print("=== K=2 compositions with NO published ATLAS bump hunt ===")
for c,n,k in uncovered:
    if k==2: print("   %-6s (appears in %d exclusive categories)" % (c,n))
print()
print("=== K=3 uncovered (first 30) ===")
print("   " + ", ".join(c for c,n,k in uncovered if k==3))
print()
print("=== K=4 uncovered ===")
print("   " + ", ".join(c for c,n,k in uncovered if k==4))
