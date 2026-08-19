from .. import io, paths
from ..core import combinatorial_budget as CB
from ..registry import stage


@stage(
    name="combinatorial-budget",
    group="scan",
    summary="the five-object combinatorial scan, one row per (category, mass group)",
    outputs=["tables/combinatorial_budget.csv"],
)
def main(options=None):
    s = CB.enumerate_scan()
    CB.report(s)
    print()
    print("top object-type multisets by #histograms:")
    for t, v in s.by_type.most_common(12):
        print(f"  {t:6s} {v:5d}")
    io.write_rows(
        paths.table("combinatorial_budget.csv"),
        ["category", "group", "charge_split", "r", "M_lo_GeV", "M_hi_GeV",
         "M_hi_fittable_GeV", "n_s"],
        [[row.cat, row.group, row.split, f"{row.r:.3f}", f"{row.lo:.0f}", f"{row.hi:.0f}",
          f"{row.hi_scan:.0f}", f"{row.n_s:.1f}"] for row in s.rows])
    print(f"\nwrote results/tables/combinatorial_budget.csv ({len(s.rows):,} group rows)")
