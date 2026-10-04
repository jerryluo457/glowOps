"""Build dispute-level escalation measures and dyadic MID history from COW MID v3.02.

Outputs (data/clean/):
  mid_bilateral.csv  one row per MID that was bilateral on day one (one originator per side),
                     with the six escalation measures and the territorial-dispute flag.
  mid_dyad_onsets.csv  every (dyad, onset year) pair from any MID (incl. multilateral),
                     used for peace-years counting.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "mid302"
OUT = ROOT / "data" / "clean"


def load():
    a = pd.read_csv(RAW / "MIDA_3.021.csv")
    b = pd.read_csv(RAW / "MIDB_3.02.csv")
    a.columns = a.columns.str.lower()
    b.columns = b.columns.str.lower()
    return a, b


def bilateral_mids(a, b):
    # "Bilateral on day one": exactly one originator on each side (footnote 2).
    o = b[b.orig == 1].groupby("dispnum").sidea.agg(n_a="sum", n="size")
    o["n_b"] = o.n - o.n_a
    keep = o[(o.n_a == 1) & (o.n_b == 1)].index

    orig = b[b.dispnum.isin(keep) & (b.orig == 1)]
    side_a = orig[orig.sidea == 1].set_index("dispnum")
    side_b = orig[orig.sidea == 0].set_index("dispnum")

    m = a[a.dispnum.isin(keep)].set_index("dispnum").copy()
    m["ccode_a"] = side_a.ccode
    m["ccode_b"] = side_b.ccode
    m["ccode1"] = m[["ccode_a", "ccode_b"]].min(axis=1)
    m["ccode2"] = m[["ccode_a", "ccode_b"]].max(axis=1)

    # Mutual use of force: each side (all participants on that side) reached hostlev >= 4.
    side_force = b.assign(force=b.hostlev >= 4).groupby(["dispnum", "sidea"]).force.max().unstack()
    m["mutual_force"] = (side_force[0] & side_force[1]).reindex(m.index).astype(int)

    m["reciprocated"] = m.recip.astype(int)
    m["use_force"] = (m.hostlev >= 4).astype(int)
    m["war"] = (m.hostlev == 5).astype(int)

    # Fatality codes: 0 none, 1 1-25, ..., 4 251-500, 5 501-999, 6 >999, -9 missing.
    # The paper's counts (536 with fatalities, 294 with >250) are only attainable if missing
    # (-9) counts as escalation, which is what Stata does once -9 is recoded to "." (missing
    # sorts above every number). We keep that coding as the primary version and the
    # substantively correct one (missing -> 0) as a robustness check.
    fat_stata = m.fatality.replace(-9, np.inf)
    m["fatal_any"] = (fat_stata > 0).astype(int)
    m["fatal_250"] = (fat_stata > 3).astype(int)
    m["fatal_any_strict"] = (m.fatality >= 1).astype(int)
    m["fatal_250_strict"] = (m.fatality >= 4).astype(int)
    m["fatality_missing"] = (m.fatality == -9).astype(int)

    # Territorial dispute: either originator seeks revision of territory (revtype 1).
    terr = orig.assign(t=(orig.revtype1 == 1) | (orig.revtype2 == 1)).groupby("dispnum").t.max()
    m["territorial_dispute"] = terr.reindex(m.index).astype(int)
    # Variant: any participant in the dispute (incl. later joiners) raises a territorial issue.
    terr_any = b.assign(t=(b.revtype1 == 1) | (b.revtype2 == 1)).groupby("dispnum").t.max()
    m["territorial_dispute_any"] = terr_any.reindex(m.index).astype(int)

    cols = ["ccode1", "ccode2", "ccode_a", "ccode_b", "styear", "stmon", "stday",
            "reciprocated", "use_force", "mutual_force", "fatal_any", "fatal_250", "war",
            "fatal_any_strict", "fatal_250_strict", "fatality_missing",
            "territorial_dispute", "territorial_dispute_any", "hostlev", "fatality", "numa", "numb"]
    return m[cols].reset_index()


def dyad_onsets(b):
    """All opposite-side participant pairs in every MID, with the year the pair's dispute began."""
    rows = []
    for disp, g in b.groupby("dispnum"):
        sa, sb = g[g.sidea == 1], g[g.sidea == 0]
        for _, x in sa.iterrows():
            for _, y in sb.iterrows():
                rows.append((disp, min(x.ccode, y.ccode), max(x.ccode, y.ccode), max(x.styear, y.styear)))
    return pd.DataFrame(rows, columns=["dispnum", "ccode1", "ccode2", "year"]).drop_duplicates()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a, b = load()
    m = bilateral_mids(a, b)
    m.to_csv(OUT / "mid_bilateral.csv", index=False)
    d = dyad_onsets(b)
    d.to_csv(OUT / "mid_dyad_onsets.csv", index=False)

    n = len(m)
    print(f"MID v3.02: {len(a)} disputes; bilateral on day one: {n} (paper: 1901 of 2331)")
    paper = {"reciprocated": 917, "use_force": 1394, "mutual_force": 576,
             "fatal_any": 536, "fatal_250": 294, "war": 104}
    for k, v in paper.items():
        print(f"  {k:13s} {m[k].sum():5d} ({m[k].mean():.0%})   paper {v} ({v / 1901:.0%})")
    print(f"  strict fatal_any {m.fatal_any_strict.sum()}, strict fatal_250 {m.fatal_250_strict.sum()}")
    print(f"dyadic MID onset pairs: {len(d)}")


if __name__ == "__main__":
    main()
