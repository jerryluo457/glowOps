"""Build the non-directed dyad-year dataset, 1919-1995, used in Braithwaite & Lemke (2011).

Data versions are the ones current around 2010 (recovered from the Wayback Machine) wherever
possible: COW system membership 2004.1, Major Powers 2004.1, NMC v4.0, Direct Contiguity v3.1,
Formal Alliances v3.03, Polity IV (p4v2009), MID v3.02, Huth & Allee (2002) territorial claims.
Rivalry uses Thompson's strategic rivalries as updated by Thompson, Sakuwa & Suhas (2021).

Output: data/clean/dyad_years.parquet (falls back to .csv.gz if parquet is unavailable).
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
ERA = RAW / "era2011"
CLEAN = ROOT / "data" / "clean"
Y0, Y1 = 1919, 1995

# Polity IV country codes that differ from COW codes over 1919-1995.
POLITY_TO_COW = {364: 365, 529: 530, 818: 816, 769: 770, 347: 345}


def dyad_key(df, a="ccode1", b="ccode2"):
    lo = df[[a, b]].min(axis=1).astype(int)
    hi = df[[a, b]].max(axis=1).astype(int)
    return lo, hi


def system_members():
    s = pd.read_csv(ERA / "system2004.csv")
    s.columns = ["abb", "ccode", "year", "version"]
    return s[["ccode", "year"]]


def frame(sysm):
    rows = []
    for y, g in sysm[sysm.year.between(Y0, Y1)].groupby("year"):
        cc = np.sort(g.ccode.unique())
        i, j = np.triu_indices(len(cc), k=1)
        rows.append(pd.DataFrame({"ccode1": cc[i], "ccode2": cc[j], "year": y}))
    return pd.concat(rows, ignore_index=True)


def add_majors(df):
    mj = pd.read_csv(ERA / "majors2004.csv")
    mj.columns = ["abb", "ccode", "sy", "sm", "sd", "ey", "em", "ed", "v"]
    years = [(r.ccode, y) for r in mj.itertuples() for y in range(r.sy, r.ey + 1)]
    major = set(years)
    m1 = np.array([(c, y) in major for c, y in zip(df.ccode1, df.year)])
    m2 = np.array([(c, y) in major for c, y in zip(df.ccode2, df.year)])
    df["both_minor"] = (~m1 & ~m2).astype(int)
    return df


def add_contiguity(df):
    c = pd.read_csv(ERA / "DirectContiguity310" / "DirectContiguity310" / "contdird.csv",
                    lineterminator="\r")
    c["ccode1"], c["ccode2"] = dyad_key(c, "state1no", "state2no")
    c = c.groupby(["ccode1", "ccode2", "year"]).conttype.min().reset_index()
    df = df.merge(c, on=["ccode1", "ccode2", "year"], how="left")
    df["conttype"] = df.conttype.fillna(6).astype(int)
    df["contig_land"] = (df.conttype == 1).astype(int)
    df["contig_24"] = (df.conttype <= 3).astype(int)   # land or <= 24 miles of water
    df["contig_150"] = (df.conttype <= 4).astype(int)  # land or <= 150 miles of water
    return df


def add_alliances(df, sysm):
    al = pd.read_csv(ERA / "Alliance_v3.03_dyadic" / "alliance_v303_annual_dyadic_single.csv")
    al["ccode1"], al["ccode2"] = dyad_key(al)
    al = al.groupby(["ccode1", "ccode2", "year"]).sstype.min().reset_index()
    df = df.merge(al, on=["ccode1", "ccode2", "year"], how="left")
    df["defense_pact"] = (df.sstype == 1).astype(int)

    # Alliance-portfolio tau-b with the system leader (UK 1919-1945, USA 1946-1995).
    # Ordinal alliance codes: defense 3, neutrality/non-aggression 2, entente 1, none 0;
    # a state is coded as having a defense pact with itself (Bueno de Mesquita 1975; EUGene).
    code = {1: 3, 2: 2, 3: 1}
    taub = {}
    for y in range(Y0, Y1 + 1):
        cc = np.sort(sysm[sysm.year == y].ccode.unique())
        idx = {c: k for k, c in enumerate(cc)}
        M = np.zeros((len(cc), len(cc)))
        np.fill_diagonal(M, 3)
        ay = al[al.year == y]
        for r in ay.itertuples():
            if r.ccode1 in idx and r.ccode2 in idx:
                M[idx[r.ccode1], idx[r.ccode2]] = M[idx[r.ccode2], idx[r.ccode1]] = code[r.sstype]
        leader = 200 if y <= 1945 else 2
        L = M[idx[leader]]
        for c in cc:
            t = kendalltau(M[idx[c]], L, variant="b").statistic
            taub[(c, y)] = 0.0 if np.isnan(t) else t
    t1 = np.array([taub.get((c, y), np.nan) for c, y in zip(df.ccode1, df.year)])
    t2 = np.array([taub.get((c, y), np.nan) for c, y in zip(df.ccode2, df.year)])
    df["taub1"], df["taub2"] = t1, t2
    df["joint_satisfaction"] = ((t1 > 0) & (t2 > 0)).astype(int)
    return df


def add_capabilities(df):
    n = pd.read_csv(RAW / "NMC_v4_0.csv")[["ccode", "year", "cinc"]]
    n.loc[n.cinc < 0, "cinc"] = np.nan
    c1 = df[["ccode1", "year"]].merge(n, left_on=["ccode1", "year"], right_on=["ccode", "year"], how="left").cinc
    c2 = df[["ccode2", "year"]].merge(n, left_on=["ccode2", "year"], right_on=["ccode", "year"], how="left").cinc
    c1, c2 = c1.to_numpy(), c2.to_numpy()
    tot = c1 + c2
    df["preponderance"] = np.where(tot > 0, np.maximum(c1, c2) / tot, np.nan)
    return df


def add_democracy(df):
    import pyreadstat
    p, _ = pyreadstat.read_sav(str(ERA / "p4v2009.sav"), usecols=["ccode", "year", "democ", "polity2"],
                               user_missing=True)
    p["ccode"] = p.ccode.replace(POLITY_TO_COW)
    # A few country-years appear twice in Polity (regime changes); keep the end-of-year value.
    p = p.dropna(subset=["ccode"]).groupby(["ccode", "year"]).last().reset_index()
    # -66 interruption, -77 interregnum, -88 transition.
    p["special"] = p.democ.isin([-66, -77, -88]).astype(int)
    p.loc[p.special == 1, "democ"] = np.nan
    for k in (1, 2):
        df = df.merge(p.rename(columns={"ccode": f"ccode{k}", "democ": f"democ{k}", "polity2": f"polity2_{k}",
                                        "special": f"polity_special{k}"}),
                      on=[f"ccode{k}", "year"], how="left")
    df["joint_democracy"] = ((df.democ1 >= 6) & (df.democ2 >= 6)).astype(int)
    df["joint_democracy_p2"] = ((df.polity2_1 >= 6) & (df.polity2_2 >= 6)).astype(int)
    df["polity_missing"] = (df.democ1.isna() | df.democ2.isna()).astype(int)
    # Stata treats missing (.) as larger than any number, so `democ1>=6 & democ2>=6` is true
    # for missing values. Two versions of that coding:
    #   _stata_special: only Polity's special codes (-66/-77/-88) recoded to missing
    #   _stata_all:     additionally, states absent from Polity are missing after the merge
    sp1, sp2 = df.polity_special1.fillna(0) == 1, df.polity_special2.fillna(0) == 1
    df["joint_democracy_stata_special"] = (((df.democ1 >= 6) | sp1) & ((df.democ2 >= 6) | sp2)).astype(int)
    df["joint_democracy_stata_all"] = (((df.democ1 >= 6) | df.democ1.isna())
                                       & ((df.democ2 >= 6) | df.democ2.isna())).astype(int)
    return df


def add_rivalry(df):
    r = pd.read_csv(RAW / "strategic_rivalries_list.csv").dropna(subset=["ccode1", "ccode2"])
    r["ccode1"], r["ccode2"] = dyad_key(r)
    r["start"] = r.start.fillna(1494)
    r["end"] = r.end.fillna(2020)
    riv = {(a, b, y) for a, b, s, e in r[["ccode1", "ccode2", "start", "end"]].itertuples(index=False)
           for y in range(int(max(s, Y0)), int(min(e, Y1)) + 1)}
    df["rivalry"] = [int((a, b, y) in riv) for a, b, y in zip(df.ccode1, df.ccode2, df.year)]
    return df


def add_territorial_claims(df):
    s = pd.read_stata(RAW / "huth_allee" / "statusquo.dta", convert_categoricals=False)
    s["year"] = pd.to_datetime(s.datetm).dt.year
    s["ccode1"], s["ccode2"] = dyad_key(s, "challenger", "target")
    # A claim is active from its first to its last recorded decision point (the data record a
    # decision at least every 12 months while the claim is alive).
    span = s.groupby(["dispno2", "ccode1", "ccode2"]).year.agg(["min", "max"]).reset_index()
    claims = {(r.ccode1, r.ccode2, y) for r in span.itertuples() for y in range(r.min, r.max + 1)}
    df["territorial_claim"] = [int((a, b, y) in claims) for a, b, y in zip(df.ccode1, df.ccode2, df.year)]
    return df


def add_mids(df, sysm):
    m = pd.read_csv(CLEAN / "mid_bilateral.csv")
    m = m[m.styear.between(Y0, Y1)].sort_values(["ccode1", "ccode2", "styear", "stmon", "stday", "dispnum"])
    m["n_bilateral_onsets"] = m.groupby(["ccode1", "ccode2", "styear"]).dispnum.transform("size")
    # One observation per dyad-year: the first bilateral MID to begin that year.
    first = m.drop_duplicates(["ccode1", "ccode2", "styear"]).rename(columns={"styear": "year"})
    keep = ["dispnum", "reciprocated", "use_force", "mutual_force", "fatal_any", "fatal_250", "war",
            "fatal_any_strict", "fatal_250_strict", "fatality_missing",
            "territorial_dispute", "territorial_dispute_any", "n_bilateral_onsets"]
    df = df.merge(first[["ccode1", "ccode2", "year"] + keep], on=["ccode1", "ccode2", "year"], how="left")
    df["mid_onset"] = df.dispnum.notna().astype(int)

    # Peace years: years since the dyad's last MID onset of any kind (counting from 1816 MIDs),
    # or since the later of the two states' first system year if they have never fought.
    on = pd.read_csv(CLEAN / "mid_dyad_onsets.csv")[["ccode1", "ccode2", "year"]].drop_duplicates()
    first_year = sysm.groupby("ccode").year.min()
    last = {}
    df = df.sort_values(["ccode1", "ccode2", "year"]).reset_index(drop=True)
    onset_years = on.groupby(["ccode1", "ccode2"]).year.apply(lambda s: np.sort(s.to_numpy())).to_dict()
    py = np.empty(len(df))
    for k, (a, b, y) in enumerate(zip(df.ccode1.to_numpy(), df.ccode2.to_numpy(), df.year.to_numpy())):
        ys = onset_years.get((a, b))
        prev = None
        if ys is not None:
            pos = np.searchsorted(ys, y)  # onsets strictly before y
            if pos > 0:
                prev = ys[pos - 1]
        start = max(first_year.get(a, 1816), first_year.get(b, 1816))
        py[k] = y - prev if prev is not None else y - start
    df["peace_years"] = py
    df["peace_years2"] = py ** 2
    df["peace_years3"] = py ** 3
    return df


def main():
    sysm = system_members()
    df = frame(sysm)
    print(f"frame: {len(df):,} non-directed dyad-years {Y0}-{Y1}")
    df = add_majors(df)
    df = add_contiguity(df)
    df = add_alliances(df, sysm)
    df = add_capabilities(df)
    df = add_democracy(df)
    df = add_rivalry(df)
    df = add_territorial_claims(df)
    df = add_mids(df, sysm)
    try:
        df.to_parquet(CLEAN / "dyad_years.parquet", index=False)
    except ImportError:
        df.to_csv(CLEAN / "dyad_years.csv.gz", index=False)

    print(f"MID onsets: {df.mid_onset.sum()}  (paper: 1352)")
    print(f"missing preponderance: {df.preponderance.isna().sum()}, missing polity: {df.polity_missing.sum()}")
    cols = ["contig_150", "contig_24", "rivalry", "both_minor", "joint_democracy", "joint_satisfaction",
            "territorial_claim", "preponderance", "defense_pact", "peace_years"]
    print(df[cols].describe().T[["mean", "min", "max"]])


if __name__ == "__main__":
    main()
