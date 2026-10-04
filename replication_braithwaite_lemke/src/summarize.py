"""Summarize agreement with the published Table 1 across all specifications."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
SPECS = ["main", "corrected", "jd_correct", "jd_stata_all", "fatal_strict", "polity2",
         "contig_24", "contig_land", "terr_any_participant"]
DVS = ["reciprocated", "use_force", "mutual_force", "fatal_any", "fatal_250", "war"]


def load(spec):
    sfx = "" if spec == "main" else f"_{spec}"
    c = pd.read_csv(OUT / f"table1_comparison{sfx}.csv", keep_default_na=False)
    c["match"] = c.match.astype(str).str.lower() == "true"
    return c


def main():
    rows, terr = [], []
    for s in SPECS:
        try:
            c = load(s)
        except FileNotFoundError:
            continue
        row = {"spec": s}
        for eq, lab in (("select", "onset"), ("outcome", "escalation"), ("rho", "rho")):
            x = c[c.equation == eq]
            row[lab] = f"{x.match.sum()}/{len(x)}"
        rows.append(row)
        t = c[(c.equation == "outcome") & (c.variable == "Territorial Dispute")].set_index("dv")
        terr.append({"spec": s, **{dv: f"{float(t.loc[dv, 'rep_coef']):.2f}{t.loc[dv, 'rep_stars']}" for dv in DVS}})
    agree = pd.DataFrame(rows)
    terr = pd.DataFrame(terr)
    paper = load("main")
    p = paper[(paper.equation == "outcome") & (paper.variable == "Territorial Dispute")].set_index("dv")
    terr.loc[len(terr)] = {"spec": "PUBLISHED", **{dv: f"{float(p.loc[dv, 'paper_coef']):.2f}{p.loc[dv, 'paper_stars']}"
                                                  for dv in DVS}}
    agree.to_csv(OUT / "agreement_by_spec.csv", index=False)
    terr.to_csv(OUT / "territorial_dispute_by_spec.csv", index=False)
    print(agree.to_string(index=False))
    print()
    print(terr.to_string(index=False))


if __name__ == "__main__":
    main()
