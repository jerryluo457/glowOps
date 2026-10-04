"""Replicate Table 1 (censored probit selection models) and Table 2 (predicted probabilities)
of Braithwaite & Lemke (2011).

Usage: python3 src/estimate.py [--spec NAME]   (specs defined in SPECS below)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from heckprob import Heckprob, bvn_cdf  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CLEAN = ROOT / "data" / "clean"
OUT = ROOT / "output"

DVS = ["reciprocated", "use_force", "mutual_force", "fatal_any", "fatal_250", "war"]
DV_LABELS = ["Reciprocated", "Use of Force", "Mutual Use of Force", "Fatalities > 0", "Fatalities > 250", "War"]

# Published Table 1 (coefficients; standard errors in parentheses in the article).
PAPER = json.loads((ROOT / "paper" / "table1_published.json").read_text())

# "main" reconstructs the coding the published table appears to rest on: Polity special codes
# (-66/-77/-88) and missing MID fatality levels both behave like Stata missing values (".")
# and therefore count as "democratic" and as "escalated". See REPORT.md for the evidence.
BASE = dict(
    contiguity="contig_150",
    joint_democracy="joint_democracy_stata_special",
    territorial_dispute="territorial_dispute",
    fatal_strict=False,
)
SPECS = {
    "main": BASE,
    "corrected": {**BASE, "joint_democracy": "joint_democracy", "fatal_strict": True},
    "jd_correct": {**BASE, "joint_democracy": "joint_democracy"},
    "jd_stata_all": {**BASE, "joint_democracy": "joint_democracy_stata_all"},
    "fatal_strict": {**BASE, "fatal_strict": True},
    "polity2": {**BASE, "joint_democracy": "joint_democracy_p2"},
    "contig_24": {**BASE, "contiguity": "contig_24"},
    "contig_land": {**BASE, "contiguity": "contig_land"},
    "terr_any_participant": {**BASE, "territorial_dispute": "territorial_dispute_any"},
}


def load():
    p = CLEAN / "dyad_years.parquet"
    df = pd.read_parquet(p) if p.exists() else pd.read_csv(CLEAN / "dyad_years.csv.gz")
    df["dyad"] = df.ccode1 * 1000 + df.ccode2
    return df


def design(df, spec):
    jd = spec["joint_democracy"]
    sel_cols = [spec["contiguity"], "rivalry", "both_minor", jd, "joint_satisfaction",
                "territorial_claim", "preponderance", "defense_pact",
                "peace_years", "peace_years2", "peace_years3"]
    sel_names = ["Contiguity", "Rivalry", "Both Minors", "Joint Democracy", "Joint Satisfaction",
                 "Territorial Claim", "Power Preponderance", "Allied",
                 "Peace Years", "Peace Years^2", "Peace Years^3"]
    out_cols = [jd, spec["territorial_dispute"], "joint_satisfaction", "preponderance", "defense_pact"]
    out_names = ["Joint Democracy", "Territorial Dispute", "Joint Satisfaction", "Power Preponderance", "Allied"]
    d = df.dropna(subset=sel_cols).copy()
    # Rescale the squared and cubic peace-years terms for numerical conditioning; their
    # coefficients are not reported and the rescaling leaves every other estimate unchanged.
    d["peace_years2"] = d.peace_years2 / 100
    d["peace_years3"] = d.peace_years3 / 10000
    Z = np.column_stack([d[sel_cols].to_numpy(float), np.ones(len(d))])
    X = np.column_stack([d[out_cols].fillna(0).to_numpy(float), np.ones(len(d))])
    return d, Z, X, sel_names + ["_cons"], out_names + ["_cons"]


def stars(p):
    return "***" if p < 0.01 else "**" if p < 0.05 else "*" if p < 0.10 else ""


def run(spec_name):
    spec = SPECS[spec_name]
    df = load()
    d, Z, X, zn, xn = design(df, spec)
    results = {}
    for dv in DVS:
        col = dv + "_strict" if spec["fatal_strict"] and dv.startswith("fatal") else dv
        y = d[col].fillna(0).to_numpy()
        m = Heckprob(y, X, d.mid_onset.to_numpy(), Z, cluster=d.dyad.to_numpy(), x_names=xn, z_names=zn)
        r = m.fit()
        results[dv] = (m, r)
        print(f"{dv:13s} converged={r.converged} N1={r.n:,} N2={r.n_selected} rho={r.rho:.3f} "
              f"chi2={r.wald_outcome()[0]:.2f}")
    return d, results, xn, zn


def table1(results, xn, zn):
    rows = []
    for eq, names in (("select", zn), ("outcome", xn)):
        for nm in names:
            if nm.startswith("Peace Years^"):
                continue
            key = f"{eq}:{nm}"
            row = {"equation": "Onset" if eq == "select" else "Escalation", "variable": nm}
            for dv in DVS:
                r = results[dv][1]
                t = r.table().loc[key]
                row[dv] = f"{t.coef:.2f}{stars(t.p)} ({t.se:.2f})"
            rows.append(row)
    for label, f in (("Rho", lambda r: f"{r.rho:.2f}{stars(r.table().loc['athrho'].p)}"),
                     ("Stage 1 N", lambda r: f"{r.n:,}"),
                     ("Stage 2 N", lambda r: f"{r.n_selected:,}"),
                     ("Model Chi2", lambda r: f"{r.wald_outcome()[0]:.2f}{stars(r.wald_outcome()[2])}")):
        rows.append({"equation": "", "variable": label, **{dv: f(results[dv][1]) for dv in DVS}})
    t = pd.DataFrame(rows)
    t.columns = ["equation", "variable"] + DV_LABELS
    return t


def comparison(results, xn, zn):
    """Long table: replicated vs published coefficient, sign/significance agreement."""
    rows = []
    for dv in DVS:
        r = results[dv][1]
        tab = r.table()
        for eq, names in (("select", zn), ("outcome", xn)):
            for nm in names:
                if nm.startswith("Peace Years^") or nm == "_cons":
                    continue
                pub = PAPER[dv].get(f"{eq}:{nm}")
                if pub is None:
                    continue
                t = tab.loc[f"{eq}:{nm}"]
                rows.append(dict(dv=dv, equation=eq, variable=nm, paper_coef=pub[0], paper_stars=pub[1],
                                 rep_coef=t.coef, rep_se=t.se, rep_p=t.p, rep_stars=stars(t.p)))
        pub = PAPER[dv]["rho"]
        rows.append(dict(dv=dv, equation="rho", variable="rho", paper_coef=pub[0], paper_stars=pub[1],
                         rep_coef=r.rho, rep_se=np.nan, rep_p=r.table().loc["athrho"].p,
                         rep_stars=stars(r.table().loc["athrho"].p)))
    c = pd.DataFrame(rows)
    sig = lambda s: s != ""  # noqa: E731
    c["same_sign"] = np.sign(c.paper_coef) == np.sign(c.rep_coef)
    c["same_significance"] = c.paper_stars.map(sig) == c.rep_stars.map(sig)
    c["match"] = np.where(c.paper_stars.map(sig), c.same_sign & c.same_significance, c.same_significance)
    return c


def predicted_probabilities(d, results, xn, zn):
    """Table 2: P(escalation | MID onset) = Phi2(xb, zg, rho) / Phi(zg) with all covariates at
    their means over the MID-onset sample. Of the variants tried (modes for binary covariates,
    full-sample means, averaged predictions, unconditional Phi(xb)), this one reproduces the
    article's "underlying probability" row most closely."""
    rows = {}
    for dv in DVS:
        m, r = results[dv]
        b, g = r.params[: r.k_out], r.params[r.k_out: r.k_out + r.k_sel]
        Zs, Xs = m.Z1, m.X1
        z0, x0 = Zs.mean(axis=0), Xs.mean(axis=0)

        def pcond(x, z):
            xb, zg = np.atleast_1d(x @ b), np.atleast_1d(z @ g)
            return float(bvn_cdf(xb, zg, r.rho)[0] / ndtr(zg)[0])

        out = {"Underlying": pcond(x0, z0)}
        # Pairs of (outcome-column, selection-column) that move together for each covariate.
        links = {"Joint Democracy": ("Joint Democracy", "Joint Democracy"),
                 "Joint Satisfaction": ("Joint Satisfaction", "Joint Satisfaction"),
                 "Territorial Dispute": ("Territorial Dispute", None),
                 "Allied": ("Allied", "Allied")}
        for lab, (xc, zc) in links.items():
            for v in (0, 1):
                x, z = x0.copy(), z0.copy()
                x[xn.index(xc)] = v
                if zc:
                    z[zn.index(zc)] = v
                out[f"{lab} = {v}"] = pcond(x, z)
        j, k = xn.index("Power Preponderance"), zn.index("Power Preponderance")
        mu, sd = Xs[:, j].mean(), Xs[:, j].std()
        for lab, v in (("Preponderance = mean - 1sd", mu - sd), ("Preponderance = mean", mu),
                       ("Preponderance = mean + 1sd", min(mu + sd, 1.0))):
            x, z = x0.copy(), z0.copy()
            x[j], z[k] = v, v
            out[lab] = pcond(x, z)
        rows[dv] = out
    t = pd.DataFrame(rows)
    t.columns = DV_LABELS
    return t


def descriptives(d):
    on = d[d.mid_onset == 1]
    rows = []
    for dv, lab in zip(DVS, DV_LABELS):
        rows.append(dict(measure=lab, n_escalated=int(on[dv].sum()), share=on[dv].mean()))
    corr = on[DVS].corr()
    return pd.DataFrame(rows), corr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", default="main", choices=list(SPECS))
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    d, results, xn, zn = run(args.spec)
    t1 = table1(results, xn, zn)
    cmp_ = comparison(results, xn, zn)
    sfx = "" if args.spec == "main" else f"_{args.spec}"
    t1.to_csv(OUT / f"table1{sfx}.csv", index=False)
    cmp_.to_csv(OUT / f"table1_comparison{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(t1.to_string(index=False))
    out_rows = cmp_[cmp_.equation != "select"]
    print(f"\nSign+significance agreement with paper: onset {cmp_[cmp_.equation == 'select'].match.mean():.0%}, "
          f"escalation {cmp_[cmp_.equation == 'outcome'].match.mean():.0%}, rho {cmp_[cmp_.equation == 'rho'].match.mean():.0%}")
    if args.spec == "main":
        t2 = predicted_probabilities(d, results, xn, zn)
        t2.to_csv(OUT / "table2.csv")
        print("\nTable 2 (predicted probabilities):")
        print(t2.round(3).to_string())
        desc, corr = descriptives(d)
        desc.to_csv(OUT / "descriptives.csv", index=False)
        corr.to_csv(OUT / "escalation_correlations.csv")
        print(desc.to_string(index=False))
        # Save parameter vectors for the report.
        allp = {dv: dict(zip(r.names, map(float, r.params))) for dv, (_, r) in results.items()}
        (OUT / "params_main.json").write_text(json.dumps(allp, indent=1))


if __name__ == "__main__":
    main()
