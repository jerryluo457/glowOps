"""Coefficient plot: escalation-equation estimates, published vs replicated (main spec)."""
import json
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output"
PAPER = json.loads((ROOT / "paper" / "table1_published.json").read_text())

DVS = ["reciprocated", "use_force", "mutual_force", "fatal_any", "fatal_250", "war"]
DV_LABELS = ["Reciprocated", "Use of Force", "Mutual Use of Force", "Fatalities > 0", "Fatalities > 250", "War"]
VARS = ["Territorial Dispute", "Joint Democracy", "Joint Satisfaction", "Power Preponderance", "Allied"]

SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
REP, PUB = "#2a78d6", "#eb6834"  # categorical slots 1 and 2 (validated pair)


def main(spec=""):
    c = pd.read_csv(OUT / f"table1_comparison{spec}.csv", keep_default_na=False)
    c = c[c.equation == "outcome"]
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "text.color": INK})
    fig, axes = plt.subplots(1, len(VARS), figsize=(13, 3.6), sharey=True, facecolor=SURFACE)
    for ax, var in zip(axes, VARS):
        ax.set_facecolor(SURFACE)
        ax.axvline(0, color=INK2, lw=1, zorder=1)
        ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
        for i, dv in enumerate(DVS):
            y = len(DVS) - 1 - i
            r = c[(c.dv == dv) & (c.variable == var)].iloc[0]
            pub_b, _, pub_se = PAPER[dv][f"outcome:{var}"]
            rb, rse = float(r.rep_coef), float(r.rep_se)
            ax.plot([pub_b - 1.96 * pub_se, pub_b + 1.96 * pub_se], [y + 0.17] * 2, color=PUB, lw=2,
                    solid_capstyle="round", zorder=2)
            ax.plot(pub_b, y + 0.17, marker="D", ms=6, mfc=SURFACE, mec=PUB, mew=2, zorder=3)
            ax.plot([rb - 1.96 * rse, rb + 1.96 * rse], [y - 0.17] * 2, color=REP, lw=2,
                    solid_capstyle="round", zorder=2)
            ax.plot(rb, y - 0.17, marker="o", ms=7, color=REP, mec=SURFACE, mew=1.5, zorder=3)
        ax.set_title(var, fontsize=10, color=INK, loc="left")
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(axis="y", length=0)
    axes[0].set_yticks(range(len(DVS)))
    axes[0].set_yticklabels(DV_LABELS[::-1])
    handles = [plt.Line2D([], [], color=PUB, marker="D", mfc=SURFACE, mew=2, lw=2, label="Published (Table 1)"),
               plt.Line2D([], [], color=REP, marker="o", lw=2, label="Replication")]
    fig.legend(handles=handles, loc="upper right", frameon=False, ncol=2)
    fig.suptitle("Escalation equation: coefficients with 95% CIs (cluster-robust)", x=0.01, ha="left",
                 fontsize=11, color=INK)
    fig.supxlabel("Probit coefficient", fontsize=9, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = OUT / f"fig_escalation_coefficients{spec}.png"
    fig.savefig(out, dpi=180, facecolor=SURFACE)
    print(out)


if __name__ == "__main__":
    main()
