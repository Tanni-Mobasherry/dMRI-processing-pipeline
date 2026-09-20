"""Paired pre/post DASS-21 depression figure in the style of aim.png.

Bar (mean +SD) + half-violin outline + stacked individual points joined by
participant, coloured by clinical responder status.
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats
from scipy.stats import gaussian_kde

RESPONDER = "#138A9A"
NONRESPONDER = "#D96B5F"

BAR_FILL = "#CBD8E8"
INK = "#1A1A1A"

CSV = "clinical_response.csv"
YLABEL = "DASS-21 depression score"
XLABELS = ["Pre-rTMS", "Post-rTMS"]

# ---------------------------------------------------------------- geometry --
BAR_W = 0.26          # bar width
VIOLIN_X = 0.24       # offset of violin baseline from bar centre
DOT_DX = 0.050        # horizontal spacing between stacked points
DOT_X0 = 0.012        # first point offset from violin baseline
CENTRES = [0.0, 1.00]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.9,
    "hatch.linewidth": 0.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def stack_positions(values, base_x, is_resp):
    """Dot-plot layout: equal values stack rightwards from the violin baseline.

    Within a row, responders are placed first so the two groups read as blocks
    rather than interleaving.
    """
    x = np.empty(len(values), dtype=float)
    seen = {}
    order = np.lexsort((np.arange(len(values)), ~is_resp))
    for i in order:
        k = round(float(values[i]), 6)
        n = seen.get(k, 0)
        x[i] = base_x + DOT_X0 + n * DOT_DX
        seen[k] = n + 1
    return x, max(seen.values())


def violin_curve(values, ylim):
    """KDE outline for one group; returned unnormalised so groups share a scale."""
    kde = gaussian_kde(values, bw_method=0.45)
    lo = max(ylim[0], values.min() - 6)
    hi = min(ylim[1] - 6, values.max() + 6)
    grid = np.linspace(lo, hi, 400)
    return grid, kde(grid)


def p_stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def main():
    df = pd.read_csv(CSV)
    pre = df["BL_MRIdepDASS"].to_numpy(float)
    post = df["FU_MRIdepDASS"].to_numpy(float)
    is_resp = (df["Responder_Status"] == "Responder").to_numpy()

    t_stat, p_val = stats.ttest_rel(pre, post)
    means = [pre.mean(), post.mean()]
    sds = [pre.std(ddof=1), post.std(ddof=1)]

    ylim = (0, 50)

    fig, ax = plt.subplots(figsize=(5.6, 4.0))

    # ---- bars: mean with upward SD whisker ---------------------------------
    hatches = [None, "///"]
    for c, m, sd, hatch in zip(CENTRES, means, sds, hatches):
        ax.bar(c, m, width=BAR_W, facecolor=BAR_FILL, edgecolor=INK,
               linewidth=0.9, hatch=hatch, zorder=1)
        ax.plot([c, c], [m, m + sd], color=INK, lw=0.9, zorder=3)
        ax.plot([c - BAR_W * 0.36, c + BAR_W * 0.36], [m + sd] * 2,
                color=INK, lw=0.9, zorder=3)

    # ---- individual points, stacked --------------------------------------
    x_pre, n_pre = stack_positions(pre, CENTRES[0] + VIOLIN_X, is_resp)
    x_post, n_post = stack_positions(post, CENTRES[1] + VIOLIN_X, is_resp)

    # ---- half violins ------------------------------------------------------
    # One shared density scale for both groups, tied to the dot spacing so the
    # outline tracks the height of the stacks rather than floating beside them.
    curves = [violin_curve(v, ylim) for v in (pre, post)]
    dmax = max(d.max() for _, d in curves)
    vscale = max(n_pre, n_post) * DOT_DX / dmax
    for c, (grid, dens) in zip(CENTRES, curves):
        ax.plot(c + VIOLIN_X + dens * vscale, grid, color=INK,
                lw=0.9, solid_capstyle="round", zorder=2, clip_on=False)
    violin_w = dmax * vscale

    # participant lines, drawn under the markers
    for i in range(len(df)):
        col = RESPONDER if is_resp[i] else NONRESPONDER
        ax.plot([x_pre[i], x_post[i]], [pre[i], post[i]],
                color=col, lw=0.6, alpha=0.45, zorder=4,
                solid_capstyle="round")

    for mask, col, marker, ms in [
        (is_resp, RESPONDER, "o", 3.6),
        (~is_resp, NONRESPONDER, "s", 3.2),
    ]:
        for xs, ys in [(x_pre, pre), (x_post, post)]:
            ax.plot(xs[mask], ys[mask], linestyle="none", marker=marker,
                    markersize=ms, markerfacecolor=col,
                    markeredgecolor="white", markeredgewidth=0.5,
                    zorder=5, clip_on=False)

    # ---- significance bracket ---------------------------------------------
    y_br = 46.5
    x0, x1 = CENTRES[0], CENTRES[1] + VIOLIN_X + violin_w * 0.55
    ax.plot([x0, x1], [y_br, y_br], color=INK, lw=1.0, zorder=6)
    ax.text((x0 + x1) / 2, y_br + 0.4, p_stars(p_val), ha="center",
            va="bottom", fontsize=11, color=INK)

    # ---- axes --------------------------------------------------------------
    ax.set_ylim(*ylim)
    ax.set_yticks(np.arange(0, 51, 10))
    ax.set_ylabel(YLABEL, fontsize=10)
    ax.set_xlim(CENTRES[0] - BAR_W * 0.95,
                CENTRES[1] + VIOLIN_X + max(violin_w, n_post * DOT_DX) + 0.06)
    ax.set_xticks(CENTRES)
    ax.set_xticklabels(XLABELS, fontsize=10)
    ax.tick_params(axis="both", direction="out", length=3.5, width=0.9,
                   colors=INK)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK)
    ax.spines["left"].set_bounds(0, 50)

    n_r, n_n = int(is_resp.sum()), int((~is_resp).sum())
    handles = [
        Line2D([], [], linestyle="none", marker="o", markersize=4.4,
               markerfacecolor=RESPONDER, markeredgecolor="white",
               markeredgewidth=0.5, label=f"Responder (n = {n_r})"),
        Line2D([], [], linestyle="none", marker="s", markersize=4.0,
               markerfacecolor=NONRESPONDER, markeredgecolor="white",
               markeredgewidth=0.5, label=f"Non-responder (n = {n_n})"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.11),
              ncol=2, frameon=False, fontsize=9, handletextpad=0.4,
              columnspacing=2.0)

    fig.tight_layout()
    fig.savefig("clinical_response_aim.png", dpi=400, bbox_inches="tight",
                facecolor="white")
    fig.savefig("clinical_response_aim.pdf", bbox_inches="tight",
                facecolor="white")
    print(f"pre  {means[0]:.2f} ± {sds[0]:.2f}")
    print(f"post {means[1]:.2f} ± {sds[1]:.2f}")
    print(f"paired t({len(df) - 1}) = {t_stat:.2f}, p = {p_val:.2e} "
          f"-> {p_stars(p_val)}")


if __name__ == "__main__":
    main()
