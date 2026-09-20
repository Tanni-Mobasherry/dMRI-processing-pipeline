from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Desktop" / "start-analysis"

INPUT = (
    BASE
    / "4-dlpfc-median-iqr"
    / "DLPFC_04_all_edge_group_statistics.csv"
)

OUT_DIR = BASE / "5-fixed-dlpfc-edges"
OUT_DIR.mkdir(exist_ok=True)

PNG_OUT = OUT_DIR / "DLPFC_05_edge_selection_plot.png"
PDF_OUT = OUT_DIR / "DLPFC_05_edge_selection_plot.pdf"


# ============================================================
# READ DATA
# ============================================================

df = pd.read_csv(INPUT)

required = [
    "Connected_Node_ID",
    "Connected_Node",
    "Median_Delta_FBC",
    "Lower_Threshold",
    "Upper_Threshold"
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise RuntimeError(
        f"Missing required columns: {missing}"
    )

if len(df) != 163:
    raise RuntimeError(
        f"Expected 163 DLPFC edges, found {len(df)}"
    )


# ============================================================
# NUMERIC VALIDATION
# ============================================================

numeric_cols = [
    "Median_Delta_FBC",
    "Lower_Threshold",
    "Upper_Threshold"
]

for col in numeric_cols:
    df[col] = pd.to_numeric(
        df[col],
        errors="raise"
    )

if not np.isfinite(
    df[numeric_cols].to_numpy()
).all():
    raise RuntimeError(
        "Non-finite values found in numerical columns."
    )


# ============================================================
# READ THRESHOLDS DIRECTLY FROM STEP 4
# ============================================================

lower_values = np.unique(
    np.round(
        df["Lower_Threshold"].to_numpy(),
        15
    )
)

upper_values = np.unique(
    np.round(
        df["Upper_Threshold"].to_numpy(),
        15
    )
)

if len(lower_values) != 1:
    raise RuntimeError(
        "More than one lower threshold found."
    )

if len(upper_values) != 1:
    raise RuntimeError(
        "More than one upper threshold found."
    )

LOWER = float(lower_values[0])
UPPER = float(upper_values[0])


# ============================================================
# READ ACTUAL MINIMUM / MAXIMUM FROM YOUR DATA
# ============================================================

DATA_MIN = float(
    df["Median_Delta_FBC"].min()
)

DATA_MAX = float(
    df["Median_Delta_FBC"].max()
)

DATA_RANGE = DATA_MAX - DATA_MIN

if DATA_RANGE <= 0:
    raise RuntimeError(
        "Median_Delta_FBC has no usable range."
    )


# ------------------------------------------------------------
# Add only a small amount of visual padding
# ------------------------------------------------------------

Y_PADDING = DATA_RANGE * 0.05

Y_MIN = DATA_MIN - Y_PADDING
Y_MAX = DATA_MAX + Y_PADDING


# ============================================================
# CLASSIFY EDGES
# ============================================================

df["Selection"] = np.select(
    [
        df["Median_Delta_FBC"] > UPPER,
        df["Median_Delta_FBC"] < LOWER
    ],
    [
        "Positive",
        "Negative"
    ],
    default="Within"
)

n_positive = int(
    (df["Selection"] == "Positive").sum()
)

n_negative = int(
    (df["Selection"] == "Negative").sum()
)

n_within = int(
    (df["Selection"] == "Within").sum()
)

n_selected = n_positive + n_negative


# ============================================================
# SORT BY ACTUAL MEDIAN ΔFBC
# ============================================================

plot_df = (
    df
    .sort_values(
        "Median_Delta_FBC",
        ascending=True
    )
    .reset_index(drop=True)
    .copy()
)

plot_df["Rank"] = np.arange(
    1,
    len(plot_df) + 1
)


# ============================================================
# VALIDATION SUMMARY
# ============================================================

print("\n" + "=" * 72)
print("DLPFC EDGE-SELECTION PLOT — DATA CHECK")
print("=" * 72)

print(
    f"\nInput file:\n{INPUT}"
)

print(
    f"\nNumber of DLPFC edges       : {len(df)}"
)

print(
    f"Actual minimum Median ΔFBC  : {DATA_MIN:.8f}"
)

print(
    f"Actual maximum Median ΔFBC  : {DATA_MAX:.8f}"
)

print(
    f"\nLower threshold             : {LOWER:.8f}"
)

print(
    f"Upper threshold             : {UPPER:.8f}"
)

print(
    f"\nBelow lower threshold       : {n_negative}"
)

print(
    f"Within thresholds           : {n_within}"
)

print(
    f"Above upper threshold       : {n_positive}"
)

print(
    f"Total selected              : {n_selected}"
)

print(
    f"\nY-axis minimum with padding : {Y_MIN:.8f}"
)

print(
    f"Y-axis maximum with padding : {Y_MAX:.8f}"
)


# ------------------------------------------------------------
# Expected Step 4 result
# ------------------------------------------------------------

if n_selected != 33:
    raise RuntimeError(
        f"Expected 33 selected edges, "
        f"but found {n_selected}."
    )

if n_positive != 29:
    raise RuntimeError(
        f"Expected 29 positive extreme edges, "
        f"but found {n_positive}."
    )

if n_negative != 4:
    raise RuntimeError(
        f"Expected 4 negative extreme edges, "
        f"but found {n_negative}."
    )


# ============================================================
# CREATE FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(12.5, 7.2)
)


# ============================================================
# MASKS
# ============================================================

within = (
    plot_df["Selection"]
    ==
    "Within"
)

positive = (
    plot_df["Selection"]
    ==
    "Positive"
)

negative = (
    plot_df["Selection"]
    ==
    "Negative"
)


# ============================================================
# SCATTER POINTS
# ============================================================

ax.scatter(
    plot_df.loc[
        within,
        "Rank"
    ],
    plot_df.loc[
        within,
        "Median_Delta_FBC"
    ],
    s=38,
    color="#9A9A9A",
    edgecolors="none",
    alpha=0.88,
    label=f"Within threshold (n = {n_within})",
    zorder=3
)

ax.scatter(
    plot_df.loc[
        positive,
        "Rank"
    ],
    plot_df.loc[
        positive,
        "Median_Delta_FBC"
    ],
    s=48,
    color="#D7192D",
    edgecolors="none",
    alpha=0.95,
    label=f"Above upper threshold (n = {n_positive})",
    zorder=4
)

ax.scatter(
    plot_df.loc[
        negative,
        "Rank"
    ],
    plot_df.loc[
        negative,
        "Median_Delta_FBC"
    ],
    s=48,
    color="#1565C0",
    edgecolors="none",
    alpha=0.95,
    label=f"Below lower threshold (n = {n_negative})",
    zorder=4
)


# ============================================================
# THRESHOLD LINES
# ============================================================

ax.axhline(
    UPPER,
    color="#D7192D",
    linestyle=(0, (5, 4)),
    linewidth=1.8,
    zorder=2
)

ax.axhline(
    LOWER,
    color="#1565C0",
    linestyle=(0, (5, 4)),
    linewidth=1.8,
    zorder=2
)


# ============================================================
# ZERO REFERENCE LINE
# ============================================================

ax.axhline(
    0,
    color="#B8B8B8",
    linewidth=1.0,
    zorder=1
)


# ============================================================
# Y AXIS — BASED ENTIRELY ON ACTUAL DATA
# ============================================================

ax.set_ylim(
    Y_MIN,
    Y_MAX
)


# ============================================================
# THRESHOLD LABELS
#
# Put them on the right, but deliberately offset them
# vertically so they remain readable even though the
# thresholds are close together relative to the full
# data range.
# ============================================================

ax.annotate(
    f"Upper threshold  {UPPER:+.6f}",
    xy=(1.0, UPPER),
    xycoords=("axes fraction", "data"),
    xytext=(10, 12),
    textcoords="offset points",
    ha="left",
    va="bottom",
    fontsize=11,
    color="#D7192D",
    fontweight="bold",
    annotation_clip=False
)

ax.annotate(
    f"Lower threshold  {LOWER:+.6f}",
    xy=(1.0, LOWER),
    xycoords=("axes fraction", "data"),
    xytext=(10, -12),
    textcoords="offset points",
    ha="left",
    va="top",
    fontsize=11,
    color="#1565C0",
    fontweight="bold",
    annotation_clip=False
)


# ============================================================
# TITLE + AXIS LABELS
# ============================================================

ax.set_title(
    "Median Longitudinal FBC Change for 163 DLPFC Connections",
    fontsize=18,
    fontweight="bold",
    pad=18
)

ax.set_xlabel(
    "DLPFC connections (n = 163)",
    fontsize=14,
    labelpad=12
)

ax.set_ylabel(
    "Median ΔFBC",
    fontsize=14,
    labelpad=12
)


# ============================================================
# X AXIS
# ============================================================

ax.set_xlim(
    0,
    166
)

ax.set_xticks(
    [
        1,
        20,
        40,
        60,
        80,
        100,
        120,
        140,
        163
    ]
)

ax.tick_params(
    axis="both",
    labelsize=11
)


# ============================================================
# GRID
# ============================================================

ax.grid(
    axis="y",
    color="#E1E1E1",
    linewidth=0.9
)

ax.grid(
    axis="x",
    visible=False
)


# ============================================================
# SPINES
# ============================================================

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.spines["left"].set_linewidth(1.2)
ax.spines["bottom"].set_linewidth(1.2)


# ============================================================
# LEGEND
# ============================================================

legend = ax.legend(
    loc="upper left",
    bbox_to_anchor=(0.02, -0.15),
    frameon=True,
    fontsize=11,
    borderpad=0.8,
    handletextpad=0.6
)

legend.get_frame().set_edgecolor(
    "#D0D0D0"
)

legend.get_frame().set_linewidth(
    1.0
)

legend.get_frame().set_facecolor(
    "white"
)


# ============================================================
# RESULT BOX
# ============================================================

percentage = (
    n_selected
    /
    len(df)
    *
    100
)

result_text = (
    f"{n_selected} of {len(df)} connections "
    f"({percentage:.1f}%)\n"
    "showed extreme longitudinal changes\n"
    f"({n_positive} positive | {n_negative} negative)"
)

ax.text(
    0.98,
    -0.16,
    result_text,
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=11.5,
    fontweight="bold",
    bbox=dict(
        boxstyle="round,pad=0.65",
        facecolor="#FFF1F2",
        edgecolor="#E7B8BD",
        linewidth=1.0
    )
)


# ============================================================
# FOOTNOTE
# ============================================================

fig.text(
    0.5,
    0.025,
    (
        "Edge selection was performed without "
        "responder/non-responder information."
    ),
    ha="center",
    va="center",
    fontsize=10.5,
    color="#3E5268"
)


# ============================================================
# LAYOUT
# ============================================================

plt.subplots_adjust(
    left=0.10,
    right=0.82,
    top=0.88,
    bottom=0.30
)


# ============================================================
# SAVE
# ============================================================

fig.savefig(
    PNG_OUT,
    dpi=400,
    bbox_inches="tight",
    facecolor="white"
)

fig.savefig(
    PDF_OUT,
    bbox_inches="tight",
    facecolor="white"
)

plt.close(fig)


# ============================================================
# DONE
# ============================================================

print("\n" + "=" * 72)
print("SAVED")
print("=" * 72)

print(
    f"\nPNG:\n{PNG_OUT}"
)

print(
    f"\nPDF:\n{PDF_OUT}"
)

print(
    "\n✓ Minimum and maximum were read directly "
    "from Median_Delta_FBC."
)

print(
    "✓ Y-axis was calculated automatically "
    "from the actual data range."
)

print(
    "✓ Thresholds were read directly from "
    "the Step 4 output."
)

print(
    "✓ No responder/non-responder information was used."
)
