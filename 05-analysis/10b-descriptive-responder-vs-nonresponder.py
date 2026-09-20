#This code descriptively inspect responder vs non-responder ΔFBC
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
    / "6-dlpfc-clinical-analysis"
    / "DLPFC_06A_fixed33_clinical_long.csv"
)

OUT_DIR = (
    BASE
    / "6-dlpfc-clinical-analysis"
    / "6B-responder-descriptives"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_OUT = OUT_DIR / "DLPFC_06B_responder_vs_nonresponder_summary.csv"
LONG_OUT = OUT_DIR / "DLPFC_06B_delta_long.csv"
RANKED_OUT = OUT_DIR / "DLPFC_06B_ranked_group_difference.csv"

FIGURE_OUT = OUT_DIR / "DLPFC_06B_delta_group_difference.png"
FIGURE_PDF = OUT_DIR / "DLPFC_06B_delta_group_difference.pdf"


# ============================================================
# SETTINGS
# ============================================================

EXPECTED_EDGES = 33
EXPECTED_RESPONDERS = 13
EXPECTED_NONRESPONDERS = 23


# ============================================================
# 1. READ DATA
# ============================================================

print("\n" + "=" * 80)
print("STEP 6B — RESPONDER VS NON-RESPONDER DESCRIPTIVE ANALYSIS")
print("=" * 80)

df = pd.read_csv(INPUT)

print(f"\nInput rows: {len(df)}")
print(f"Participants: {df['ID'].nunique()}")
print(f"Unique edges: {df['Connected_Node_ID'].nunique()}")


# ============================================================
# 2. REQUIRED COLUMNS
# ============================================================

required = [
    "ID",
    "Connected_Node_ID",
    "Connected_Node",
    "BL_FBC",
    "FU_FBC",
    "Delta_FBC",
    "Responder_Status"
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise RuntimeError(
        f"Missing required columns: {missing}"
    )


# ============================================================
# 3. KEEP ELIGIBLE PARTICIPANTS
# ============================================================

if "Eligible_for_Responder_Analysis" in df.columns:

    eligible = df["Eligible_for_Responder_Analysis"]

    # Handle bool or text representation safely
    if eligible.dtype == object:
        eligible = (
            eligible.astype(str)
            .str.strip()
            .str.lower()
            .map({
                "true": True,
                "false": False
            })
        )

    df = df.loc[eligible == True].copy()


# ============================================================
# 4. BASIC VALIDATION
# ============================================================

valid_groups = {
    "Responder",
    "Non-responder"
}

found_groups = set(
    df["Responder_Status"]
    .dropna()
    .unique()
)

if found_groups != valid_groups:
    raise RuntimeError(
        f"Unexpected responder groups: {found_groups}"
    )


# Participant counts by group
subject_groups = (
    df[
        ["ID", "Responder_Status"]
    ]
    .drop_duplicates()
)

# Ensure one group per subject
group_count_per_subject = (
    subject_groups
    .groupby("ID")
    ["Responder_Status"]
    .nunique()
)

if (group_count_per_subject != 1).any():
    raise RuntimeError(
        "At least one participant has more than one "
        "responder classification."
    )


n_resp = (
    subject_groups.loc[
        subject_groups["Responder_Status"] == "Responder",
        "ID"
    ]
    .nunique()
)

n_nonresp = (
    subject_groups.loc[
        subject_groups["Responder_Status"] == "Non-responder",
        "ID"
    ]
    .nunique()
)


print("\nGroup counts:")
print(f"  Responders     : {n_resp}")
print(f"  Non-responders : {n_nonresp}")


if n_resp != EXPECTED_RESPONDERS:
    raise RuntimeError(
        f"Expected {EXPECTED_RESPONDERS} responders, "
        f"found {n_resp}."
    )

if n_nonresp != EXPECTED_NONRESPONDERS:
    raise RuntimeError(
        f"Expected {EXPECTED_NONRESPONDERS} non-responders, "
        f"found {n_nonresp}."
    )


# ============================================================
# 5. VALIDATE FIXED 33 EDGES
# ============================================================

if df["Connected_Node_ID"].nunique() != EXPECTED_EDGES:
    raise RuntimeError(
        f"Expected {EXPECTED_EDGES} unique edges, "
        f"found {df['Connected_Node_ID'].nunique()}."
    )


edge_counts = (
    df.groupby("ID")
    ["Connected_Node_ID"]
    .nunique()
)

bad = edge_counts[
    edge_counts != EXPECTED_EDGES
]

if len(bad):
    raise RuntimeError(
        "Some participants do not have exactly "
        f"{EXPECTED_EDGES} edges:\n{bad}"
    )


# Check exact same edge set across participants
reference_id = sorted(df["ID"].unique())[0]

reference_edges = set(
    df.loc[
        df["ID"] == reference_id,
        "Connected_Node_ID"
    ].astype(int)
)

for sid in sorted(df["ID"].unique()):

    current_edges = set(
        df.loc[
            df["ID"] == sid,
            "Connected_Node_ID"
        ].astype(int)
    )

    if current_edges != reference_edges:
        raise RuntimeError(
            f"{sid}: edge set differs from fixed 33."
        )


print("✓ Same fixed 33 edges retained for all participants.")


# ============================================================
# 6. VALIDATE DELTA
# ============================================================

calculated_delta = (
    pd.to_numeric(df["FU_FBC"], errors="raise")
    -
    pd.to_numeric(df["BL_FBC"], errors="raise")
)

stored_delta = pd.to_numeric(
    df["Delta_FBC"],
    errors="raise"
)

if not np.allclose(
    calculated_delta,
    stored_delta,
    rtol=1e-12,
    atol=1e-15
):
    raise RuntimeError(
        "Delta validation failed."
    )

print("✓ Delta_FBC = FU_FBC - BL_FBC validated.")


# ============================================================
# 7. CHECK MISSING / NON-FINITE VALUES
# ============================================================

for col in [
    "BL_FBC",
    "FU_FBC",
    "Delta_FBC"
]:
    values = pd.to_numeric(
        df[col],
        errors="coerce"
    )

    if values.isna().any():
        raise RuntimeError(
            f"Missing/non-numeric values found in {col}."
        )

    if not np.isfinite(values).all():
        raise RuntimeError(
            f"Non-finite values found in {col}."
        )


# ============================================================
# 8. SAVE CLEAN LONG-FORM DELTA DATA
# ============================================================

long_cols = [
    "ID",
    "Responder_Status",
    "Connected_Node_ID",
    "Connected_Node",
    "BL_FBC",
    "FU_FBC",
    "Delta_FBC"
]

for optional in [
    "Age",
    "Sex",
    "Gender",
    "BL_MRIdepDASS",
    "FU_MRIdepDASS",
    "DASS_Reduction_Calculated",
    "DASS_Improvement_Percent"
]:
    if optional in df.columns:
        long_cols.append(optional)


delta_long = (
    df[long_cols]
    .sort_values(
        [
            "Connected_Node_ID",
            "Responder_Status",
            "ID"
        ]
    )
    .reset_index(drop=True)
)

delta_long.to_csv(
    LONG_OUT,
    index=False
)


# ============================================================
# 9. DESCRIPTIVE STATISTICS FOR EACH EDGE AND GROUP
# ============================================================

rows = []

for node_id in sorted(
    df["Connected_Node_ID"].astype(int).unique()
):

    edge = df.loc[
        df["Connected_Node_ID"].astype(int) == node_id
    ].copy()

    names = edge["Connected_Node"].dropna().unique()

    if len(names) != 1:
        raise RuntimeError(
            f"Node {node_id}: inconsistent anatomical names: "
            f"{names}"
        )

    node_name = names[0]

    result = {
        "Connected_Node_ID": node_id,
        "Connected_Node": node_name,
        "Connected_Node_ID_Name":
            f"{node_id}|{node_name}"
    }

    for group, prefix in [
        ("Responder", "Responder"),
        ("Non-responder", "NonResponder")
    ]:

        g = edge.loc[
            edge["Responder_Status"] == group
        ].copy()

        bl = g["BL_FBC"].astype(float)
        fu = g["FU_FBC"].astype(float)
        delta = g["Delta_FBC"].astype(float)

        result[f"{prefix}_N"] = len(g)

        # -------------------------------
        # Baseline
        # -------------------------------

        result[f"{prefix}_BL_Mean"] = bl.mean()
        result[f"{prefix}_BL_SD"] = bl.std(ddof=1)
        result[f"{prefix}_BL_Median"] = bl.median()

        # -------------------------------
        # Follow-up
        # -------------------------------

        result[f"{prefix}_FU_Mean"] = fu.mean()
        result[f"{prefix}_FU_SD"] = fu.std(ddof=1)
        result[f"{prefix}_FU_Median"] = fu.median()

        # -------------------------------
        # Change: FU - BL
        # -------------------------------

        result[f"{prefix}_Delta_Mean"] = delta.mean()
        result[f"{prefix}_Delta_SD"] = delta.std(ddof=1)

        result[f"{prefix}_Delta_Median"] = delta.median()

        q1 = delta.quantile(0.25)
        q3 = delta.quantile(0.75)

        result[f"{prefix}_Delta_Q1"] = q1
        result[f"{prefix}_Delta_Q3"] = q3
        result[f"{prefix}_Delta_IQR"] = q3 - q1

        result[f"{prefix}_N_Increase"] = int(
            (delta > 0).sum()
        )

        result[f"{prefix}_N_Decrease"] = int(
            (delta < 0).sum()
        )

        result[f"{prefix}_N_NoChange"] = int(
            np.isclose(
                delta,
                0,
                rtol=0,
                atol=1e-15
            ).sum()
        )


    # ========================================================
    # 10. DESCRIPTIVE BETWEEN-GROUP DIFFERENCE IN CHANGE
    #
    # Positive:
    # responders changed more positively than non-responders
    #
    # Negative:
    # responders changed more negatively than non-responders
    #
    # This is descriptive only — NOT a statistical test.
    # ========================================================

    result[
        "Median_Delta_Difference_Resp_minus_NonResp"
    ] = (
        result["Responder_Delta_Median"]
        -
        result["NonResponder_Delta_Median"]
    )

    result[
        "Mean_Delta_Difference_Resp_minus_NonResp"
    ] = (
        result["Responder_Delta_Mean"]
        -
        result["NonResponder_Delta_Mean"]
    )


    # Directional description
    med_diff = result[
        "Median_Delta_Difference_Resp_minus_NonResp"
    ]

    if med_diff > 0:
        direction = "Responder_more_positive"

    elif med_diff < 0:
        direction = "Responder_more_negative"

    else:
        direction = "Equal_median_change"

    result[
        "Descriptive_Direction"
    ] = direction

    rows.append(result)


summary = pd.DataFrame(rows)


# ============================================================
# 11. VALIDATE GROUP N FOR EVERY EDGE
# ============================================================

if not (
    summary["Responder_N"] == EXPECTED_RESPONDERS
).all():
    raise RuntimeError(
        "Not every edge contains all responders."
    )

if not (
    summary["NonResponder_N"] == EXPECTED_NONRESPONDERS
).all():
    raise RuntimeError(
        "Not every edge contains all non-responders."
    )


print(
    "✓ Every edge contains 13 responders "
    "and 23 non-responders."
)


# ============================================================
# 12. SAVE FULL SUMMARY
# ============================================================

summary = summary.sort_values(
    "Connected_Node_ID"
).reset_index(drop=True)

summary.to_csv(
    SUMMARY_OUT,
    index=False
)


# ============================================================
# 13. RANK BY ABSOLUTE MEDIAN GROUP DIFFERENCE
#
# Ranking is for DESCRIPTIVE INSPECTION ONLY.
# It is NOT significance ranking.
# ============================================================

ranked = summary.copy()

ranked[
    "Absolute_Median_Delta_Difference"
] = ranked[
    "Median_Delta_Difference_Resp_minus_NonResp"
].abs()

ranked = ranked.sort_values(
    "Absolute_Median_Delta_Difference",
    ascending=False
).reset_index(drop=True)

ranked.insert(
    0,
    "Descriptive_Rank",
    np.arange(1, len(ranked) + 1)
)

ranked.to_csv(
    RANKED_OUT,
    index=False
)


# ============================================================
# 14. TERMINAL SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("DESCRIPTIVE GROUP DIFFERENCE IN ΔFBC")
print("=" * 80)

terminal_cols = [
    "Connected_Node_ID",
    "Connected_Node",
    "Responder_Delta_Median",
    "NonResponder_Delta_Median",
    "Median_Delta_Difference_Resp_minus_NonResp",
    "Descriptive_Direction"
]

print(
    ranked[
        terminal_cols
    ].to_string(
        index=False,
        formatters={
            "Responder_Delta_Median":
                lambda x: f"{x:.8f}",

            "NonResponder_Delta_Median":
                lambda x: f"{x:.8f}",

            "Median_Delta_Difference_Resp_minus_NonResp":
                lambda x: f"{x:.8f}"
        }
    )
)


# ============================================================
# 15. FIGURE
#
# Ranked by:
# responder median ΔFBC - non-responder median ΔFBC
#
# No p-values.
# ============================================================

plot_df = summary.sort_values(
    "Median_Delta_Difference_Resp_minus_NonResp"
).reset_index(drop=True)

x = np.arange(len(plot_df))

resp = plot_df[
    "Responder_Delta_Median"
].to_numpy()

nonresp = plot_df[
    "NonResponder_Delta_Median"
].to_numpy()

difference = plot_df[
    "Median_Delta_Difference_Resp_minus_NonResp"
].to_numpy()


fig, ax = plt.subplots(
    figsize=(13, 7)
)

# Median ΔFBC for each group
ax.scatter(
    x,
    resp,
    marker="o",
    label=f"Responders (n={EXPECTED_RESPONDERS})"
)

ax.scatter(
    x,
    nonresp,
    marker="s",
    label=f"Non-responders (n={EXPECTED_NONRESPONDERS})"
)

# Connect group medians for each edge
for i in range(len(plot_df)):
    ax.plot(
        [x[i], x[i]],
        [nonresp[i], resp[i]],
        linewidth=0.8,
        alpha=0.5
    )


ax.axhline(
    0,
    linewidth=1,
    linestyle="--"
)

ax.set_xlabel(
    "Fixed DLPFC connections ranked by responder–non-responder difference"
)

ax.set_ylabel(
    "Median ΔFBC (Follow-up − Baseline)"
)

ax.set_title(
    "Longitudinal FBC Change in Responders vs Non-responders\n"
    "33 Fixed Left-DLPFC Connections"
)

ax.set_xticks(x)

ax.set_xticklabels(
    plot_df["Connected_Node_ID"].astype(str),
    rotation=90,
    fontsize=8
)

ax.legend(
    frameon=False
)

fig.tight_layout()

fig.savefig(
    FIGURE_OUT,
    dpi=300,
    bbox_inches="tight"
)

fig.savefig(
    FIGURE_PDF,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# 16. FINAL QC SUMMARY
# ============================================================

n_resp_more_positive = int(
    (
        summary[
            "Median_Delta_Difference_Resp_minus_NonResp"
        ]
        >
        0
    ).sum()
)

n_resp_more_negative = int(
    (
        summary[
            "Median_Delta_Difference_Resp_minus_NonResp"
        ]
        <
        0
    ).sum()
)

n_equal = int(
    np.isclose(
        summary[
            "Median_Delta_Difference_Resp_minus_NonResp"
        ],
        0,
        rtol=0,
        atol=1e-15
    ).sum()
)


print("\n" + "=" * 80)
print("STEP 6B SUMMARY")
print("=" * 80)

print(f"\nParticipants        : {df['ID'].nunique()}")
print(f"Responders          : {n_resp}")
print(f"Non-responders      : {n_nonresp}")
print(f"Fixed edges         : {summary.shape[0]}")

print(
    "\nEdges where responder median ΔFBC "
    "is more positive than non-responder:"
)
print(f"  {n_resp_more_positive}")

print(
    "\nEdges where responder median ΔFBC "
    "is more negative than non-responder:"
)
print(f"  {n_resp_more_negative}")

print(
    "\nEdges with equal median ΔFBC:"
)
print(f"  {n_equal}")


print("\n" + "=" * 80)
print("OUTPUT FILES")
print("=" * 80)

print(f"\nFull descriptive summary:\n{SUMMARY_OUT}")
print(f"\nClean long-form data:\n{LONG_OUT}")
print(f"\nDescriptively ranked edges:\n{RANKED_OUT}")
print(f"\nFigure PNG:\n{FIGURE_OUT}")
print(f"\nFigure PDF:\n{FIGURE_PDF}")


print("\n" + "=" * 80)
print("DONE")
print("=" * 80)

print(
    "\n✓ Responders and non-responders analysed separately."
)

print(
    "✓ Longitudinal change defined as FU_FBC - BL_FBC."
)

print(
    "✓ Same fixed 33 DLPFC edges retained."
)

print(
    "✓ No inferential statistical tests performed."
)

print(
    "✓ No p-values or multiple-comparison correction performed."
)
