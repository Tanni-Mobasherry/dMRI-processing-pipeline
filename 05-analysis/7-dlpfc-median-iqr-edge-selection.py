from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Desktop" / "start-analysis"

CLINICAL = BASE / "clinical-output.xlsx"

OUTPUT_DIR = BASE / "4-dlpfc-median-iqr"
OUTPUT_DIR.mkdir(exist_ok=True)

ALL_EDGES_OUT = (
    OUTPUT_DIR /
    "DLPFC_04_all_edge_group_statistics.csv"
)

SELECTED_OUT = (
    OUTPUT_DIR /
    "DLPFC_04_selected_extreme_edges.csv"
)

SUMMARY_OUT = (
    OUTPUT_DIR /
    "DLPFC_04_selection_summary.csv"
)

LONG_OUT = (
    OUTPUT_DIR /
    "DLPFC_04_all_subject_edge_deltas_long.csv"
)


# ============================================================
# SETTINGS
# ============================================================

EXPECTED_DLPFC_ID = 15
EXPECTED_DLPFC_NAME = "ctx_lh_G_front_middle"

N_POSSIBLE_EDGES = 163


# ============================================================
# 1. READ SUBJECT IDs
# ============================================================

clinical = pd.read_excel(
    CLINICAL,
    sheet_name="Sheet1"
)

subjects = (
    clinical["ID"]
    .dropna()
    .astype(str)
    .str.strip()
)

subjects = [
    s for s in subjects
    if s.startswith("YTH")
]


print("\n" + "=" * 78)
print("DLPFC GROUP-LEVEL MEDIAN ± 3×IQR EDGE SELECTION")
print("=" * 78)

print(f"\nSubjects expected: {len(subjects)}")

print(
    "\nIMPORTANT:"
    "\n  • Uses signed Delta_FBC = FU - BL"
    "\n  • Uses NON-ZERO participant-level connections only"
    "\n  • Missing/absent connections are NOT replaced with zero"
    "\n  • Responder/non-responder labels are NOT used"
)


# ============================================================
# 2. READ EACH PARTICIPANT'S DLPFC DELTA FILE
# ============================================================

all_subject_edges = []

subject_audit = []

reference_node_map = {}


for sid in subjects:

    infile = (
        BASE /
        sid /
        "DLPFC_delta_edges.csv"
    )

    print(f"\n{sid}")

    # --------------------------------------------------------
    # FILE EXISTS
    # --------------------------------------------------------

    if not infile.exists():

        raise FileNotFoundError(
            f"{sid}: missing file:\n{infile}"
        )


    df = pd.read_csv(infile)


    # --------------------------------------------------------
    # REQUIRED COLUMNS
    # --------------------------------------------------------

    required_columns = [
        "DLPFC_Node_ID",
        "DLPFC_Node",
        "Connected_Node_ID",
        "Connected_Node",
        "BL_FBC",
        "FU_FBC",
        "Delta_FBC"
    ]

    missing_columns = [
        col
        for col in required_columns
        if col not in df.columns
    ]

    if missing_columns:

        raise RuntimeError(
            f"{sid}: missing columns: "
            f"{missing_columns}"
        )


    # --------------------------------------------------------
    # DLPFC VALIDATION
    # --------------------------------------------------------

    if not (
        df["DLPFC_Node_ID"]
        .astype(int)
        ==
        EXPECTED_DLPFC_ID
    ).all():

        raise RuntimeError(
            f"{sid}: DLPFC node ID is not "
            f"{EXPECTED_DLPFC_ID}."
        )


    # Allow either:
    #
    # ctx_lh_G_front_middle
    #
    # OR
    #
    # 15|ctx_lh_G_front_middle
    #
    # depending on which version of the earlier script
    # generated the participant-level delta file.

    dlpfc_names = (
        df["DLPFC_Node"]
        .astype(str)
        .str.split("|", regex=False)
        .str[-1]
    )

    if not (
        dlpfc_names
        ==
        EXPECTED_DLPFC_NAME
    ).all():

        raise RuntimeError(
            f"{sid}: unexpected DLPFC anatomical label."
        )


    # --------------------------------------------------------
    # CLEAN CONNECTED NODE NAME
    #
    # If labels are:
    #
    # 23|ctx_lh_...
    #
    # retain ID separately and anatomical name separately.
    # --------------------------------------------------------

    df["Connected_Node_ID"] = (
        df["Connected_Node_ID"]
        .astype(int)
    )

    df["Connected_Node"] = (
        df["Connected_Node"]
        .astype(str)
        .str.split("|", regex=False)
        .str[-1]
    )


    # --------------------------------------------------------
    # CHECK NODE IDs
    # --------------------------------------------------------

    if df["Connected_Node_ID"].duplicated().any():

        duplicates = (
            df.loc[
                df["Connected_Node_ID"].duplicated(),
                "Connected_Node_ID"
            ]
            .tolist()
        )

        raise RuntimeError(
            f"{sid}: duplicate connected node IDs: "
            f"{duplicates}"
        )


    # DLPFC itself must not be present.
    if (
        df["Connected_Node_ID"]
        ==
        EXPECTED_DLPFC_ID
    ).any():

        raise RuntimeError(
            f"{sid}: DLPFC self-connection "
            "is still present."
        )


    # Valid node range
    if not (
        df["Connected_Node_ID"]
        .between(1, 164)
    ).all():

        raise RuntimeError(
            f"{sid}: connected node ID "
            "outside range 1–164."
        )


    # --------------------------------------------------------
    # NUMERIC VALIDATION
    # --------------------------------------------------------

    for col in [
        "BL_FBC",
        "FU_FBC",
        "Delta_FBC"
    ]:

        df[col] = pd.to_numeric(
            df[col],
            errors="raise"
        )


    if not np.isfinite(
        df[
            [
                "BL_FBC",
                "FU_FBC",
                "Delta_FBC"
            ]
        ].to_numpy()
    ).all():

        raise RuntimeError(
            f"{sid}: NaN or infinite FBC value found."
        )


    # --------------------------------------------------------
    # VERIFY DELTA = FU - BL
    # --------------------------------------------------------

    expected_delta = (
        df["FU_FBC"]
        -
        df["BL_FBC"]
    )

    if not np.allclose(
        df["Delta_FBC"],
        expected_delta,
        rtol=1e-12,
        atol=1e-15
    ):

        max_diff = np.max(
            np.abs(
                df["Delta_FBC"]
                -
                expected_delta
            )
        )

        raise RuntimeError(
            f"{sid}: Delta_FBC does not equal "
            f"FU_FBC - BL_FBC. "
            f"Max difference={max_diff:.3e}"
        )


    # --------------------------------------------------------
    # VERIFY ZERO-ZERO EDGES WERE REMOVED
    # --------------------------------------------------------

    zero_zero = (
        (df["BL_FBC"] == 0)
        &
        (df["FU_FBC"] == 0)
    )

    if zero_zero.any():

        raise RuntimeError(
            f"{sid}: found "
            f"{int(zero_zero.sum())} "
            "BL=0/FU=0 edges. "
            "Step 3 should have removed them."
        )


    # --------------------------------------------------------
    # NUMBER OF RETAINED CONNECTIONS
    # --------------------------------------------------------

    n_edges = len(df)

    if n_edges > N_POSSIBLE_EDGES:

        raise RuntimeError(
            f"{sid}: found {n_edges} edges; "
            f"maximum possible is "
            f"{N_POSSIBLE_EDGES}."
        )


    # --------------------------------------------------------
    # VERIFY NODE ID -> NAME CONSISTENCY ACROSS SUBJECTS
    # --------------------------------------------------------

    for _, row in df.iterrows():

        node_id = int(
            row["Connected_Node_ID"]
        )

        node_name = str(
            row["Connected_Node"]
        )

        if node_id in reference_node_map:

            if (
                reference_node_map[node_id]
                !=
                node_name
            ):

                raise RuntimeError(
                    f"{sid}: Node {node_id} has "
                    "inconsistent anatomical names: "
                    f"{reference_node_map[node_id]} vs "
                    f"{node_name}"
                )

        else:

            reference_node_map[
                node_id
            ] = node_name


    # --------------------------------------------------------
    # ADD SUBJECT ID
    # --------------------------------------------------------

    subject_data = df[
        [
            "Connected_Node_ID",
            "Connected_Node",
            "BL_FBC",
            "FU_FBC",
            "Delta_FBC"
        ]
    ].copy()

    subject_data.insert(
        0,
        "ID",
        sid
    )

    all_subject_edges.append(
        subject_data
    )


    # --------------------------------------------------------
    # SUBJECT AUDIT
    # --------------------------------------------------------

    n_increase = int(
        (
            df["Delta_FBC"] > 0
        ).sum()
    )

    n_decrease = int(
        (
            df["Delta_FBC"] < 0
        ).sum()
    )

    n_no_change = int(
        (
            df["Delta_FBC"] == 0
        ).sum()
    )


    subject_audit.append({
        "ID": sid,
        "N_edges": n_edges,
        "N_increase": n_increase,
        "N_decrease": n_decrease,
        "N_no_change": n_no_change
    })


    print(
        f"  ✓ {n_edges} non-zero "
        "DLPFC connections"
    )


# ============================================================
# 3. COMBINE ALL SUBJECT-LEVEL DATA
# ============================================================

long_df = pd.concat(
    all_subject_edges,
    ignore_index=True
)


print("\n" + "=" * 78)
print("COMBINED DATA")
print("=" * 78)

print(
    f"Subjects loaded: "
    f"{long_df['ID'].nunique()}"
)

print(
    f"Total participant-edge observations: "
    f"{len(long_df)}"
)

print(
    f"Unique anatomical DLPFC edges: "
    f"{long_df['Connected_Node_ID'].nunique()}"
)


# ============================================================
# 4. VERIFY ALL 163 POSSIBLE DLPFC EDGES ARE REPRESENTED
#    SOMEWHERE IN THE SAMPLE
# ============================================================

expected_connected_ids = set(
    range(1, 165)
)

expected_connected_ids.remove(
    EXPECTED_DLPFC_ID
)

observed_connected_ids = set(
    long_df["Connected_Node_ID"]
    .astype(int)
    .unique()
)

missing_from_entire_sample = sorted(
    expected_connected_ids
    -
    observed_connected_ids
)


if missing_from_entire_sample:

    print(
        "\nWARNING:"
        "\nThe following DLPFC edges were absent "
        "in every participant:"
    )

    print(
        missing_from_entire_sample
    )

else:

    print(
        "\n✓ All 163 possible DLPFC edges are "
        "represented in at least one participant."
    )


# ============================================================
# 5. SAVE LONG-FORM DATA FOR AUDIT
# ============================================================

long_df = long_df.sort_values(
    [
        "Connected_Node_ID",
        "ID"
    ]
).reset_index(
    drop=True
)

long_df.to_csv(
    LONG_OUT,
    index=False
)


# ============================================================
# 6. CALCULATE GROUP STATISTICS FOR EACH ANATOMICAL EDGE
# ============================================================
#
# IMPORTANT:
#
# For each DLPFC -> connected-region edge:
#
# Median_Delta_FBC =
# median signed FU-BL change across participants
# in whom that connection exists.
#
# An absent connection is NOT replaced with zero.
#
# ============================================================

edge_rows = []


for node_id in sorted(
    observed_connected_ids
):

    edge = long_df[
        long_df["Connected_Node_ID"]
        ==
        node_id
    ].copy()


    # --------------------------------------------------------
    # Node name
    # --------------------------------------------------------

    names = (
        edge["Connected_Node"]
        .dropna()
        .unique()
    )

    if len(names) != 1:

        raise RuntimeError(
            f"Node {node_id}: expected exactly "
            f"one anatomical name, found "
            f"{list(names)}"
        )

    node_name = names[0]


    # --------------------------------------------------------
    # Delta values
    # --------------------------------------------------------

    delta = edge[
        "Delta_FBC"
    ].astype(float)


    n = len(delta)

    median_delta = float(
        delta.median()
    )

    q1_delta = float(
        delta.quantile(0.25)
    )

    q3_delta = float(
        delta.quantile(0.75)
    )

    iqr_delta = (
        q3_delta
        -
        q1_delta
    )

    mean_delta = float(
        delta.mean()
    )

    sd_delta = float(
        delta.std(ddof=1)
    ) if n > 1 else np.nan


    n_increase = int(
        (delta > 0).sum()
    )

    n_decrease = int(
        (delta < 0).sum()
    )

    n_no_change = int(
        (delta == 0).sum()
    )


    edge_rows.append({

        "Connected_Node_ID":
            int(node_id),

        "Connected_Node":
            node_name,

        "N":
            int(n),

        "N_missing_or_absent":
            int(
                len(subjects)
                -
                n
            ),

        "Median_Delta_FBC":
            median_delta,

        "Q1_Delta_FBC":
            q1_delta,

        "Q3_Delta_FBC":
            q3_delta,

        "IQR_Delta_FBC":
            iqr_delta,

        "Mean_Delta_FBC":
            mean_delta,

        "SD_Delta_FBC":
            sd_delta,

        "N_increase":
            n_increase,

        "N_decrease":
            n_decrease,

        "N_no_change":
            n_no_change
    })


edge_stats = pd.DataFrame(
    edge_rows
)


# ============================================================
# 7. VERIFY EDGE COUNTS
# ============================================================

if len(edge_stats) != 163:

    print(
        f"\nWARNING: group table contains "
        f"{len(edge_stats)} edges rather than 163."
    )

    print(
        "This means at least one anatomical edge "
        "was absent in every participant."
    )

else:

    print(
        "\n✓ Group table contains "
        "163 anatomical DLPFC edges."
    )


# ============================================================
# 8. CALCULATE DISTRIBUTION OF THE EDGE-LEVEL MEDIANS
# ============================================================
#
# We now have one group-level median Delta_FBC
# for each anatomical DLPFC edge.
#
# Sjoerd's criterion is applied ACROSS these edge-level
# difference values.
#
# ============================================================

edge_medians = (
    edge_stats[
        "Median_Delta_FBC"
    ]
)


global_median = float(
    edge_medians.median()
)

global_q1 = float(
    edge_medians.quantile(0.25)
)

global_q3 = float(
    edge_medians.quantile(0.75)
)

global_iqr = (
    global_q3
    -
    global_q1
)


# ============================================================
# 9. MEDIAN ± 3 × IQR THRESHOLDS
# ============================================================

lower_threshold = (
    global_median
    -
    3 * global_iqr
)

upper_threshold = (
    global_median
    +
    3 * global_iqr
)


print("\n" + "=" * 78)
print("GROUP-LEVEL EDGE-MEDIAN DISTRIBUTION")
print("=" * 78)

print(
    f"Number of edge medians : "
    f"{len(edge_medians)}"
)

print(
    f"Q1                    : "
    f"{global_q1:.10f}"
)

print(
    f"Median                : "
    f"{global_median:.10f}"
)

print(
    f"Q3                    : "
    f"{global_q3:.10f}"
)

print(
    f"IQR                   : "
    f"{global_iqr:.10f}"
)

print(
    f"Lower threshold       : "
    f"{lower_threshold:.10f}"
)

print(
    f"Upper threshold       : "
    f"{upper_threshold:.10f}"
)


# ============================================================
# 10. APPLY EXTREME-EDGE SELECTION
# ============================================================

edge_stats[
    "Group_Median_of_Edge_Medians"
] = global_median

edge_stats[
    "Group_Q1_of_Edge_Medians"
] = global_q1

edge_stats[
    "Group_Q3_of_Edge_Medians"
] = global_q3

edge_stats[
    "Group_IQR_of_Edge_Medians"
] = global_iqr

edge_stats[
    "Lower_Threshold"
] = lower_threshold

edge_stats[
    "Upper_Threshold"
] = upper_threshold


edge_stats[
    "Extreme"
] = (
    (
        edge_stats[
            "Median_Delta_FBC"
        ]
        <
        lower_threshold
    )
    |
    (
        edge_stats[
            "Median_Delta_FBC"
        ]
        >
        upper_threshold
    )
)


edge_stats[
    "Extreme_Direction"
] = np.where(

    edge_stats[
        "Median_Delta_FBC"
    ]
    <
    lower_threshold,

    "Extreme_Decrease",

    np.where(

        edge_stats[
            "Median_Delta_FBC"
        ]
        >
        upper_threshold,

        "Extreme_Increase",

        "Not_Extreme"
    )
)


# ============================================================
# 11. ADD COMBINED NODE LABEL
# ============================================================

edge_stats.insert(
    2,
    "Connected_Node_ID_Name",
    (
        edge_stats[
            "Connected_Node_ID"
        ].astype(str)
        +
        "|"
        +
        edge_stats[
            "Connected_Node"
        ].astype(str)
    )
)


# ============================================================
# 12. SORT COMPLETE TABLE BY NODE ID
# ============================================================

edge_stats = edge_stats.sort_values(
    "Connected_Node_ID"
).reset_index(
    drop=True
)


# ============================================================
# 13. SELECT EXTREME EDGES
# ============================================================

selected = edge_stats[
    edge_stats["Extreme"]
].copy()


# Sort selected edges from most negative to most positive.
selected = selected.sort_values(
    "Median_Delta_FBC"
).reset_index(
    drop=True
)


n_extreme_decrease = int(
    (
        selected[
            "Extreme_Direction"
        ]
        ==
        "Extreme_Decrease"
    ).sum()
)

n_extreme_increase = int(
    (
        selected[
            "Extreme_Direction"
        ]
        ==
        "Extreme_Increase"
    ).sum()
)

n_selected = len(
    selected
)


# ============================================================
# 14. SAVE COMPLETE AND SELECTED TABLES
# ============================================================

edge_stats.to_csv(
    ALL_EDGES_OUT,
    index=False
)

selected.to_csv(
    SELECTED_OUT,
    index=False
)


# ============================================================
# 15. SAVE SUMMARY
# ============================================================

summary = pd.DataFrame(
    [{
        "N_subjects":
            len(subjects),

        "N_unique_DLPFC_edges":
            len(edge_stats),

        "Global_Q1_Edge_Median_Delta_FBC":
            global_q1,

        "Global_Median_Edge_Median_Delta_FBC":
            global_median,

        "Global_Q3_Edge_Median_Delta_FBC":
            global_q3,

        "Global_IQR_Edge_Median_Delta_FBC":
            global_iqr,

        "Lower_Threshold_Median_minus_3IQR":
            lower_threshold,

        "Upper_Threshold_Median_plus_3IQR":
            upper_threshold,

        "N_Extreme_Decrease":
            n_extreme_decrease,

        "N_Extreme_Increase":
            n_extreme_increase,

        "N_Selected_Total":
            n_selected
    }]
)

summary.to_csv(
    SUMMARY_OUT,
    index=False
)


# ============================================================
# 16. FINAL RESULTS
# ============================================================

print("\n" + "=" * 78)
print("EDGE SELECTION RESULT")
print("=" * 78)

print(
    f"\nTotal anatomical DLPFC edges : "
    f"{len(edge_stats)}"
)

print(
    f"Extreme decreases            : "
    f"{n_extreme_decrease}"
)

print(
    f"Extreme increases            : "
    f"{n_extreme_increase}"
)

print(
    f"TOTAL SELECTED               : "
    f"{n_selected}"
)


# ============================================================
# 17. DISPLAY SELECTED EDGES
# ============================================================

if n_selected > 0:

    print("\nSELECTED EDGES:\n")

    display_columns = [
        "Connected_Node_ID",
        "Connected_Node",
        "N",
        "Median_Delta_FBC",
        "Extreme_Direction"
    ]

    print(
        selected[
            display_columns
        ].to_string(
            index=False
        )
    )

else:

    print(
        "\nNo edges exceeded "
        "median ± 3×IQR."
    )


# ============================================================
# 18. COVERAGE INFORMATION
# ============================================================

print("\n" + "=" * 78)
print("EDGE COVERAGE")
print("=" * 78)

print(
    f"Minimum participants contributing "
    f"to an edge : {edge_stats['N'].min()}"
)

print(
    f"Maximum participants contributing "
    f"to an edge : {edge_stats['N'].max()}"
)

print(
    f"Median participants contributing "
    f"to an edge : {edge_stats['N'].median():.1f}"
)


# ============================================================
# 19. OUTPUT FILES
# ============================================================

print("\n" + "=" * 78)
print("OUTPUT FILES")
print("=" * 78)

print(
    f"\nAll 163 edge statistics:\n"
    f"{ALL_EDGES_OUT}"
)

print(
    f"\nSelected extreme edges only:\n"
    f"{SELECTED_OUT}"
)

print(
    f"\nSelection summary:\n"
    f"{SUMMARY_OUT}"
)

print(
    f"\nLong-form participant data:\n"
    f"{LONG_OUT}"
)


print("\n" + "=" * 78)
print("DONE")
print("=" * 78)

print(
    "\nNo participant files were modified."
)

print(
    "No clinical response-group labels "
    "were used for edge selection."
)

print(
    "Absent/non-connected edges were "
    "not replaced with zero."
)

