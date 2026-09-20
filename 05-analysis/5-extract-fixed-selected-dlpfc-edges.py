from pathlib import Path
import pandas as pd
import numpy as np


# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Desktop" / "start-analysis"

CLINICAL = BASE / "clinical-output.xlsx"

SELECTED_FILE = (
    BASE /
    "4-dlpfc-median-iqr" /
    "DLPFC_04_selected_extreme_edges.csv"
)

OUTPUT_DIR = BASE / "5-fixed-dlpfc-edges"
OUTPUT_DIR.mkdir(exist_ok=True)

LONG_OUT = (
    OUTPUT_DIR /
    "DLPFC_05_fixed_33_edges_long.csv"
)

WIDE_BL_OUT = (
    OUTPUT_DIR /
    "DLPFC_05_fixed_33_edges_BL_wide.csv"
)

WIDE_FU_OUT = (
    OUTPUT_DIR /
    "DLPFC_05_fixed_33_edges_FU_wide.csv"
)

WIDE_DELTA_OUT = (
    OUTPUT_DIR /
    "DLPFC_05_fixed_33_edges_Delta_wide.csv"
)

SELECTED_NODES_OUT = (
    OUTPUT_DIR /
    "DLPFC_05_locked_selected_nodes.csv"
)

AUDIT_OUT = (
    OUTPUT_DIR /
    "DLPFC_05_extraction_audit.csv"
)


# ============================================================
# SETTINGS
# ============================================================

DLPFC_ID = 15
DLPFC_NAME = "ctx_lh_G_front_middle"

EXPECTED_MATRIX_SIZE = 164


# ============================================================
# HELPER: PARSE LABEL
# ============================================================

def parse_label(label):
    """
    Expected format:
        15|ctx_lh_G_front_middle

    Returns:
        node_id, node_name
    """

    label = str(label).strip()

    if "|" not in label:
        raise RuntimeError(
            f"Expected Node_ID|Node_Name format, "
            f"but found: {label}"
        )

    node_id_text, node_name = label.split("|", 1)

    return int(node_id_text), node_name


# ============================================================
# HELPER: READ + VALIDATE MATRIX
# ============================================================

def read_matrix(path, sid, timepoint):

    if not path.exists():
        raise FileNotFoundError(
            f"{sid}: missing {timepoint} matrix:\n"
            f"{path}"
        )

    df = pd.read_csv(
        path,
        index_col=0
    )

    if df.shape != (
        EXPECTED_MATRIX_SIZE,
        EXPECTED_MATRIX_SIZE
    ):
        raise RuntimeError(
            f"{sid} {timepoint}: expected "
            f"164x164 matrix, got {df.shape}"
        )

    # --------------------------------------------------------
    # Row/column labels must match exactly
    # --------------------------------------------------------

    row_labels = [
        str(x).strip()
        for x in df.index
    ]

    col_labels = [
        str(x).strip()
        for x in df.columns
    ]

    if row_labels != col_labels:
        raise RuntimeError(
            f"{sid} {timepoint}: "
            "row and column labels do not match."
        )

    # --------------------------------------------------------
    # Parse node IDs/names
    # --------------------------------------------------------

    parsed = [
        parse_label(x)
        for x in row_labels
    ]

    node_ids = [
        x[0]
        for x in parsed
    ]

    node_names = [
        x[1]
        for x in parsed
    ]

    if node_ids != list(
        range(1, 165)
    ):
        raise RuntimeError(
            f"{sid} {timepoint}: "
            "node IDs are not exactly 1–164."
        )

    if (
        node_ids[DLPFC_ID - 1]
        !=
        DLPFC_ID
    ):
        raise RuntimeError(
            f"{sid} {timepoint}: "
            "Node 15 position is incorrect."
        )

    if (
        node_names[DLPFC_ID - 1]
        !=
        DLPFC_NAME
    ):
        raise RuntimeError(
            f"{sid} {timepoint}: "
            "Node 15 anatomical label is incorrect."
        )

    # --------------------------------------------------------
    # Numerical validation
    # --------------------------------------------------------

    values = df.to_numpy(
        dtype=float
    )

    if not np.isfinite(
        values
    ).all():

        raise RuntimeError(
            f"{sid} {timepoint}: "
            "matrix contains NaN or infinity."
        )

    return df


# ============================================================
# HELPER: GET EDGE VALUE
# ============================================================

def get_edge_value(
    matrix,
    node_a_id,
    node_b_id
):
    """
    Matrices are stored upper triangular.

    Therefore retrieve the value from whichever
    triangle contains the anatomical edge.

    Node IDs are 1-based.
    Pandas/NumPy positions are 0-based.
    """

    i = node_a_id - 1
    j = node_b_id - 1

    a = float(
        matrix.iloc[i, j]
    )

    b = float(
        matrix.iloc[j, i]
    )

    # Diagonal is irrelevant here because
    # node_a != node_b.

    # --------------------------------------------------------
    # Expected structure:
    # one triangle contains the edge,
    # the opposite triangle is zero.
    # --------------------------------------------------------

    if a != 0 and b != 0:

        # If both triangles contain values, they should
        # represent the same edge. This is not expected
        # for the current stored matrices, so stop rather
        # than silently double-count.
        if not np.isclose(
            a,
            b,
            rtol=1e-12,
            atol=1e-15
        ):
            raise RuntimeError(
                f"Edge {node_a_id}-{node_b_id}: "
                f"both matrix triangles are non-zero "
                f"but differ ({a} vs {b})."
            )

        return a

    if a != 0:
        return a

    if b != 0:
        return b

    # Real absent/zero connection
    return 0.0


# ============================================================
# 1. READ THE 33 SELECTED NODES
# ============================================================

if not SELECTED_FILE.exists():
    raise FileNotFoundError(
        f"Selected-edge file not found:\n"
        f"{SELECTED_FILE}"
    )

selected = pd.read_csv(
    SELECTED_FILE
)

required_selected_columns = [
    "Connected_Node_ID",
    "Connected_Node"
]

missing = [
    x
    for x in required_selected_columns
    if x not in selected.columns
]

if missing:
    raise RuntimeError(
        f"Selected-edge file missing columns: "
        f"{missing}"
    )


selected["Connected_Node_ID"] = (
    selected["Connected_Node_ID"]
    .astype(int)
)

selected["Connected_Node"] = (
    selected["Connected_Node"]
    .astype(str)
)


# ------------------------------------------------------------
# Clean name if combined format somehow appears
# ------------------------------------------------------------

selected["Connected_Node"] = (
    selected["Connected_Node"]
    .str.split("|", regex=False)
    .str[-1]
)


# ------------------------------------------------------------
# Validate selected node set
# ------------------------------------------------------------

if selected[
    "Connected_Node_ID"
].duplicated().any():

    raise RuntimeError(
        "Selected node file contains duplicate node IDs."
    )


if (
    selected[
        "Connected_Node_ID"
    ]
    ==
    DLPFC_ID
).any():

    raise RuntimeError(
        "DLPFC Node 15 itself appears "
        "in selected connections."
    )


N_SELECTED = len(
    selected
)


print("\n" + "=" * 78)
print("STEP 5 — EXTRACT FIXED SELECTED DLPFC EDGES")
print("=" * 78)

print(
    f"\nLocked selected edges: "
    f"{N_SELECTED}"
)

if N_SELECTED != 33:

    print(
        f"\nWARNING: expected 33 selected edges "
        f"from Step 4, but found {N_SELECTED}."
    )


# ------------------------------------------------------------
# Sort by node ID for fixed consistent order
# ------------------------------------------------------------

selected = selected.sort_values(
    "Connected_Node_ID"
).reset_index(
    drop=True
)


# ------------------------------------------------------------
# Save LOCKED node list
# ------------------------------------------------------------

locked_nodes = selected[
    [
        "Connected_Node_ID",
        "Connected_Node"
    ]
].copy()

locked_nodes.insert(
    0,
    "DLPFC_Node_ID",
    DLPFC_ID
)

locked_nodes.insert(
    1,
    "DLPFC_Node",
    DLPFC_NAME
)

locked_nodes.to_csv(
    SELECTED_NODES_OUT,
    index=False
)


print(
    "\n✓ Selected node set locked."
)

print(
    "✓ The SAME node IDs will be extracted "
    "for every participant."
)


# ============================================================
# 2. READ SUBJECT IDs
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
    sid
    for sid in subjects
    if sid.startswith("YTH")
]


print(
    f"\nParticipants expected: "
    f"{len(subjects)}"
)


# ============================================================
# 3. EXTRACT SAME FIXED EDGES FOR EVERY SUBJECT
# ============================================================

all_rows = []
audit_rows = []


for sid in subjects:

    print(f"\n{sid}")

    subject_dir = (
        BASE /
        sid
    )

    BL_FILE = (
        subject_dir /
        "BL_fbc_labelled.csv"
    )

    FU_FILE = (
        subject_dir /
        "FU_fbc_labelled.csv"
    )


    # --------------------------------------------------------
    # Read matrices
    # --------------------------------------------------------

    BL = read_matrix(
        BL_FILE,
        sid,
        "BL"
    )

    FU = read_matrix(
        FU_FILE,
        sid,
        "FU"
    )


    # --------------------------------------------------------
    # BL/FU labels must be identical
    # --------------------------------------------------------

    if list(BL.index) != list(FU.index):

        raise RuntimeError(
            f"{sid}: BL/FU row labels differ."
        )

    if list(BL.columns) != list(FU.columns):

        raise RuntimeError(
            f"{sid}: BL/FU column labels differ."
        )


    subject_rows = []


    # --------------------------------------------------------
    # Extract EXACTLY the locked selected nodes
    # --------------------------------------------------------

    for _, edge in selected.iterrows():

        connected_id = int(
            edge[
                "Connected_Node_ID"
            ]
        )

        expected_name = str(
            edge[
                "Connected_Node"
            ]
        )


        # ----------------------------------------------------
        # Validate anatomical label against matrix
        # ----------------------------------------------------

        matrix_label = str(
            BL.index[
                connected_id - 1
            ]
        )

        parsed_id, parsed_name = (
            parse_label(
                matrix_label
            )
        )

        if parsed_id != connected_id:

            raise RuntimeError(
                f"{sid}: expected Node "
                f"{connected_id}, found "
                f"{parsed_id}."
            )

        if parsed_name != expected_name:

            raise RuntimeError(
                f"{sid}: Node {connected_id} "
                "name mismatch: "
                f"{parsed_name} vs "
                f"{expected_name}"
            )


        # ----------------------------------------------------
        # Extract BL and FU
        # ----------------------------------------------------

        bl_fbc = get_edge_value(
            BL,
            DLPFC_ID,
            connected_id
        )

        fu_fbc = get_edge_value(
            FU,
            DLPFC_ID,
            connected_id
        )

        delta_fbc = (
            fu_fbc
            -
            bl_fbc
        )


        subject_rows.append({

            "ID":
                sid,

            "DLPFC_Node_ID":
                DLPFC_ID,

            "DLPFC_Node":
                DLPFC_NAME,

            "Connected_Node_ID":
                connected_id,

            "Connected_Node":
                expected_name,

            "BL_FBC":
                bl_fbc,

            "FU_FBC":
                fu_fbc,

            "Delta_FBC":
                delta_fbc,

            "BL_nonzero":
                bool(bl_fbc != 0),

            "FU_nonzero":
                bool(fu_fbc != 0)
        })


    subject_df = pd.DataFrame(
        subject_rows
    )


    # --------------------------------------------------------
    # STRICT VALIDATION:
    # every subject gets exactly same selected nodes
    # --------------------------------------------------------

    observed_ids = (
        subject_df[
            "Connected_Node_ID"
        ]
        .astype(int)
        .tolist()
    )

    expected_ids = (
        selected[
            "Connected_Node_ID"
        ]
        .astype(int)
        .tolist()
    )

    if observed_ids != expected_ids:

        raise RuntimeError(
            f"{sid}: extracted node set "
            "does not match locked node set."
        )


    if len(subject_df) != N_SELECTED:

        raise RuntimeError(
            f"{sid}: expected "
            f"{N_SELECTED} extracted edges, "
            f"found {len(subject_df)}."
        )


    # --------------------------------------------------------
    # Count zero connections
    # --------------------------------------------------------

    n_bl_zero = int(
        (
            subject_df[
                "BL_FBC"
            ]
            ==
            0
        ).sum()
    )

    n_fu_zero = int(
        (
            subject_df[
                "FU_FBC"
            ]
            ==
            0
        ).sum()
    )

    n_zero_zero = int(
        (
            (
                subject_df[
                    "BL_FBC"
                ]
                ==
                0
            )
            &
            (
                subject_df[
                    "FU_FBC"
                ]
                ==
                0
            )
        ).sum()
    )


    # --------------------------------------------------------
    # Save participant-level extraction
    # --------------------------------------------------------

    subject_out = (
        subject_dir /
        "DLPFC_fixed_selected_edges.csv"
    )

    subject_df.to_csv(
        subject_out,
        index=False
    )


    all_rows.append(
        subject_df
    )


    audit_rows.append({

        "ID":
            sid,

        "N_expected_edges":
            N_SELECTED,

        "N_extracted_edges":
            len(subject_df),

        "N_BL_zero":
            n_bl_zero,

        "N_FU_zero":
            n_fu_zero,

        "N_BL_and_FU_zero":
            n_zero_zero,

        "Node_set_matches":
            True
    })


    print(
        f"  ✓ {len(subject_df)} / "
        f"{N_SELECTED} fixed edges extracted"
    )

    print(
        f"    BL zero: {n_bl_zero}"
    )

    print(
        f"    FU zero: {n_fu_zero}"
    )

    print(
        f"    BL+FU both zero: "
        f"{n_zero_zero}"
    )


# ============================================================
# 4. COMBINE ALL PARTICIPANTS
# ============================================================

long_df = pd.concat(
    all_rows,
    ignore_index=True
)


expected_total_rows = (
    len(subjects)
    *
    N_SELECTED
)


if len(long_df) != expected_total_rows:

    raise RuntimeError(
        "Combined long table has incorrect "
        f"number of rows: {len(long_df)} "
        f"instead of {expected_total_rows}."
    )


# ============================================================
# 5. FINAL CROSS-SUBJECT NODE-SET CHECK
# ============================================================

expected_set = set(
    selected[
        "Connected_Node_ID"
    ].astype(int)
)


for sid, sub in long_df.groupby(
    "ID"
):

    observed_set = set(
        sub[
            "Connected_Node_ID"
        ].astype(int)
    )

    if observed_set != expected_set:

        raise RuntimeError(
            f"{sid}: final node-set validation failed."
        )


print("\n" + "=" * 78)
print("CROSS-SUBJECT VALIDATION")
print("=" * 78)

print(
    f"\n✓ {len(subjects)} participants"
)

print(
    f"✓ {N_SELECTED} identical fixed "
    "DLPFC edges per participant"
)

print(
    f"✓ {len(long_df)} total rows "
    f"({len(subjects)} × {N_SELECTED})"
)


# ============================================================
# 6. SAVE LONG-FORM MASTER DATASET
# ============================================================

long_df = long_df.sort_values(
    [
        "ID",
        "Connected_Node_ID"
    ]
).reset_index(
    drop=True
)

long_df.to_csv(
    LONG_OUT,
    index=False
)


# ============================================================
# 7. CREATE WIDE MATRICES
# ============================================================

def make_wide(
    df,
    value_column
):

    wide = df.pivot(
        index="ID",
        columns="Connected_Node_ID",
        values=value_column
    )

    # Keep locked node order
    wide = wide.reindex(
        columns=selected[
            "Connected_Node_ID"
        ].tolist()
    )

    # Rename columns with ID + name
    name_map = dict(
        zip(
            selected[
                "Connected_Node_ID"
            ],
            selected[
                "Connected_Node"
            ]
        )
    )

    wide.columns = [
        f"{node_id}|{name_map[node_id]}"
        for node_id in wide.columns
    ]

    wide = wide.reset_index()

    return wide


BL_wide = make_wide(
    long_df,
    "BL_FBC"
)

FU_wide = make_wide(
    long_df,
    "FU_FBC"
)

DELTA_wide = make_wide(
    long_df,
    "Delta_FBC"
)


BL_wide.to_csv(
    WIDE_BL_OUT,
    index=False
)

FU_wide.to_csv(
    WIDE_FU_OUT,
    index=False
)

DELTA_wide.to_csv(
    WIDE_DELTA_OUT,
    index=False
)


# ============================================================
# 8. SAVE AUDIT
# ============================================================

audit = pd.DataFrame(
    audit_rows
)

audit.to_csv(
    AUDIT_OUT,
    index=False
)


# ============================================================
# 9. FINAL SUMMARY
# ============================================================

print("\n" + "=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

print(
    f"\nParticipants              : "
    f"{len(subjects)}"
)

print(
    f"Locked DLPFC edges        : "
    f"{N_SELECTED}"
)

print(
    f"Edges per participant     : "
    f"{N_SELECTED}"
)

print(
    f"Total observations        : "
    f"{len(long_df)}"
)

print(
    f"Expected observations     : "
    f"{expected_total_rows}"
)

print(
    f"\nBL=0 observations         : "
    f"{int((long_df['BL_FBC'] == 0).sum())}"
)

print(
    f"FU=0 observations         : "
    f"{int((long_df['FU_FBC'] == 0).sum())}"
)

print(
    f"BL=FU=0 observations      : "
    f"{int(((long_df['BL_FBC'] == 0) & (long_df['FU_FBC'] == 0)).sum())}"
)


print("\n" + "=" * 78)
print("OUTPUT FILES")
print("=" * 78)

print(
    f"\nLocked 33-node list:\n"
    f"{SELECTED_NODES_OUT}"
)

print(
    f"\nMaster long-format dataset:\n"
    f"{LONG_OUT}"
)

print(
    f"\nBL wide matrix:\n"
    f"{WIDE_BL_OUT}"
)

print(
    f"\nFU wide matrix:\n"
    f"{WIDE_FU_OUT}"
)

print(
    f"\nDelta wide matrix:\n"
    f"{WIDE_DELTA_OUT}"
)

print(
    f"\nExtraction audit:\n"
    f"{AUDIT_OUT}"
)


print("\n" + "=" * 78)
print("DONE")
print("=" * 78)

print(
    "\n✓ Same selected nodes used for every participant."
)

print(
    "✓ FBC values extracted directly from "
    "BL/FU labelled source matrices."
)

print(
    "✓ Delta_FBC recalculated as FU - BL."
)

print(
    "✓ No responder/non-responder information used."
)

print(
    "✓ No source matrix was modified."
)

