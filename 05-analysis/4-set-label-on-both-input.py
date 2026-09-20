from pathlib import Path
import pandas as pd
import numpy as np
import csv


# ============================================================
# PATHS
# ============================================================

START = Path.home() / "Desktop" / "start-analysis"

CLINICAL = START / "clinical-output.xlsx"
LUT = START / "fs_a2009s.txt"

HARD = Path("/Volumes/Toshiba-Ext/raw-data")

AUDIT = START / "FBC_labelled_export_audit.csv"


# ============================================================
# 1. READ CLINICAL SUBJECT IDs
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
    .tolist()
)

subjects = [
    s for s in subjects
    if s.startswith("YTH")
]

print("\n========================================")
print("CLINICAL SUBJECTS")
print("========================================")
print(f"N = {len(subjects)}")
print(subjects)


# ============================================================
# 2. READ fs_a2009s.txt
# ============================================================
#
# Expected:
#
# 0  Unknown
# 1  ctx_lh_G_and_S_frontomargin
# ...
#
# column 1 = node ID
# column 2 = node name
#
# Node 0 is background and is excluded.
# ============================================================

lut = {}

with open(LUT, "r") as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        parts = line.split()

        if len(parts) < 2:
            continue

        try:
            node_id = int(parts[0])

        except ValueError:
            continue

        node_name = parts[1]

        if node_id == 0:
            continue

        lut[node_id] = node_name


print("\n========================================")
print("LOOKUP TABLE")
print("========================================")

print(
    f"Number of non-zero LUT entries found: "
    f"{len(lut)}"
)


# ============================================================
# 3. VALIDATE NODES 1–164
# ============================================================

missing_nodes = [
    node
    for node in range(1, 165)
    if node not in lut
]

if missing_nodes:

    raise RuntimeError(
        f"Missing node IDs in fs_a2009s.txt: "
        f"{missing_nodes}"
    )


# ============================================================
# 4. CREATE NODE IDS + NAMES + COMBINED LABELS
# ============================================================

node_ids = list(range(1, 165))

node_names = [
    lut[node]
    for node in node_ids
]

# Preserve BOTH node number and anatomical name.
#
# Example:
#
# 15|ctx_lh_G_front_middle

node_labels = [
    f"{node_id}|{node_name}"
    for node_id, node_name
    in zip(node_ids, node_names)
]


print(
    "\nLabels for nodes 1–164 successfully found."
)

print("\nFirst 5 nodes:")

for i in range(5):

    print(
        f"Node {node_ids[i]} -> "
        f"{node_names[i]} -> "
        f"{node_labels[i]}"
    )


# ============================================================
# DLPFC VALIDATION
# ============================================================

print("\nDLPFC CHECK:")

print(
    f"Matrix index 14 -> "
    f"Node {node_ids[14]} -> "
    f"{node_names[14]}"
)

print(
    f"Combined label: "
    f"{node_labels[14]}"
)

if node_ids[14] != 15:

    raise RuntimeError(
        "Expected matrix index 14 to correspond "
        "to anatomical Node 15."
    )

if node_names[14] != "ctx_lh_G_front_middle":

    raise RuntimeError(
        "Expected Node 15 to be "
        "ctx_lh_G_front_middle, "
        f"but found {node_names[14]}"
    )

print(
    "✓ Node 15 = ctx_lh_G_front_middle"
)


print("\nLast 5 nodes:")

for i in range(159, 164):

    print(
        f"Node {node_ids[i]} -> "
        f"{node_names[i]}"
    )


# ============================================================
# 5. FUNCTION TO SAVE LABELLED MATRIX
# ============================================================

def save_labelled_matrix(matrix, output_file):

    if matrix.shape != (164, 164):

        raise ValueError(
            f"Expected 164x164 matrix, "
            f"got {matrix.shape}"
        )

    with open(
        output_file,
        "w",
        newline=""
    ) as f:

        writer = csv.writer(f)

        # ----------------------------------------------------
        # HEADER
        #
        # First CSV column = row identifier
        # Remaining 164 columns = matrix Nodes 1–164
        #
        # Example:
        # 15|ctx_lh_G_front_middle
        # ----------------------------------------------------

        writer.writerow(
            ["Node_ID|Node_Name"]
            +
            node_labels
        )

        # ----------------------------------------------------
        # MATRIX ROWS
        # ----------------------------------------------------

        for i in range(164):

            writer.writerow(
                [node_labels[i]]
                +
                matrix[i, :].tolist()
            )


# ============================================================
# 6. PROCESS ALL SUBJECTS
# ============================================================

audit = []

print("\n========================================")
print("PROCESSING SUBJECTS")
print("========================================")


for subject in subjects:

    connectome_dir = (
        HARD
        / subject
        / "longi_FT"
        / "template"
        / "connectomes"
    )

    bl_file = (
        connectome_dir
        / "BL_fbc.csv"
    )

    fu_file = (
        connectome_dir
        / "FU_fbc.csv"
    )

    bl_out = (
        connectome_dir
        / "BL_fbc_labelled.csv"
    )

    fu_out = (
        connectome_dir
        / "FU_fbc_labelled.csv"
    )


    print(f"\n{subject}")


    # --------------------------------------------------------
    # CHECK CONNECTOME FOLDER
    # --------------------------------------------------------

    if not connectome_dir.exists():

        print(
            "  ✗ connectomes folder not found"
        )

        audit.append({
            "ID": subject,
            "status": "CONNECTOME_FOLDER_NOT_FOUND",
            "BL_found": False,
            "FU_found": False,
            "BL_labelled_saved": False,
            "FU_labelled_saved": False,
            "folder": str(connectome_dir)
        })

        continue


    # --------------------------------------------------------
    # CHECK ORIGINAL MATRICES
    # --------------------------------------------------------

    bl_exists = bl_file.exists()
    fu_exists = fu_file.exists()

    if not bl_exists or not fu_exists:

        print(
            f"  ✗ missing matrix:"
            f" BL={bl_exists},"
            f" FU={fu_exists}"
        )

        audit.append({
            "ID": subject,
            "status": "MISSING_FBC_MATRIX",
            "BL_found": bl_exists,
            "FU_found": fu_exists,
            "BL_labelled_saved": False,
            "FU_labelled_saved": False,
            "folder": str(connectome_dir)
        })

        continue


    # --------------------------------------------------------
    # LOAD ORIGINAL NUMERICAL MATRICES
    # --------------------------------------------------------

    try:

        BL = np.loadtxt(
            bl_file,
            delimiter=","
        )

        FU = np.loadtxt(
            fu_file,
            delimiter=","
        )

    except Exception as e:

        print(
            f"  ✗ could not read matrix: {e}"
        )

        audit.append({
            "ID": subject,
            "status": f"READ_ERROR: {e}",
            "BL_found": True,
            "FU_found": True,
            "BL_labelled_saved": False,
            "FU_labelled_saved": False,
            "folder": str(connectome_dir)
        })

        continue


    print(
        f"  BL shape: {BL.shape}"
    )

    print(
        f"  FU shape: {FU.shape}"
    )


    # --------------------------------------------------------
    # VALIDATE DIMENSIONS
    # --------------------------------------------------------

    if (
        BL.shape != (164, 164)
        or
        FU.shape != (164, 164)
    ):

        print(
            "  ✗ matrix is not 164 × 164"
        )

        audit.append({
            "ID": subject,
            "status": (
                f"WRONG_MATRIX_SIZE "
                f"BL={BL.shape}, "
                f"FU={FU.shape}"
            ),
            "BL_found": True,
            "FU_found": True,
            "BL_labelled_saved": False,
            "FU_labelled_saved": False,
            "folder": str(connectome_dir)
        })

        continue


    # --------------------------------------------------------
    # VALIDATE FINITE VALUES
    # --------------------------------------------------------

    if not np.isfinite(BL).all():

        raise RuntimeError(
            f"{subject}: BL contains "
            "NaN or infinite values."
        )

    if not np.isfinite(FU).all():

        raise RuntimeError(
            f"{subject}: FU contains "
            "NaN or infinite values."
        )


    # --------------------------------------------------------
    # SAVE LABELLED MATRICES
    # --------------------------------------------------------

    save_labelled_matrix(
        BL,
        bl_out
    )

    save_labelled_matrix(
        FU,
        fu_out
    )


    # --------------------------------------------------------
    # VERIFY OUTPUT EXISTS
    # --------------------------------------------------------

    bl_saved = bl_out.exists()
    fu_saved = fu_out.exists()

    if not (
        bl_saved
        and
        fu_saved
    ):

        print(
            "  ✗ output verification failed"
        )

        audit.append({
            "ID": subject,
            "status": "OUTPUT_ERROR",
            "BL_found": True,
            "FU_found": True,
            "BL_labelled_saved": bl_saved,
            "FU_labelled_saved": fu_saved,
            "folder": str(connectome_dir)
        })

        continue


    print(
        "  ✓ BL_fbc_labelled.csv"
    )

    print(
        "  ✓ FU_fbc_labelled.csv"
    )


    # ========================================================
    # 7. READ LABELLED MATRICES BACK AND VALIDATE
    # ========================================================

    BL_check = pd.read_csv(
        bl_out,
        index_col=0
    )

    FU_check = pd.read_csv(
        fu_out,
        index_col=0
    )


    # --------------------------------------------------------
    # VALIDATE LABELLED DIMENSIONS
    # --------------------------------------------------------

    if BL_check.shape != (164, 164):

        raise RuntimeError(
            f"{subject}: labelled BL has "
            f"shape {BL_check.shape}; "
            "expected (164, 164)."
        )

    if FU_check.shape != (164, 164):

        raise RuntimeError(
            f"{subject}: labelled FU has "
            f"shape {FU_check.shape}; "
            "expected (164, 164)."
        )


    # --------------------------------------------------------
    # VALIDATE ROW LABELS
    # --------------------------------------------------------

    if BL_check.index.tolist() != node_labels:

        raise RuntimeError(
            f"{subject}: BL row labels "
            "do not match Nodes 1–164."
        )

    if FU_check.index.tolist() != node_labels:

        raise RuntimeError(
            f"{subject}: FU row labels "
            "do not match Nodes 1–164."
        )


    # --------------------------------------------------------
    # VALIDATE COLUMN LABELS
    # --------------------------------------------------------

    if BL_check.columns.tolist() != node_labels:

        raise RuntimeError(
            f"{subject}: BL column labels "
            "do not match Nodes 1–164."
        )

    if FU_check.columns.tolist() != node_labels:

        raise RuntimeError(
            f"{subject}: FU column labels "
            "do not match Nodes 1–164."
        )


    # --------------------------------------------------------
    # VALIDATE DLPFC POSITION
    # --------------------------------------------------------

    expected_dlpfc = (
        "15|ctx_lh_G_front_middle"
    )

    if BL_check.index[14] != expected_dlpfc:

        raise RuntimeError(
            f"{subject}: BL DLPFC row "
            "is not anatomical Node 15."
        )

    if FU_check.index[14] != expected_dlpfc:

        raise RuntimeError(
            f"{subject}: FU DLPFC row "
            "is not anatomical Node 15."
        )

    if BL_check.columns[14] != expected_dlpfc:

        raise RuntimeError(
            f"{subject}: BL DLPFC column "
            "is not anatomical Node 15."
        )

    if FU_check.columns[14] != expected_dlpfc:

        raise RuntimeError(
            f"{subject}: FU DLPFC column "
            "is not anatomical Node 15."
        )


    # --------------------------------------------------------
    # VALIDATE NUMERICAL VALUES
    #
    # CSV write/read can cause tiny floating-point
    # representation differences.
    #
    # Therefore:
    #
    # 1. calculate the maximum absolute difference
    # 2. use a very tight numerical tolerance
    #
    # --------------------------------------------------------

    BL_readback = (
        BL_check
        .to_numpy(dtype=float)
    )

    FU_readback = (
        FU_check
        .to_numpy(dtype=float)
    )


    bl_abs_diff = np.abs(
        BL_readback - BL
    )

    fu_abs_diff = np.abs(
        FU_readback - FU
    )


    bl_max_diff = float(
        np.max(bl_abs_diff)
    )

    fu_max_diff = float(
        np.max(fu_abs_diff)
    )


    if not np.allclose(
        BL_readback,
        BL,
        rtol=1e-12,
        atol=1e-15
    ):

        raise RuntimeError(
            f"{subject}: BL numerical values "
            f"differ after labelling. "
            f"Max absolute difference = "
            f"{bl_max_diff:.3e}"
        )


    if not np.allclose(
        FU_readback,
        FU,
        rtol=1e-12,
        atol=1e-15
    ):

        raise RuntimeError(
            f"{subject}: FU numerical values "
            f"differ after labelling. "
            f"Max absolute difference = "
            f"{fu_max_diff:.3e}"
        )


    # --------------------------------------------------------
    # SUBJECT PASSED ALL CHECKS
    # --------------------------------------------------------

    print(
        "  ✓ row labels validated"
    )

    print(
        "  ✓ column labels validated"
    )

    print(
        "  ✓ Node 15 = "
        "ctx_lh_G_front_middle"
    )

    print(
        "  ✓ numerical FBC values preserved"
    )

    print(
        f"    max |difference| BL = "
        f"{bl_max_diff:.3e}"
    )

    print(
        f"    max |difference| FU = "
        f"{fu_max_diff:.3e}"
    )


    # --------------------------------------------------------
    # AUDIT
    # --------------------------------------------------------

    audit.append({
        "ID": subject,
        "status": "OK",
        "BL_found": True,
        "FU_found": True,
        "BL_labelled_saved": True,
        "FU_labelled_saved": True,
        "DLPFC_Node_ID": 15,
        "DLPFC_Node_Name":
            "ctx_lh_G_front_middle",
        "DLPFC_Combined_Label":
            "15|ctx_lh_G_front_middle",
        "BL_max_abs_numeric_diff":
            bl_max_diff,
        "FU_max_abs_numeric_diff":
            fu_max_diff,
        "folder":
            str(connectome_dir)
    })


# ============================================================
# 8. SAVE AUDIT
# ============================================================

audit_df = pd.DataFrame(
    audit
)

audit_df.to_csv(
    AUDIT,
    index=False
)


# ============================================================
# 9. FINAL SUMMARY
# ============================================================

print("\n\n========================================")
print("FINISHED")
print("========================================")

print(
    f"Clinical subjects: "
    f"{len(subjects)}"
)


if len(audit_df):

    n_ok = int(
        (
            audit_df["status"]
            ==
            "OK"
        ).sum()
    )

    failed = audit_df[
        audit_df["status"]
        !=
        "OK"
    ]

    print(
        f"Successfully labelled: "
        f"{n_ok}"
    )

    print(
        f"Not completed: "
        f"{len(failed)}"
    )


    if len(failed):

        print(
            "\nSubjects needing attention:"
        )

        print(
            failed[
                [
                    "ID",
                    "status"
                ]
            ].to_string(
                index=False
            )
        )


# ============================================================
# FINAL GLOBAL CHECK
# ============================================================

if (
    len(audit_df) == len(subjects)
    and
    (audit_df["status"] == "OK").all()
):

    print(
        "\n✓ ALL SUBJECTS PASSED"
    )

    print(
        f"✓ {len(subjects)} / "
        f"{len(subjects)} subjects successfully labelled"
    )

    print(
        "✓ Node IDs retained"
    )

    print(
        "✓ Anatomical names retained"
    )

    print(
        "✓ Numerical FBC matrices preserved"
    )

else:

    print(
        "\n⚠ NOT ALL SUBJECTS PASSED."
    )

    print(
        "Check the audit before continuing."
    )


print("\n========================================")
print("LABEL FORMAT")
print("========================================")

print(
    "Node_ID|Node_Name"
)

print(
    "Example:"
)

print(
    "15|ctx_lh_G_front_middle"
)


print("\nAudit saved to:")
print(AUDIT)


print(
    "\nOriginal BL_fbc.csv and "
    "FU_fbc.csv files were NOT modified."
)

