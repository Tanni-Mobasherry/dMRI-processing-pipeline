python - <<'PY'
from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path.home() / "Desktop" / "start-analysis"

DLPFC_NODE = 15
DLPFC_NAME = "ctx_lh_G_front_middle"

# ============================================================
# SUBJECTS FROM CLINICAL FILE
# ============================================================

clinical = pd.read_excel(
    BASE / "clinical-output.xlsx",
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


# ============================================================
# SYMMETRISE FUNCTION
# ============================================================

def symmetrise(df):

    A = df.to_numpy(dtype=float)

    if A.shape != (164, 164):
        raise ValueError(
            f"Expected 164x164 matrix, got {A.shape}"
        )

    # These FBC matrices should contain data only
    # in the upper triangle.
    lower = np.tril(A, k=-1)

    if not np.allclose(lower, 0):
        raise ValueError(
            "Lower triangle contains non-zero values."
        )

    # Mirror upper triangle into lower triangle.
    # Keep diagonal unchanged.
    A_sym = A + A.T - np.diag(np.diag(A))

    return pd.DataFrame(
        A_sym,
        index=df.index,
        columns=df.columns
    )


# ============================================================
# PROCESS EACH SUBJECT
# ============================================================

summary = []

print("=" * 70)
print("DLPFC EDGE EXTRACTION")
print(f"Node {DLPFC_NODE}: {DLPFC_NAME}")
print("=" * 70)


for sid in subjects:

    subject_dir = BASE / sid

    print(f"\n{sid}")

    subject_result = {
        "ID": sid,
        "BL_nonzero": np.nan,
        "FU_nonzero": np.nan,
        "status": "OK"
    }

    for tp in ["BL", "FU"]:

        input_file = (
            subject_dir /
            f"{tp}_fbc_labelled.csv"
        )

        output_file = (
            subject_dir /
            f"DLPFC_{tp}_edges.csv"
        )

        # ----------------------------------------------------
        # CHECK INPUT
        # ----------------------------------------------------

        if not input_file.exists():

            print(f"  {tp}: INPUT FILE MISSING")

            subject_result["status"] = (
                f"MISSING_{tp}"
            )

            continue


        # ----------------------------------------------------
        # READ
        # ----------------------------------------------------

        df = pd.read_csv(
            input_file,
            index_col=0
        )

        # ----------------------------------------------------
        # VALIDATE
        # ----------------------------------------------------

        if df.shape != (164, 164):
            raise RuntimeError(
                f"{sid} {tp}: "
                f"wrong matrix shape {df.shape}"
            )

        if list(df.index) != list(df.columns):
            raise RuntimeError(
                f"{sid} {tp}: "
                "row/column labels differ"
            )

        if DLPFC_NAME not in df.index:
            raise RuntimeError(
                f"{sid} {tp}: "
                f"DLPFC label not found"
            )


        # ----------------------------------------------------
        # SYMMETRISE IN MEMORY
        # ----------------------------------------------------

        sym = symmetrise(df)


        # ----------------------------------------------------
        # EXTRACT DLPFC ROW
        # ----------------------------------------------------

        row = (
            sym
            .loc[DLPFC_NAME]
            .astype(float)
            .copy()
        )

        # Remove DLPFC → DLPFC
        row = row.drop(DLPFC_NAME)

        assert len(row) == 163


        # ----------------------------------------------------
        # CREATE OUTPUT TABLE
        # ----------------------------------------------------

        output = pd.DataFrame({

            "DLPFC_Node_ID":
                DLPFC_NODE,

            "DLPFC_Node":
                DLPFC_NAME,

            "Connected_Node_ID":
                [
                    df.index.get_loc(node) + 1
                    for node in row.index
                ],

            "Connected_Node":
                row.index,

            "FBC":
                row.values,

            "Nonzero":
                row.values != 0
        })


        # ----------------------------------------------------
        # SAVE IN SUBJECT FOLDER
        # ----------------------------------------------------

        output.to_csv(
            output_file,
            index=False
        )


        # ----------------------------------------------------
        # COUNT
        # ----------------------------------------------------

        nonzero_count = int(
            output["Nonzero"].sum()
        )

        subject_result[
            f"{tp}_nonzero"
        ] = nonzero_count


        print(
            f"  {tp}: "
            f"{nonzero_count} / 163 non-zero"
        )

        print(
            f"      saved: {output_file.name}"
        )


    summary.append(subject_result)


# ============================================================
# SAVE OVERALL SUMMARY
# ============================================================

summary_df = pd.DataFrame(summary)

summary_file = (
    BASE /
    "DLPFC_01_nonzero_summary.csv"
)

summary_df.to_csv(
    summary_file,
    index=False
)


# ============================================================
# DISPLAY SUMMARY
# ============================================================

print("\n")
print("=" * 70)
print("SUMMARY")
print("=" * 70)

print(
    summary_df[
        ["ID", "BL_nonzero", "FU_nonzero", "status"]
    ].to_string(index=False)
)

print("\nOverall summary saved:")
print(summary_file)

print(
    "\nOriginal labelled matrices were NOT modified."
)

PY
