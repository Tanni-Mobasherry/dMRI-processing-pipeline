from pathlib import Path
import pandas as pd
import numpy as np

BASE = Path.home() / "Desktop" / "start-analysis"

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

audit = []

print("=" * 78)
print("DLPFC LONGITUDINAL CHANGE — NON-ZERO CONNECTIONS ONLY")
print("Delta FBC = FU - BL")
print("=" * 78)


# ============================================================
# PROCESS EACH SUBJECT
# ============================================================

for sid in subjects:

    subject_dir = BASE / sid

    bl_file = subject_dir / "DLPFC_BL_edges.csv"
    fu_file = subject_dir / "DLPFC_FU_edges.csv"

    out_file = subject_dir / "DLPFC_delta_edges.csv"

    print(f"\n{sid}")

    # --------------------------------------------------------
    # CHECK INPUT FILES
    # --------------------------------------------------------

    if not bl_file.exists() or not fu_file.exists():

        print("  SKIPPED: BL/FU DLPFC edge file missing")

        audit.append({
            "ID": sid,
            "status": "MISSING_FILE"
        })

        continue


    # --------------------------------------------------------
    # READ BL AND FU
    # --------------------------------------------------------

    BL = pd.read_csv(bl_file)
    FU = pd.read_csv(fu_file)


    # --------------------------------------------------------
    # BASIC VALIDATION
    # --------------------------------------------------------

    if len(BL) != 163 or len(FU) != 163:

        raise RuntimeError(
            f"{sid}: expected 163 possible DLPFC edges, "
            f"got BL={len(BL)}, FU={len(FU)}"
        )


    # Check DLPFC anatomical label
    if not (BL["DLPFC_Node"] == DLPFC_NAME).all():

        raise RuntimeError(
            f"{sid}: unexpected DLPFC label in BL"
        )

    if not (FU["DLPFC_Node"] == DLPFC_NAME).all():

        raise RuntimeError(
            f"{sid}: unexpected DLPFC label in FU"
        )


    # --------------------------------------------------------
    # CHECK BL/FU NODE ALIGNMENT
    # --------------------------------------------------------

    if not BL["Connected_Node_ID"].equals(
        FU["Connected_Node_ID"]
    ):

        raise RuntimeError(
            f"{sid}: BL/FU Connected_Node_ID mismatch"
        )


    if not BL["Connected_Node"].equals(
        FU["Connected_Node"]
    ):

        raise RuntimeError(
            f"{sid}: BL/FU Connected_Node mismatch"
        )


    # --------------------------------------------------------
    # GET FBC VALUES
    # --------------------------------------------------------

    bl_fbc = BL["FBC"].astype(float)
    fu_fbc = FU["FBC"].astype(float)


    # --------------------------------------------------------
    # CHECK NON-ZERO SUPPORT
    # --------------------------------------------------------

    bl_nonzero = bl_fbc != 0
    fu_nonzero = fu_fbc != 0

    n_bl_nonzero = int(bl_nonzero.sum())
    n_fu_nonzero = int(fu_nonzero.sum())

    support_mismatch = bl_nonzero != fu_nonzero
    n_support_mismatch = int(support_mismatch.sum())


    print(
        f"  BL non-zero : {n_bl_nonzero}"
    )

    print(
        f"  FU non-zero : {n_fu_nonzero}"
    )


    if n_support_mismatch == 0:

        print(
            "  ✓ BL/FU non-zero support identical"
        )

    else:

        print(
            f"  WARNING: {n_support_mismatch} edges "
            "have different BL/FU support"
        )


    # --------------------------------------------------------
    # KEEP ACTUAL CONNECTIONS ONLY
    #
    # Include if the edge exists at BL OR FU.
    #
    # Therefore:
    # BL=0, FU=0 --> EXCLUDED
    #
    # BL>0, FU>0 --> INCLUDED
    #
    # If an edge ever appears/disappears longitudinally,
    # it is retained rather than silently discarded.
    # --------------------------------------------------------

    keep = bl_nonzero | fu_nonzero

    n_kept = int(keep.sum())
    n_excluded_zero_zero = int((~keep).sum())


    # --------------------------------------------------------
    # CREATE FILTERED OUTPUT
    # --------------------------------------------------------

    result = pd.DataFrame({

        "DLPFC_Node_ID":
            BL.loc[keep, "DLPFC_Node_ID"].values,

        "DLPFC_Node":
            BL.loc[keep, "DLPFC_Node"].values,

        "Connected_Node_ID":
            BL.loc[keep, "Connected_Node_ID"].values,

        "Connected_Node":
            BL.loc[keep, "Connected_Node"].values,

        "BL_FBC":
            bl_fbc.loc[keep].values,

        "FU_FBC":
            fu_fbc.loc[keep].values

    })


    # --------------------------------------------------------
    # SIGNED LONGITUDINAL DIFFERENCE
    # --------------------------------------------------------

    result["Delta_FBC"] = (
        result["FU_FBC"]
        -
        result["BL_FBC"]
    )


    # --------------------------------------------------------
    # DIRECTION
    # --------------------------------------------------------

    result["Direction"] = np.where(
        result["Delta_FBC"] > 0,
        "Increase",
        np.where(
            result["Delta_FBC"] < 0,
            "Decrease",
            "No_change"
        )
    )


    # --------------------------------------------------------
    # SAFETY CHECKS
    # --------------------------------------------------------

    # There must be no zero-zero edges left.
    zero_zero = (
        (result["BL_FBC"] == 0)
        &
        (result["FU_FBC"] == 0)
    )

    if zero_zero.any():

        raise RuntimeError(
            f"{sid}: zero-zero edge remained after filtering"
        )


    # Delta must equal FU - BL.
    if not np.allclose(
        result["Delta_FBC"],
        result["FU_FBC"] - result["BL_FBC"],
        rtol=1e-12,
        atol=1e-12
    ):

        raise RuntimeError(
            f"{sid}: Delta FBC validation failed"
        )


    # --------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------

    n_inc = int(
        (result["Delta_FBC"] > 0).sum()
    )

    n_dec = int(
        (result["Delta_FBC"] < 0).sum()
    )

    n_no_change = int(
        (result["Delta_FBC"] == 0).sum()
    )

    median_delta = float(
        result["Delta_FBC"].median()
    )


    # --------------------------------------------------------
    # SAVE — OVERWRITES OLD DLPFC_delta_edges.csv
    # --------------------------------------------------------

    result.to_csv(
        out_file,
        index=False
    )


    print(f"  Possible edges       : 163")
    print(f"  Zero-zero excluded   : {n_excluded_zero_zero}")
    print(f"  Connections retained : {n_kept}")
    print(f"  Increases            : {n_inc}")
    print(f"  Decreases            : {n_dec}")
    print(f"  No change            : {n_no_change}")
    print(f"  Median ΔFBC          : {median_delta:.8f}")
    print(f"  ✓ saved: {out_file.name}")


    audit.append({

        "ID":
            sid,

        "status":
            "OK",

        "Possible_edges":
            163,

        "BL_nonzero":
            n_bl_nonzero,

        "FU_nonzero":
            n_fu_nonzero,

        "Support_mismatch":
            n_support_mismatch,

        "Zero_zero_excluded":
            n_excluded_zero_zero,

        "N_edges_retained":
            n_kept,

        "N_increase":
            n_inc,

        "N_decrease":
            n_dec,

        "N_no_change":
            n_no_change,

        "Median_Delta_FBC":
            median_delta
    })


# ============================================================
# SAVE SUMMARY
# ============================================================

audit_df = pd.DataFrame(audit)

summary_file = (
    BASE /
    "DLPFC_03_nonzero_delta_summary.csv"
)

audit_df.to_csv(
    summary_file,
    index=False
)


# ============================================================
# FINAL DISPLAY
# ============================================================

print("\n")
print("=" * 78)
print("FINAL SUMMARY")
print("=" * 78)

cols = [
    "ID",
    "BL_nonzero",
    "FU_nonzero",
    "Zero_zero_excluded",
    "N_edges_retained",
    "N_increase",
    "N_decrease",
    "N_no_change"
]

print(
    audit_df[cols].to_string(index=False)
)


# ============================================================
# GLOBAL VALIDATION
# ============================================================

if (
    audit_df["N_edges_retained"]
    ==
    audit_df["BL_nonzero"]
).all() and (
    audit_df["Support_mismatch"] == 0
).all():

    print(
        "\n✓ VALIDATION PASSED:"
        "\n  retained edge count matches the "
        "previous BL non-zero count for every subject."
    )

else:

    print(
        "\n⚠ CHECK REQUIRED:"
        "\n  At least one subject has a support/count mismatch."
    )


print("\nSummary saved:")
print(summary_file)

print(
    "\nExisting DLPFC_delta_edges.csv files were replaced "
    "with the non-zero-only versions."
)

print(
    "Original DLPFC_BL_edges.csv and "
    "DLPFC_FU_edges.csv were NOT modified."
)

