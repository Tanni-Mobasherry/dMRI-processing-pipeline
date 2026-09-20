from pathlib import Path
import pandas as pd
import numpy as np

# ============================================================
# PATHS
# ============================================================

BASE = Path.home() / "Desktop" / "start-analysis"

CLINICAL_FILE = BASE / "clinical-output.xlsx"

FBC_FILE = (
    BASE
    / "5-fixed-dlpfc-edges"
    / "DLPFC_05_fixed_33_edges_long.csv"
)

OUT_DIR = BASE / "6-dlpfc-clinical-analysis"
OUT_DIR.mkdir(exist_ok=True)

CLINICAL_QC_OUT = (
    OUT_DIR
    / "DLPFC_06A_clinical_responder_QC.csv"
)

MERGED_OUT = (
    OUT_DIR
    / "DLPFC_06A_fixed33_clinical_long.csv"
)

SUBJECT_QC_OUT = (
    OUT_DIR
    / "DLPFC_06A_subject_merge_QC.csv"
)


# ============================================================
# SETTINGS
# ============================================================

EXPECTED_EDGES = 33
RESPONDER_THRESHOLD = 0.50

# Use MRI-day DASS scores for responder classification
BL_DASS_COL = "BL_MRIdepDASS"
FU_DASS_COL = "FU_MRIdepDASS"

# Existing spreadsheet reduction column
EXISTING_REDUCTION_COL = "Reduction"


# ============================================================
# 1. READ CLINICAL DATA
# ============================================================

print("\n" + "=" * 78)
print("STEP 6A — CLINICAL MERGE + RESPONDER QC")
print("=" * 78)

print("\nReading clinical data...")

clinical = pd.read_excel(
    CLINICAL_FILE,
    sheet_name="Sheet1"
)

if "ID" not in clinical.columns:
    raise RuntimeError(
        "Clinical file does not contain an ID column."
    )

# Remove completely blank / non-subject rows
clinical["ID"] = (
    clinical["ID"]
    .astype("string")
    .str.strip()
)

clinical = clinical[
    clinical["ID"].notna()
].copy()

clinical = clinical[
    clinical["ID"].str.startswith("YTH", na=False)
].copy()

clinical["ID"] = clinical["ID"].astype(str)


# ============================================================
# 2. VALIDATE CLINICAL COLUMNS
# ============================================================

required_clinical = [
    "ID",
    BL_DASS_COL,
    FU_DASS_COL
]

missing = [
    c for c in required_clinical
    if c not in clinical.columns
]

if missing:
    raise RuntimeError(
        f"Missing required clinical columns: {missing}"
    )


# ------------------------------------------------------------
# IDs must be unique
# ------------------------------------------------------------

if clinical["ID"].duplicated().any():

    duplicates = clinical.loc[
        clinical["ID"].duplicated(keep=False),
        "ID"
    ].tolist()

    raise RuntimeError(
        f"Duplicate clinical IDs found: {duplicates}"
    )


print(
    f"✓ Clinical participants found: {len(clinical)}"
)


# ============================================================
# 3. CONVERT DASS TO NUMERIC
# ============================================================

clinical[BL_DASS_COL] = pd.to_numeric(
    clinical[BL_DASS_COL],
    errors="coerce"
)

clinical[FU_DASS_COL] = pd.to_numeric(
    clinical[FU_DASS_COL],
    errors="coerce"
)


# ============================================================
# 4. IDENTIFY MISSING DASS
# ============================================================

clinical["DASS_complete"] = (
    clinical[BL_DASS_COL].notna()
    &
    clinical[FU_DASS_COL].notna()
)

missing_dass = clinical.loc[
    ~clinical["DASS_complete"],
    [
        "ID",
        BL_DASS_COL,
        FU_DASS_COL
    ]
].copy()


# ============================================================
# 5. CHECK FOR INVALID BASELINE = 0
# ============================================================

zero_baseline = clinical.loc[
    clinical["DASS_complete"]
    &
    (clinical[BL_DASS_COL] == 0),
    "ID"
].tolist()

if zero_baseline:
    raise RuntimeError(
        "Cannot calculate proportional DASS reduction "
        "because baseline DASS = 0 for: "
        f"{zero_baseline}"
    )


# ============================================================
# 6. INDEPENDENTLY RECALCULATE DASS REDUCTION
#
# Reduction = (BL - FU) / BL
#
# 0.50 = 50% improvement
# ============================================================

clinical["DASS_Reduction_Calculated"] = np.nan

mask = clinical["DASS_complete"]

clinical.loc[
    mask,
    "DASS_Reduction_Calculated"
] = (
    (
        clinical.loc[mask, BL_DASS_COL]
        -
        clinical.loc[mask, FU_DASS_COL]
    )
    /
    clinical.loc[mask, BL_DASS_COL]
)


clinical["DASS_Improvement_Percent"] = (
    clinical["DASS_Reduction_Calculated"]
    * 100
)


# ============================================================
# 7. COMPARE AGAINST EXISTING REDUCTION COLUMN
# ============================================================

if EXISTING_REDUCTION_COL in clinical.columns:

    clinical[
        "Reduction_Excel"
    ] = pd.to_numeric(
        clinical[EXISTING_REDUCTION_COL],
        errors="coerce"
    )

    clinical[
        "Reduction_Difference"
    ] = (
        clinical[
            "DASS_Reduction_Calculated"
        ]
        -
        clinical[
            "Reduction_Excel"
        ]
    )

    comparison_mask = (
        clinical[
            "DASS_Reduction_Calculated"
        ].notna()
        &
        clinical[
            "Reduction_Excel"
        ].notna()
    )

    clinical[
        "Reduction_Matches_Excel"
    ] = pd.NA

    clinical.loc[
        comparison_mask,
        "Reduction_Matches_Excel"
    ] = np.isclose(
        clinical.loc[
            comparison_mask,
            "DASS_Reduction_Calculated"
        ],
        clinical.loc[
            comparison_mask,
            "Reduction_Excel"
        ],
        rtol=1e-9,
        atol=1e-9
    )

else:

    clinical["Reduction_Excel"] = np.nan
    clinical["Reduction_Difference"] = np.nan
    clinical["Reduction_Matches_Excel"] = pd.NA


# ============================================================
# 8. RESPONDER CLASSIFICATION
#
# Responder = >= 50% DASS improvement
# ============================================================

clinical["Responder_Status"] = pd.NA

clinical.loc[
    clinical["DASS_complete"]
    &
    (
        clinical[
            "DASS_Reduction_Calculated"
        ]
        >=
        RESPONDER_THRESHOLD
    ),
    "Responder_Status"
] = "Responder"

clinical.loc[
    clinical["DASS_complete"]
    &
    (
        clinical[
            "DASS_Reduction_Calculated"
        ]
        <
        RESPONDER_THRESHOLD
    ),
    "Responder_Status"
] = "Non-responder"


# Binary version for later modelling
clinical["Responder_Binary"] = np.nan

clinical.loc[
    clinical["Responder_Status"]
    ==
    "Responder",
    "Responder_Binary"
] = 1

clinical.loc[
    clinical["Responder_Status"]
    ==
    "Non-responder",
    "Responder_Binary"
] = 0


# ============================================================
# 9. CLINICAL QC SUMMARY
# ============================================================

n_complete = int(
    clinical["DASS_complete"].sum()
)

n_missing = int(
    (~clinical["DASS_complete"]).sum()
)

n_responders = int(
    (
        clinical["Responder_Status"]
        ==
        "Responder"
    ).sum()
)

n_nonresponders = int(
    (
        clinical["Responder_Status"]
        ==
        "Non-responder"
    ).sum()
)


print("\n" + "-" * 78)
print("RESPONDER CLASSIFICATION")
print("-" * 78)

print(
    f"\nTotal clinical participants : {len(clinical)}"
)

print(
    f"Complete MRI-DASS data      : {n_complete}"
)

print(
    f"Missing MRI-DASS data       : {n_missing}"
)

print(
    f"\nResponders (>=50%)          : {n_responders}"
)

print(
    f"Non-responders (<50%)       : {n_nonresponders}"
)


if n_missing > 0:

    print("\nParticipants with missing MRI-DASS:")

    print(
        missing_dass.to_string(
            index=False
        )
    )


# ============================================================
# 10. CHECK EXCEL REDUCTION AGREEMENT
# ============================================================

if EXISTING_REDUCTION_COL in clinical.columns:

    comparable = clinical[
        clinical[
            "Reduction_Matches_Excel"
        ].notna()
    ]

    n_compared = len(comparable)

    n_match = int(
        (
            comparable[
                "Reduction_Matches_Excel"
            ]
            ==
            True
        ).sum()
    )

    n_mismatch = int(
        (
            comparable[
                "Reduction_Matches_Excel"
            ]
            ==
            False
        ).sum()
    )

    print("\n" + "-" * 78)
    print("EXISTING REDUCTION COLUMN CHECK")
    print("-" * 78)

    print(
        f"\nParticipants compared : {n_compared}"
    )

    print(
        f"Exact/tolerance match  : {n_match}"
    )

    print(
        f"Mismatch               : {n_mismatch}"
    )

    if n_mismatch > 0:

        print("\nMISMATCHES:")

        print(
            comparable.loc[
                comparable[
                    "Reduction_Matches_Excel"
                ]
                ==
                False,
                [
                    "ID",
                    BL_DASS_COL,
                    FU_DASS_COL,
                    "Reduction_Excel",
                    "DASS_Reduction_Calculated",
                    "Reduction_Difference"
                ]
            ].to_string(
                index=False
            )
        )


# ============================================================
# 11. SHOW PARTICIPANT CLASSIFICATION
# ============================================================

display_cols = [
    "ID",
    BL_DASS_COL,
    FU_DASS_COL,
    "DASS_Reduction_Calculated",
    "DASS_Improvement_Percent",
    "Responder_Status"
]

print("\n" + "-" * 78)
print("PARTICIPANT RESPONDER CLASSIFICATION")
print("-" * 78)

print(
    clinical[
        display_cols
    ].to_string(
        index=False,
        formatters={
            "DASS_Reduction_Calculated":
                lambda x: (
                    ""
                    if pd.isna(x)
                    else f"{x:.4f}"
                ),

            "DASS_Improvement_Percent":
                lambda x: (
                    ""
                    if pd.isna(x)
                    else f"{x:.1f}%"
                )
        }
    )
)


# ============================================================
# 12. SAVE CLINICAL QC BEFORE MERGING
# ============================================================

qc_columns = [
    "ID",
    BL_DASS_COL,
    FU_DASS_COL,
    "DASS_Reduction_Calculated",
    "DASS_Improvement_Percent",
    "Reduction_Excel",
    "Reduction_Difference",
    "Reduction_Matches_Excel",
    "DASS_complete",
    "Responder_Status",
    "Responder_Binary"
]

# Add useful participant covariates if present
for optional_col in [
    "Age",
    "Sex",
    "Gender"
]:
    if optional_col in clinical.columns:
        qc_columns.append(
            optional_col
        )

clinical[
    qc_columns
].to_csv(
    CLINICAL_QC_OUT,
    index=False
)


# ============================================================
# 13. READ FIXED-33 FBC DATA
# ============================================================

print("\n" + "=" * 78)
print("READING FIXED-33 FBC DATA")
print("=" * 78)

fbc = pd.read_csv(
    FBC_FILE
)

required_fbc = [
    "ID",
    "DLPFC_Node_ID",
    "DLPFC_Node",
    "Connected_Node_ID",
    "Connected_Node",
    "BL_FBC",
    "FU_FBC",
    "Delta_FBC"
]

missing_fbc = [
    c for c in required_fbc
    if c not in fbc.columns
]

if missing_fbc:
    raise RuntimeError(
        f"Missing FBC columns: {missing_fbc}"
    )

fbc["ID"] = (
    fbc["ID"]
    .astype(str)
    .str.strip()
)


# ============================================================
# 14. VALIDATE FBC DATA BEFORE MERGE
# ============================================================

fbc_subjects = sorted(
    fbc["ID"].unique()
)

print(
    f"\nFBC participants : {len(fbc_subjects)}"
)

print(
    f"FBC rows         : {len(fbc)}"
)


# Every subject should have exactly 33 rows
counts = (
    fbc.groupby("ID")
    .size()
)

bad_counts = counts[
    counts != EXPECTED_EDGES
]

if len(bad_counts) > 0:

    raise RuntimeError(
        "Some participants do not have exactly "
        f"{EXPECTED_EDGES} FBC edges:\n"
        f"{bad_counts}"
    )


print(
    f"✓ Every FBC participant has exactly "
    f"{EXPECTED_EDGES} edges."
)


# ============================================================
# 15. VERIFY SAME 33 NODE IDs ACROSS PARTICIPANTS
# ============================================================

reference_subject = fbc_subjects[0]

reference_nodes = set(
    fbc.loc[
        fbc["ID"] == reference_subject,
        "Connected_Node_ID"
    ].astype(int)
)

for sid in fbc_subjects:

    nodes = set(
        fbc.loc[
            fbc["ID"] == sid,
            "Connected_Node_ID"
        ].astype(int)
    )

    if nodes != reference_nodes:

        raise RuntimeError(
            f"{sid}: selected node set differs "
            "from the fixed 33-node set."
        )


print(
    "✓ Same fixed 33 node IDs present "
    "for every FBC participant."
)


# ============================================================
# 16. VALIDATE DELTA AGAIN
# ============================================================

delta_check = (
    pd.to_numeric(
        fbc["FU_FBC"],
        errors="raise"
    )
    -
    pd.to_numeric(
        fbc["BL_FBC"],
        errors="raise"
    )
)

if not np.allclose(
    delta_check,
    pd.to_numeric(
        fbc["Delta_FBC"],
        errors="raise"
    ),
    rtol=1e-12,
    atol=1e-15
):

    raise RuntimeError(
        "Delta_FBC validation failed: "
        "Delta_FBC != FU_FBC - BL_FBC."
    )


print(
    "✓ Delta_FBC = FU_FBC - BL_FBC validated."
)


# ============================================================
# 17. CHECK ID OVERLAP BEFORE MERGING
# ============================================================

clinical_ids = set(
    clinical["ID"]
)

fbc_ids = set(
    fbc["ID"]
)

only_clinical = sorted(
    clinical_ids - fbc_ids
)

only_fbc = sorted(
    fbc_ids - clinical_ids
)

print("\n" + "-" * 78)
print("ID MATCH CHECK")
print("-" * 78)

print(
    f"\nClinical IDs : {len(clinical_ids)}"
)

print(
    f"FBC IDs      : {len(fbc_ids)}"
)

print(
    f"Shared IDs   : "
    f"{len(clinical_ids & fbc_ids)}"
)

print(
    f"\nClinical only: "
    f"{only_clinical if only_clinical else 'None'}"
)

print(
    f"FBC only     : "
    f"{only_fbc if only_fbc else 'None'}"
)


# ============================================================
# 18. PREPARE CLINICAL COLUMNS FOR MERGE
# ============================================================

merge_columns = [
    "ID",
    BL_DASS_COL,
    FU_DASS_COL,
    "DASS_Reduction_Calculated",
    "DASS_Improvement_Percent",
    "DASS_complete",
    "Responder_Status",
    "Responder_Binary"
]

# Add covariates/descriptive variables if present
for optional_col in [
    "Age",
    "Sex",
    "Gender"
]:
    if optional_col in clinical.columns:
        merge_columns.append(
            optional_col
        )

clinical_for_merge = clinical[
    merge_columns
].copy()


# ============================================================
# 19. MERGE
#
# many_to_one means:
# many FBC edge rows -> one clinical row per participant
# ============================================================

merged = fbc.merge(
    clinical_for_merge,
    on="ID",
    how="left",
    validate="many_to_one",
    indicator=True
)


# ============================================================
# 20. MERGE VALIDATION
# ============================================================

merge_counts = (
    merged["_merge"]
    .value_counts()
)

print("\n" + "-" * 78)
print("MERGE VALIDATION")
print("-" * 78)

print(
    f"\n{merge_counts.to_string()}"
)

if (
    merged["_merge"]
    !=
    "both"
).any():

    failed_ids = sorted(
        merged.loc[
            merged["_merge"] != "both",
            "ID"
        ].unique()
    )

    raise RuntimeError(
        "Some FBC rows did not match clinical data: "
        f"{failed_ids}"
    )

merged = merged.drop(
    columns="_merge"
)


# ============================================================
# 21. SUBJECT-LEVEL MERGE QC
# ============================================================

subject_qc_rows = []

for sid in fbc_subjects:

    sub = merged[
        merged["ID"] == sid
    ]

    status_values = (
        sub["Responder_Status"]
        .dropna()
        .unique()
    )

    if len(status_values) == 1:
        status = status_values[0]

    elif len(status_values) == 0:
        status = "Missing"

    else:
        raise RuntimeError(
            f"{sid}: more than one responder "
            "status after merge."
        )

    subject_qc_rows.append({
        "ID":
            sid,

        "N_edges":
            len(sub),

        "DASS_complete":
            bool(
                sub[
                    "DASS_complete"
                ].iloc[0]
            ),

        "Responder_Status":
            status,

        "BL_DASS":
            sub[
                BL_DASS_COL
            ].iloc[0],

        "FU_DASS":
            sub[
                FU_DASS_COL
            ].iloc[0],

        "DASS_Improvement_Percent":
            sub[
                "DASS_Improvement_Percent"
            ].iloc[0]
    })


subject_qc = pd.DataFrame(
    subject_qc_rows
)


# ============================================================
# 22. FINAL ANALYSIS ELIGIBILITY
# ============================================================

merged[
    "Eligible_for_Responder_Analysis"
] = (
    merged["DASS_complete"]
    &
    merged["Responder_Status"].notna()
)

eligible_subjects = sorted(
    merged.loc[
        merged[
            "Eligible_for_Responder_Analysis"
        ],
        "ID"
    ].unique()
)

excluded_subjects = sorted(
    set(fbc_subjects)
    -
    set(eligible_subjects)
)


# ============================================================
# 23. FINAL COUNTS
# ============================================================

eligible_data = merged[
    merged[
        "Eligible_for_Responder_Analysis"
    ]
].copy()

n_eligible = len(
    eligible_subjects
)

n_responder_final = (
    eligible_data.loc[
        eligible_data[
            "Responder_Status"
        ]
        ==
        "Responder",
        "ID"
    ]
    .nunique()
)

n_nonresponder_final = (
    eligible_data.loc[
        eligible_data[
            "Responder_Status"
        ]
        ==
        "Non-responder",
        "ID"
    ]
    .nunique()
)


expected_eligible_rows = (
    n_eligible
    *
    EXPECTED_EDGES
)

if len(eligible_data) != expected_eligible_rows:

    raise RuntimeError(
        "Eligible analysis dataset row count "
        "does not equal participants × 33."
    )


# ============================================================
# 24. SAVE OUTPUTS
# ============================================================

merged = merged.sort_values(
    [
        "ID",
        "Connected_Node_ID"
    ]
).reset_index(
    drop=True
)

merged.to_csv(
    MERGED_OUT,
    index=False
)

subject_qc.to_csv(
    SUBJECT_QC_OUT,
    index=False
)


# ============================================================
# 25. FINAL SUMMARY
# ============================================================

print("\n" + "=" * 78)
print("FINAL STEP 6A SUMMARY")
print("=" * 78)

print(
    f"\nFBC participants                    : "
    f"{len(fbc_subjects)}"
)

print(
    f"Fixed edges per participant         : "
    f"{EXPECTED_EDGES}"
)

print(
    f"Total merged rows                   : "
    f"{len(merged)}"
)

print(
    f"\nEligible for responder analysis     : "
    f"{n_eligible}"
)

print(
    f"Responders                          : "
    f"{n_responder_final}"
)

print(
    f"Non-responders                      : "
    f"{n_nonresponder_final}"
)

print(
    f"Excluded for missing clinical data  : "
    f"{len(excluded_subjects)}"
)

print(
    f"Excluded IDs                        : "
    f"{excluded_subjects if excluded_subjects else 'None'}"
)

print(
    f"\nEligible edge observations          : "
    f"{len(eligible_data)}"
)

print(
    f"Expected eligible observations      : "
    f"{expected_eligible_rows}"
)


print("\n" + "=" * 78)
print("OUTPUT FILES")
print("=" * 78)

print(
    f"\nClinical responder QC:\n"
    f"{CLINICAL_QC_OUT}"
)

print(
    f"\nMerged fixed-33 + clinical dataset:\n"
    f"{MERGED_OUT}"
)

print(
    f"\nSubject-level merge QC:\n"
    f"{SUBJECT_QC_OUT}"
)


print("\n" + "=" * 78)
print("DONE")
print("=" * 78)

print(
    "\n✓ Responder status independently calculated "
    "from MRI-day BL/FU DASS."
)

print(
    "✓ >=50% improvement classified as Responder."
)

print(
    "✓ Fixed 33-edge set was not changed."
)

print(
    "✓ Clinical data merged by participant ID."
)

print(
    "✓ No statistical testing performed."
)
