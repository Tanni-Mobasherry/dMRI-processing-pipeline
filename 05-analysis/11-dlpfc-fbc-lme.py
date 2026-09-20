'''
 The code refits the LME for each of the 33 edges:
The model is:
(FBC ~ Time * ResponseGroup + Age + Sex + (1|Participant)\)

\]Then it calculates/checks:
- Model convergence and optimiser used.
- Random-intercept variance and residual variance.
- ICC = proportion of variance attributable to between-participant differences.
- Standardized residuals, including counts with \(|r|>2\) and \(|r|>3\).
- Shapiro–Wilk test of residual normality — used only diagnostically, not as a pass/fail rule.
- Residual–fitted correlation.
- Residual-vs-fitted plots to inspect model fit/heteroscedasticity or patterns.
- Q–Q plots to visually assess residual normality.
- Checks that every participant has exactly one BL and one FU measurement for every edge.
'''
from pathlib import Path
import warnings
import numpy as np
import pandas as pd

import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


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
    / "7-dlpfc-fbc-lme"
)

OUT_DIR.mkdir(parents=True, exist_ok=True)

RESULTS_OUT = OUT_DIR / "DLPFC_07_LME_all_33_edges.csv"
SIGNIFICANT_OUT = OUT_DIR / "DLPFC_07_LME_FDR_significant_edges.csv"
MODEL_LONG_OUT = OUT_DIR / "DLPFC_07_model_input_long.csv"
QC_OUT = OUT_DIR / "DLPFC_07_model_QC.csv"


# ============================================================
# SETTINGS
# ============================================================

EXPECTED_SUBJECTS = 36
EXPECTED_RESPONDERS = 13
EXPECTED_NONRESPONDERS = 23
EXPECTED_EDGES = 33

ALPHA = 0.05


# ============================================================
# 1. READ DATA
# ============================================================

print("\n" + "=" * 82)
print("STEP 7 — EDGE-WISE LINEAR MIXED-EFFECTS MODELS")
print("=" * 82)

df = pd.read_csv(INPUT)

required = [
    "ID",
    "Connected_Node_ID",
    "Connected_Node",
    "BL_FBC",
    "FU_FBC",
    "Responder_Status",
    "Age",
    "Sex"
]

missing = [c for c in required if c not in df.columns]

if missing:
    raise RuntimeError(
        f"Missing required columns: {missing}"
    )


# ============================================================
# 2. BASIC CLEANING
# ============================================================

df["ID"] = df["ID"].astype(str).str.strip()

df["Responder_Status"] = (
    df["Responder_Status"]
    .astype(str)
    .str.strip()
)

df["Connected_Node_ID"] = pd.to_numeric(
    df["Connected_Node_ID"],
    errors="raise"
).astype(int)

df["Age"] = pd.to_numeric(
    df["Age"],
    errors="coerce"
)

df["BL_FBC"] = pd.to_numeric(
    df["BL_FBC"],
    errors="coerce"
)

df["FU_FBC"] = pd.to_numeric(
    df["FU_FBC"],
    errors="coerce"
)


# ============================================================
# 3. ELIGIBILITY
# ============================================================

if "Eligible_for_Responder_Analysis" in df.columns:

    eligibility = (
        df["Eligible_for_Responder_Analysis"]
        .astype(str)
        .str.lower()
        .map({
            "true": True,
            "false": False
        })
    )

    df = df.loc[
        eligibility == True
    ].copy()


# ============================================================
# 4. VALIDATE SUBJECT/GROUP COUNTS
# ============================================================

subject_info = (
    df[
        [
            "ID",
            "Responder_Status",
            "Age",
            "Sex"
        ]
    ]
    .drop_duplicates()
)

if subject_info["ID"].duplicated().any():

    print("\nChecking repeated participant covariates...")

    for sid, g in df.groupby("ID"):

        if g["Responder_Status"].nunique() != 1:
            raise RuntimeError(
                f"{sid}: inconsistent responder status."
            )

        if g["Age"].nunique(dropna=False) != 1:
            raise RuntimeError(
                f"{sid}: inconsistent Age."
            )

        if g["Sex"].nunique(dropna=False) != 1:
            raise RuntimeError(
                f"{sid}: inconsistent Sex."
            )

    subject_info = (
        df.groupby("ID", as_index=False)
        .first()[
            [
                "ID",
                "Responder_Status",
                "Age",
                "Sex"
            ]
        ]
    )


n_subjects = subject_info["ID"].nunique()

n_resp = (
    subject_info[
        "Responder_Status"
    ]
    .eq("Responder")
    .sum()
)

n_nonresp = (
    subject_info[
        "Responder_Status"
    ]
    .eq("Non-responder")
    .sum()
)

n_edges = df["Connected_Node_ID"].nunique()


print(f"\nParticipants     : {n_subjects}")
print(f"Responders       : {n_resp}")
print(f"Non-responders   : {n_nonresp}")
print(f"Fixed edges      : {n_edges}")


if n_subjects != EXPECTED_SUBJECTS:
    raise RuntimeError(
        f"Expected {EXPECTED_SUBJECTS} participants, "
        f"found {n_subjects}."
    )

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

if n_edges != EXPECTED_EDGES:
    raise RuntimeError(
        f"Expected {EXPECTED_EDGES} edges, "
        f"found {n_edges}."
    )


# ============================================================
# 5. CHECK MISSING AGE / SEX
# ============================================================

print("\nAge missing:")
print(subject_info["Age"].isna().sum())

print("\nSex values:")
print(
    subject_info["Sex"]
    .value_counts(dropna=False)
    .to_string()
)

if subject_info["Age"].isna().any():
    raise RuntimeError(
        "Missing Age values found."
    )

if subject_info["Sex"].isna().any():
    raise RuntimeError(
        "Missing Sex values found."
    )


# ============================================================
# 6. RESHAPE BL/FU INTO LONG FORMAT
# ============================================================

base_cols = [
    "ID",
    "Connected_Node_ID",
    "Connected_Node",
    "Responder_Status",
    "Age",
    "Sex"
]

bl = df[
    base_cols + ["BL_FBC"]
].copy()

bl = bl.rename(
    columns={
        "BL_FBC": "FBC"
    }
)

bl["Time"] = "BL"


fu = df[
    base_cols + ["FU_FBC"]
].copy()

fu = fu.rename(
    columns={
        "FU_FBC": "FBC"
    }
)

fu["Time"] = "FU"


long = pd.concat(
    [bl, fu],
    ignore_index=True
)


# ============================================================
# 7. EXPLICIT CATEGORY CODING
#
# Reference Time  = BL
# Reference Group = Non-responder
# ============================================================

long["Time"] = pd.Categorical(
    long["Time"],
    categories=[
        "BL",
        "FU"
    ],
    ordered=True
)

long["Responder_Status"] = pd.Categorical(
    long["Responder_Status"],
    categories=[
        "Non-responder",
        "Responder"
    ],
    ordered=True
)


# ============================================================
# 8. VALIDATE MODEL DATA
# ============================================================

if long["FBC"].isna().any():
    raise RuntimeError(
        "Missing FBC values found."
    )

if not np.isfinite(
    long["FBC"].astype(float)
).all():
    raise RuntimeError(
        "Non-finite FBC values found."
    )


expected_rows = (
    EXPECTED_SUBJECTS
    * EXPECTED_EDGES
    * 2
)

if len(long) != expected_rows:
    raise RuntimeError(
        f"Expected {expected_rows} long-format rows, "
        f"found {len(long)}."
    )


# Each subject-edge should have BL and FU
counts = (
    long.groupby(
        [
            "ID",
            "Connected_Node_ID"
        ],
        observed=True
    )
    .size()
)

if not (counts == 2).all():
    raise RuntimeError(
        "Not every subject-edge combination "
        "contains exactly BL and FU."
    )


print(
    f"\n✓ Model long-format rows: {len(long)}"
)

print(
    "✓ Every participant-edge has BL and FU."
)

print(
    "✓ Reference group = Non-responder."
)

print(
    "✓ Reference time = Baseline."
)


long.to_csv(
    MODEL_LONG_OUT,
    index=False
)


# ============================================================
# 9. MODEL FORMULA
#
# FBC ~ Time * ResponseGroup + Age + Sex
# random intercept = Participant
#
# Explicit Treatment coding ensures:
# BL = Time reference
# Non-responder = Group reference
# ============================================================

formula = (
    'FBC ~ '
    'C(Time, Treatment(reference="BL")) * '
    'C(Responder_Status, Treatment(reference="Non-responder")) '
    '+ Age + C(Sex)'
)


# ============================================================
# 10. FIT ONE MODEL PER EDGE
# ============================================================

results = []
qc_rows = []

edge_ids = sorted(
    long["Connected_Node_ID"].unique()
)

print("\n" + "=" * 82)
print("FITTING 33 EDGE-WISE MIXED MODELS")
print("=" * 82)


for number, node_id in enumerate(
    edge_ids,
    start=1
):

    edge = long.loc[
        long["Connected_Node_ID"] == node_id
    ].copy()

    names = edge[
        "Connected_Node"
    ].dropna().unique()

    if len(names) != 1:
        raise RuntimeError(
            f"Node {node_id}: inconsistent node names."
        )

    node_name = names[0]

    print(
        f"\n[{number:02d}/{EXPECTED_EDGES}] "
        f"Node {node_id}: {node_name}"
    )

    # --------------------------------------------------------
    # Check expected rows
    # --------------------------------------------------------

    expected_edge_rows = (
        EXPECTED_SUBJECTS * 2
    )

    if len(edge) != expected_edge_rows:
        raise RuntimeError(
            f"Node {node_id}: expected "
            f"{expected_edge_rows} observations, "
            f"found {len(edge)}."
        )


    # --------------------------------------------------------
    # Fit MixedLM
    # --------------------------------------------------------

    fit_success = False
    fit = None
    used_method = None
    error_messages = []

    # Try several optimizers if necessary.
    for method in [
        "lbfgs",
        "powell",
        "cg"
    ]:

        try:

            with warnings.catch_warnings(
                record=True
            ) as caught:

                warnings.simplefilter(
                    "always"
                )

                model = smf.mixedlm(
                    formula,
                    data=edge,
                    groups=edge["ID"],
                    re_formula="1"
                )

                fit = model.fit(
                    reml=True,
                    method=method,
                    maxiter=2000,
                    disp=False
                )

            fit_success = True
            used_method = method

            warning_text = " | ".join(
                str(w.message)
                for w in caught
            )

            break

        except Exception as exc:

            error_messages.append(
                f"{method}: {exc}"
            )


    if not fit_success:

        qc_rows.append({
            "Connected_Node_ID":
                node_id,

            "Connected_Node":
                node_name,

            "Fit_Success":
                False,

            "Optimizer":
                None,

            "Converged":
                False,

            "Warnings":
                "",

            "Errors":
                " || ".join(
                    error_messages
                )
        })

        print("   ✗ MODEL FAILED")

        continue


    # --------------------------------------------------------
    # Identify interaction coefficient
    # --------------------------------------------------------

    interaction_candidates = [
        term
        for term in fit.params.index
        if (
            "Time" in term
            and
            "Responder_Status" in term
            and
            ":" in term
        )
    ]

    if len(interaction_candidates) != 1:

        raise RuntimeError(
            f"Node {node_id}: expected exactly one "
            f"Time × Group interaction term; found "
            f"{interaction_candidates}"
        )

    interaction = interaction_candidates[0]


    # --------------------------------------------------------
    # Extract interaction statistics
    # --------------------------------------------------------

    beta = float(
        fit.params[
            interaction
        ]
    )

    se = float(
        fit.bse[
            interaction
        ]
    )

    p = float(
        fit.pvalues[
            interaction
        ]
    )

    ci = fit.conf_int().loc[
        interaction
    ]

    ci_low = float(ci.iloc[0])
    ci_high = float(ci.iloc[1])


    # --------------------------------------------------------
    # Also extract fixed effects for interpretation
    # --------------------------------------------------------

    fixed_names = list(
        fit.fe_params.index
    )

    time_candidates = [
        x for x in fixed_names
        if (
            "Time" in x
            and
            "Responder_Status" not in x
        )
    ]

    group_candidates = [
        x for x in fixed_names
        if (
            "Responder_Status" in x
            and
            "Time" not in x
        )
    ]

    time_term = (
        time_candidates[0]
        if len(time_candidates) == 1
        else None
    )

    group_term = (
        group_candidates[0]
        if len(group_candidates) == 1
        else None
    )


    time_beta = (
        float(fit.params[time_term])
        if time_term is not None
        else np.nan
    )

    group_beta = (
        float(fit.params[group_term])
        if group_term is not None
        else np.nan
    )


    # --------------------------------------------------------
    # Descriptive actual changes
    # --------------------------------------------------------

    original_edge = df.loc[
        df["Connected_Node_ID"] == node_id
    ].copy()

    resp_delta = (
        original_edge.loc[
            original_edge[
                "Responder_Status"
            ]
            ==
            "Responder",
            "FU_FBC"
        ].to_numpy()
        -
        original_edge.loc[
            original_edge[
                "Responder_Status"
            ]
            ==
            "Responder",
            "BL_FBC"
        ].to_numpy()
    )

    nonresp_delta = (
        original_edge.loc[
            original_edge[
                "Responder_Status"
            ]
            ==
            "Non-responder",
            "FU_FBC"
        ].to_numpy()
        -
        original_edge.loc[
            original_edge[
                "Responder_Status"
            ]
            ==
            "Non-responder",
            "BL_FBC"
        ].to_numpy()
    )


    results.append({

        "Connected_Node_ID":
            node_id,

        "Connected_Node":
            node_name,

        "Connected_Node_ID_Name":
            f"{node_id}|{node_name}",

        "N_Total":
            EXPECTED_SUBJECTS,

        "N_Responder":
            EXPECTED_RESPONDERS,

        "N_NonResponder":
            EXPECTED_NONRESPONDERS,

        # Descriptive change
        "Responder_Mean_Delta_FBC":
            np.mean(resp_delta),

        "NonResponder_Mean_Delta_FBC":
            np.mean(nonresp_delta),

        "Responder_Median_Delta_FBC":
            np.median(resp_delta),

        "NonResponder_Median_Delta_FBC":
            np.median(nonresp_delta),

        # Fixed effects
        "Time_FU_Beta":
            time_beta,

        "Group_Responder_Beta":
            group_beta,

        # Main result
        "Time_x_Group_Beta":
            beta,

        "Time_x_Group_SE":
            se,

        "Time_x_Group_CI95_Lower":
            ci_low,

        "Time_x_Group_CI95_Upper":
            ci_high,

        "Time_x_Group_P":
            p,

        "Converged":
            bool(fit.converged),

        "Optimizer":
            used_method
    })


    qc_rows.append({

        "Connected_Node_ID":
            node_id,

        "Connected_Node":
            node_name,

        "Fit_Success":
            True,

        "Optimizer":
            used_method,

        "Converged":
            bool(fit.converged),

        "Warnings":
            warning_text,

        "Errors":
            ""
    })


    print(
        f"   β interaction = {beta:+.8f}"
    )

    print(
        f"   95% CI        = "
        f"[{ci_low:+.8f}, {ci_high:+.8f}]"
    )

    print(
        f"   raw p         = {p:.6g}"
    )

    print(
        f"   converged     = {fit.converged}"
    )


# ============================================================
# 11. SAVE QC
# ============================================================

qc = pd.DataFrame(qc_rows)

qc.to_csv(
    QC_OUT,
    index=False
)


# ============================================================
# 12. REQUIRE ALL 33 MODELS
# ============================================================

results = pd.DataFrame(results)

if len(results) != EXPECTED_EDGES:

    raise RuntimeError(
        f"Only {len(results)} of "
        f"{EXPECTED_EDGES} models were successfully fitted. "
        f"Check:\n{QC_OUT}"
    )


if not results["Converged"].all():

    bad = results.loc[
        ~results["Converged"],
        [
            "Connected_Node_ID",
            "Connected_Node",
            "Optimizer"
        ]
    ]

    print(
        "\nWARNING: Some models report "
        "non-convergence:"
    )

    print(
        bad.to_string(index=False)
    )


# ============================================================
# 13. FDR CORRECTION
#
# Benjamini-Hochberg across exactly the 33
# Time × ResponseGroup interaction tests.
# ============================================================

reject, p_fdr, _, _ = multipletests(
    results[
        "Time_x_Group_P"
    ].to_numpy(),
    alpha=ALPHA,
    method="fdr_bh"
)

results[
    "Time_x_Group_P_FDR"
] = p_fdr

results[
    "FDR_Significant"
] = reject


# ============================================================
# 14. ADD INTERPRETATION OF BETA DIRECTION
# ============================================================

results[
    "Interaction_Direction"
] = np.where(
    results[
        "Time_x_Group_Beta"
    ] > 0,

    "Responder_change_more_positive",

    np.where(
        results[
            "Time_x_Group_Beta"
        ] < 0,

        "Responder_change_more_negative",

        "No_difference"
    )
)


# ============================================================
# 15. SORT BY FDR P
# ============================================================

results = results.sort_values(
    [
        "Time_x_Group_P_FDR",
        "Time_x_Group_P"
    ]
).reset_index(drop=True)


results.insert(
    0,
    "FDR_Rank",
    np.arange(
        1,
        len(results) + 1
    )
)


# ============================================================
# 16. SAVE RESULTS
# ============================================================

results.to_csv(
    RESULTS_OUT,
    index=False
)

significant = results.loc[
    results[
        "FDR_Significant"
    ]
].copy()

significant.to_csv(
    SIGNIFICANT_OUT,
    index=False
)


# ============================================================
# 17. TERMINAL RESULTS
# ============================================================

print("\n" + "=" * 82)
print("TIME × RESPONSE GROUP RESULTS — SORTED BY FDR")
print("=" * 82)

show_cols = [
    "FDR_Rank",
    "Connected_Node_ID",
    "Connected_Node",
    "Time_x_Group_Beta",
    "Time_x_Group_CI95_Lower",
    "Time_x_Group_CI95_Upper",
    "Time_x_Group_P",
    "Time_x_Group_P_FDR",
    "FDR_Significant"
]

print(
    results[
        show_cols
    ].to_string(
        index=False,
        formatters={
            "Time_x_Group_Beta":
                lambda x: f"{x:+.8f}",

            "Time_x_Group_CI95_Lower":
                lambda x: f"{x:+.8f}",

            "Time_x_Group_CI95_Upper":
                lambda x: f"{x:+.8f}",

            "Time_x_Group_P":
                lambda x: f"{x:.6g}",

            "Time_x_Group_P_FDR":
                lambda x: f"{x:.6g}"
        }
    )
)


# ============================================================
# 18. FINAL SUMMARY
# ============================================================

n_raw = int(
    (
        results[
            "Time_x_Group_P"
        ]
        <
        ALPHA
    ).sum()
)

n_fdr = int(
    results[
        "FDR_Significant"
    ].sum()
)


print("\n" + "=" * 82)
print("FINAL SUMMARY")
print("=" * 82)

print(
    f"\nModels fitted             : "
    f"{len(results)} / {EXPECTED_EDGES}"
)

print(
    f"Raw p < .05              : "
    f"{n_raw}"
)

print(
    f"FDR q < .05              : "
    f"{n_fdr}"
)

print(
    "\nPrimary effect tested:"
)

print(
    "Time × Response Group interaction"
)

print(
    "\nInterpretation:"
)

print(
    "Positive β = responders show a more positive "
    "BL→FU FBC change than non-responders."
)

print(
    "Negative β = responders show a more negative "
    "BL→FU FBC change than non-responders."
)


print("\n" + "=" * 82)
print("OUTPUT FILES")
print("=" * 82)

print(
    f"\nAll 33 LME results:\n"
    f"{RESULTS_OUT}"
)

print(
    f"\nFDR-significant edges:\n"
    f"{SIGNIFICANT_OUT}"
)

print(
    f"\nModel long-format input:\n"
    f"{MODEL_LONG_OUT}"
)

print(
    f"\nModel QC:\n"
    f"{QC_OUT}"
)


print("\n" + "=" * 82)
print("DONE")
print("=" * 82)

print(
    "\n✓ 33 separate mixed-effects models fitted."
)

print(
    "✓ BL is the Time reference."
)

print(
    "✓ Non-responder is the Group reference."
)

print(
    "✓ Participant random intercept included."
)

print(
    "✓ Age and Sex included as covariates."
)

print(
    "✓ FDR correction applied only to the "
    "33 Time × Group interaction p-values."
)
