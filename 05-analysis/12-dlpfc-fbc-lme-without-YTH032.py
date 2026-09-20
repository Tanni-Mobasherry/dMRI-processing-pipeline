from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# ============================================================
# PATHS
# ============================================================

ROOT = Path.home() / "Desktop" / "start-analysis"

INPUT = (
    ROOT
    / "6-dlpfc-clinical-analysis"
    / "DLPFC_06A_fixed33_clinical_long.csv"
)

ORIGINAL_RESULTS = (
    ROOT
    / "7-dlpfc-fbc-lme"
    / "DLPFC_07_LME_all_33_edges.csv"
)

OUTDIR = ROOT / "8-dlpfc-fbc-lme-without-YTH032"
OUTDIR.mkdir(parents=True, exist_ok=True)

EXCLUDE_ID = "YTH032"

# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(INPUT)

print("=" * 80)
print("STEP 8 — LME SENSITIVITY ANALYSIS EXCLUDING YTH032")
print("=" * 80)

print(f"\nInput: {INPUT}")
print(f"Original rows: {len(df)}")
print(f"Original participants: {df['ID'].nunique()}")

# ============================================================
# VALIDATE YTH032 BEFORE EXCLUSION
# ============================================================

if EXCLUDE_ID not in set(df["ID"].astype(str)):
    raise ValueError(f"{EXCLUDE_ID} is not present in the input dataset.")

yth = df.loc[df["ID"].astype(str) == EXCLUDE_ID].copy()

print(f"\n{EXCLUDE_ID} rows before exclusion: {len(yth)}")
print(
    f"{EXCLUDE_ID} edges before exclusion: "
    f"{yth['Connected_Node_ID'].nunique()}"
)

if yth["Connected_Node_ID"].nunique() != 33:
    raise ValueError(
        f"{EXCLUDE_ID} does not have exactly 33 edges."
    )

# ============================================================
# EXCLUDE YTH032
# ============================================================

df = df.loc[df["ID"].astype(str) != EXCLUDE_ID].copy()

print("\nAfter exclusion:")
print(f"Participants: {df['ID'].nunique()}")
print(f"Edges       : {df['Connected_Node_ID'].nunique()}")
print(f"Wide rows   : {len(df)}")

if df["ID"].nunique() != 35:
    raise ValueError("Expected 35 participants after excluding YTH032.")

if df["Connected_Node_ID"].nunique() != 33:
    raise ValueError("Expected the same fixed 33 edges.")

# ============================================================
# GROUP COUNTS
# ============================================================

subject_groups = (
    df[["ID", "Responder_Status"]]
    .drop_duplicates()
)

group_counts = subject_groups["Responder_Status"].value_counts()

print("\nResponse groups:")
print(group_counts.to_string())

n_resp = int(group_counts.get("Responder", 0))
n_nonresp = int(group_counts.get("Non-responder", 0))

if n_resp != 13 or n_nonresp != 22:
    raise ValueError(
        f"Expected 13 responders and 22 non-responders; "
        f"found {n_resp} and {n_nonresp}."
    )

# ============================================================
# VALIDATE ONE ROW PER SUBJECT × EDGE
# ============================================================

dupes = df.duplicated(["ID", "Connected_Node_ID"]).sum()

if dupes != 0:
    raise ValueError(
        f"Found {dupes} duplicated ID-edge rows."
    )

counts = (
    df.groupby("ID")["Connected_Node_ID"]
    .nunique()
)

if not (counts == 33).all():
    raise ValueError(
        "Not every remaining participant has all 33 fixed edges."
    )

# ============================================================
# VALIDATE DELTA = FU - BL
# ============================================================

if "Delta_FBC" in df.columns:
    calc_delta = df["FU_FBC"] - df["BL_FBC"]

    if not np.allclose(
        df["Delta_FBC"],
        calc_delta,
        rtol=1e-12,
        atol=1e-15
    ):
        raise ValueError(
            "Delta_FBC is not equal to FU_FBC - BL_FBC."
        )

    print("\n✓ Delta_FBC validated as FU − BL.")

# ============================================================
# RESHAPE TO LONG BL/FU FORMAT
# ============================================================

id_vars = [
    "ID",
    "Connected_Node_ID",
    "Connected_Node",
    "Responder_Status",
    "Age",
    "Sex"
]

long = df.melt(
    id_vars=id_vars,
    value_vars=["BL_FBC", "FU_FBC"],
    var_name="Time",
    value_name="FBC"
)

long["Time"] = long["Time"].map(
    {
        "BL_FBC": "BL",
        "FU_FBC": "FU"
    }
)

long["Time"] = pd.Categorical(
    long["Time"],
    categories=["BL", "FU"],
    ordered=True
)

long["Responder_Status"] = pd.Categorical(
    long["Responder_Status"],
    categories=["Non-responder", "Responder"],
    ordered=True
)

# ============================================================
# EXACT PAIR VALIDATION
# ============================================================

pair_check = (
    long.groupby(
        ["ID", "Connected_Node_ID"],
        observed=True
    )["Time"]
    .agg(
        n_rows="size",
        n_times="nunique"
    )
)

invalid = pair_check[
    (pair_check["n_rows"] != 2)
    | (pair_check["n_times"] != 2)
]

if len(invalid) != 0:
    raise ValueError(
        f"Found {len(invalid)} invalid BL/FU participant-edge pairs."
    )

duplicate_long = long.duplicated(
    ["ID", "Connected_Node_ID", "Time"]
).sum()

if duplicate_long != 0:
    raise ValueError(
        f"Found {duplicate_long} duplicate ID-edge-Time rows."
    )

expected_rows = 35 * 33 * 2

if len(long) != expected_rows:
    raise ValueError(
        f"Expected {expected_rows} long rows, found {len(long)}."
    )

print("✓ Every remaining participant-edge has exactly one BL and one FU.")
print("✓ No duplicate ID-edge-Time observations.")
print(f"✓ Long rows = {len(long)} = 35 × 33 × 2")

# ============================================================
# SAVE MODEL INPUT
# ============================================================

long.to_csv(
    OUTDIR / "DLPFC_08_model_input_without_YTH032.csv",
    index=False
)

# ============================================================
# MODEL
# ============================================================

formula = (
    'FBC ~ '
    'C(Time, Treatment(reference="BL")) * '
    'C(Responder_Status, Treatment(reference="Non-responder")) '
    '+ Age + C(Sex)'
)

interaction_term = (
    'C(Time, Treatment(reference="BL"))[T.FU]:'
    'C(Responder_Status, Treatment(reference="Non-responder"))[T.Responder]'
)

edge_ids = sorted(long["Connected_Node_ID"].unique())

results = []

print("\n" + "=" * 80)
print("FITTING SAME 33 EDGE-WISE MODELS")
print("=" * 80)

for i, node_id in enumerate(edge_ids, start=1):

    edge = long.loc[
        long["Connected_Node_ID"] == node_id
    ].copy()

    node_name = edge["Connected_Node"].iloc[0]

    # Actual observed FU-BL changes
    wide_edge = df.loc[
        df["Connected_Node_ID"] == node_id
    ].copy()

    wide_edge["Observed_Delta"] = (
        wide_edge["FU_FBC"] - wide_edge["BL_FBC"]
    )

    r_delta = wide_edge.loc[
        wide_edge["Responder_Status"].astype(str) == "Responder",
        "Observed_Delta"
    ]

    nr_delta = wide_edge.loc[
        wide_edge["Responder_Status"].astype(str) == "Non-responder",
        "Observed_Delta"
    ]

    fitted = None
    used_optimizer = None
    warning_messages = []

    for optimizer in ["lbfgs", "powell", "cg"]:

        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")

                model = smf.mixedlm(
                    formula,
                    data=edge,
                    groups=edge["ID"],
                    re_formula="1"
                )

                fit = model.fit(
                    reml=True,
                    method=optimizer,
                    disp=False
                )

                warning_messages = [
                    str(w.message)
                    for w in caught
                ]

            fitted = fit
            used_optimizer = optimizer

            if fit.converged:
                break

        except Exception as exc:
            warning_messages.append(
                f"{optimizer}: {exc}"
            )

    if fitted is None:
        raise RuntimeError(
            f"Model failed completely for Node {node_id} {node_name}"
        )

    if interaction_term not in fitted.params.index:
        raise RuntimeError(
            f"Interaction term not found for Node {node_id}."
        )

    beta = float(fitted.params[interaction_term])
    se = float(fitted.bse[interaction_term])
    p = float(fitted.pvalues[interaction_term])

    ci = fitted.conf_int().loc[interaction_term]
    ci_low = float(ci.iloc[0])
    ci_high = float(ci.iloc[1])

    # Direct observed difference in mean FU-BL change
    observed_difference = (
        float(r_delta.mean())
        - float(nr_delta.mean())
    )

    # This should match interaction beta very closely
    beta_difference = beta - observed_difference

    results.append(
        {
            "Connected_Node_ID": int(node_id),
            "Connected_Node": node_name,

            "N_total": wide_edge["ID"].nunique(),
            "N_Responder": len(r_delta),
            "N_NonResponder": len(nr_delta),

            "Mean_Delta_Responder": float(r_delta.mean()),
            "Mean_Delta_NonResponder": float(nr_delta.mean()),
            "Median_Delta_Responder": float(r_delta.median()),
            "Median_Delta_NonResponder": float(nr_delta.median()),

            "Observed_Mean_Delta_Difference_R_minus_NR":
                observed_difference,

            "Interaction_Beta": beta,
            "Interaction_SE": se,
            "CI95_Lower": ci_low,
            "CI95_Upper": ci_high,
            "Raw_P": p,

            "Beta_minus_direct_delta_difference":
                beta_difference,

            "Converged": bool(fitted.converged),
            "Optimizer": used_optimizer,
            "Model_Warnings": " | ".join(warning_messages)
        }
    )

    print(
        f"[{i:02d}/33] Node {node_id:3d} "
        f"{node_name:<38} "
        f"beta={beta:+.6f} "
        f"p={p:.4f} "
        f"converged={fitted.converged}"
    )

# ============================================================
# RESULTS + BH-FDR
# ============================================================

res = pd.DataFrame(results)

reject, qvals, _, _ = multipletests(
    res["Raw_P"].values,
    alpha=0.05,
    method="fdr_bh"
)

res["FDR_Q"] = qvals
res["FDR_Significant"] = reject

res = res.sort_values(
    "Raw_P"
).reset_index(drop=True)

res.to_csv(
    OUTDIR / "DLPFC_08_LME_without_YTH032_all_33_edges.csv",
    index=False
)

res.loc[
    res["FDR_Significant"]
].to_csv(
    OUTDIR / "DLPFC_08_LME_without_YTH032_FDR_significant.csv",
    index=False
)

# ============================================================
# COMPARE AGAINST ORIGINAL N=36 ANALYSIS
# ============================================================

comparison_saved = False

if ORIGINAL_RESULTS.exists():

    original = pd.read_csv(ORIGINAL_RESULTS)

    # Try to identify the relevant columns from Step 7
    possible_beta = [
        "Interaction_Beta",
        "Beta_Time_x_Group",
        "Interaction_Estimate"
    ]

    possible_p = [
        "Raw_P",
        "P_Time_x_Group",
        "Interaction_P"
    ]

    possible_q = [
        "FDR_Q",
        "FDR_P",
        "Interaction_FDR"
    ]

    def first_existing(columns, candidates):
        for c in candidates:
            if c in columns:
                return c
        return None

    beta_col = first_existing(
        original.columns,
        possible_beta
    )
    p_col = first_existing(
        original.columns,
        possible_p
    )
    q_col = first_existing(
        original.columns,
        possible_q
    )

    keep = [
        "Connected_Node_ID",
        "Connected_Node"
    ]

    rename = {}

    if beta_col:
        keep.append(beta_col)
        rename[beta_col] = "Interaction_Beta_N36"

    if p_col:
        keep.append(p_col)
        rename[p_col] = "Raw_P_N36"

    if q_col:
        keep.append(q_col)
        rename[q_col] = "FDR_Q_N36"

    original_small = (
        original[keep]
        .rename(columns=rename)
        .copy()
    )

    sensitivity_small = res[
        [
            "Connected_Node_ID",
            "Connected_Node",
            "Interaction_Beta",
            "Raw_P",
            "FDR_Q"
        ]
    ].rename(
        columns={
            "Interaction_Beta": "Interaction_Beta_N35",
            "Raw_P": "Raw_P_N35",
            "FDR_Q": "FDR_Q_N35"
        }
    )

    comp = original_small.merge(
        sensitivity_small,
        on=[
            "Connected_Node_ID",
            "Connected_Node"
        ],
        how="outer",
        validate="one_to_one"
    )

    if (
        "Interaction_Beta_N36" in comp.columns
        and "Interaction_Beta_N35" in comp.columns
    ):
        comp["Beta_Change_after_excluding_YTH032"] = (
            comp["Interaction_Beta_N35"]
            - comp["Interaction_Beta_N36"]
        )

        comp["Beta_Direction_N36"] = np.sign(
            comp["Interaction_Beta_N36"]
        )

        comp["Beta_Direction_N35"] = np.sign(
            comp["Interaction_Beta_N35"]
        )

        comp["Direction_Flipped"] = (
            comp["Beta_Direction_N36"]
            != comp["Beta_Direction_N35"]
        )

    comp.to_csv(
        OUTDIR / "DLPFC_08_compare_N36_vs_N35.csv",
        index=False
    )

    comparison_saved = True

# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("FINAL SENSITIVITY SUMMARY")
print("=" * 80)

print(f"Participants analysed : {df['ID'].nunique()}")
print(f"Responders            : {n_resp}")
print(f"Non-responders        : {n_nonresp}")
print(f"Edges analysed        : {len(res)}")
print(
    "Models converged      : "
    f"{int(res['Converged'].sum())} / {len(res)}"
)

print(
    "Raw p < .05           : "
    f"{int((res['Raw_P'] < 0.05).sum())} / {len(res)}"
)

print(
    "FDR q < .05           : "
    f"{int((res['FDR_Q'] < 0.05).sum())} / {len(res)}"
)

best = res.iloc[0]

print("\nSmallest raw p:")
print(
    f"Node {int(best['Connected_Node_ID'])} "
    f"{best['Connected_Node']}"
)
print(
    f"beta = {best['Interaction_Beta']:+.6f}"
)
print(
    f"95% CI = "
    f"[{best['CI95_Lower']:+.6f}, "
    f"{best['CI95_Upper']:+.6f}]"
)
print(
    f"p = {best['Raw_P']:.6f}"
)
print(
    f"q = {best['FDR_Q']:.6f}"
)

max_beta_check = (
    res["Beta_minus_direct_delta_difference"]
    .abs()
    .max()
)

print(
    "\nMaximum |LME beta - direct difference "
    "in mean FU-BL change|:"
)
print(f"{max_beta_check:.12g}")

print("\nOutputs:")
print(
    OUTDIR
    / "DLPFC_08_LME_without_YTH032_all_33_edges.csv"
)

if comparison_saved:
    print(
        OUTDIR
        / "DLPFC_08_compare_N36_vs_N35.csv"
    )

print("\nIMPORTANT:")
print("This is a sensitivity analysis.")
print("The original N=36 model remains the primary analysis.")
print("YTH032 has not been deleted from any source dataset.")

print("\nDONE")
