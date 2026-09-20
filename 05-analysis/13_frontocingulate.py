#!/usr/bin/env python3
"""
Step 10 - Independent analysis of the 6 fronto-cingulate DLPFC edges.

Built directly from the 36 per-subject files
    fronto-cingulate/<ID>_frontocingulate_edges.csv
and validated against each participant's raw labelled connectivity matrices.

Edges: left-DLPFC (node 15) to the anterior, mid-anterior and mid-posterior
cingulate in BOTH hemispheres - a 3 (region) x 2 (hemisphere) homologous set.

THREE QUESTIONS, each its own multiplicity family:

  A (primary, 6 tests)  Is there ANY longitudinal BL->FU FBC change on each edge?
                        Within-participant one-sample contrast on dFBC.
  B (primary, 6 tests)  Does that change differ between responders and
                        non-responders? Same interaction model as the 33-edge
                        analysis, so the two are directly comparable.
  C (exploratory, 3)    Does the change differ between homologous right and left
                        edges? Paired within participant; no edge is averaged away.

Inference does not rest on normal theory: every test is accompanied by an exact
or Monte-Carlo permutation test and a rank-based test. FDR is Benjamini-Hochberg
within each family.
"""
from pathlib import Path
import glob
import warnings
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

ROOT = Path.home() / "Desktop" / "start-analysis"
OUT = ROOT / "independent-fbc-analysis" / "frontocingulate"
SRC = ROOT / "fronto-cingulate"
SEED = 20260920
N_PERM = 20000
DLPFC_NODE = 15

# homologous left/right pairs: (region label, left node, right node)
HOMOLOGUES = [("Anterior cingulate", 6, 95),
              ("Mid-anterior cingulate", 7, 96),
              ("Mid-posterior cingulate", 8, 97)]


def load_matrix(path):
    m = pd.read_csv(path, index_col=0)
    m.index = [int(str(i).split("|")[0]) for i in m.index]
    m.columns = [int(str(c).split("|")[0]) for c in m.columns]
    return m.astype(float)


def edge_value(m, a, b):
    """Matrices are upper-triangular: the value lives at [min, max]."""
    return m.loc[min(a, b), max(a, b)]


def sign_flip_p(x, n_perm, rng):
    """Exact-if-small, else Monte-Carlo sign-flip test of H0: median change = 0."""
    obs = abs(x.mean())
    n = len(x)
    if n <= 20:                                   # exact enumeration
        signs = np.array(np.meshgrid(*[[-1, 1]] * n)).T.reshape(-1, n)
    else:
        signs = rng.choice([-1, 1], size=(n_perm, n))
    null = np.abs((signs * x).mean(axis=1))
    return (1.0 + np.sum(null >= obs - 1e-15)) / (len(null) + 1.0)


def label_perm_p(delta, is_resp, n_perm, rng):
    obs = abs(delta[is_resp].mean() - delta[~is_resp].mean())
    idx = np.array([rng.permutation(len(delta)) for _ in range(n_perm)])
    dp = delta[idx]
    null = np.abs(dp[:, is_resp].mean(axis=1) - dp[:, ~is_resp].mean(axis=1))
    return (1.0 + np.sum(null >= obs - 1e-15)) / (n_perm + 1.0)


def build():
    """Read all 36 per-subject files, QC them, verify against raw matrices."""
    files = sorted(glob.glob(str(SRC / "*_frontocingulate_edges.csv")))
    frames, checks = [], []
    for f in files:
        df = pd.read_csv(f)
        stem_id = Path(f).name.split("_")[0]
        assert df["ID"].nunique() == 1 and df["ID"].iloc[0] == stem_id, f
        frames.append(df)
    d = pd.concat(frames, ignore_index=True)

    def rec(cid, desc, exp, obs, ok):
        checks.append(dict(check_id=cid, description=desc, expected=str(exp),
                           observed=str(obs), status="PASS" if ok else "FAIL"))

    rec(1, "Per-subject files read", 36, len(files), len(files) == 36)
    rec(2, "Total participant-edge rows", 216, len(d), len(d) == 216)
    rec(3, "Participants", 36, d["ID"].nunique(), d["ID"].nunique() == 36)
    rec(4, "Distinct fronto-cingulate edges", 6, d["Target_Node_ID"].nunique(),
        d["Target_Node_ID"].nunique() == 6)
    same = d.groupby("ID")["Target_Node_ID"].apply(lambda x: tuple(sorted(x))).nunique() == 1
    rec(5, "Identical 6-edge set for every participant", "1 unique set",
        f"identical={bool(same)}", bool(same))
    rec(6, "All edges seed from DLPFC node 15", "single seed",
        sorted(d["DLPFC_Node_ID"].unique()), d["DLPFC_Node_ID"].nunique() == 1)
    ndup = int(d.duplicated(["ID", "Target_Node_ID"]).sum())
    rec(7, "No duplicated participant-edge rows", 0, ndup, ndup == 0)
    nna = int(d[["BL_FBC", "FU_FBC"]].isna().sum().sum())
    rec(8, "No missing BL/FU values", 0, nna, nna == 0)
    pos = bool((d["BL_FBC"] > 0).all() and (d["FU_FBC"] > 0).all())
    rec(9, "All FBC strictly positive",
        "all > 0", f"BL min={d.BL_FBC.min():.5f}, FU min={d.FU_FBC.min():.5f}", pos)
    dok = bool(np.allclose(d["Delta_FBC"], d["FU_FBC"] - d["BL_FBC"]))
    rec(10, "Stored Delta_FBC equals FU_FBC - BL_FBC", "exact", f"match={dok}", dok)
    dirok = bool(((d["Delta_FBC"] > 0) == (d["Direction"] == "Increase")).all())
    rec(11, "Direction column consistent with the sign of Delta", "consistent",
        f"consistent={dirok}", dirok)

    # provenance against the raw per-subject matrices
    maxdev = 0.0
    for pid, g in d.groupby("ID"):
        bl = load_matrix(ROOT / pid / "BL_fbc_labelled.csv")
        fu = load_matrix(ROOT / pid / "FU_fbc_labelled.csv")
        for _, row in g.iterrows():
            t = int(row["Target_Node_ID"])
            maxdev = max(maxdev,
                         abs(edge_value(bl, DLPFC_NODE, t) - row["BL_FBC"]),
                         abs(edge_value(fu, DLPFC_NODE, t) - row["FU_FBC"]))
    rec(12, "BL/FU traced to raw <ID>/BL_fbc_labelled.csv and FU_fbc_labelled.csv",
        "max deviation < 1e-9", f"{maxdev:.3e}", maxdev < 1e-9)

    # clinical merge, reusing the already-validated participant-level table
    clin = (pd.read_csv(ROOT / "independent-fbc-analysis" / "qc"
                        / "ANALYSIS_dataset_wide_participant_edge.csv")
            .drop_duplicates("ID")[["ID", "Responder_Status", "Age", "Sex",
                                    "DASS_Improvement_Percent"]])
    d = d.merge(clin, on="ID", how="left", validate="many_to_one")
    rec(13, "Clinical data merged for every participant", "0 unmatched",
        int(d["Responder_Status"].isna().sum()), int(d["Responder_Status"].isna().sum()) == 0)
    vc = d.drop_duplicates("ID")["Responder_Status"].value_counts()
    rec(14, "Responder / Non-responder counts", "13 / 23",
        f"{vc.get('Responder',0)} / {vc.get('Non-responder',0)}",
        vc.get("Responder", 0) == 13 and vc.get("Non-responder", 0) == 23)

    # Was this edge already selected by the Step-04 median +/- 3 x IQR screen,
    # which selected on the edge-level MEDIAN of participant Delta_FBC? If so,
    # testing "is the longitudinal change non-zero" on that edge is circular.
    sel = pd.read_csv(ROOT / "4-dlpfc-median-iqr"
                      / "DLPFC_04_all_edge_group_statistics.csv")
    sel = sel.set_index("Connected_Node_ID")[["Extreme", "Extreme_Direction",
                                              "Median_Delta_FBC"]]
    d["Preselected_by_Step04"] = d["Target_Node_ID"].map(
        sel["Extreme"].astype(bool)).fillna(False)
    d["Step04_selection_direction"] = d["Target_Node_ID"].map(sel["Extreme_Direction"])
    n_pre = int(d.drop_duplicates("Target_Node_ID")["Preselected_by_Step04"].sum())
    rec(15, "Edges also selected by the Step-04 extreme-change screen "
            "(circularity flag for question A)",
        "reported, not a pass/fail gate",
        f"{n_pre}/6 pre-selected: "
        f"{sorted(d.loc[d.Preselected_by_Step04,'Target_Node_ID'].unique())} "
        f"(both Extreme_Increase)", True)

    d["Delta"] = d["FU_FBC"] - d["BL_FBC"]
    d["Pct_Change"] = 100 * d["Delta"] / d["BL_FBC"]
    d["Hemisphere"] = np.where(d["Target_Node"].str.startswith("ctx_rh"), "Right", "Left")
    d["Region"] = d["Target_Node"].str.replace("ctx_lh_G_and_S_cingul-", "", regex=False) \
                                  .str.replace("ctx_rh_G_and_S_cingul-", "", regex=False)
    d = d.sort_values(["Target_Node_ID", "ID"]).reset_index(drop=True)

    qc = pd.DataFrame(checks)
    qc.to_csv(OUT / "FC_QC_summary.csv", index=False)
    d.to_csv(OUT / "FC_analysis_dataset.csv", index=False)
    return d, qc


def analysis_a(d, rng):
    """A: is there any within-participant BL->FU change on each edge?"""
    rows = []
    for tid, g in d.groupby("Target_Node_ID"):
        x = g["Delta"].to_numpy(float)
        n = len(x)
        t = stats.ttest_1samp(x, 0)
        ci = t.confidence_interval()
        w = stats.wilcoxon(x)
        rows.append(dict(
            Target_Node_ID=tid, Target_Node=g["Target_Node"].iloc[0],
            Hemisphere=g["Hemisphere"].iloc[0], Region=g["Region"].iloc[0],
            Preselected_by_Step04=bool(g["Preselected_by_Step04"].iloc[0]),
            Step04_selection_direction=g["Step04_selection_direction"].iloc[0],
            N=n,
            BL_mean=g["BL_FBC"].mean(), BL_sd=g["BL_FBC"].std(ddof=1),
            FU_mean=g["FU_FBC"].mean(), FU_sd=g["FU_FBC"].std(ddof=1),
            Change_mean=x.mean(), Change_sd=x.std(ddof=1), Change_median=np.median(x),
            Change_SE=x.std(ddof=1) / np.sqrt(n),
            CI95_lower=float(ci.low), CI95_upper=float(ci.high),
            Median_pct_change=np.median(g["Pct_Change"]),
            N_increased=int((x > 0).sum()),
            Cohen_dz=x.mean() / x.std(ddof=1),
            Direction="Increase" if x.mean() > 0 else "Decrease",
            Shapiro_p=float(stats.shapiro(x).pvalue),
            Skew=float(stats.skew(x, bias=False)),
            p_ttest=float(t.pvalue),
            p_wilcoxon=float(w.pvalue),
            p_signflip=sign_flip_p(x, N_PERM, rng),
        ))
    a = pd.DataFrame(rows).sort_values("p_signflip").reset_index(drop=True)
    for c in ["p_ttest", "p_wilcoxon", "p_signflip"]:
        a[c.replace("p_", "q_FDR_")] = multipletests(a[c], method="fdr_bh")[1]
    a["p_Bonferroni_signflip"] = multipletests(a["p_signflip"], method="bonferroni")[1]
    # Selection-independent family: the 4 edges NOT pre-selected by Step 04.
    # These are the only ones on which question A is inferentially valid.
    clean = ~a["Preselected_by_Step04"]
    a["q_FDR_signflip_clean4"] = np.nan
    a.loc[clean, "q_FDR_signflip_clean4"] = multipletests(
        a.loc[clean, "p_signflip"], method="fdr_bh")[1]
    return a


def analysis_b(d, rng):
    """B: does the BL->FU change differ between responders and non-responders?"""
    rows = []
    for tid, g in d.groupby("Target_Node_ID"):
        g = g.sort_values("ID")
        lng = pd.melt(g, id_vars=["ID", "Responder_Status", "Age", "Sex"],
                      value_vars=["BL_FBC", "FU_FBC"],
                      var_name="Time", value_name="FBC")
        lng["Time"] = pd.Categorical(lng["Time"].map({"BL_FBC": "BL", "FU_FBC": "FU"}),
                                     categories=["BL", "FU"])
        lng["Group"] = pd.Categorical(lng["Responder_Status"],
                                      categories=["Non-responder", "Responder"])
        lng["Sex"] = pd.Categorical(lng["Sex"], categories=["F", "M"])
        with warnings.catch_warnings(record=True) as wl:
            warnings.simplefilter("always")
            res = smf.mixedlm("FBC ~ C(Time)*C(Group) + Age + C(Sex)",
                              lng, groups=lng["ID"]).fit(reml=True, method="lbfgs")
            wmsg = "; ".join(sorted({str(x.message).split(chr(10))[0][:100] for x in wl}))
        term = "C(Time)[T.FU]:C(Group)[T.Responder]"
        lo, hi = res.conf_int().loc[term].astype(float)

        delta = g["Delta"].to_numpy(float)
        is_r = (g["Responder_Status"].to_numpy() == "Responder")
        nr, nn = int(is_r.sum()), int((~is_r).sum())
        rows.append(dict(
            Target_Node_ID=tid, Target_Node=g["Target_Node"].iloc[0],
            Hemisphere=g["Hemisphere"].iloc[0], Region=g["Region"].iloc[0],
            N_responder=nr, N_nonresponder=nn,
            BL_mean_responder=g.loc[is_r, "BL_FBC"].mean(),
            BL_sd_responder=g.loc[is_r, "BL_FBC"].std(ddof=1),
            BL_mean_nonresponder=g.loc[~is_r, "BL_FBC"].mean(),
            BL_sd_nonresponder=g.loc[~is_r, "BL_FBC"].std(ddof=1),
            FU_mean_responder=g.loc[is_r, "FU_FBC"].mean(),
            FU_sd_responder=g.loc[is_r, "FU_FBC"].std(ddof=1),
            FU_mean_nonresponder=g.loc[~is_r, "FU_FBC"].mean(),
            FU_sd_nonresponder=g.loc[~is_r, "FU_FBC"].std(ddof=1),
            Change_mean_responder=delta[is_r].mean(),
            Change_mean_nonresponder=delta[~is_r].mean(),
            Interaction_estimate=float(res.params[term]),
            Interaction_SE=float(res.bse[term]),
            CI95_lower=float(lo), CI95_upper=float(hi),
            p_LME=float(res.pvalues[term]),
            p_permutation=label_perm_p(delta, is_r, N_PERM, rng),
            p_MannWhitney=float(stats.mannwhitneyu(delta[is_r], delta[~is_r]).pvalue),
            Cohen_d=(delta[is_r].mean() - delta[~is_r].mean()) / np.sqrt(
                ((nr - 1) * delta[is_r].var(ddof=1) + (nn - 1) * delta[~is_r].var(ddof=1))
                / (nr + nn - 2)),
            Direction=("Responders more positive change"
                       if float(res.params[term]) > 0 else "Responders more negative change"),
            p_OLS_changescore=float(sm.OLS(
                delta, sm.add_constant(is_r.astype(float))).fit().pvalues[1]),
            Converged=bool(res.converged), Model_warnings=wmsg,
            RandomIntercept_var=float(np.asarray(res.cov_re)[0, 0]),
            Residual_var=float(res.scale),
            ICC=float(np.asarray(res.cov_re)[0, 0])
            / (float(np.asarray(res.cov_re)[0, 0]) + float(res.scale)),
            Resid_shapiro_p=float(stats.shapiro(np.asarray(res.resid)).pvalue),
        ))
    b = pd.DataFrame(rows).sort_values("p_LME").reset_index(drop=True)
    for c in ["p_LME", "p_permutation", "p_MannWhitney"]:
        b[c.replace("p_", "q_FDR_")] = multipletests(b[c], method="fdr_bh")[1]
    return b


def analysis_c(d, rng):
    """C (exploratory): right vs left change, within homologous cingulate pairs."""
    piv = d.pivot(index="ID", columns="Target_Node_ID", values="Delta")
    rows = []
    for label, lnode, rnode in HOMOLOGUES:
        diff = (piv[rnode] - piv[lnode]).to_numpy(float)
        t = stats.ttest_1samp(diff, 0)
        ci = t.confidence_interval()
        pre_l = bool(d.loc[d.Target_Node_ID == lnode, "Preselected_by_Step04"].iloc[0])
        pre_r = bool(d.loc[d.Target_Node_ID == rnode, "Preselected_by_Step04"].iloc[0])
        rows.append(dict(
            Region=label, Left_node=lnode, Right_node=rnode,
            Left_preselected=pre_l, Right_preselected=pre_r,
            Contrast_confounded_by_selection=(pre_l != pre_r),
            Left_change_mean=piv[lnode].mean(), Right_change_mean=piv[rnode].mean(),
            RightMinusLeft_mean=diff.mean(),
            RightMinusLeft_SE=diff.std(ddof=1) / np.sqrt(len(diff)),
            CI95_lower=float(ci.low), CI95_upper=float(ci.high),
            Cohen_dz=diff.mean() / diff.std(ddof=1),
            N_right_greater=int((diff > 0).sum()), N=len(diff),
            p_ttest=float(t.pvalue),
            p_wilcoxon=float(stats.wilcoxon(diff).pvalue),
            p_signflip=sign_flip_p(diff, N_PERM, rng)))
    c = pd.DataFrame(rows)
    c["q_FDR_signflip"] = multipletests(c["p_signflip"], method="fdr_bh")[1]
    return c


def sensitivity(d, a, rng):
    """Influence and robustness for question A - the family with signal."""
    rows = []
    for tid, g in d.groupby("Target_Node_ID"):
        g = g.sort_values("ID")
        x = g["Delta"].to_numpy(float)
        ids = g["ID"].tolist()
        full_p = float(stats.ttest_1samp(x, 0).pvalue)

        # leave-one-participant-out
        lo_p, lo_est = [], []
        for i in range(len(x)):
            keep = np.ones(len(x), bool)
            keep[i] = False
            lo_p.append(float(stats.ttest_1samp(x[keep], 0).pvalue))
            lo_est.append(x[keep].mean())
        lo_p, lo_est = np.array(lo_p), np.array(lo_est)
        worst = int(np.argmax(lo_p))

        # exclude YTH032
        k32 = np.array([i != "YTH032" for i in ids])
        x32 = x[k32]

        # log-ratio (proportional change)
        lr = np.log(g["FU_FBC"].to_numpy()) - np.log(g["BL_FBC"].to_numpy())

        # trimmed mean / robust location
        rob = sm.RLM(x, np.ones((len(x), 1)), M=sm.robust.norms.HuberT()).fit()

        rows.append(dict(
            Target_Node_ID=tid, Target_Node=g["Target_Node"].iloc[0],
            Primary_mean=x.mean(), Primary_p=full_p,
            LOPO_max_p=float(lo_p.max()),
            LOPO_worst_participant=ids[worst],
            LOPO_min_p=float(lo_p.min()),
            LOPO_flips_at_05=bool((full_p < 0.05) != (lo_p.max() < 0.05)),
            Excl_YTH032_mean=x32.mean(),
            Excl_YTH032_p=float(stats.ttest_1samp(x32, 0).pvalue),
            Excl_YTH032_wilcoxon_p=float(stats.wilcoxon(x32).pvalue),
            LogRatio_mean=lr.mean(),
            LogRatio_p=float(stats.ttest_1samp(lr, 0).pvalue),
            LogRatio_wilcoxon_p=float(stats.wilcoxon(lr).pvalue),
            Robust_Huber_est=float(rob.params[0]), Robust_Huber_p=float(rob.pvalues[0]),
            Trimmed20_mean=float(stats.trim_mean(x, 0.2)),
        ))
    return pd.DataFrame(rows).sort_values("Primary_p").reset_index(drop=True)


def main():
    rng = np.random.default_rng(SEED)
    d, qc = build()
    a = analysis_a(d, rng)
    b = analysis_b(d, rng)
    c = analysis_c(d, rng)
    s = sensitivity(d, a, rng)

    a.to_csv(OUT / "FC_A_within_participant_change.csv", index=False)
    b.to_csv(OUT / "FC_B_group_difference_in_change.csv", index=False)
    c.to_csv(OUT / "FC_C_laterality_exploratory.csv", index=False)
    s.to_csv(OUT / "FC_sensitivity.csv", index=False)

    pd.set_option("display.width", 250)
    print(f"QC: {(qc.status=='PASS').sum()}/{len(qc)} PASS, {(qc.status=='FAIL').sum()} FAIL")
    print(f"Provenance max deviation vs raw matrices: "
          f"{qc.loc[qc.check_id==12,'observed'].iloc[0]}")

    print("\n=== A  Within-participant BL->FU change (6 tests) ===")
    print(a[["Target_Node_ID", "Target_Node", "Preselected_by_Step04", "Change_mean",
             "CI95_lower", "CI95_upper", "Cohen_dz", "N_increased", "p_ttest",
             "p_wilcoxon", "p_signflip", "q_FDR_signflip",
             "q_FDR_signflip_clean4"]].round(5).to_string(index=False))
    print("  NOTE: nodes 95 and 96 were pre-selected by the Step-04 median +/- 3 x IQR")
    print("  screen AS Extreme_Increase, i.e. selected on the very statistic question A")
    print("  tests. Their p-values are circular and are NOT valid evidence of change.")
    print("  q_FDR_signflip_clean4 corrects only across the 4 selection-independent edges.")

    print("\n=== B  Responder vs non-responder difference in change (6 tests) ===")
    print(b[["Target_Node_ID", "Target_Node", "Interaction_estimate", "CI95_lower",
             "CI95_upper", "Cohen_d", "p_LME", "p_permutation", "q_FDR_LME"]]
          .round(5).to_string(index=False))

    print("\n=== C  Right vs left homologous change (exploratory, 3 tests) ===")
    print(c[["Region", "Left_change_mean", "Right_change_mean", "RightMinusLeft_mean",
             "CI95_lower", "CI95_upper", "Cohen_dz", "N_right_greater",
             "p_signflip", "q_FDR_signflip",
             "Contrast_confounded_by_selection"]].round(5).to_string(index=False))

    print("\n=== Sensitivity for question A ===")
    print(s[["Target_Node_ID", "Target_Node", "Primary_p", "LOPO_max_p",
             "LOPO_worst_participant", "Excl_YTH032_p", "Excl_YTH032_wilcoxon_p",
             "LogRatio_p", "Robust_Huber_p"]].round(4).to_string(index=False))

    print("\nB: boundary-warning edges vs the equivalent change-score OLS "
          "(no boundary issue):")
    print(b[["Target_Node_ID", "Target_Node", "p_LME", "p_OLS_changescore",
             "Model_warnings"]].round(4).to_string(index=False))
    print("\nB diagnostics: converged "
          f"{int(b.Converged.sum())}/6, warnings "
          f"{int((b.Model_warnings.fillna('').str.len()>0).sum())}, "
          f"ICC {b.ICC.min():.3f}-{b.ICC.max():.3f}")


if __name__ == "__main__":
    main()
