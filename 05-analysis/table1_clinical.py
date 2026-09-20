"""Table 1 - demographics and clinical outcomes by treatment response.

Journal-style table rendered in the layout of aim-tabel.png.
Sources: clinical-output.xlsx (Sheet1, demographics) + clinical_response.csv.
Outputs: table1_clinical.{png,pdf,csv,md}
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats

XLSX, CSV = "clinical-output.xlsx", "clinical_response.csv"
RESPONDER_THRESHOLD = 0.50

SERIF = ["Times New Roman", "DejaVu Serif"]
plt.rcParams.update({
    "font.family": "serif", "font.serif": SERIF,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

# ------------------------------------------------------------ formatting --
def ms(v):
    return f"{np.mean(v):.2f} ± {np.std(v, ddof=1):.2f}"


def ms_pct(v):
    return f"{np.mean(v) * 100:.2f}% ± {np.std(v, ddof=1) * 100:.2f}%"


def n_pct(k, n):
    return f"{k} ({100 * k / n:.2f}%)"


def fmt_p(p):
    return "<0.001" if p < 0.001 else f"{p:.4f}"


def sup(letter, text):
    return rf"$^{{\mathrm{{{letter}}}}}$ {text}"


# ------------------------------------------------------------------ data --
def build():
    dx = pd.read_excel(XLSX, "Sheet1").dropna(subset=["ID"])
    dc = pd.read_csv(CSV)
    assert set(dx.ID) == set(dc.ID), "ID mismatch between workbook and CSV"
    m = dx.merge(dc[["ID", "DASS_Reduction", "Responder_Status"]], on="ID")
    assert np.allclose(m.BL_MRIdepDASS_x if "BL_MRIdepDASS_x" in m
                       else m.BL_MRIdepDASS,
                       dc.set_index("ID").loc[m.ID, "BL_MRIdepDASS"].values)

    r = m[m.Responder_Status == "Responder"]
    nr = m[m.Responder_Status == "Non-responder"]
    nR, nN, nT = len(r), len(nr), len(m)

    # --- tests ---------------------------------------------------------
    p_age = stats.ttest_ind(r.Age, nr.Age, equal_var=False).pvalue
    sex_ct = [[(r.Sex == "M").sum(), (r.Sex == "F").sum()],
              [(nr.Sex == "M").sum(), (nr.Sex == "F").sum()]]
    p_sex = stats.fisher_exact(sex_ct)[1]
    p_pre = stats.ttest_ind(r.BL_MRIdepDASS, nr.BL_MRIdepDASS,
                            equal_var=False).pvalue
    p_post = stats.ttest_ind(r.FU_MRIdepDASS, nr.FU_MRIdepDASS,
                             equal_var=False).pvalue
    dR, dN = r.BL_MRIdepDASS - r.FU_MRIdepDASS, nr.BL_MRIdepDASS - nr.FU_MRIdepDASS
    p_chg = stats.ttest_ind(dR, dN, equal_var=False).pvalue
    p_red = stats.ttest_ind(r.DASS_Reduction, nr.DASS_Reduction,
                            equal_var=False).pvalue

    # within-group paired change, for footnote d
    paired = {}
    for lab, g in [("Responder", r), ("Non-responder", nr), ("Total", m)]:
        dd = g.BL_MRIdepDASS - g.FU_MRIdepDASS
        t, p = stats.ttest_rel(g.BL_MRIdepDASS, g.FU_MRIdepDASS)
        paired[lab] = (dd.mean(), dd.std(ddof=1), len(g) - 1, t, p,
                       dd.mean() / dd.std(ddof=1))

    # --- rows ----------------------------------------------------------
    rows = [
        ("row", sup("a", "Age (years)"), ms(r.Age), ms(nr.Age), fmt_p(p_age)),
        ("row", sup("b", "Sex (male/female)"),
         f"{sex_ct[0][0]}/{sex_ct[0][1]}", f"{sex_ct[1][0]}/{sex_ct[1][1]}",
         fmt_p(p_sex)),
        ("row", sup("c", "DASS-21 depression: pre-rTMS"),
         ms(r.BL_MRIdepDASS), ms(nr.BL_MRIdepDASS), fmt_p(p_pre)),
        ("row", "Post-rTMS", ms(r.FU_MRIdepDASS), ms(nr.FU_MRIdepDASS),
         fmt_p(p_post)),
        ("row", sup("d", "Change (pre − post)"), ms(dR), ms(dN), fmt_p(p_chg)),
        ("row", sup("d", "Reductive rate"), ms_pct(r.DASS_Reduction),
         ms_pct(nr.DASS_Reduction), fmt_p(p_red)),
    ]

    pr, pn, pt = paired["Responder"], paired["Non-responder"], paired["Total"]
    notes = [
        ("n", "Note: Bold values indicate a highly significant result at p < 0.001. "
              "Groups are defined by a "
              f"{RESPONDER_THRESHOLD:.0%} reduction in DASS-21 depression score."),
        ("n", f"Whole cohort (N = {nT}): age {ms(m.Age)} years; sex "
              f"{(m.Sex == 'M').sum()}/{(m.Sex == 'F').sum()} male/female; "
              f"pre-rTMS {ms(m.BL_MRIdepDASS)}; post-rTMS {ms(m.FU_MRIdepDASS)}; "
              f"reductive rate {ms_pct(m.DASS_Reduction)}."),
        ("n", "Abbreviations: DASS-21, Depression Anxiety Stress Scales-21 "
              "(depression subscale, doubled scoring, range 0-42); rTMS, "
              "repetitive transcranial magnetic stimulation; N/A, not applicable."),
        ("f", sup("a", "p-Value was obtained from a Welch two-sample t-test;")),
        ("f", sup("b", "p-Value was obtained from Fisher's exact test "
                       "(expected cell counts < 5 preclude a chi-square test);")),
        ("f", sup("c", "p-Values in this block compare the two groups. The "
                       "within-group change was significant in both: responders "
                       f"{pr[0]:.2f} ± {pr[1]:.2f}, paired t({pr[2]}) = {pr[3]:.2f}, "
                       f"p < 0.001, d$_z$ = {pr[5]:.2f}; non-responders "
                       f"{pn[0]:.2f} ± {pn[1]:.2f}, paired t({pn[2]}) = {pn[3]:.2f}, "
                       f"p = {pn[4]:.3f}, d$_z$ = {pn[5]:.2f}; whole cohort "
                       f"{pt[0]:.2f} ± {pt[1]:.2f}, paired t({pt[2]}) = {pt[3]:.2f}, "
                       f"p < 0.001, d$_z$ = {pt[5]:.2f};")),
        ("f", sup("d", "Circular - group membership is derived from the change "
                       "score, so these contrasts are significant by construction "
                       "and are descriptive only, not findings.")),
    ]
    return rows, notes, nR, nN


# ---------------------------------------------------------------- render --
def render(rows, notes, nR, nN):
    COLS = [0.000, 0.500, 0.735, 1.000]   # char | resp(ctr) | nonresp(ctr) | p(right)
    H_ROW, H_SUB, H_NOTE, FS = 1.00, 0.92, 0.78, 10.0

    total = 2.6 + 1.5 + sum(H_ROW if k != "sub" else H_SUB for k, *_ in rows)
    total += 0.8 + sum(H_NOTE for _ in notes)
    fig_h = total * 0.205
    fig, ax = plt.subplots(figsize=(9.6, fig_h))
    ax.set_xlim(0, 1)
    ax.set_ylim(total, 0)
    ax.axis("off")

    def rule(y, lw):
        ax.add_line(Line2D([0, 1], [y, y], color="black", lw=lw,
                           solid_capstyle="butt"))

    y = 0.9
    ax.text(0, y, "TABLE 1", fontsize=FS - 0.5, fontweight="bold", va="center")
    ax.text(0.063, y, "|", fontsize=FS - 0.5, va="center", color="#555555")
    ax.text(0.080, y, "Demographics and clinical outcomes by treatment response.",
            fontsize=FS - 0.5, va="center")

    y = 2.0
    rule(y, 1.4)
    y += 1.0
    hdr = ["Characteristics", f"Responder ($n$ = {nR})",
           f"Non-responder ($n$ = {nN})", "$p$-Value"]
    ax.text(COLS[0], y, hdr[0], fontsize=FS, fontweight="bold", va="center")
    ax.text(COLS[1], y, hdr[1], fontsize=FS, fontweight="bold", va="center",
            ha="center")
    ax.text(COLS[2], y, hdr[2], fontsize=FS, fontweight="bold", va="center",
            ha="center")
    ax.text(COLS[3], y, hdr[3], fontsize=FS, fontweight="bold", va="center",
            ha="right")
    y += 0.55
    rule(y, 0.8)

    for kind, label, c1, c2, p in rows:
        y += H_SUB if kind == "sub" else H_ROW
        indent = 0.022 if kind == "sub" else 0.0
        ax.text(COLS[0] + indent, y, label, fontsize=FS, va="center")
        for x, val in [(COLS[1], c1), (COLS[2], c2)]:
            if val:
                ax.text(x, y, val, fontsize=FS, va="center", ha="center")
        if p:
            ax.text(COLS[3], y, p, fontsize=FS, va="center", ha="right",
                    fontweight="bold" if p == "<0.001" else "normal")
    y += 0.55
    rule(y, 1.4)

    y += 0.45
    for kind, text in notes:
        y += H_NOTE
        if kind == "i":
            continue
        ax.text(0, y, text, fontsize=FS - 1.8, va="center",
                style="italic" if kind == "n" else "normal")

    fig.subplots_adjust(left=0.012, right=0.988, top=0.995, bottom=0.005)
    fig.savefig("table1_clinical.png", dpi=400, facecolor="white",
                bbox_inches="tight", pad_inches=0.06)
    fig.savefig("table1_clinical.pdf", facecolor="white",
                bbox_inches="tight", pad_inches=0.06)


def export_text(rows, notes, nR, nN):
    """Plain-text CSV + markdown versions, mathtext stripped."""
    import re

    def clean(s):
        s = re.sub(r"\$\^\{\\mathrm\{(\w)\}\}\$\s*", r"[\1] ", s)
        return s.replace("$_z$", "z").replace("$", "")

    hdr = ["Characteristics", f"Responder (n={nR})",
           f"Non-responder (n={nN})", "p-Value"]
    body = [[clean(l), c1, c2, p] for _, l, c1, c2, p in rows]
    pd.DataFrame(body, columns=hdr).to_csv("table1_clinical.csv", index=False)

    with open("table1_clinical.md", "w") as f:
        f.write("**TABLE 1** | Demographics and clinical outcomes by "
                "treatment response.\n\n")
        f.write("| " + " | ".join(hdr) + " |\n")
        f.write("|---|---:|---:|---:|\n")
        for kind, l, c1, c2, p in rows:
            lab = ("&nbsp;&nbsp;&nbsp;&nbsp;" if kind == "sub" else "") + clean(l)
            pv = f"**{p}**" if p == "<0.001" else p
            f.write(f"| {lab} | {c1} | {c2} | {pv} |\n")
        f.write("\n")
        for kind, t in notes:
            if kind != "i":
                f.write(f"{clean(t)}\n\n")


if __name__ == "__main__":
    rows, notes, nR, nN = build()
    render(rows, notes, nR, nN)
    export_text(rows, notes, nR, nN)
    print(f"Table 1 written: {nR} responders, {nN} non-responders")
    for _, l, c1, c2, p in rows:
        import re
        print(f"  {re.sub(r'[$^{}]|\\\\mathrm', '', l)[:34]:<36} {c1:>18} "
              f"{c2:>18} {p:>8}")
