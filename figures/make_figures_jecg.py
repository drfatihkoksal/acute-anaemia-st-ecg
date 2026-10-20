"""Journal of Electrocardiology version: Figures 1–2 (vector PDF, plus PNG previews) and the supplementary
data workbook, from the frozen result files in results/. Outputs are written to figures/jecg/.

Fig. 2 plots the prespecified lead-level secondary outcomes already in results_*.json; no new analysis.
The equivalence margin and the follow-up window are read from the frozen analysis code, not typed here.
"""
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIG = HERE / "jecg"
FIG.mkdir(exist_ok=True)
RES = {p.stem: json.loads(p.read_text()) for p in (ROOT / "results").glob("*.json")}


def get(path):
    """'results_mgh.secondary.d_st60_V1.ci90.0' -> value in results/results_mgh.json"""
    stem, *keys = path.split(".")
    v = RES[stem]
    for k in keys:
        v = v[int(k)] if isinstance(v, list) else v[k]
    return v


def frozen_const(script, pattern):
    m = re.search(pattern, (ROOT / "analysis" / script).read_text(), re.M)
    return float(m.group(1))


LEADS = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]
CENTRES = [("mgh", "MGH"), ("mimic", "BIDMC")]
UV = 1000
MARGIN = frozen_const("06_analysis.py", r"^MARGIN = ([\d.]+)") * UV
FOLLOWUP = int(frozen_const("10_build_tolerance_mgh.py", r"^T\[\"fu30_complete\"\] = \(T\.post_time \+ pd\.Timedelta\(days=(\d+)\)"))

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42, "ps.fonttype": 42})
BLUE, ORANGE, GREY, INK2 = "#2a78d6", "#eb6834", "#8c8c8c", "#444444"
COL = {"mgh": BLUE, "mimic": ORANGE}
MARK = {"mgh": "o", "mimic": "s"}


def mdash(x, fmt):
    return format(x, fmt).replace("-", "−")


def rng(lo, hi, fmt):
    sep = " to " if (lo < 0 or hi < 0) else "–"
    return f"{mdash(lo, fmt)}{sep}{mdash(hi, fmt)}"


def strip(ax):
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)


# ------------------------------------------------------------------ Figure 1
panel_a = pd.DataFrame([{
    "Centre": lab, "Patients": get(f"results_{c}.H1_st60_min.n"),
    "Estimate_uV_per_gdL": get(f"results_{c}.H1_st60_min.beta") * UV,
    "CI90_low": get(f"results_{c}.H1_st60_min.ci90.0") * UV,
    "CI90_high": get(f"results_{c}.H1_st60_min.ci90.1") * UV} for c, lab in CENTRES])


def orow(label, kind, key, est="OR"):
    return {"Row": label, "Group": kind, "Odds_ratio": get(f"{key}.{est}"),
            "CI95_low": get(f"{key}.ci95.0"), "CI95_high": get(f"{key}.ci95.1")}


panel_b = pd.DataFrame([
    orow("Acute anaemia, MGH", "centre", "tolerance_results_mgh.primary_OR_new_std"),
    orow("Acute anaemia, BIDMC", "centre", "tolerance_results_mimic.primary_OR_new_std"),
    orow("Acute anaemia, pooled", "pooled", "tolerance_pooled.primary"),
    orow("Stable normal Hb, MGH", "contrast", "tolerance_results_mgh.specificity.OR_contrast"),
    orow("Stable normal Hb, BIDMC", "contrast", "tolerance_results_mimic.specificity.OR_contrast"),
])

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(174 / 25.4, 120 / 25.4), gridspec_kw={"height_ratios": [2, 5]})
fig.subplots_adjust(left=0.49, right=0.97, top=0.93, bottom=0.10, hspace=0.75)
lim = 1.5 * MARGIN
ax1.axvspan(-MARGIN, MARGIN, ymax=0.83, color="#e6e6e6", zorder=0)
ax1.axvline(0, ymax=0.83, color="black", lw=0.6)
for i, r in panel_a.iterrows():
    y = len(panel_a) - 1 - i
    c = CENTRES[i][0]
    ax1.plot([r.CI90_low, r.CI90_high], [y, y], color=COL[c], lw=1.5)
    ax1.plot(r.Estimate_uV_per_gdL, y, MARK[c], color=COL[c], ms=4.5)
    ax1.text(-0.94, y, f"{r.Centre} (n={r.Patients:,})", transform=ax1.get_yaxis_transform(), va="center")
    ax1.text(-0.40, y, f"{mdash(r.Estimate_uV_per_gdL, '.2f')} ({rng(r.CI90_low, r.CI90_high, '.2f')})",
             transform=ax1.get_yaxis_transform(), va="center", color=INK2)
ax1.set_xlim(-lim, lim); ax1.set_ylim(-0.7, len(panel_a) + 0.05); ax1.set_yticks([])
ax1.set_xlabel("µV per g/dL (positive = ST higher after the fall)")
strip(ax1)
ax1.set_title("A   Change in most negative ST60 per g/dL Hb fall (90% CI)", loc="left", x=-0.94, fontweight="bold")
ax1.text(0, len(panel_a) - 0.28, f"grey band = equivalence margin ±{MARGIN:.0f} µV per g/dL",
         fontsize=7, color="#555555", ha="center")

ypos, y = [], 0
for i, r in panel_b.iterrows():
    if r.Group == "contrast" and (i == 0 or panel_b.Group[i - 1] != "contrast"):
        y += 0.6
    ypos.append(y); y += 1
ypos = [max(ypos) - v for v in ypos]
ax2.axvline(1, color=GREY, lw=0.8, ls="--")
for (i, r), yy in zip(panel_b.iterrows(), ypos):
    col = "black" if r.Group == "pooled" else (BLUE if "MGH" in r.Row else ORANGE)
    ax2.plot([r.CI95_low, r.CI95_high], [yy, yy], color=col, lw=1.5)
    ax2.plot(r.Odds_ratio, yy, "D" if r.Group == "pooled" else ("o" if "MGH" in r.Row else "s"), color=col,
             ms=6 if r.Group == "pooled" else 4.5)
    ax2.text(-0.94, yy, r.Row, transform=ax2.get_yaxis_transform(), va="center",
             fontweight="bold" if r.Group == "pooled" else "normal")
    ax2.text(-0.40, yy, f"{r.Odds_ratio:.2f} ({rng(r.CI95_low, r.CI95_high, '.2f')})",
             transform=ax2.get_yaxis_transform(), va="center", color=INK2)
ax2.set_xscale("log"); ax2.set_xlim(0.5, 2.5)
ticks = [0.5, 0.75, 1, 1.5, 2, 2.5]
ax2.set_xticks(ticks); ax2.set_xticklabels([str(t) for t in ticks]); ax2.minorticks_off()
ax2.set_yticks([]); ax2.set_ylim(-0.6, max(ypos) + 0.6)
ax2.set_xlabel(f"Adjusted odds ratio for death within {FOLLOWUP} days (log scale)")
strip(ax2)
ax2.set_title(f"B   New vs no new ST depression and death within {FOLLOWUP} days (95% CI)",
              loc="left", x=-0.94, fontweight="bold")
fig.savefig(FIG / "Figure1.pdf")
fig.savefig(FIG / "Figure1_preview.png", dpi=150)
plt.close(fig)

# ------------------------------------------------------------------ Figure 2
rows = []
for c, lab in CENTRES:
    for L in LEADS:
        k = f"results_{c}.secondary.d_st60_{L}"
        rows.append({"Centre": lab, "Lead": L, "Patients": get(f"{k}.n"), "Estimate_uV_per_gdL": get(f"{k}.beta") * UV,
                     "CI90_low": get(f"{k}.ci90.0") * UV, "CI90_high": get(f"{k}.ci90.1") * UV})
fig2 = pd.DataFrame(rows)

fig, ax = plt.subplots(figsize=(84 / 25.4, 110 / 25.4))
fig.subplots_adjust(left=0.16, right=0.97, top=0.90, bottom=0.12)
ax.axvspan(-MARGIN, MARGIN, color="#e6e6e6", zorder=0)
ax.axvline(0, color="black", lw=0.6)
off = {"mgh": 0.17, "mimic": -0.17}
for c, lab in CENTRES:
    d = fig2[fig2.Centre == lab]
    for j, r in enumerate(d.itertuples()):
        y = len(LEADS) - 1 - j + off[c]
        ax.plot([r.CI90_low, r.CI90_high], [y, y], color=COL[c], lw=1.2)
        ax.plot(r.Estimate_uV_per_gdL, y, MARK[c], color=COL[c], ms=4, mec="white", mew=0.6,
                label=f"{lab} (n={r.Patients:,})" if j == 0 else None)
ax.set_yticks(range(len(LEADS))); ax.set_yticklabels(LEADS[::-1])
ax.tick_params(axis="y", length=0)
ax.set_xlim(-lim, lim); ax.set_ylim(-0.7, len(LEADS) - 0.3)
ax.set_xlabel("Change in ST60, µV per g/dL Hb fall (90% CI)")
strip(ax)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=7, handletextpad=0.3)
fig.savefig(FIG / "Figure2.pdf")
fig.savefig(FIG / "Figure2_preview.png", dpi=200)
plt.close(fig)

# ------------------------------------------------------------------ supplementary data
with pd.ExcelWriter(FIG / "Supplementary_data.xlsx") as xw:
    pd.DataFrame({"Item": ["Fig. 1A", "Fig. 1B", "Fig. 2", "Source"], "Description": [
        "Adjusted within-patient change in the most negative ST60 per g/dL haemoglobin fall (µV), 90% CI; one ECG pair per patient.",
        f"Adjusted odds ratios (95% CI) for death within {FOLLOWUP} days, new vs no new ST depression; acute anaemia cohorts, pooled fixed-effect estimate, and stable-normal-Hb contrast cohorts.",
        "Lead-level adjusted change in ST60 per g/dL haemoglobin fall (µV), 90% CI; complete cases; prespecified secondary outcomes.",
        "Aggregate estimates from the prespecified, frozen analysis code (06_analysis.py sha1 "
        f"{get('results_mgh.code_sha1')[:8]}; 11_tolerance_analysis.py sha1 {get('tolerance_results_mgh.code_sha1')[:8]}). "
        "No individual-level data are included, as required by the data-use agreements."]}
                 ).to_excel(xw, sheet_name="README", index=False)
    panel_a.to_excel(xw, sheet_name="Fig1A", index=False)
    panel_b.to_excel(xw, sheet_name="Fig1B", index=False)
    fig2.to_excel(xw, sheet_name="Fig2", index=False)
print("figures and supplementary data written to", FIG)
