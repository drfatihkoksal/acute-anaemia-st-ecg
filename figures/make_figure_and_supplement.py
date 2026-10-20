"""Figure1.tif (600 dpi, LZW) and the Online Data Supplement (xlsx) from the frozen result files."""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent          # outputs written next to this script
OUT = HERE.parent / "results"                    # aggregate result files of this repository
A_mgh = json.load(open(OUT / "results_mgh.json"))["H1_st60_min"]
A_mim = json.load(open(OUT / "results_mimic.json"))["H1_st60_min"]
T_mgh = json.load(open(OUT / "tolerance_results_mgh.json"))
T_mim = json.load(open(OUT / "tolerance_results_mimic.json"))
pooled = json.load(open(OUT / "tolerance_pooled.json"))["primary"]

panel_a = pd.DataFrame([
    {"Centre": "MGH", "Patients": 154187, "Estimate_uV_per_gdL": A_mgh["beta"] * 1000,
     "CI90_low": A_mgh["ci90"][0] * 1000, "CI90_high": A_mgh["ci90"][1] * 1000},
    {"Centre": "BIDMC", "Patients": 24741, "Estimate_uV_per_gdL": A_mim["beta"] * 1000,
     "CI90_low": A_mim["ci90"][0] * 1000, "CI90_high": A_mim["ci90"][1] * 1000},
])


def orow(label, kind, r):
    return {"Row": label, "Group": kind, "Odds_ratio": r["OR"], "CI95_low": r["ci95"][0], "CI95_high": r["ci95"][1]}


panel_b = pd.DataFrame([
    orow("Acute anaemia, MGH", "centre", T_mgh["primary_OR_new_std"]),
    orow("Acute anaemia, BIDMC", "centre", T_mim["primary_OR_new_std"]),
    orow("Acute anaemia, pooled", "pooled", pooled),
    orow("Stable normal Hb, MGH", "contrast", T_mgh["specificity"]["OR_contrast"]),
    orow("Stable normal Hb, BIDMC", "contrast", T_mim["specificity"]["OR_contrast"]),
])

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
BLUE, GREY = "#1f5fbf", "#8c8c8c"
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(174 / 25.4, 120 / 25.4), gridspec_kw={"height_ratios": [2, 5]})
fig.subplots_adjust(left=0.49, right=0.97, top=0.93, bottom=0.10, hspace=0.75)

# Panel A
ax1.axvspan(-10, 10, ymax=0.83, color="#e6e6e6", zorder=0)
ax1.axvline(0, ymax=0.83, color="black", lw=0.6)
for v in (-10, 10):
    ax1.axvline(v, ymax=0.83, color=GREY, lw=0.6, ls="--")
for i, r in panel_a.iterrows():
    y = len(panel_a) - 1 - i
    ax1.plot([r.CI90_low, r.CI90_high], [y, y], color=BLUE, lw=1.5)
    ax1.plot(r.Estimate_uV_per_gdL, y, "o", color=BLUE, ms=4)
    ax1.text(-0.94, y, f"{r.Centre} (n={r.Patients:,})", transform=ax1.get_yaxis_transform(), va="center", ha="left")
    ax1.text(-0.40, y, f"{r.Estimate_uV_per_gdL:.2f} ({r.CI90_low:.2f} to {r.CI90_high:.2f})",
             transform=ax1.get_yaxis_transform(), va="center", ha="left", color="#444444")
ax1.set_xlim(-15, 15); ax1.set_ylim(-0.7, len(panel_a) + 0.05); ax1.set_yticks([])
ax1.set_xlabel("µV per g/dL (positive = ST higher after the fall)")
for s in ("top", "right", "left"):
    ax1.spines[s].set_visible(False)
ax1.set_title("A   Change in most negative ST60 per 1 g/dL Hb fall (90% CI)", loc="left", x=-0.94, fontweight="bold")
ax1.text(0, len(panel_a) - 0.28, "grey band = equivalence margin ±10 µV per g/dL", fontsize=7, color="#555555", ha="center")

# Panel B
ypos, y = [], 0
for i, r in panel_b.iterrows():
    if r.Group == "contrast" and (i == 0 or panel_b.Group[i - 1] != "contrast"):
        y += 0.6
    ypos.append(y); y += 1
ypos = [max(ypos) - v for v in ypos]
ax2.axvline(1, color=GREY, lw=0.8, ls="--")
for (i, r), yy in zip(panel_b.iterrows(), ypos):
    col = BLUE if r.Group == "pooled" else GREY
    ax2.plot([r.CI95_low, r.CI95_high], [yy, yy], color=col, lw=1.5)
    ax2.plot(r.Odds_ratio, yy, "D" if r.Group == "pooled" else "s", color=col, ms=6 if r.Group == "pooled" else 4.5)
    ax2.text(-0.94, yy, r.Row, transform=ax2.get_yaxis_transform(), va="center", ha="left",
             fontweight="bold" if r.Group == "pooled" else "normal")
    ax2.text(-0.40, yy, f"{r.Odds_ratio:.2f} ({r.CI95_low:.2f}–{r.CI95_high:.2f})",
             transform=ax2.get_yaxis_transform(), va="center", ha="left", color="#444444")
ax2.set_xscale("log"); ax2.set_xlim(0.5, 2.5)
ticks = [0.5, 0.75, 1, 1.5, 2, 2.5]
ax2.set_xticks(ticks); ax2.set_xticklabels([str(t) for t in ticks]); ax2.minorticks_off()
ax2.set_yticks([]); ax2.set_ylim(-0.6, max(ypos) + 0.6)
ax2.set_xlabel("Adjusted odds ratio for death within 30 days (log scale)")
for s in ("top", "right", "left"):
    ax2.spines[s].set_visible(False)
ax2.set_title("B   New vs no new ST depression and death within 30 days (95% CI)", loc="left", x=-0.94, fontweight="bold")

fig.savefig(HERE / "Figure1.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(HERE / "Figure1_preview.png", dpi=150)

with pd.ExcelWriter(HERE / "Online_Data_Supplement.xlsx") as xw:
    pd.DataFrame({"Item": ["Figure 1A", "Figure 1B", "Source"],
                  "Description": [
                      "Adjusted within-patient change in the most negative ST60 per 1 g/dL haemoglobin fall (µV), 90% CI; one ECG pair per patient.",
                      "Adjusted odds ratios (95% CI) for death within 30 days, new vs no new ST depression; acute anaemia cohorts, pooled fixed-effect estimate, and stable-normal-Hb contrast cohorts.",
                      "Aggregate estimates from the protocol analysis code (06_analysis.py sha1 788c5b4d; 11_tolerance_analysis.py sha1 5954a4ad). No individual-level data are included, as required by the data-use agreements."]}
                 ).to_excel(xw, sheet_name="README", index=False)
    panel_a.to_excel(xw, sheet_name="Figure1A", index=False)
    panel_b.to_excel(xw, sheet_name="Figure1B", index=False)
print(panel_a.round(3).to_string()); print(panel_b.round(3).to_string())
