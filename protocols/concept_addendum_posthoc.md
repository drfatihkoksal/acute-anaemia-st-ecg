# Addendum P: post hoc analyses anticipating peer review

**Status:** locked 2026-10-26, before any of the analyses below was run. The SHA-256 hashes of this file and of `heedb/12_posthoc_reviewer.py` are recorded in `concept_addendum_posthoc.sha256` before the script is run. Keeping them in a separate file means recording them does not change this file's hash.

**Context.** Protocols v3.1 (concept.md §11) and the tolerance protocol (concept_tolerance.md v1.0) were locked and analysed on 2026-10-02. All of their results are known. The analyses in this addendum were therefore defined **after** those results were seen. Their purpose is to answer questions a reviewer is expected to ask (internal review, 2026-10-26).

**Reporting rules.**
- Every analysis here is labelled *post hoc* in the manuscript and repository, whatever its result.
- No result is omitted.
- None of these analyses replaces or redefines a primary or secondary outcome of the locked protocols.
- Only aggregate results are written to disk. Any count below 10 is suppressed.

## 1. Data and unchanged elements

- **Inputs:** `heedb/out/analysis_{mgh,mimic}.parquet` (Analysis 1) and `heedb/out/tolerance_{mgh,mimic}.parquet` (Analysis 2). For MGH, J-point amplitudes of all leads are joined from `heedb/out/measure_mgh/part_*.parquet` by file name. MIMIC already has them.
- **Reused as is (imported, not modified):**
  - Analysis 1 pair selection and model: `06_analysis.prep`, `fit` and `tost`; adjustment set `ADJ` plus baseline ST60; HC3; MI m=5; margin ±0.01 mV per g/dL.
  - Analysis 2 cohort, covariates, imputation and model: `11_tolerance_analysis.imputations` and `primary`; full adjustment set incl. heart rate; MI m=10.
- **Pooling:** fixed-effect inverse-variance on the log scale, with Cochran's Q p for heterogeneity.

## 2. Analysis 1 (Hb decline and ST60)

**P1. Adults only.** The primary model and equivalence test, restricted to primary pairs with age ≥18 years.

**P2. Categories of Hb change.** In primary pairs, ΔHb (earlier minus later) is grouped as:
- rise ≤−0.5
- stable (−0.5, 0.5) — reference
- 0.5 to <2
- 2 to <3
- ≥3 g/dL

OLS of ΔST60(min) on category indicators plus `ADJ` and baseline ST60; complete cases; HC3. Reported: adjusted mean difference versus stable, with 90% CI, for each category. For the ≥3 g/dL category, equivalence (two one-sided tests) against ±0.03 mV, i.e. the prespecified ±0.01 mV per g/dL times 3 g/dL.

**P3. Very low nadir.** The primary slope model in primary pairs with later Hb ≤7 g/dL. The ≤8 g/dL and coded-bleeding subgroups are already in the frozen results and are reported from there, not re-run.

**P4. Frequency of new ST depression by Hb change.** In adult primary pairs, the share with new ST depression in each category of P2. New ST depression has the Analysis 2 definition: ≥2 of the 11 non-aVR leads with later ST60 ≤−0.05 mV and a fall of ≥0.05 mV. Reported as crude percentages with Wilson 95% CIs; counts below 10 are suppressed.

**P5. Measurement noise.** The SD of ΔST60(min) in pairs with |ΔHb| <0.5 g/dL. This is reported from the frozen results (`noise_floor_sd_stable_hb`), not recomputed.

## 3. Analysis 2 (new ST depression and 30-day mortality)

**P6. J-point definition.** New ST depression at the J point: ≥2 non-aVR leads with later J amplitude ≤−0.05 mV and a fall of ≥0.05 mV from the earlier ECG. Same cohort, model and imputation as the primary analysis; OR per centre and pooled.

**P7. Morphology.** Among patients with new ST depression (primary definition), a qualifying lead is classified as *horizontal or downsloping* when its later ST60 is no more than 0.01 mV above its later J amplitude; otherwise it is *upsloping*.
- The patient is *horizontal/downsloping* if ≥2 qualifying leads are horizontal or downsloping, and *upsloping-only* otherwise.
- Model: the primary model with two indicators (horizontal/downsloping; upsloping-only) against no new ST depression. OR for each, per centre and pooled.

**P8. Without heart-rate adjustment, pooled.** The per-centre estimates of the locked sensitivity analysis `no_hr_adjust` are pooled in the same way as the primary estimate. This is presentation only; no new model is fitted.

**P9. E-values.** For the pooled primary OR and the pooled OR without heart-rate adjustment (P8):
- Each OR is converted to an approximate risk ratio, RR ≈ √OR, because 30-day mortality is above 15%.
- E-value = RR + √(RR(RR−1)), for the point estimate and for the confidence limit closer to 1.

## 4. Output

`heedb/out/posthoc_reviewer.json` (aggregates only), plus the code SHA-1 inside the file.

## 5. Log

| Date | Event |
| --- | --- |
| 2026-10-26 | Addendum written and locked; hashes in `concept_addendum_posthoc.sha256`. Any later change to this file or to the script is logged in that file with a reason. |
