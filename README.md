# Acute haemoglobin decline and the ST segment — analysis code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23110835.svg)](https://doi.org/10.5281/zenodo.23110835)

Code for the study *"Acute haemoglobin decline and the ST segment: within-patient analysis of serial ECGs in two hospitals"* (Koksal F; manuscript under review).

The study asks two questions, each analysed under a written protocol, in serial 12-lead ECGs from two hospitals: Massachusetts General Hospital (HEEDB, discovery) and Beth Israel Deaconess Medical Center (MIMIC-IV / MIMIC-IV-ECG, external validation).

1. Does an acute within-patient haemoglobin (Hb) fall shift the ST segment? This is tested by equivalence, with a margin of ±0.01 mV per g/dL.
2. When Hb falls by ≥1 g/dL to ≤8 g/dL, is new ST depression associated with 30-day mortality, and is that association specific to anaemia?

**This repository contains no patient data.** It holds the analysis code, the study protocols and aggregate result files only. HEEDB and MIMIC-IV are available to credentialed users under their data-use agreements (see *Data access*). Depositing this code fulfils the code-sharing clauses of the PhysioNet Credentialed Health Data License 1.5.0 (§8) and the BDSP Credentialed Health Data License 1.5.0 (§9).

## Repository layout

| Path | Contents |
| --- | --- |
| `analysis/` | The full pipeline: ECG measurement, cohort building, protocol analyses and exploratory analyses (scripts `01`–`11`, `ecgmeasure_gpu.py`). |
| `mimic_manifest/` | Builds the MIMIC ECG-pair manifest (`manifests_aligned/pairs_all.parquet`) used by `analysis/05b`. `ecganemia/` is the minimal package it needs. |
| `protocols/` | The study protocols with their dated change logs, plus literature notes. Written in Turkish; the analysis decisions, thresholds and hashes are language-independent. |
| `results/` | Aggregate outputs (JSON) behind every number in the manuscript. |
| `figures/` | Scripts that draw the figures from `results/`. `make_figures_jecg.py` draws Figures 1–2 and the supplementary data of the current manuscript (written to `figures/jecg/`); `make_figure_and_supplement.py` draws Figure 1 and the Online Data Supplement of the first submission. |

## Frozen code

**Registration status.** The protocols were kept as local documents with dated change logs and were not registered on a public platform before analysis; the first public timestamp of this repository (2 October 2026) is after the analyses were run. The dates and hashes in the protocol change logs are therefore the author's own record and cannot be verified independently.

Order of events, as recorded in the change logs: the Analysis 1 protocol was written before any Hb–ST analysis; `06_analysis.py` was run once on MGH and its hash recorded, then run unchanged once on MIMIC. Analysis 2 was formulated after the Analysis 1 results and exploratory analyses; crude 30-day mortality by new ST depression was tabulated in both centres (`09_feasibility_tolerance.py`) before its protocol was locked; `10_build_tolerance_mgh.py` and `11_tolerance_analysis.py` were then frozen and hashed before any model was run. `SHA1SUMS` lists every file in `analysis/`. The files are byte-identical to the versions that were run.

The docstrings of `06_analysis.py` and `11_tolerance_analysis.py` still say "pre-registered". That wording predates the correction in v1.0.1 and means prespecified in the internal protocol. The files were left unchanged so that their hashes still match.

| Step | Script | sha1 (first 8) |
| --- | --- | --- |
| ECG measurement (GPU median beat + Kors VCG) | `analysis/ecgmeasure_gpu.py` | f47a1270 |
| Batch measurement | `analysis/04_measure.py` | de33fc0b |
| Analysis 1, MGH table | `analysis/05_build_analysis_mgh.py` | 1e05bd4f |
| Analysis 1, MIMIC table | `analysis/05b_build_analysis_mimic.py` | ffff2916 |
| Analysis 1, protocol models | `analysis/06_analysis.py` | 788c5b4d |
| Analysis 2, MGH table | `analysis/10_build_tolerance_mgh.py` | df7fae37 |
| Analysis 2, MIMIC table | `analysis/10b_build_tolerance_mimic.py` | 3bfac8b4 |
| Analysis 2, protocol models | `analysis/11_tolerance_analysis.py` | 5954a4ad |

Exploratory analyses, labelled as such in the paper:

- `07_exploratory.py`: centre heterogeneity, definite bleeding
- `08_tail_st.py`: new ST depression by nadir Hb, troponin
- `09_feasibility_tolerance.py`: feasibility counts for Analysis 2

Feasibility and validation scripts:

- `01`, `02`: HEEDB cohort feasibility
- `03_validate_ludb.py`: QRS onset/offset tuning on odd LUDB records, testing on even records

## Running the pipeline

The scripts were run on Linux with Python 3.13, a CUDA GPU (RTX 5090) and the package versions in `requirements.txt`. Data paths are module-level constants at the top of each script (for example `HEEDB`, `HOSP`, `ROOTS`). They were left unchanged so that the files keep their recorded hashes; edit them to point to your local copies. Intermediate files (`out/*.parquet`) contain patient-level data and are not part of this repository.

1. **MIMIC manifest.**
   - Set the paths in `mimic_manifest/ecganemia/config.py`.
   - Run `python -m ecganemia.data.index`, then `python scripts/30_build_ecg_aligned_longitudinal_manifest.py` from `mimic_manifest/`.
   - Place the resulting `manifests_aligned/` folder next to `analysis/`.
2. **HEEDB pairs:** `python analysis/01_heedb_feasibility.py`. This writes `out/pairs_mgh.parquet`.
3. **Measurement:**
   - `python analysis/03_validate_ludb.py` (LUDB check).
   - `python analysis/04_measure.py --cohort mgh` and `--cohort mimic`.
4. **Analysis 1:**
   - Build the tables with `05_build_analysis_mgh.py` and `05b_build_analysis_mimic.py`.
   - Run `06_analysis.py --cohort mgh`, then `--cohort mimic`.
5. **Analysis 2:**
   - Build the tables with `10_build_tolerance_mgh.py` and `10b_build_tolerance_mimic.py`.
   - Run `11_tolerance_analysis.py --cohort mgh`, then `--cohort mimic`.
6. **Figures and supplementary data:** `python figures/make_figures_jecg.py` (current manuscript). The first submission's version is `python figures/make_figure_and_supplement.py`.

## Data access

- **HEEDB v5.0**, Brain Data Science Platform. Koscova Z, Moura Junior V, Reyna M, et al. Harvard-Emory ECG Database (version 5.0). 2026. doi:10.60508/rv6h-7d10. Publication of HEEDB is supported by NHLBI grant R01HL161253.
- **MIMIC-IV v3.1**, PhysioNet. doi:10.13026/kpb9-mt58.
- **MIMIC-IV-ECG v1.0**, PhysioNet. doi:10.13026/4nqg-sb35.
- **LUDB**, PhysioNet. Kalyakulina et al., IEEE Access 2020.

Please also cite Goldberger et al., *Circulation* 2000;101:e215–e220 (PhysioNet).

## Citation

Koksal F. Acute haemoglobin decline and the ST segment: analysis code. Zenodo; 2026. doi:10.5281/zenodo.23110835 (concept DOI; always resolves to the latest version). v1.0.0: doi:10.5281/zenodo.23110836. See also `CITATION.cff`.

## Version history

- **v1.0.1** (October 2026): documentation only. Removed the description of the study as "pre-registered" and of the protocols as "time-stamped", added the registration status and order of events, and updated the study title. Added `figures/make_figures_jecg.py`, which draws the current manuscript's figures, including the lead-level Fig. 2, from the existing result files. No file in `analysis/`, `protocols/` or `results/` changed; in `make_figure_and_supplement.py`, only one descriptive string was reworded.
- **v1.0.0** (2 October 2026): first public release.

## Licence

MIT (code). The data are not covered by this licence; their use is governed by the respective data-use agreements.
