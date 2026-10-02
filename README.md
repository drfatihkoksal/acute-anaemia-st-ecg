# Acute haemoglobin decline and ST-segment depression — analysis code

Code for the study *"Acute haemoglobin decline and ST-segment depression in hospitalised patients: a pre-registered, two-centre, within-patient electrocardiographic study"* (Koksal F; manuscript submitted to the *International Journal of Cardiology*).

The study asks two pre-specified questions in serial 12-lead ECGs from two hospitals: Massachusetts General Hospital (HEEDB, discovery) and Beth Israel Deaconess Medical Center (MIMIC-IV / MIMIC-IV-ECG, external validation).

1. Does an acute within-patient haemoglobin (Hb) fall shift the ST segment? This is tested by equivalence, with a margin of ±0.01 mV per g/dL.
2. When Hb falls by ≥1 g/dL to ≤8 g/dL, is new ST depression associated with 30-day mortality, and is that association specific to anaemia?

**This repository contains no patient data.** It holds the analysis code, the locked protocols and aggregate result files only. HEEDB and MIMIC-IV are available to credentialed users under their data-use agreements (see *Data access*). Depositing this code fulfils the code-sharing clauses of the PhysioNet Credentialed Health Data License 1.5.0 (§8) and the BDSP Credentialed Health Data License 1.5.0 (§9).

## Repository layout

| Path | Contents |
| --- | --- |
| `analysis/` | The full pipeline: ECG measurement, cohort building, pre-registered analyses and exploratory analyses (scripts `01`–`11`, `ecgmeasure_gpu.py`). |
| `mimic_manifest/` | Builds the MIMIC ECG-pair manifest (`manifests_aligned/pairs_all.parquet`) used by `analysis/05b`. `ecganemia/` is the minimal package it needs. |
| `protocols/` | The locked, time-stamped protocols with their change logs, plus literature notes. Written in Turkish; the analysis decisions, thresholds and hashes are language-independent. |
| `results/` | Aggregate outputs (JSON) behind every number in the manuscript. |
| `figures/` | Script that draws Figure 1 and writes the Online Data Supplement from `results/`. |

## Frozen code

Hashes were recorded in the protocols before each outcome analysis. `SHA1SUMS` lists every file in `analysis/`. The files are byte-identical to the versions that were run.

| Step | Script | sha1 (first 8) |
| --- | --- | --- |
| ECG measurement (GPU median beat + Kors VCG) | `analysis/ecgmeasure_gpu.py` | f47a1270 |
| Batch measurement | `analysis/04_measure.py` | de33fc0b |
| Analysis 1, MGH table | `analysis/05_build_analysis_mgh.py` | 1e05bd4f |
| Analysis 1, MIMIC table | `analysis/05b_build_analysis_mimic.py` | ffff2916 |
| Analysis 1, pre-registered models | `analysis/06_analysis.py` | 788c5b4d |
| Analysis 2, MGH table | `analysis/10_build_tolerance_mgh.py` | df7fae37 |
| Analysis 2, MIMIC table | `analysis/10b_build_tolerance_mimic.py` | 3bfac8b4 |
| Analysis 2, pre-registered models | `analysis/11_tolerance_analysis.py` | 5954a4ad |

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
6. **Figure and supplement:** `python figures/make_figure_and_supplement.py`.

## Data access

- **HEEDB v5.0**, Brain Data Science Platform. Koscova Z, Moura Junior V, Reyna M, et al. Harvard-Emory ECG Database (version 5.0). 2026. doi:10.60508/rv6h-7d10. Publication of HEEDB is supported by NHLBI grant R01HL161253.
- **MIMIC-IV v3.1**, PhysioNet. doi:10.13026/kpb9-mt58.
- **MIMIC-IV-ECG v1.0**, PhysioNet. doi:10.13026/4nqg-sb35.
- **LUDB**, PhysioNet. Kalyakulina et al., IEEE Access 2020.

Please also cite Goldberger et al., *Circulation* 2000;101:e215–e220 (PhysioNet).

## Citation

See `CITATION.cff`. A Zenodo DOI will be added at the first release.

## Licence

MIT (code). The data are not covered by this licence; their use is governed by the respective data-use agreements.
