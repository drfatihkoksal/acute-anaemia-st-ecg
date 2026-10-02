"""EXPLORATORY (not pre-registered; concept.md §11.9 2026-10-02 entry):
 (1) between-center heterogeneity of the H2a QRS effects
     a. MGH restricted to pairs whose two Hb draws share one encounter (proxy for MIMIC's same-admission rule)
     b. MGH reweighted to the MIMIC covariate distribution (inverse-odds weights)
     c. descriptive comparison
 (2) replication in definite bleeding: D62 / ICD-9 285.1 (acute posthemorrhagic anaemia), and major bleeding codes.
Model functions are imported unchanged from the frozen 06_analysis.py.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.linear_model import LogisticRegression

HERE = Path(__file__).resolve().parent
HEEDB = Path("/home/fatih/datasets/heedb")
HOSP = Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimiciv/3.1/hosp")
spec = importlib.util.spec_from_file_location("a06", HERE / "06_analysis.py")
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)

OUTC = ["st60_min", "m_early", "m_late", "sum_qrs_pp", "sum_p_abs"]
D62_10, D62_9 = ["D62"], ["2851"]
MAJOR_10 = ["K920", "K921", "K922", "K250", "K252", "K254", "K256", "K260", "K262", "K264", "K266", "K270", "K272",
            "K274", "K276", "K280", "K282", "K284", "K286", "I8501", "I8511", "K625", "K226", "T810", "D62"]
MAJOR_9 = ["578", "5310", "5312", "5314", "5316", "5320", "5322", "5324", "5326", "5330", "5332", "5334", "5336",
           "5340", "5342", "5344", "5346", "4560", "45620", "5693", "53082", "9981", "2851"]


def like(col, codes):
    return " OR ".join(f"{col} LIKE '{c}%'" for c in codes)


con = duckdb.connect()
con.execute(f"SET temp_directory='{HEEDB / 'duckdb_tmp'}'; SET threads=24; SET memory_limit='40GB'")

# ------------------------------------------------------------------ MGH extras
mgh = pd.read_parquet(HERE / "out/analysis_mgh.parquet")
con.register("p", mgh[["pair_id", "pid", "post_time", "base_hb_time", "post_hb_time"]])
enc = con.execute(f"""
WITH hb AS (
  SELECT BDSPPatientID pid, CAST(COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) AS TIMESTAMP) t,
         min(BDSPEncounterID) enc
  FROM read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet')
  WHERE LoincTXT IN ('718-7','30313-1','30350-3','20509-6') AND BDSPPatientID IN (SELECT DISTINCT pid FROM p)
  GROUP BY 1, 2)
SELECT p.pair_id, b.enc base_enc, q.enc post_enc FROM p
LEFT JOIN hb b ON b.pid = p.pid AND b.t = p.base_hb_time
LEFT JOIN hb q ON q.pid = p.pid AND q.t = p.post_hb_time""").df()
mgh = mgh.merge(enc, on="pair_id", how="left")
mgh["same_enc"] = (mgh.base_enc == mgh.post_enc) & mgh.base_enc.notna()


def mgh_flag(codes10, codes9):
    q = f"""
    WITH b AS (
      SELECT CAST(BDSPPatientID AS INTEGER) pid, CAST(RECORDED_DT AS TIMESTAMP) t
      FROM read_csv('{HEEDB}/ECG/I0001/ICD_codes/icd10_codes.csv', header=true, all_varchar=true)
      WHERE {like("replace(DIAGNOSIS_ICD10_CD, '.', '')", codes10)}
      UNION ALL
      SELECT CAST(BDSPPatientID AS INTEGER) pid, CAST(RECORDED_DT AS TIMESTAMP) t
      FROM read_csv('{HEEDB}/ECG/I0001/ICD_codes/icd9_codes.csv', header=true, all_varchar=true)
      WHERE {like("replace(DIAGNOSIS_ICD9_CD, '.', '')", codes9)})
    SELECT DISTINCT p.pair_id FROM p JOIN b ON p.pid = b.pid
     AND b.t BETWEEN p.post_time - INTERVAL 3 DAY AND p.post_time + INTERVAL 3 DAY"""
    return set(con.execute(q).df().pair_id)


mgh["bleed_d62"] = mgh.pair_id.isin(mgh_flag(D62_10, D62_9))
mgh["bleed_major"] = mgh.pair_id.isin(mgh_flag(MAJOR_10, MAJOR_9))

# ------------------------------------------------------------------ MIMIC extras
mim = pd.read_parquet(HERE / "out/analysis_mimic.parquet")
dx = con.execute(f"""SELECT hadm_id, icd_version, replace(icd_code,'.','') code, seq_num
  FROM read_csv_auto('{HOSP / "diagnoses_icd.csv.gz"}', types={{'icd_code':'VARCHAR'}})""").df()
con.register("dx", dx)


def mim_flag(codes10, codes9, principal=False):
    extra = "AND seq_num = 1" if principal else ""
    q = f"""SELECT DISTINCT hadm_id FROM dx WHERE ((icd_version=10 AND ({like('code', codes10)}))
            OR (icd_version=9 AND ({like('code', codes9)}))) {extra}"""
    return set(con.execute(q).df().hadm_id)


mim["bleed_d62"] = mim.hadm_id.isin(mim_flag(D62_10, D62_9))
mim["bleed_major"] = mim.hadm_id.isin(mim_flag(MAJOR_10, MAJOR_9))
mim["bleed_major_principal"] = mim.hadm_id.isin(mim_flag(MAJOR_10, MAJOR_9, principal=True))
mim["same_enc"] = True

# ------------------------------------------------------------------ analyses
mgh["post_key"] = mgh.post_fname  # as in 06 main()
Pm, flow_m = A.prep(mgh); Pm = Pm[Pm.primary]
Pi, flow_i = A.prep(mim); Pi = Pi[Pi.primary]
fit = lambda d, y: A.fit(d, f"d_{y}", impute=False)
R = {"note": "EXPLORATORY, not pre-registered"}


def block(d):
    return {"n": int(len(d)), **{y: fit(d, y) for y in OUTC}}


R["mgh_all"] = block(Pm)
R["mimic_all"] = block(Pi)
R["mgh_same_encounter"] = block(Pm[Pm.same_enc])
R["mgh_diff_encounter"] = block(Pm[~Pm.same_enc])
R["mgh_same_enc_frac"] = float(Pm.same_enc.mean())

# inverse-odds reweighting of MGH to MIMIC
cov = ["age", "male", "gap_h", "base_hb", "delta_hb", "base_hr", "base_m_qrs", "base_qrs_dur_ms", "base_k", "base_cr"]
both = pd.concat([Pm.assign(_c=0), Pi.assign(_c=1)])[cov + ["_c"]].dropna()
lr = LogisticRegression(max_iter=2000).fit((both[cov] - both[cov].mean()) / both[cov].std(), both._c)
Xm = Pm[cov].dropna()
ps = lr.predict_proba((Xm - both[cov].mean()) / both[cov].std())[:, 1]
w = pd.Series(ps / (1 - ps), index=Xm.index)
w = (w.clip(upper=w.quantile(0.99)) / w.mean())


def wfit(d, y, w):
    y_ = f"d_{y}"
    cols = ["delta_hb"] + A.ADJ + [f"base_{y}"]
    dd = d.loc[w.index, [y_] + cols].dropna()
    r = sm.WLS(dd[y_], sm.add_constant(dd[cols]), weights=w.loc[dd.index]).fit(cov_type="HC3")
    b, se = r.params["delta_hb"], r.bse["delta_hb"]
    return {"beta": float(b), "se": float(se), "ci95": [float(b - 1.96 * se), float(b + 1.96 * se)],
            "p": float(r.pvalues["delta_hb"]), "n": int(r.nobs)}


R["mgh_reweighted_to_mimic"] = {"ess": float(w.sum() ** 2 / (w ** 2).sum()), **{y: wfit(Pm, y, w) for y in OUTC}}
R["descriptives"] = {
    n: d[["age", "male", "gap_h", "base_hb", "post_hb", "delta_hb", "base_hr", "base_m_qrs", "base_qrs_dur_ms",
          "base_alb", "base_cr"]].describe().loc[["mean", "50%"]].round(3).to_dict()
    for n, d in [("mgh", Pm), ("mimic", Pi)]}

# bleeding
for nm, P in [("mgh", Pm), ("mimic", Pi)]:
    R[f"{nm}_bleed_d62"] = block(P[P.bleed_d62])
    R[f"{nm}_bleed_major"] = block(P[P.bleed_major])
R["mimic_bleed_major_principal"] = block(Pi[Pi.bleed_major_principal])

# pooled bleeding estimate (fixed-effect inverse variance) for D62
R["pooled_bleed_d62"] = {}
for y in OUTC:
    a, b = R["mgh_bleed_d62"][y], R["mimic_bleed_d62"][y]
    wa, wb = 1 / a["se"] ** 2, 1 / b["se"] ** 2
    est = (a["beta"] * wa + b["beta"] * wb) / (wa + wb); se = (wa + wb) ** -0.5
    q = wa * (a["beta"] - est) ** 2 + wb * (b["beta"] - est) ** 2
    R["pooled_bleed_d62"][y] = {"beta": est, "ci95": [est - 1.96 * se, est + 1.96 * se],
                                "p": float(2 * stats.norm.sf(abs(est / se))), "het_p": float(stats.chi2.sf(q, 1))}
(HERE / "out/exploratory_heterogeneity_bleeding.json").write_text(json.dumps(R, indent=1, default=float))
print("done")
