"""EXPLORATORY: does acute Hb fall produce ST depression in a SUBSET (tail), especially with low nadir Hb,
known CAD, and myocardial injury (type 2 MI)? Mean-shift analysis (H1) cannot see a minority effect.

Per pair (primary one-per-patient set from frozen 06.prep):
  new_std     : >=2 non-aVR leads with post ST60 <= -0.05 mV AND ΔST60 <= -0.05 mV   (0.5 mm, new)
  new_std_1mm : same with 0.10 mV
  trop_rise   : troponin in [base ECG, post ECG + 24 h] with max > URL and max >= 1.2 x min (>=2 values)
  cad         : CAD ICD (I20-I25 / 410-414) before the post ECG (MGH) or in the same admission (MIMIC)
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

HERE = Path(__file__).resolve().parent
HEEDB = Path("/home/fatih/datasets/heedb")
HOSP = Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimiciv/3.1/hosp")
spec = importlib.util.spec_from_file_location("a06", HERE / "06_analysis.py")
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
LEADS = ["I", "II", "III", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]

con = duckdb.connect()
con.execute(f"SET temp_directory='{HEEDB / 'duckdb_tmp'}'; SET threads=24; SET memory_limit='40GB'")


def mgh_extras(df):
    con.register("p", df[["pair_id", "pid", "base_time", "post_time"]])
    trop = con.execute(f"""
    WITH t AS (
      SELECT BDSPPatientID pid, CAST(COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) AS TIMESTAMP) t,
             TRY_CAST(regexp_replace(ResultTXT, '[<>]', '', 'g') AS DOUBLE)
               / CASE LoincTXT WHEN '67151-1' THEN 19.0 ELSE 0.01 END v
      FROM read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet')
      WHERE LoincTXT IN ('67151-1','6598-7') AND BDSPPatientID IN (SELECT DISTINCT pid FROM p))
    SELECT p.pair_id, count(t.v) n, max(t.v) mx, min(t.v) mn FROM p JOIN t ON t.pid = p.pid
     AND t.t BETWEEN p.base_time AND p.post_time + INTERVAL 24 HOUR WHERE t.v IS NOT NULL GROUP BY 1""").df()
    cad = con.execute(f"""
    WITH c AS (
      SELECT CAST(BDSPPatientID AS INTEGER) pid, CAST(RECORDED_DT AS TIMESTAMP) t
      FROM read_csv('{HEEDB}/ECG/I0001/ICD_codes/icd10_codes.csv', header=true, all_varchar=true)
      WHERE regexp_matches(DIAGNOSIS_ICD10_CD, '^I2[0-5]')
      UNION ALL
      SELECT CAST(BDSPPatientID AS INTEGER) pid, CAST(RECORDED_DT AS TIMESTAMP) t
      FROM read_csv('{HEEDB}/ECG/I0001/ICD_codes/icd9_codes.csv', header=true, all_varchar=true)
      WHERE regexp_matches(DIAGNOSIS_ICD9_CD, '^41[0-4]'))
    SELECT DISTINCT p.pair_id FROM p JOIN c ON c.pid = p.pid AND c.t <= p.post_time""").df()
    return trop, set(cad.pair_id)


def mimic_extras(df):
    con.register("p", df[["pair_id", "pid", "base_time", "post_time", "hadm_id"]])
    trop = con.execute(f"""
    WITH t AS (SELECT subject_id pid, charttime t, valuenum / 0.01 v
               FROM read_csv_auto('{HOSP / "labevents.csv.gz"}', types={{'valuenum':'DOUBLE','value':'VARCHAR','comments':'VARCHAR'}})
               WHERE itemid = 51003 AND valuenum IS NOT NULL AND subject_id IN (SELECT DISTINCT pid FROM p))
    SELECT p.pair_id, count(t.v) n, max(t.v) mx, min(t.v) mn FROM p JOIN t ON t.pid = p.pid
     AND t.t BETWEEN p.base_time AND p.post_time + INTERVAL 24 HOUR GROUP BY 1""").df()
    cad = con.execute(f"""SELECT DISTINCT p.pair_id FROM p JOIN read_csv_auto('{HOSP / "diagnoses_icd.csv.gz"}', types={{'icd_code':'VARCHAR'}}) d
      ON d.hadm_id = p.hadm_id WHERE (d.icd_version = 10 AND regexp_matches(d.icd_code, '^I2[0-5]'))
                               OR (d.icd_version = 9 AND regexp_matches(d.icd_code, '^41[0-4]'))""").df()
    return trop, set(cad.pair_id)


def run(cohort):
    raw = pd.read_parquet(HERE / f"out/analysis_{cohort}.parquet")
    if "post_key" not in raw:
        raw["post_key"] = raw.post_fname
    df, _ = A.prep(raw)
    P = df[df.primary].copy()
    trop, cad = (mgh_extras if cohort == "mgh" else mimic_extras)(P)
    P = P.merge(trop, on="pair_id", how="left")
    P["trop_measured"] = P.n.fillna(0) >= 2
    P["trop_rise"] = P.trop_measured & (P.mx > 1) & (P.mx >= 1.2 * P.mn)
    P["cad"] = P.pair_id.isin(cad)
    for thr, nm in [(0.05, "new_std"), (0.10, "new_std_1mm")]:
        hit = sum(((P[f"post_st60_{l}"] <= -thr) & ((P[f"post_st60_{l}"] - P[f"base_st60_{l}"]) <= -thr)).astype(int)
                  for l in LEADS)
        P[nm] = (hit >= 2).astype(int)
    P["dhb_cat"] = pd.cut(P.delta_hb, [-np.inf, -0.5, 0.5, 1.5, 2.5, np.inf],
                          labels=["rise", "stable", "1-1.5", "1.5-2.5", ">=2.5"])
    P["nadir_cat"] = pd.cut(P.post_hb, [0, 7, 8, 10, np.inf], labels=["<=7", "7-8", "8-10", ">10"])
    P["d_hr"] = P.post_hr - P.base_hr
    R = {"n": len(P)}
    for nm in ["new_std", "new_std_1mm"]:
        R[f"{nm}_by_dhb"] = P.groupby("dhb_cat", observed=True)[nm].agg(["mean", "size"]).round(4).to_dict("index")
        R[f"{nm}_by_nadir_drop>=1.5"] = P[P.delta_hb >= 1.5].groupby("nadir_cat", observed=True)[nm].agg(["mean", "size"]).round(4).to_dict("index")
        R[f"{nm}_by_nadir_stable"] = P[P.delta_hb.abs() < 0.5].groupby("nadir_cat", observed=True)[nm].agg(["mean", "size"]).round(4).to_dict("index")
        R[f"{nm}_cad_x_drop"] = P.assign(drop=P.delta_hb >= 1.5).groupby(["cad", "drop"])[nm].agg(["mean", "size"]).round(4).reset_index().to_dict("records")
        d = P.dropna(subset=["delta_hb", "post_hb", "age", "male"])
        forms = {"no_hr": f"{nm} ~ delta_hb + I(post_hb<=8) + cad + age + male",
                 "with_dhr": f"{nm} ~ delta_hb + I(post_hb<=8) + cad + age + male + d_hr",
                 "interaction": f"{nm} ~ delta_hb * cad + delta_hb * I(post_hb<=8) + age + male"}
        R[f"{nm}_logit"] = {}
        for k, f in forms.items():
            m = smf.logit(f, d.dropna(subset=["d_hr"]) if "d_hr" in f else d).fit(disp=0)
            ci = m.conf_int()
            R[f"{nm}_logit"][k] = {t: {"OR": float(np.exp(m.params[t])), "ci": [float(np.exp(ci.loc[t, 0])), float(np.exp(ci.loc[t, 1]))],
                                       "p": float(m.pvalues[t])} for t in m.params.index if t != "Intercept"}
    T = P[P.trop_measured]
    R["trop"] = {"measured_frac": float(P.trop_measured.mean()),
                 "rise_by_dhb": T.groupby("dhb_cat", observed=True).trop_rise.agg(["mean", "size"]).round(4).to_dict("index"),
                 "new_std_by_trop_rise": T.groupby("trop_rise").new_std.agg(["mean", "size"]).round(4).to_dict("index"),
                 "new_std_by_trop_rise_drop>=1.5": T[T.delta_hb >= 1.5].groupby("trop_rise").new_std.agg(["mean", "size"]).round(4).to_dict("index")}
    return R


out = {"note": "EXPLORATORY", "mgh": run("mgh"), "mimic": run("mimic")}
(HERE / "out/exploratory_tail_st.json").write_text(json.dumps(out, indent=1, default=str))
print("done")
