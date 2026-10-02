"""Feasibility: ECG 'anaemia intolerance' study.

Index event per patient = earliest analysable pair (frozen 06.prep exclusions) whose post ECG has
Hb <= 8 g/dL after an acute fall (ΔHb >= 1 g/dL). Counts: cohort size, new ST depression, troponin
measurement/rise, deaths at 30 d / 1 y with follow-up-complete denominators.
MGH death: Massachusetts registry (DateOfDeathMARegistryData; complete 2011-01..2023-04).
MIMIC death: patients.dod (hospital + state records, censored 1 y after last discharge).
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
HEEDB = Path("/home/fatih/datasets/heedb")
HOSP = Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimiciv/3.1/hosp")
spec = importlib.util.spec_from_file_location("a06", HERE / "06_analysis.py")
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
spec8 = importlib.util.spec_from_file_location("a08", HERE / "08_tail_st.py")
LEADS = ["I", "II", "III", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]
MGH_REG_END = pd.Timestamp("2023-04-30")


def index_cohort(cohort):
    raw = pd.read_parquet(HERE / f"out/analysis_{cohort}.parquet")
    if "post_key" not in raw:
        raw["post_key"] = raw.post_fname
    df, _ = A.prep(raw)
    c = df[(df.post_hb <= 8) & (df.delta_hb >= 1)].sort_values("post_time")
    idx = c.groupby("pid").head(1).copy()
    hit = sum(((idx[f"post_st60_{l}"] <= -0.05) & ((idx[f"post_st60_{l}"] - idx[f"base_st60_{l}"]) <= -0.05)).astype(int)
              for l in LEADS)
    idx["new_std"] = hit >= 2
    idx["d_hr"] = idx.post_hr - idx.base_hr
    return idx


con = duckdb.connect()
con.execute(f"SET temp_directory='{HEEDB / 'duckdb_tmp'}'; SET threads=24; SET memory_limit='40GB'")
out = {}
for cohort in ["mgh", "mimic"]:
    idx = index_cohort(cohort)
    con.register("p", idx[["pair_id", "pid", "base_time", "post_time"]])
    if cohort == "mgh":
        death = con.execute(f"""SELECT CAST(BDSPPatientID AS INTEGER) pid, max(CAST(DateOfDeathMARegistryData AS TIMESTAMP)) dod
            FROM read_csv('{HEEDB}/ECG/I0001/metadata/metadata.csv', header=true, all_varchar=true)
            WHERE CAST(BDSPPatientID AS INTEGER) IN (SELECT pid FROM p) GROUP BY 1""").df()
        end = MGH_REG_END
        trop = con.execute(f"""
        WITH t AS (SELECT BDSPPatientID pid, CAST(COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) AS TIMESTAMP) t,
                   TRY_CAST(regexp_replace(ResultTXT,'[<>]','','g') AS DOUBLE)/CASE LoincTXT WHEN '67151-1' THEN 19.0 ELSE 0.01 END v
                   FROM read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet')
                   WHERE LoincTXT IN ('67151-1','6598-7') AND BDSPPatientID IN (SELECT pid FROM p))
        SELECT p.pair_id, count(t.v) n, max(t.v) mx, min(t.v) mn FROM p JOIN t ON t.pid=p.pid
          AND t.t BETWEEN p.post_time - INTERVAL 24 HOUR AND p.post_time + INTERVAL 48 HOUR WHERE t.v IS NOT NULL GROUP BY 1""").df()
    else:
        death = pd.read_csv(HOSP / "patients.csv.gz", usecols=["subject_id", "dod"]).rename(columns={"subject_id": "pid"})
        death["dod"] = pd.to_datetime(death.dod)
        end = None
        trop = con.execute(f"""
        WITH t AS (SELECT subject_id pid, charttime t, valuenum/0.01 v FROM read_csv_auto('{HOSP / "labevents.csv.gz"}',
                   types={{'valuenum':'DOUBLE','value':'VARCHAR','comments':'VARCHAR'}})
                   WHERE itemid=51003 AND valuenum IS NOT NULL AND subject_id IN (SELECT pid FROM p))
        SELECT p.pair_id, count(t.v) n, max(t.v) mx, min(t.v) mn FROM p JOIN t ON t.pid=p.pid
          AND t.t BETWEEN p.post_time - INTERVAL 24 HOUR AND p.post_time + INTERVAL 48 HOUR GROUP BY 1""").df()
    idx = idx.merge(death, on="pid", how="left").merge(trop, on="pair_id", how="left")
    idx["days_to_death"] = (idx.dod.dt.normalize() - idx.post_time.dt.normalize()).dt.days
    r = {"n_patients": int(len(idx)), "new_std_rate": float(idx.new_std.mean()), "n_new_std": int(idx.new_std.sum()),
         "post_hb_median": float(idx.post_hb.median()), "delta_hb_median": float(idx.delta_hb.median()),
         "post_hb_le7": int((idx.post_hb <= 7).sum()),
         "trop_measured_frac": float((idx.n.fillna(0) >= 1).mean()),
         "trop_above_url_frac_of_measured": float((idx.mx > 1)[idx.n.fillna(0) >= 1].mean())}
    for h in [30, 365]:
        if end is not None:
            ok = idx.post_time + pd.Timedelta(days=h) <= end
        else:  # MIMIC: dod complete to 1 y after last discharge; assume complete for index-based 30 d/1 y
            ok = pd.Series(True, index=idx.index)
        d = idx[ok]
        dead = d.days_to_death.between(0, h)
        r[f"death_{h}d"] = {"n_followup_complete": int(len(d)), "deaths": int(dead.sum()), "rate": float(dead.mean()),
                            "rate_new_std": float(dead[d.new_std].mean()), "rate_no_std": float(dead[~d.new_std].mean()),
                            "deaths_new_std": int(dead[d.new_std].sum())}
    r["by_year"] = idx.post_time.dt.year.value_counts().sort_index().to_dict()
    out[cohort] = r
(HERE / "out/feasibility_tolerance.json").write_text(json.dumps(out, indent=1, default=str))
print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "by_year"} for k, v in out.items()}, indent=1))
