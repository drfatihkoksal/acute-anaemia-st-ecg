"""HEEDB-MGH feasibility: replicate the MIMIC aligned longitudinal ECG/Hb pair cohort.

MIMIC definition (scripts/30): each ECG matched to the closest Hb within +/-12 h;
consecutive distinct ECGs with distinct Hb measurements, 0 < gap <= 7 d, form a pair.
MIMIC additionally required both ECGs inside one admission; HEEDB ECG metadata has no
encounter link, so here the Hb-within-12h requirement is the only "in-care" proxy.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import duckdb

HEEDB = Path("/home/fatih/datasets/heedb")
OUT = Path(__file__).resolve().parent / "out"
OUT.mkdir(exist_ok=True)

con = duckdb.connect()
con.execute(f"SET temp_directory='{HEEDB / 'duckdb_tmp'}'")
con.execute("SET threads=24")
con.execute("SET memory_limit='40GB'")
t0 = time.time()

con.execute(f"""
CREATE TABLE hb AS
SELECT BDSPPatientID AS pid,
       CAST(COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) AS TIMESTAMP) AS hb_time,
       TRY_CAST(ResultTXT AS DOUBLE) AS hb,
       LoincTXT AS loinc
FROM read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet')
WHERE LoincTXT IN ('718-7','30313-1','30350-3','20509-6')
  AND coalesce(ResultStatusDSC,'Final') IN ('Final','Corrected')
  AND COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) IS NOT NULL
  AND TRY_CAST(ResultTXT AS DOUBLE) BETWEEN 2 AND 25
""")
# collapse duplicate draws (same patient, same timestamp)
con.execute("CREATE TABLE hb1 AS SELECT pid, hb_time, avg(hb) hb FROM hb GROUP BY pid, hb_time")
print("hb rows", con.execute("SELECT count(*), count(DISTINCT pid) FROM hb1").fetchone(), time.time() - t0)

con.execute(f"""
CREATE TABLE ecg AS
SELECT BDSPPatientID AS pid, FileName AS fname, CAST(ECGAcquisitionTime AS TIMESTAMP) AS ecg_time,
       AgeAtAcquisition AS age_days, SexDSC AS sex
FROM read_csv('{HEEDB}/ECG/I0001/metadata/metadata.csv', header=true, all_varchar=true)
WHERE BDSPPatientID IS NOT NULL AND ECGAcquisitionTime IS NOT NULL
""")
con.execute("UPDATE ecg SET pid = trim(pid)")
print("ecg rows", con.execute("SELECT count(*), count(DISTINCT pid) FROM ecg").fetchone(), time.time() - t0)

con.execute("""
CREATE TABLE m AS
SELECT * EXCLUDE (rn) FROM (
  SELECT e.*, h.hb_time, h.hb,
         abs(epoch(e.ecg_time) - epoch(h.hb_time)) / 3600.0 AS hb_gap_h,
         row_number() OVER (PARTITION BY e.fname
                            ORDER BY abs(epoch(e.ecg_time) - epoch(h.hb_time)), h.hb_time) rn
  FROM ecg e JOIN hb1 h
    ON CAST(e.pid AS INTEGER) = h.pid
   AND h.hb_time BETWEEN e.ecg_time - INTERVAL 12 HOUR AND e.ecg_time + INTERVAL 12 HOUR
) WHERE rn = 1
""")
print("hb-matched ecg", con.execute("SELECT count(*), count(DISTINCT pid) FROM m").fetchone(), time.time() - t0)

con.execute("""
CREATE TABLE pairs AS
SELECT * FROM (
  SELECT pid,
         lag(fname) OVER w AS base_fname, fname AS post_fname,
         lag(ecg_time) OVER w AS base_time, ecg_time AS post_time,
         lag(hb_time) OVER w AS base_hb_time, hb_time AS post_hb_time,
         lag(hb) OVER w AS base_hb, hb AS post_hb,
         age_days, sex
  FROM m WINDOW w AS (PARTITION BY pid ORDER BY ecg_time, fname)
)
WHERE base_fname IS NOT NULL
  AND post_time > base_time
  AND post_time - base_time <= INTERVAL 7 DAY
  AND post_hb_time <> base_hb_time
""")
con.execute("ALTER TABLE pairs ADD COLUMN delta_hb DOUBLE")
con.execute("UPDATE pairs SET delta_hb = base_hb - post_hb")  # positive = Hb fell (MIMIC sign)
con.execute(f"COPY pairs TO '{OUT}/pairs_mgh.parquet' (FORMAT parquet)")

q = lambda s: con.execute(s).df().to_dict("records")
summary = {
    "pairs": q("SELECT count(*) n, count(DISTINCT pid) subjects FROM pairs")[0],
    "delta_hb": q("SELECT avg(delta_hb) mean, stddev(delta_hb) sd, quantile_cont(delta_hb,[0.05,0.25,0.5,0.75,0.95]) q FROM pairs")[0],
    "acute_drop_ge2": q("SELECT count(*) n, count(DISTINCT pid) subjects FROM pairs WHERE delta_hb >= 2")[0],
    "acute_drop_ge3": q("SELECT count(*) n, count(DISTINCT pid) subjects FROM pairs WHERE delta_hb >= 3")[0],
    "by_year": q("SELECT year(post_time) y, count(*) n FROM pairs GROUP BY 1 ORDER BY 1"),
    "elapsed_s": time.time() - t0,
}
summary["delta_hb"]["q"] = list(summary["delta_hb"]["q"])
(OUT / "feasibility_mgh.json").write_text(json.dumps(summary, indent=1, default=str))
print(json.dumps(summary, indent=1, default=str))
