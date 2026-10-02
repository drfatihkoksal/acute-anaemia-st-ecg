"""HEEDB-MGH feasibility for the clinical question MIMIC could not power (n_true=40):
acute Hb drop + ischemic-looking post ECG -> true ischemia vs pseudo-ischemia.

ECG ischemic-looking: 12SL v24 statement of ST depression / subendocardial injury,
injury/ST elevation, acute MI/ACS, or "T wave abnormality, consider ischemia".
True ischemia (4th UDMI-like): hs-TnT gen5 (LOINC 67151-1) peak > 19 ng/L (MGH 99th
pct URL used at MGH) with >=20% rise within [-24h, +48h] of the post ECG, OR ICD-10
I21/I22 recorded within [-3d, +14d].  Pseudo: >=2 hs-TnT values in window, all <= 19.
Older conventional TnT (6598-7, ng/mL, URL 0.01) treated analogously.
"""
from __future__ import annotations

import json
from pathlib import Path

import duckdb

HEEDB = Path("/home/fatih/datasets/heedb")
OUT = Path(__file__).resolve().parent / "out"
ISCH_CODES = ([1002, 1022, 1024, 1040, 1050, 1060, 1070, 1071, 1080, 1081, 1082]
              + list(range(920, 969)) + [821, 822, 826, 827, 828, 829]
              + [1145, 1150, 1151, 1160, 1161, 1170, 1171, 1172, 1180, 1181])

con = duckdb.connect()
con.execute(f"SET temp_directory='{HEEDB / 'duckdb_tmp'}'")
con.execute("SET threads=24; SET memory_limit='40GB'")
con.execute(f"CREATE TABLE pairs AS SELECT * FROM '{OUT}/pairs_mgh.parquet'")

con.execute(f"""
CREATE TABLE dx AS
SELECT regexp_replace(regexp_replace(regexp_replace(FileName, '\\s', '', 'g'), '^\\.', ''), '\\.hea$', '') AS fname,
       list_transform(string_split(regexp_replace(codes, '\\s', '', 'g'), ','), x -> TRY_CAST(x AS INTEGER)) AS codes
FROM read_csv('{HEEDB}/ECG/I0001/12SL_diagnoses/diagnoses_v24.csv', header=true, all_varchar=true)
""")
con.execute(f"""
CREATE TABLE p AS
SELECT p.*,
       list_has_any(d.codes, {ISCH_CODES}) AS post_isch,
       list_has_any(db.codes, {ISCH_CODES}) AS base_isch,
       d.codes IS NOT NULL AS has_dx
FROM pairs p
LEFT JOIN dx d ON regexp_replace(p.post_fname, '\\s', '', 'g') = d.fname
LEFT JOIN dx db ON regexp_replace(p.base_fname, '\\s', '', 'g') = db.fname
""")

con.execute(f"""
CREATE TABLE trop AS
SELECT BDSPPatientID pid,
       CAST(COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) AS TIMESTAMP) t,
       TRY_CAST(regexp_replace(ResultTXT, '[<>]', '', 'g') AS DOUBLE) v, LoincTXT loinc
FROM read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet')
WHERE LoincTXT IN ('67151-1', '6598-7')
  AND COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) IS NOT NULL
""")
con.execute("DELETE FROM trop WHERE v IS NULL")
con.execute("UPDATE trop SET v = v / 19.0 WHERE loinc = '67151-1'")   # multiples of URL
con.execute("UPDATE trop SET v = v / 0.01 WHERE loinc = '6598-7'")

con.execute(f"""
CREATE TABLE mi AS
SELECT CAST(BDSPPatientID AS INTEGER) pid, CAST(RECORDED_DT AS TIMESTAMP) t
FROM read_csv('{HEEDB}/ECG/I0001/ICD_codes/icd10_codes.csv', header=true, all_varchar=true)
WHERE DIAGNOSIS_ICD10_CD LIKE 'I21%' OR DIAGNOSIS_ICD10_CD LIKE 'I22%'
""")

con.execute("""
CREATE TABLE cand AS
SELECT p.*, row_number() OVER () rid FROM p WHERE delta_hb >= 2 AND post_isch
""")
con.execute("""
CREATE TABLE lab AS
SELECT c.rid, count(t.v) n_trop, max(t.v) peak, min(t.v) lo
FROM cand c LEFT JOIN trop t
  ON CAST(c.pid AS INTEGER) = t.pid
 AND t.t BETWEEN c.post_time - INTERVAL 24 HOUR AND c.post_time + INTERVAL 48 HOUR
GROUP BY c.rid
""")
con.execute("""
CREATE TABLE icd AS
SELECT c.rid, count(m.t) > 0 AS mi_icd
FROM cand c LEFT JOIN mi m
  ON CAST(c.pid AS INTEGER) = m.pid
 AND m.t BETWEEN c.post_time - INTERVAL 3 DAY AND c.post_time + INTERVAL 14 DAY
GROUP BY c.rid
""")
con.execute("""
CREATE TABLE lab_c AS
SELECT c.*, l.n_trop, l.peak, l.lo, i.mi_icd,
  CASE WHEN (l.peak > 1 AND l.peak >= 1.2 * l.lo AND l.n_trop >= 2) OR i.mi_icd THEN 'true'
       WHEN l.n_trop >= 2 AND l.peak <= 1 THEN 'pseudo'
       ELSE 'indeterminate' END AS grp
FROM cand c JOIN lab l USING (rid) JOIN icd i USING (rid)
""")
con.execute(f"COPY lab_c TO '{OUT}/pseudo_true_candidates_mgh.parquet' (FORMAT parquet)")

q = lambda s: con.execute(s).df().to_dict("records")
res = {
    "pairs_with_12sl": q("SELECT count(*) n, avg(has_dx::INT) frac FROM p")[0],
    "post_isch_rate_by_dhb": q("""SELECT CASE WHEN delta_hb>=2 THEN 'drop>=2' WHEN delta_hb<=-1 THEN 'rise>=1'
         WHEN abs(delta_hb)<0.5 THEN 'stable' ELSE 'other' END g, count(*) n,
         avg(post_isch::INT) post_isch, avg((post_isch AND NOT base_isch)::INT) new_isch
         FROM p WHERE has_dx GROUP BY 1 ORDER BY 1"""),
    "acute_drop_ischemic_ecg": q("SELECT grp, count(*) n, count(DISTINCT pid) subjects FROM lab_c GROUP BY 1"),
    "acute_drop_NEW_ischemic_ecg": q("SELECT grp, count(*) n, count(DISTINCT pid) subjects FROM lab_c WHERE NOT coalesce(base_isch,false) GROUP BY 1"),
    "true_by_source": q("SELECT mi_icd, (peak>1 AND n_trop>=2) trop_pos, count(*) n FROM lab_c WHERE grp='true' GROUP BY ALL"),
}
(OUT / "pseudo_true_feasibility_mgh.json").write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str))
