"""Build the MGH pair-level analysis table (protocol v3.1, concept.md §11.4-11.7).

Per pair: ΔHb, base/post ECG features (ecgmeasure_gpu), nearest labs within +-12 h of each ECG
(K, Ca, Mg, Na, BUN, albumin, creatinine), 12SL paced flag, bleeding ICD within +-3 d of the post ECG.
Output: out/analysis_mgh.parquet
"""
from __future__ import annotations

import glob
from pathlib import Path

import duckdb
import pandas as pd

HERE = Path(__file__).resolve().parent
HEEDB = Path("/home/fatih/datasets/heedb")
OUT = HERE / "out"

LABS = {  # LOINC (MGH); total calcium only (17863-2 is mislabelled ionized in mg/dL -> excluded)
    "k": ["2823-3", "6298-4"], "ca": ["17861-6"], "mg": ["2601-3", "19123-9"], "na": ["2951-2", "2947-0"],
    "bun": ["3094-0", "6299-2"], "alb": ["1751-7", "61151-7"], "cr": ["2160-0", "38483-4"],
}
RANGES = {"k": (1.5, 9), "ca": (4, 16), "mg": (0.5, 6), "na": (100, 180), "bun": (1, 250), "alb": (0.5, 7), "cr": (0.1, 25)}
PACED = [183, 184, 185, 186, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 326, 1669]
AF = [161, 162, 273, 288]
BLEED_ICD10 = ["K92.0", "K92.1", "K92.2", "K25.0", "K25.2", "K25.4", "K25.6", "K26.0", "K26.2", "K26.4", "K26.6",
               "K27", "K28", "K62.5", "R31", "N93", "I85.01", "D62", "K22.6", "K22.8", "R58", "T81.0"]
FEATS = ["ok", "hr", "rr_ms", "rr_cv", "qrs_dur_ms", "m_early", "m_late", "m_qrs", "sum_r", "sum_qrs_pp",
         "sum_p_abs", "st60_min", "s_min_V6", "n_beats"] + [f"st60_{l}" for l in
         ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]] + \
        [f"t_amp_{l}" for l in ["I", "II", "V2", "V5", "V6"]] + [f"j_amp_{l}" for l in ["II", "V5"]]

con = duckdb.connect()
con.execute(f"SET temp_directory='{HEEDB / 'duckdb_tmp'}'; SET threads=24; SET memory_limit='40GB'")
pairs = pd.read_parquet(OUT / "pairs_mgh.parquet")
pairs["base_fname"] = pairs.base_fname.str.strip()
pairs["post_fname"] = pairs.post_fname.str.strip()
pairs["pid"] = pairs.pid.astype(int)
pairs["pair_id"] = range(len(pairs))
con.register("pairs", pairs)

# --- labs: nearest value within +-12 h of each ECG -----------------------------------------
ecg = pd.concat([pairs[["pid", "base_fname", "base_time"]].set_axis(["pid", "fname", "t"], axis=1),
                 pairs[["pid", "post_fname", "post_time"]].set_axis(["pid", "fname", "t"], axis=1)]).drop_duplicates("fname")
con.register("ecg", ecg)
case = " ".join(f"WHEN LoincTXT IN ({','.join(repr(c) for c in v)}) THEN '{k}'" for k, v in LABS.items())
allc = ",".join(repr(c) for v in LABS.values() for c in v)
con.execute(f"""
CREATE TABLE lab AS
SELECT BDSPPatientID pid, CASE {case} END var,
       CAST(COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) AS TIMESTAMP) t,
       TRY_CAST(regexp_replace(ResultTXT, '[<>]', '', 'g') AS DOUBLE) v
FROM read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet')
WHERE LoincTXT IN ({allc}) AND coalesce(ResultStatusDSC,'Final') IN ('Final','Corrected')
  AND COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) IS NOT NULL
  AND BDSPPatientID IN (SELECT DISTINCT pid FROM pairs)
""")
rng = " OR ".join(f"(var='{k}' AND v BETWEEN {a} AND {b})" for k, (a, b) in RANGES.items())
con.execute(f"DELETE FROM lab WHERE NOT ({rng}) OR v IS NULL")
near = con.execute("""
SELECT fname, var, v FROM (
  SELECT e.fname, l.var, l.v,
         row_number() OVER (PARTITION BY e.fname, l.var ORDER BY abs(epoch(e.t)-epoch(l.t)), l.t) rn
  FROM ecg e JOIN lab l ON e.pid = l.pid AND l.t BETWEEN e.t - INTERVAL 12 HOUR AND e.t + INTERVAL 12 HOUR
) WHERE rn = 1""").df().pivot(index="fname", columns="var", values="v")
print("labs matched", near.notna().mean().round(3).to_dict(), flush=True)

# --- 12SL paced / AF flags ----------------------------------------------------------------
dx = con.execute(f"""
SELECT regexp_replace(regexp_replace(regexp_replace(FileName, '\\s', '', 'g'), '^\\.', ''), '\\.hea$', '') AS fname,
       list_has_any(list_transform(string_split(regexp_replace(codes, '\\s', '', 'g'), ','), x -> TRY_CAST(x AS INTEGER)), {PACED}) paced,
       list_has_any(list_transform(string_split(regexp_replace(codes, '\\s', '', 'g'), ','), x -> TRY_CAST(x AS INTEGER)), {AF}) af
FROM read_csv('{HEEDB}/ECG/I0001/12SL_diagnoses/diagnoses_v24.csv', header=true, all_varchar=true)
WHERE regexp_replace(regexp_replace(regexp_replace(FileName, '\\s', '', 'g'), '^\\.', ''), '\\.hea$', '') IN (SELECT fname FROM ecg)
""").df().drop_duplicates("fname").set_index("fname")

# --- bleeding ICD within +-3 d of the post ECG --------------------------------------------
like = " OR ".join(f"DIAGNOSIS_ICD10_CD LIKE '{c}%'" for c in BLEED_ICD10)
bleed = con.execute(f"""
WITH b AS (SELECT CAST(BDSPPatientID AS INTEGER) pid, CAST(RECORDED_DT AS TIMESTAMP) t
           FROM read_csv('{HEEDB}/ECG/I0001/ICD_codes/icd10_codes.csv', header=true, all_varchar=true) WHERE {like})
SELECT DISTINCT p.pair_id FROM pairs p JOIN b ON p.pid = b.pid
 AND b.t BETWEEN p.post_time - INTERVAL 3 DAY AND p.post_time + INTERVAL 3 DAY""").df()

# --- ECG features -------------------------------------------------------------------------
meas = pd.concat([pd.read_parquet(f, columns=["key"] + FEATS) for f in sorted(glob.glob(str(OUT / "measure_mgh/part_*.parquet")))])
meas = meas.set_index("key")

df = pairs.copy()
for side, col in (("base", "base_fname"), ("post", "post_fname")):
    df = df.join(meas.add_prefix(f"{side}_"), on=col)
    df = df.join(near.add_prefix(f"{side}_"), on=col)
    df = df.join(dx.add_prefix(f"{side}_"), on=col)
df["bleed_icd"] = df.pair_id.isin(set(bleed.pair_id))
df["gap_h"] = (df.post_time - df.base_time).dt.total_seconds() / 3600
df["age"] = df.age_days.astype(float) / 365.25
df["male"] = (df.sex == "Male").astype(float)
df.to_parquet(OUT / "analysis_mgh.parquet", index=False)
print(df.shape, "bleed", df.bleed_icd.mean().round(4), "paced", df.post_paced.mean().round(4), flush=True)
