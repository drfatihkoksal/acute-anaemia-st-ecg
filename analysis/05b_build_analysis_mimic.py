"""Build the MIMIC (BIDMC) pair-level analysis table with the same schema as 05_build_analysis_mgh.py.

Cohort: manifests_aligned/pairs_all.parquet (scripts/30, original definition incl. same admission).
MIMIC-specific substitutions (logged in concept.md §11.9):
  - bleeding ICD: MIMIC diagnoses are admission-level (no date) -> "bleeding code in the same admission"
  - paced flag: machine-report text (MIMIC-IV-ECG machine_measurements) instead of 12SL codes
  - labs: MIMIC labevents itemids equivalent to the MGH LOINC set, nearest within +-12 h
  - age: anchor_age + (ECG year - anchor_year)
"""
from __future__ import annotations

import glob
from pathlib import Path

import duckdb
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
HOSP = Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimiciv/3.1/hosp")
ECGDIR = Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimic-iv-ecg/1.0")

LABS = {"k": [50971, 52610, 50822], "ca": [50893], "mg": [50960], "na": [50983, 52623, 50824],
        "bun": [51006, 52647], "alb": [50862, 53085], "cr": [50912, 52546]}
RANGES = {"k": (1.5, 9), "ca": (4, 16), "mg": (0.5, 6), "na": (100, 180), "bun": (1, 250), "alb": (0.5, 7), "cr": (0.1, 25)}
BLEED_ICD10 = ["K920", "K921", "K922", "K250", "K252", "K254", "K256", "K260", "K262", "K264", "K266",
               "K27", "K28", "K625", "R31", "N93", "I8501", "D62", "K226", "K228", "R58", "T810"]
BLEED_ICD9 = ["578", "5307", "5310", "5312", "5314", "5316", "5320", "5322", "5324", "5326", "5340", "5342",
              "5344", "5346", "5693", "5997", "7847", "8531", "9218", "2851", "4590", "9981"]
PACED_RE = r"pac(ed|ing|emaker)"
NOT_PACED_RE = r"wandering|junctional pacemaker|ectopic atrial pacemaker"

con = duckdb.connect()
con.execute("SET threads=24; SET memory_limit='40GB'")
p = pd.read_parquet(HERE.parent / "manifests_aligned/pairs_all.parquet")
p = p.rename(columns={"subject_id": "pid", "baseline_path": "base_fname", "post_path": "post_fname",
                      "baseline_ecg_time": "base_time", "post_ecg_time": "post_time", "baseline_hb": "base_hb"})
p = p.drop(columns=["baseline_hr", "post_hr", "delta_hr", "baseline_qrs", "post_qrs", "delta_qrs"])  # machine values; ours used
p["pair_id"] = range(len(p))
con.register("pairs", p)

ecg = pd.concat([p[["pid", "base_fname", "base_time"]].set_axis(["pid", "fname", "t"], axis=1),
                 p[["pid", "post_fname", "post_time"]].set_axis(["pid", "fname", "t"], axis=1)]).drop_duplicates("fname")
con.register("ecg", ecg)
case = " ".join(f"WHEN itemid IN ({','.join(map(str, v))}) THEN '{k}'" for k, v in LABS.items())
allc = ",".join(str(c) for v in LABS.values() for c in v)
con.execute(f"""
CREATE TABLE lab AS
SELECT subject_id pid, CASE {case} END var, charttime t, valuenum v
FROM read_csv_auto('{HOSP / "labevents.csv.gz"}', types={{'valuenum': 'DOUBLE', 'value': 'VARCHAR', 'comments': 'VARCHAR'}})
WHERE itemid IN ({allc}) AND valuenum IS NOT NULL AND subject_id IN (SELECT DISTINCT pid FROM pairs)
""")
rng = " OR ".join(f"(var='{k}' AND v BETWEEN {a} AND {b})" for k, (a, b) in RANGES.items())
con.execute(f"DELETE FROM lab WHERE NOT ({rng})")
near = con.execute("""
SELECT fname, var, v FROM (
  SELECT e.fname, l.var, l.v,
         row_number() OVER (PARTITION BY e.fname, l.var ORDER BY abs(epoch(e.t)-epoch(l.t)), l.t) rn
  FROM ecg e JOIN lab l ON e.pid = l.pid AND l.t BETWEEN e.t - INTERVAL 12 HOUR AND e.t + INTERVAL 12 HOUR
) WHERE rn = 1""").df().pivot(index="fname", columns="var", values="v")
print("labs matched", near.notna().mean().round(3).to_dict(), flush=True)

# paced flag from machine reports
rep = [f"report_{i}" for i in range(18)]
mm = pd.read_csv(ECGDIR / "machine_measurements.csv", usecols=["study_id"] + rep, dtype=str)
txt = mm[rep].fillna("").agg(" | ".join, axis=1)
mm["paced"] = txt.str.contains(PACED_RE, case=False, regex=True) & ~txt.str.contains(NOT_PACED_RE, case=False, regex=True)
paced = dict(zip(mm.study_id.astype(int), mm.paced))
rl = pd.read_csv(ECGDIR / "record_list.csv", usecols=["study_id", "path"])
path2study = dict(zip(rl.path, rl.study_id))

# bleeding ICD in the same admission
like10 = " OR ".join(f"icd_code LIKE '{c}%'" for c in BLEED_ICD10)
like9 = " OR ".join(f"icd_code LIKE '{c}%'" for c in BLEED_ICD9)
bleed = set(con.execute(f"""SELECT DISTINCT hadm_id FROM read_csv_auto('{HOSP / "diagnoses_icd.csv.gz"}', types={{'icd_code':'VARCHAR'}})
  WHERE (icd_version = 10 AND ({like10})) OR (icd_version = 9 AND ({like9}))""").df().hadm_id)

pat = pd.read_csv(HOSP / "patients.csv.gz", usecols=["subject_id", "gender", "anchor_age", "anchor_year"])

FEATS = pd.read_parquet(sorted(glob.glob(str(OUT / "measure_mgh/part_*.parquet")))[0]).columns  # same schema
meas = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(str(OUT / "measure_mimic/part_*.parquet")))]).set_index("key")
meas = meas.drop(columns=["read_err"])

df = p.copy()
for side, col in (("base", "base_fname"), ("post", "post_fname")):
    df = df.join(meas.add_prefix(f"{side}_"), on=col)
    df = df.join(near.add_prefix(f"{side}_"), on=col)
    df[f"{side}_paced"] = df[col].map(path2study).map(paced).fillna(False)
df["bleed_icd"] = df.hadm_id.isin(bleed)
df["gap_h"] = (df.post_time - df.base_time).dt.total_seconds() / 3600
df = df.merge(pat, left_on="pid", right_on="subject_id", how="left")
df["age"] = df.anchor_age + (df.post_time.dt.year - df.anchor_year)
df["male"] = (df.gender == "M").astype(float)
df["post_key"] = df.post_fname
df.to_parquet(OUT / "analysis_mimic.parquet", index=False)
print(df.shape, "bleed", round(df.bleed_icd.mean(), 4), "paced", round(df.post_paced.mean(), 4),
      "age", df.age.describe()[["mean", "min", "max"]].round(1).to_dict(), flush=True)
