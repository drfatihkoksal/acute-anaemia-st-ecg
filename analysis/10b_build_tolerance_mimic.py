"""concept_tolerance.md v1.0 — MIMIC (BIDMC) analysis table, same schema as 10_build_tolerance_mgh.py.

MIMIC-specific substitutions (logged in concept_tolerance.md §10):
  - labs: itemid equivalents of the MGH LOINC set; troponin T 51003 (URL 0.01 ng/mL)
  - comorbidities: ICD codes of PRIOR admissions only (same-admission codes overlap the outcome)
  - bleeding: D62 / major-bleeding codes in the SAME admission (MIMIC diagnoses are undated)
  - year: dates are shifted per patient -> approximate calendar year = midpoint(anchor_year_group) + (ECG year - anchor_year)
  - death: patients.dod; 30-d follow-up treated as complete; 1-y outcome subject to MIMIC's 1-y-after-discharge censoring
"""
from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import polars as pl
from comorbidipy import comorbidity
from comorbidipy.codemaps.weights import weights as CW

HERE = Path(__file__).resolve().parent
HOSP = Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimiciv/3.1/hosp")
spec = importlib.util.spec_from_file_location("a06", HERE / "06_analysis.py")
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
spec10 = importlib.util.spec_from_file_location("b10", HERE / "10_build_tolerance_mgh.py")
EXTRA = {"wbc": [51301], "plt": [51265], "bili": [50885], "inr": [51237], "lactate": [50813],
         "glucose": [50931, 50809], "hco3": [50882, 50803]}
RANGES = {"wbc": (0.01, 500), "plt": (1, 2000), "bili": (0.05, 60), "inr": (0.5, 15), "lactate": (0.1, 30),
          "glucose": (10, 2000), "hco3": (3, 60)}
BNP = [50963]
LEADS = ["I", "II", "III", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]
D62_10, D62_9 = ["D62"], ["2851"]
MAJOR_10 = ["K920", "K921", "K922", "K250", "K252", "K254", "K256", "K260", "K262", "K264", "K266", "K270", "K272",
            "K274", "K276", "K280", "K282", "K284", "K286", "I8501", "I8511", "K625", "K226", "T810", "D62"]
MAJOR_9 = ["578", "5310", "5312", "5314", "5316", "5320", "5322", "5324", "5326", "5330", "5332", "5334", "5336",
           "5340", "5342", "5344", "5346", "4560", "45620", "5693", "53082", "9981", "2851"]
FLAGS = {"cad": (r"^I2[0-5]", r"^41[0-4]"), "hf": (r"^I50", r"^428"), "ckd": (r"^N18", r"^585"),
         "dm": (r"^E(0[8-9]|1[0-3])", r"^250")}


def st_flags(d):  # identical definition to the MGH builder
    def hits(thr):
        return sum(((d[f"post_st60_{l}"] <= -thr) & ((d[f"post_st60_{l}"] - d[f"base_st60_{l}"]) <= -thr)).astype(int)
                   for l in LEADS)
    d["new_std"] = (hits(0.05) >= 2).astype(int)
    d["new_std_1mm"] = (hits(0.10) >= 2).astype(int)
    burden = 0
    for l in LEADS:
        dl = d[f"post_st60_{l}"] - d[f"base_st60_{l}"]
        q = (d[f"post_st60_{l}"] <= -0.05) & (dl <= -0.05)
        burden = burden + np.where(q, -dl, 0.0)
    d["std_burden"] = burden / 0.1
    d["prior_std"] = (sum((d[f"base_st60_{l}"] <= -0.05).astype(int) for l in LEADS) >= 2).astype(int)
    return d


raw = pd.read_parquet(HERE / "out/analysis_mimic.parquet")
df, flow = A.prep(raw)
df = df[df.age >= 18]
an = df[(df.post_hb <= 8) & (df.delta_hb >= 1)].copy()
an["anemic"] = 1
an = an.sort_values(["pid", "post_time"])
an["first"] = ~an.pid.duplicated()
an["rh"] = [hashlib.md5(f"tol:{p}:{k}".encode()).hexdigest() for p, k in zip(an.pid, an.post_key)]
an["random"] = an.rh == an.groupby("pid").rh.transform("min")
ct = df[(df.delta_hb.abs() < 0.5) & (df.post_hb > 10)].sort_values(["pid", "post_time"])
ct = ct[~ct.pid.duplicated()].copy()
ct["anemic"], ct["first"], ct["random"] = 0, True, False
T = pd.concat([an, ct], ignore_index=True)
T["row"] = range(len(T))
T = st_flags(T)
T["hr_rise10"] = ((T.post_hr - T.base_hr) >= 10).astype(int)

pat = pd.read_csv(HOSP / "patients.csv.gz", usecols=["subject_id", "anchor_year_group", "dod"])
pat["ayg_mid"] = pat.anchor_year_group.str[:4].astype(int) + 1
T = T.merge(pat[["subject_id", "ayg_mid", "dod"]], left_on="pid", right_on="subject_id", how="left", suffixes=("", "_p"))
T["year"] = T.ayg_mid + (T.post_time.dt.year - T.anchor_year)

con = duckdb.connect(); con.execute("SET threads=24; SET memory_limit='40GB'")
con.register("p", T[["row", "pid", "hadm_id", "post_time"]])
case = " ".join(f"WHEN itemid IN ({','.join(map(str, v))}) THEN '{k}'" for k, v in EXTRA.items())
allc = ",".join(str(c) for v in EXTRA.values() for c in v)
con.execute(f"""CREATE TABLE lab AS SELECT subject_id pid, itemid, CASE {case} END var, charttime t, valuenum v
  FROM read_csv_auto('{HOSP / "labevents.csv.gz"}', types={{'valuenum':'DOUBLE','value':'VARCHAR','comments':'VARCHAR'}})
  WHERE (itemid IN ({allc}) OR itemid IN ({','.join(map(str, BNP))}) OR itemid = 51003)
    AND valuenum IS NOT NULL AND subject_id IN (SELECT DISTINCT pid FROM p)""")
rng = " OR ".join(f"(var='{k}' AND v BETWEEN {a} AND {b})" for k, (a, b) in RANGES.items())
near = con.execute(f"""SELECT row, var, v FROM (SELECT p.row, l.var, l.v,
  row_number() OVER (PARTITION BY p.row, l.var ORDER BY abs(epoch(p.post_time)-epoch(l.t))) rn
  FROM p JOIN (SELECT * FROM lab WHERE {rng}) l ON l.pid=p.pid
  AND l.t BETWEEN p.post_time - INTERVAL 12 HOUR AND p.post_time + INTERVAL 12 HOUR) WHERE rn=1""").df()
near = near.pivot(index="row", columns="var", values="v").add_prefix("post_")
bnp = set(con.execute(f"""SELECT DISTINCT p.row FROM p JOIN lab l ON l.pid=p.pid AND l.itemid IN ({','.join(map(str, BNP))})
  AND l.t BETWEEN p.post_time - INTERVAL 48 HOUR AND p.post_time + INTERVAL 48 HOUR""").df().row)
T = T.join(near, on="row")
T["lactate_measured"] = T.post_lactate.notna().astype(int)
T["inr_measured"] = T.post_inr.notna().astype(int)
T["bnp_measured"] = T.row.isin(bnp).astype(int)
trop = con.execute("""SELECT p.row, count(l.v) trop_n, max(l.v/0.01) trop_max, min(l.v/0.01) trop_min FROM p JOIN lab l
  ON l.pid=p.pid AND l.itemid=51003 AND l.t BETWEEN p.post_time - INTERVAL 24 HOUR AND p.post_time + INTERVAL 48 HOUR GROUP BY 1""").df()
T = T.merge(trop, on="row", how="left")
T["trop_measured"] = (T.trop_n.fillna(0) >= 2).astype(int)
T["mi_injury"] = np.where(T.trop_measured == 1, ((T.trop_max > 1) & (T.trop_max >= 1.2 * T.trop_min)).astype(float), np.nan)

# comorbidities from PRIOR admissions; bleeding from the SAME admission
adm = pd.read_csv(HOSP / "admissions.csv.gz", usecols=["subject_id", "hadm_id", "admittime"], parse_dates=["admittime"])
dx = pd.read_csv(HOSP / "diagnoses_icd.csv.gz", usecols=["subject_id", "hadm_id", "icd_code", "icd_version"], dtype={"icd_code": str})
dx = dx.merge(adm, on=["subject_id", "hadm_id"])
cur = T[["row", "pid", "hadm_id"]].merge(adm.rename(columns={"subject_id": "pid", "admittime": "cur_adm"}), on=["pid", "hadm_id"], how="left")
x = cur.merge(dx.rename(columns={"subject_id": "pid", "hadm_id": "dx_hadm"}), on="pid")
prior = x[x.admittime < x.cur_adm]
for k, (r10, r9) in FLAGS.items():
    m = ((prior.icd_version == 10) & prior.icd_code.str.match(r10)) | ((prior.icd_version == 9) & prior.icd_code.str.match(r9))
    T[k] = T.row.isin(set(prior.row[m])).astype(int)
flags = []
for ver in (9, 10):
    sub = prior[prior.icd_version == ver][["row", "icd_code"]].drop_duplicates()
    if len(sub):
        r = comorbidity(pl.DataFrame({"id": sub.row.to_numpy(), "code": sub.icd_code.to_numpy()}), id_col="id", code_col="code",
                        score="elixhauser", icd=f"icd{ver}", variant="quan", weighting="van_walraven").to_pandas()
        flags.append(r.drop(columns="comorbidity_score").set_index("id"))
allf = pd.concat(flags).groupby(level=0).max()
w = CW["elixhauser_icd10_quan"]["van_walraven"]
allf["elix_vw"] = sum(allf[c] * wt for c, wt in w.items() if c in allf)
T = T.join(allf[["elix_vw"]], on="row"); T["elix_vw"] = T.elix_vw.fillna(0.0)
same = x[x.dx_hadm == x.hadm_id]


def has(c10, c9):
    m = ((same.icd_version == 10) & same.icd_code.str.startswith(tuple(c10))) | ((same.icd_version == 9) & same.icd_code.str.startswith(tuple(c9)))
    return T.row.isin(set(same.row[m])).astype(int)


T["bleed_d62"] = has(D62_10, D62_9)
T["bleed_major"] = has(MAJOR_10, MAJOR_9)

dod = pd.to_datetime(T.dod)
dd = (dod.dt.normalize() - T.post_time.dt.normalize()).dt.days
T["fu30_complete"] = 1
T["died30"] = dd.between(0, 30).astype(int)
T["time365"] = np.where(dd.between(0, 365), dd, 365).clip(min=0.5)
T["event365"] = dd.between(0, 365).astype(int)
T["early_death_before_index"] = (dd < 0).astype(int)
T["base_fs"], T["post_fs"] = 500, 500
T.drop(columns=["dod"]).to_parquet(HERE / "out/tolerance_mimic.parquet", index=False)
print("saved", T.shape, "anemic first", int(((T.anemic == 1) & T["first"]).sum()), "contrast", int((T.anemic == 0).sum()), flush=True)
