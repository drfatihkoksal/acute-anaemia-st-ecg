"""concept_tolerance.md v1.0, §3-6: build the MGH analysis table (anaemic eligible pairs + contrast cohort).

Rows: all eligible anaemic pairs (post Hb <= 8, ΔHb >= 1; flags `first` and `random` per patient)
      + first eligible contrast pair per patient (|ΔHb| < 0.5, post Hb > 10).
No outcome analysis here.
"""
from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import polars as pl
import wfdb
from comorbidipy import comorbidity
from comorbidipy.codemaps.weights import weights as CW

HERE = Path(__file__).resolve().parent
HEEDB = Path("/home/fatih/datasets/heedb")
WF = HEEDB / "ECG/I0001/WFDB"
REG_END = pd.Timestamp("2023-04-30")
spec = importlib.util.spec_from_file_location("a06", HERE / "06_analysis.py")
A = importlib.util.module_from_spec(spec); spec.loader.exec_module(A)
LEADS = ["I", "II", "III", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]
EXTRA = {"wbc": ["6690-2", "26464-8"], "plt": ["777-3", "26515-7"], "bili": ["1975-2"], "inr": ["6301-6", "34714-6"],
         "lactate": ["2524-7", "32693-4"], "glucose": ["2345-7"], "hco3": ["2028-9", "1963-8", "1959-6"]}
RANGES = {"wbc": (0.01, 500), "plt": (1, 2000), "bili": (0.05, 60), "inr": (0.5, 15), "lactate": (0.1, 30),
          "glucose": (10, 2000), "hco3": (3, 60)}
BNP = ["33762-6", "30934-4", "83107-3"]
D62_10, D62_9 = ["D62"], ["2851"]
MAJOR_10 = ["K920", "K921", "K922", "K250", "K252", "K254", "K256", "K260", "K262", "K264", "K266", "K270", "K272",
            "K274", "K276", "K280", "K282", "K284", "K286", "I8501", "I8511", "K625", "K226", "T810", "D62"]
MAJOR_9 = ["578", "5310", "5312", "5314", "5316", "5320", "5322", "5324", "5326", "5330", "5332", "5334", "5336",
           "5340", "5342", "5344", "5346", "4560", "45620", "5693", "53082", "9981", "2851"]
FLAGS = {"cad": (r"^I2[0-5]", r"^41[0-4]"), "hf": (r"^I50", r"^428"), "ckd": (r"^N18", r"^585"),
         "dm": (r"^E(0[8-9]|1[0-3])", r"^250")}


def st_flags(d):
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


raw = pd.read_parquet(HERE / "out/analysis_mgh.parquet")
raw["post_key"] = raw.post_fname
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
T["year"] = T.post_time.dt.year
T["hr_rise10"] = ((T.post_hr - T.base_hr) >= 10).astype(int)
print("rows", len(T), "anemic first", int((T.anemic == 1)[T["first"]].sum()), "contrast", int((T.anemic == 0).sum()), flush=True)

con = duckdb.connect()
con.execute(f"SET temp_directory='{HEEDB / 'duckdb_tmp'}'; SET threads=24; SET memory_limit='40GB'")
con.register("p", T[["row", "pid", "post_time"]])

# extra labs at post ECG (+-12 h, nearest), BNP measured (+-48 h)
case = " ".join(f"WHEN LoincTXT IN ({','.join(repr(c) for c in v)}) THEN '{k}'" for k, v in EXTRA.items())
allc = ",".join(repr(c) for v in EXTRA.values() for c in v)
con.execute(f"""CREATE TABLE lab AS SELECT BDSPPatientID pid, CASE {case} END var,
  CAST(COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) AS TIMESTAMP) t,
  TRY_CAST(regexp_replace(ResultTXT,'[<>]','','g') AS DOUBLE) v
  FROM read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet')
  WHERE LoincTXT IN ({allc}) AND coalesce(ResultStatusDSC,'Final') IN ('Final','Corrected')
    AND COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) IS NOT NULL AND BDSPPatientID IN (SELECT pid FROM p)""")
rng = " OR ".join(f"(var='{k}' AND v BETWEEN {a} AND {b})" for k, (a, b) in RANGES.items())
con.execute(f"DELETE FROM lab WHERE v IS NULL OR NOT ({rng})")
near = con.execute("""SELECT row, var, v FROM (SELECT p.row, l.var, l.v,
  row_number() OVER (PARTITION BY p.row, l.var ORDER BY abs(epoch(p.post_time)-epoch(l.t))) rn
  FROM p JOIN lab l ON l.pid=p.pid AND l.t BETWEEN p.post_time - INTERVAL 12 HOUR AND p.post_time + INTERVAL 12 HOUR) WHERE rn=1
""").df().pivot(index="row", columns="var", values="v").add_prefix("post_")
bnp = set(con.execute(f"""SELECT DISTINCT p.row FROM p JOIN read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet') l
  ON l.BDSPPatientID=p.pid AND l.LoincTXT IN ({','.join(repr(c) for c in BNP)})
  AND CAST(COALESCE(l.SpecimenTakenTimeDTS, l.ComponentObservedDTS) AS TIMESTAMP)
      BETWEEN p.post_time - INTERVAL 48 HOUR AND p.post_time + INTERVAL 48 HOUR""").df().row)
T = T.join(near, on="row")
T["lactate_measured"] = T.post_lactate.notna().astype(int)
T["inr_measured"] = T.post_inr.notna().astype(int)
T["bnp_measured"] = T.row.isin(bnp).astype(int)

# troponin [-24 h, +48 h] in URL multiples
trop = con.execute(f"""WITH t AS (SELECT BDSPPatientID pid, CAST(COALESCE(SpecimenTakenTimeDTS, ComponentObservedDTS) AS TIMESTAMP) t,
   TRY_CAST(regexp_replace(ResultTXT,'[<>]','','g') AS DOUBLE)/CASE LoincTXT WHEN '67151-1' THEN 19.0 ELSE 0.01 END v
   FROM read_parquet('{HEEDB}/EHR/I0001/Labs/*.parquet') WHERE LoincTXT IN ('67151-1','6598-7') AND BDSPPatientID IN (SELECT pid FROM p))
  SELECT p.row, count(t.v) trop_n, max(t.v) trop_max, min(t.v) trop_min FROM p JOIN t ON t.pid=p.pid
  AND t.t BETWEEN p.post_time - INTERVAL 24 HOUR AND p.post_time + INTERVAL 48 HOUR WHERE t.v IS NOT NULL GROUP BY 1""").df()
T = T.merge(trop, on="row", how="left")
T["trop_measured"] = (T.trop_n.fillna(0) >= 2).astype(int)
T["mi_injury"] = np.where(T.trop_measured == 1, ((T.trop_max > 1) & (T.trop_max >= 1.2 * T.trop_min)).astype(float), np.nan)

# ICD before index (comorbidities) and bleeding +-3 d
icd = con.execute(f"""
  SELECT CAST(BDSPPatientID AS INTEGER) pid, CAST(RECORDED_DT AS TIMESTAMP) t, replace(DIAGNOSIS_ICD10_CD,'.','') code, 10 ver
  FROM read_csv('{HEEDB}/ECG/I0001/ICD_codes/icd10_codes.csv', header=true, all_varchar=true)
  WHERE CAST(BDSPPatientID AS INTEGER) IN (SELECT pid FROM p)
  UNION ALL
  SELECT CAST(BDSPPatientID AS INTEGER), CAST(RECORDED_DT AS TIMESTAMP), replace(DIAGNOSIS_ICD9_CD,'.',''), 9
  FROM read_csv('{HEEDB}/ECG/I0001/ICD_codes/icd9_codes.csv', header=true, all_varchar=true)
  WHERE CAST(BDSPPatientID AS INTEGER) IN (SELECT pid FROM p)""").df()
icd = icd.merge(T[["row", "pid", "post_time"]], on="pid")
prior = icd[icd.t < icd.post_time.dt.normalize()]
for k, (r10, r9) in FLAGS.items():
    m = ((prior.ver == 10) & prior.code.str.match(r10)) | ((prior.ver == 9) & prior.code.str.match(r9))
    T[k] = T.row.isin(set(prior.row[m])).astype(int)
flags = []
for ver in (9, 10):
    sub = prior[prior.ver == ver][["row", "code"]].drop_duplicates()
    if len(sub):
        r = comorbidity(pl.DataFrame({"id": sub.row.to_numpy(), "code": sub.code.to_numpy()}), id_col="id", code_col="code",
                        score="elixhauser", icd=f"icd{ver}", variant="quan", weighting="van_walraven").to_pandas()
        flags.append(r.drop(columns="comorbidity_score").set_index("id"))
allf = pd.concat(flags).groupby(level=0).max()
w = CW["elixhauser_icd10_quan"]["van_walraven"]
allf["elix_vw"] = sum(allf[c] * wt for c, wt in w.items() if c in allf)
T = T.join(allf[["elix_vw"]], on="row")
T["elix_vw"] = T.elix_vw.fillna(0.0)
win = icd[(icd.t >= icd.post_time - pd.Timedelta(days=3)) & (icd.t <= icd.post_time + pd.Timedelta(days=3))]


def has(codes10, codes9):
    m = ((win.ver == 10) & win.code.str.startswith(tuple(codes10))) | ((win.ver == 9) & win.code.str.startswith(tuple(codes9)))
    return T.row.isin(set(win.row[m])).astype(int)


T["bleed_d62"] = has(D62_10, D62_9)
T["bleed_major"] = has(MAJOR_10, MAJOR_9)

# death (MA registry)
death = con.execute(f"""SELECT CAST(BDSPPatientID AS INTEGER) pid, max(CAST(DateOfDeathMARegistryData AS TIMESTAMP)) dod
  FROM read_csv('{HEEDB}/ECG/I0001/metadata/metadata.csv', header=true, all_varchar=true)
  WHERE CAST(BDSPPatientID AS INTEGER) IN (SELECT pid FROM p) GROUP BY 1""").df()
T = T.merge(death, on="pid", how="left")
dd = (T.dod.dt.normalize() - T.post_time.dt.normalize()).dt.days
T["fu30_complete"] = (T.post_time + pd.Timedelta(days=30) <= REG_END).astype(int)
T["died30"] = dd.between(0, 30).astype(int)
end = np.minimum((REG_END - T.post_time.dt.normalize()).dt.days, 365)
T["time365"] = np.where(dd.between(0, 365) & (dd <= end), dd, end).clip(min=0.5)
T["event365"] = (dd.between(0, 365) & (dd <= end)).astype(int)
T["early_death_before_index"] = (dd < 0).astype(int)

# sampling frequency of both ECGs
def fs(k):
    try:
        return wfdb.rdheader(str(WF) + k).fs
    except Exception:  # noqa: BLE001
        return np.nan
T["base_fs"] = [fs(k) for k in T.base_fname]
T["post_fs"] = [fs(k) for k in T.post_fname]

T.to_parquet(HERE / "out/tolerance_mgh.parquet", index=False)
print("saved", T.shape, flush=True)
