"""Build an outcome-agnostic longitudinal ECG/Hb pair cohort.

Unlike the original acute-drop cohort, eligibility does not depend on a future
haemoglobin trajectory.  Every ECG recorded during an admission is matched to
the closest Hb measurement within +/- 12 hours.  Consecutive, distinct ECGs
with distinct Hb measurements and a positive gap of at most seven days form a
pair.  Splits are subject-disjoint.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from ecganemia import config as C
from ecganemia.utils import get_logger, timer

log = get_logger("aligned_longitudinal_manifest")


def split_of(subject_id: int, salt: str = "aligned-v1") -> str:
    value = int(hashlib.md5(f"{salt}:{subject_id}".encode()).hexdigest(), 16) % 100
    return "train" if value < 70 else "val" if value < 80 else "test"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hb-window-hours", type=float, default=12.0)
    ap.add_argument("--max-pair-gap-days", type=float, default=7.0)
    ap.add_argument("--out-dir", type=Path, default=C.ROOT / "manifests_aligned")
    args = ap.parse_args()

    idx = pd.read_parquet(
        C.INDEX_PARQUET,
        columns=["subject_id", "study_id", "ecg_time", "path", "hr_bpm",
                 "qrs_dur", "pr_int", "rr_interval"],
    )
    idx["ecg_time"] = pd.to_datetime(idx["ecg_time"])
    con = duckdb.connect()
    con.register("ecg", idx)

    # Limit the 18-GB laboratory scan to admissions that contain at least two ECGs.
    sql = f"""
    WITH ecg_adm AS (
      SELECT e.*, a.hadm_id
      FROM ecg e
      JOIN read_csv_auto('{C.hosp("admissions.csv")}') a
        ON e.subject_id = a.subject_id
       AND e.ecg_time BETWEEN a.admittime AND a.dischtime
    ), eligible_adm AS (
      SELECT hadm_id FROM ecg_adm GROUP BY hadm_id HAVING count(DISTINCT study_id) >= 2
    ), hb AS (
      SELECT l.hadm_id, l.charttime AS hb_time, l.valuenum AS hb
      FROM read_csv_auto('{C.hosp("labevents.csv")}',
                         types={{'valuenum': 'DOUBLE'}}) l
      JOIN eligible_adm a USING (hadm_id)
      WHERE l.itemid IN ({','.join(map(str, C.HGB_LAB))})
        AND l.valuenum BETWEEN 2 AND 25
    ), candidates AS (
      SELECT e.*, h.hb_time, h.hb,
             abs(epoch(e.ecg_time)-epoch(h.hb_time))/3600.0 AS hb_gap_hours,
             row_number() OVER (
               PARTITION BY e.study_id
               ORDER BY abs(epoch(e.ecg_time)-epoch(h.hb_time)), h.hb_time
             ) AS rn
      FROM ecg_adm e
      JOIN eligible_adm a USING (hadm_id)
      JOIN hb h ON e.hadm_id = h.hadm_id
               AND h.hb_time BETWEEN e.ecg_time - INTERVAL '{args.hb_window_hours} hours'
                                  AND e.ecg_time + INTERVAL '{args.hb_window_hours} hours'
    )
    SELECT * EXCLUDE (rn) FROM candidates WHERE rn=1
    """
    with timer(log, "ECG-admission-Hb temporal matching"):
        matched = con.sql(sql).df()
    log.info("Hb-aligned ECGs: %d ECG / %d admissions / %d subjects",
             len(matched), matched.hadm_id.nunique(), matched.subject_id.nunique())

    matched = matched.sort_values(["subject_id", "hadm_id", "ecg_time", "study_id"])
    rows: list[dict] = []
    for (_subject, _hadm), g in matched.groupby(["subject_id", "hadm_id"], sort=False):
        g = g.drop_duplicates("study_id").reset_index(drop=True)
        if len(g) < 2:
            continue
        for i in range(len(g) - 1):
            b, p = g.iloc[i], g.iloc[i + 1]
            gap_days = (p.ecg_time - b.ecg_time).total_seconds() / 86400.0
            if not (0 < gap_days <= args.max_pair_gap_days):
                continue
            # A single laboratory measurement assigned to both ECGs supplies no
            # observed longitudinal Hb change and is therefore not a valid target.
            if pd.Timestamp(b.hb_time) == pd.Timestamp(p.hb_time):
                continue
            rows.append({
                "subject_id": int(b.subject_id), "hadm_id": int(b.hadm_id),
                "baseline_study_id": int(b.study_id), "post_study_id": int(p.study_id),
                "baseline_ecg_time": b.ecg_time, "post_ecg_time": p.ecg_time,
                "baseline_path": b.path, "post_path": p.path,
                "baseline_hb_time": b.hb_time, "post_hb_time": p.hb_time,
                "baseline_hb": float(b.hb), "post_hb": float(p.hb),
                "delta_hb": float(b.hb - p.hb),
                "baseline_hb_gap_hours": float(b.hb_gap_hours),
                "post_hb_gap_hours": float(p.hb_gap_hours),
                "pair_gap_days": float(gap_days),
                "baseline_hr": float(b.hr_bpm) if pd.notna(b.hr_bpm) else np.nan,
                "post_hr": float(p.hr_bpm) if pd.notna(p.hr_bpm) else np.nan,
                "delta_hr": float(p.hr_bpm - b.hr_bpm)
                if pd.notna(b.hr_bpm) and pd.notna(p.hr_bpm) else np.nan,
                "baseline_qrs": float(b.qrs_dur) if pd.notna(b.qrs_dur) else np.nan,
                "post_qrs": float(p.qrs_dur) if pd.notna(p.qrs_dur) else np.nan,
                "delta_qrs": float(p.qrs_dur - b.qrs_dur)
                if pd.notna(b.qrs_dur) and pd.notna(p.qrs_dur) else np.nan,
            })

    pairs = pd.DataFrame(rows)
    if pairs.empty:
        raise RuntimeError("No eligible pairs were produced")
    pairs["split"] = pairs.subject_id.map(split_of)
    assert (pairs.baseline_study_id != pairs.post_study_id).all()
    assert (pairs.post_ecg_time > pairs.baseline_ecg_time).all()
    assert not (pairs.groupby("subject_id").split.nunique() > 1).any()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for split in ("train", "val", "test"):
        part = pairs[pairs.split == split].reset_index(drop=True)
        part.to_parquet(args.out_dir / f"pairs_{split}.parquet", index=False)
    pairs.to_parquet(args.out_dir / "pairs_all.parquet", index=False)

    summary = pairs.groupby("split").agg(
        pairs=("study_id" if "study_id" in pairs else "baseline_study_id", "size"),
        subjects=("subject_id", "nunique"), admissions=("hadm_id", "nunique"),
        mean_delta_hb=("delta_hb", "mean"), sd_delta_hb=("delta_hb", "std"),
    )
    log.info("\n%s", summary.to_string())
    log.info("delta Hb distribution:\n%s", pairs.delta_hb.describe(percentiles=[.01,.1,.25,.5,.75,.9,.99]))
    log.info("wrote %s", args.out_dir)


if __name__ == "__main__":
    main()
