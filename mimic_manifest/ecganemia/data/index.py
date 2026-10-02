"""Master EKG indeksi: record_list.csv ⨝ machine_measurements.csv.

identity/ecgid/data/index.py'den uyarlandı (değişmeden reusable) — kayıt başına bir
satır, dalga formu yolu, rapor metni, interval ölçümleri, kalp hızı proxy'si (HR
karıştırıcı kontrolü için, bkz. concept.md §6.3) ve aynı-hasta-içi göreli zaman
(rel_days) türetilir. Bu proje ayrıca klinik/lab kohortuyla subject_id üzerinden
join edilecek (bkz. cohort.py, Faz 1).
"""
from __future__ import annotations

import pandas as pd

from .. import config as C
from ..utils import get_logger, timer

log = get_logger("index")


def build_index() -> pd.DataFrame:
    with timer(log, "record_list + machine_measurements okuma"):
        rec = pd.read_csv(
            C.RECORD_LIST_CSV,
            dtype={"subject_id": "int64", "study_id": "int64",
                   "file_name": "string", "path": "string"},
            parse_dates=["ecg_time"],
        )
        mm_cols = (["subject_id", "study_id", "cart_id"] + C.REPORT_COLS +
                   ["rr_interval", "p_onset", "p_end", "qrs_onset", "qrs_end",
                    "t_end", "p_axis", "qrs_axis", "t_axis"])
        mm = pd.read_csv(C.MACHINE_MEAS_CSV, usecols=mm_cols,
                         dtype={"subject_id": "int64", "study_id": "int64"})

    with timer(log, "join (subject_id, study_id)"):
        df = rec.merge(mm, on=["subject_id", "study_id"], how="left",
                       validate="one_to_one")

    log.info("birleşik satır: %d (record_list: %d, machine_meas: %d)",
             len(df), len(rec), len(mm))

    # Rapor metni: 18 alanı tek stringe birleştir (ischemia/MI regex etiketleme + LLM
    # çapraz-kontrol buradan besleniyor, bkz. concept.md §6.1, §6.7.6).
    reports = df[C.REPORT_COLS].fillna("").astype(str)
    df["report_text"] = reports.agg(" | ".join, axis=1).str.strip(" |")

    # Interval türetmeleri (ms).
    df["qrs_dur"] = df["qrs_end"] - df["qrs_onset"]
    df["pr_int"] = df["qrs_onset"] - df["p_onset"]

    # Kalp hızı proxy'si: rr_interval (ms) -> bpm. HR karıştırıcı kontrolü (§6.3) burada
    # başlar; vital-signs HR (chartevents) ile Faz 1'de çapraz-doğrulanacak.
    df["hr_bpm"] = 60000.0 / df["rr_interval"].where(df["rr_interval"] > 0)

    # Aynı-hasta-içi göreli zaman (gün). MIMIC mutlak tarihi kaydırır ama hasta-içi
    # göreli zaman korunur -> olay-öncesi/sonrası pencerelemeyi mümkün kılar (§4.1).
    df = df.sort_values(["subject_id", "ecg_time"]).reset_index(drop=True)
    first = df.groupby("subject_id")["ecg_time"].transform("min")
    df["rel_days"] = (df["ecg_time"] - first).dt.total_seconds() / 86400.0
    df["n_records_subject"] = df.groupby("subject_id")["study_id"].transform("count")

    return df


def main() -> None:
    C.ARTIFACTS.mkdir(parents=True, exist_ok=True)
    df = build_index()
    with timer(log, "ecg_index.parquet yazma"):
        df.to_parquet(C.INDEX_PARQUET, index=False)
    log.info("Yazıldı: %s", C.INDEX_PARQUET)
    log.info("hasta sayısı: %d | ≥2 kayıtlı hasta: %d",
             df["subject_id"].nunique(),
             (df.groupby("subject_id").size() >= 2).sum())
    hr = df["hr_bpm"].dropna()
    log.info("hr_bpm — medyan %.0f, %%5-%%95: %.0f-%.0f",
             hr.median(), hr.quantile(0.05), hr.quantile(0.95))


if __name__ == "__main__":
    main()
