"""Tek kaynak yapılandırma: yollar, itemid/ICD sabitleri, pencere parametreleri.

concept.md ile birebir hizalı. Önceki JBHI projesinin (`/data/identity/ecgid/config.py`)
yapısal paternini izler, ama bu projeye özgü: MIMIC-IV klinik/lab bağlantısı (Hb,
troponin, transfüzyon, kanama ICD) + within-subject anemi/iskemi pencereleme.
"""
from __future__ import annotations

from pathlib import Path

# --------------------------------------------------------------------------- #
# Yollar
# --------------------------------------------------------------------------- #
ROOT = Path("/data/mimicAnemia")

# MIMIC-IV klinik/lab (v3.1, düz CSV). NOT: pci_spo2_nephropathy/scripts/config.py'deki
# "/data/archive (1)/mimic-iv-3.1" yolu bayat (o dizin artık yok) — burada doğru yol.
MIMIC_CLINICAL_DIR = Path("/data/mimic 4 3 1/mimic-iv-3.1")

def hosp(fname: str) -> Path:
    return MIMIC_CLINICAL_DIR / "hosp" / fname

def icu(fname: str) -> Path:
    return MIMIC_CLINICAL_DIR / "icu" / fname

# MIMIC-IV-ECG (800K WFDB kaydı, 500 Hz, 12-lead) — /data/identity altında barınıyor
# (önceki JBHI projesiyle paylaşılan salt-okunur kaynak; bu proje sadece okur).
MIMIC_ECG_DIR = (
    Path("/data/identity")
    / "mimic-iv-ecg-diagnostic-electrocardiogram-matched-subset-1.0"
    / "mimic-iv-ecg-diagnostic-electrocardiogram-matched-subset-1.0"
)
RECORD_LIST_CSV = MIMIC_ECG_DIR / "record_list.csv"
MACHINE_MEAS_CSV = MIMIC_ECG_DIR / "machine_measurements.csv"
REPORT_COLS = [f"report_{i}" for i in range(18)]

ARTIFACTS = ROOT / "data"
INDEX_PARQUET = ARTIFACTS / "ecg_index.parquet"
COHORT_PARQUET = ARTIFACTS / "cohort.parquet"
LABELS_PARQUET = ARTIFACTS / "labels.parquet"
MANIFEST_DIR = ROOT / "manifests"
CACHE_DIR = ROOT / "cache"          # float16 .npy dalga formu cache
CKPT_DIR = ROOT / "checkpoints"
RESULTS_DIR = ROOT / "outputs"
FIG_DIR = ROOT / "figures"
PAPER_DIR = ROOT / "paper"

# --------------------------------------------------------------------------- #
# Sinyal parametreleri (MIMIC-IV-ECG: 12-lead, 500 Hz, 10 s) — identity/ecgid ile aynı
# --------------------------------------------------------------------------- #
FS = 500
SIG_LEN = 5000
LEADS = ["I", "II", "III", "aVR", "aVF", "aVL", "V1", "V2", "V3", "V4", "V5", "V6"]
N_LEADS = 12
BANDPASS = (0.5, 40.0)

# --------------------------------------------------------------------------- #
# labevents itemid'leri (hosp/d_labitems.csv üzerinden 2026-07-12'de doğrulandı)
# --------------------------------------------------------------------------- #
# Hemoglobin — pci_spo2_nephropathy/scripts/config.py ile aynı (Hematology + Blood Gas)
HGB_LAB = [51222, 50811]
# Troponin — 51002/52642 Troponin I (52642 daha yeni assay), 51003 Troponin T.
# Not: I ve T farklı assay/eşik; Grup A/B ayrımında ikisi ayrı sütun olarak tutulmalı,
# birleştirilmemeli (bkz. concept.md §6.1 azaltma: delta-troponin explicit kriter).
TROPONIN_I_LAB = [51002, 52642]
TROPONIN_T_LAB = [51003]

# --------------------------------------------------------------------------- #
# Transfüzyon proxy'si — icu/inputevents.csv itemid'leri (icu/d_items.csv üzerinden
# 2026-07-12'de doğrulandı). ÖNEMLİ KISIT: bu yalnızca ICU kalışlarını kapsar;
# ward-düzeyi (ICU-dışı) transfüzyonlar pharmacy.csv/emar.csv'de ayrı bir kan bankası
# kaydı olarak GÖRÜNMÜYOR (medication alanında blood-product string'i yok — 2026-07-12
# taramasında doğrulandı, sadece "Plasma-Lyte 148" gibi yanlış-pozitifler çıktı).
# Bu, Faz 1 feasibility'de kohort daralmasına katkı yapabilecek bilinen bir kısıt
# (concept.md §6.4) — birincil kriter Hb-düşüşü olacak, transfüzyon ikincil/doğrulayıcı.
TRANSFUSION_ITEMIDS = {
    "packed_red_cells": [225168, 220996],   # PRBC (225168 güncel, 220996 eski/not-in-use)
    "ffp": [220970],
    "platelets": [225170, 226369, 227071],  # + OR/PACU intraop varyantları
    "cryoprecipitate": [225171, 226371],
}
TRANSFUSION_ITEMIDS_ALL = sorted(
    {i for v in TRANSFUSION_ITEMIDS.values() for i in v}
)

# --------------------------------------------------------------------------- #
# Kanama ICD-9/10 kodları (taslak — Faz 1'de d_icd_diagnoses.csv üzerinden
# genişletilip doğrulanacak). SQL'de LIKE/regex ile prefix eşleşmesi kullanılır.
# --------------------------------------------------------------------------- #
BLEEDING_ICD9_PREFIXES = [
    "578",    # GI hemorrhage
    "5307",   # esophageal hemorrhage
    "5310", "5312", "5314", "5316",  # gastric/duodenal ulcer with hemorrhage
    "5320", "5322", "5324", "5326",  # peptic ulcer with hemorrhage
    "5340", "5342", "5344", "5346",  # gastrojejunal ulcer with hemorrhage
    "5693",  # hemorrhage of rectum/anus
    "5997",  # hematuria
    "7847",  # epistaxis (destek/duyarlılık)
    "8531",  # subdural hemorrhage (travma-ilişkili, ayrı kova)
    "9218",  # traumatic hemorrhage (destek)
]
BLEEDING_ICD10_PREFIXES = [
    "K92.0", "K92.1", "K92.2",  # hematemesis, melena, GI hemorrhage NOS
    "K25.0", "K25.2", "K25.4", "K25.6",  # gastric ulcer with hemorrhage
    "K26.0", "K26.2", "K26.4", "K26.6",  # duodenal ulcer with hemorrhage
    "K27", "K28",  # peptic/gastrojejunal ulcer (subkodlar Faz1'de daraltılacak)
    "K62.5",  # hemorrhage of anus/rectum
    "R31",    # hematuria
    "N93",    # abnormal uterine/vaginal bleeding
    "I85.01", # esophageal varices with bleeding
]

# --------------------------------------------------------------------------- #
# Anemia ICD (yalnızca referans/duyarlılık — birincil kriter labevents Hb düşüşü)
# pci_spo2_nephropathy/scripts/03_covariates.py'deki tanımla aynı.
# --------------------------------------------------------------------------- #
ANEMIA_ICD9_REGEX = r"^28[0-5]"
ANEMIA_ICD10_REGEX = r"^D5|^D6[0-4]"

# --------------------------------------------------------------------------- #
# Pencereleme parametreleri (concept.md §5.1, §6.2 — Faz 2'de kesinleşecek taslak)
# --------------------------------------------------------------------------- #
ACUTE_HB_DROP_G_DL = 2.0
ACUTE_HB_DROP_WINDOW_H = 48
BASELINE_WINDOW_DAYS = 7      # olay öncesi 0-7 gün
POST_WINDOW_H = 24            # nadir Hb'ye en yakın ±24s

# --------------------------------------------------------------------------- #
# Eğitim / model hyperparametreleri (Faz 4+ — identity/ecgid ile tutarlı varsayılanlar)
# --------------------------------------------------------------------------- #
SEED = 1234
EMBED_DIM = 256
BATCH_SIZE = 256
NUM_WORKERS = 16
LR = 1e-3
WEIGHT_DECAY = 1e-4
EPOCHS = 40
