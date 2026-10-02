"""Measure every ECG in the pair cohorts with ecgmeasure_gpu (WFDB reading in CPU workers, measurement on GPU).

  python3 04_measure.py --cohort mgh   [--limit N]
  python3 04_measure.py --cohort mimic [--limit N]

Writes out/measure_<cohort>/part_XXXXX.parquet (one row per ECG, keyed by `key`); finished parts are
skipped on restart. Records that cannot be read/mapped get ok=False with a `read_err` reason.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import wfdb
from scipy.signal import resample_poly
from torch.utils.data import DataLoader, Dataset

import ecgmeasure_gpu as G

HERE = Path(__file__).resolve().parent
ROOTS = {
    "mgh": Path("/home/fatih/datasets/heedb/ECG/I0001/WFDB"),
    "mimic": Path("/media/fatih/48bd54fc-58ee-46af-8646-96be97c8d33c/mimic-iv-ecg/1.0"),
}
LEAD_KEY = [l.lower() for l in G.LEADS]


def ecg_list(cohort: str) -> list[str]:
    if cohort == "mgh":
        p = pd.read_parquet(HERE / "out/pairs_mgh.parquet", columns=["base_fname", "post_fname"])
        keys = pd.concat([p.base_fname, p.post_fname]).str.strip()
    else:
        p = pd.read_parquet(HERE.parent / "manifests_aligned/pairs_all.parquet", columns=["baseline_path", "post_path"])
        keys = pd.concat([p.baseline_path, p.post_path])
    return sorted(pd.unique(keys))


class ECGs(Dataset):
    def __init__(self, keys, root):
        self.keys, self.root = keys, root

    def __len__(self):
        return len(self.keys)

    def __getitem__(self, i):
        k = self.keys[i]
        x = np.zeros((12, G.N), np.float32)
        err = ""
        try:
            r = wfdb.rdrecord(str(self.root) + ("/" if not k.startswith("/") else "") + k)
            names = [n.lower() for n in r.sig_name]
            sig = r.p_signal
            if r.fs == 250:
                sig = resample_poly(sig, 2, 1, axis=0)
            elif r.fs != 500:
                raise ValueError(f"fs={r.fs}")
            if sig.shape[0] < G.N:
                raise ValueError(f"len={sig.shape[0]}")
            for li, ln in enumerate(LEAD_KEY):
                x[li] = sig[:G.N, names.index(ln)]
            if not np.isfinite(x).all():
                bad = ~np.isfinite(x)
                if bad.mean() > 0.01:
                    raise ValueError("nan")
                x[bad] = 0
        except Exception as e:  # noqa: BLE001 - record-level failure is data, not a crash
            err = f"{type(e).__name__}: {e}"[:120]
        return k, torch.from_numpy(x), err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", choices=list(ROOTS), required=True)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--part", type=int, default=50_000)
    ap.add_argument("--batch", type=int, default=512)
    ap.add_argument("--workers", type=int, default=24)
    args = ap.parse_args()

    keys = ecg_list(args.cohort)
    if args.limit:
        keys = keys[:: max(1, len(keys) // args.limit)][: args.limit]
    out_dir = HERE / "out" / (f"measure_{args.cohort}" + ("_test" if args.limit else ""))
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"{len(keys)} ECGs -> {out_dir}", flush=True)

    t0 = time.time()
    for pi, s in enumerate(range(0, len(keys), args.part)):
        f = out_dir / f"part_{pi:05d}.parquet"
        if f.exists():
            continue
        dl = DataLoader(ECGs(keys[s:s + args.part], ROOTS[args.cohort]), batch_size=args.batch,
                        num_workers=args.workers, pin_memory=True, persistent_workers=False)
        frames = []
        for k, x, err in dl:
            res = G.to_rows(G.measure(x.cuda(non_blocking=True)))
            df = pd.DataFrame(res)
            df.insert(0, "key", list(k))
            df["read_err"] = list(err)
            df.loc[df.read_err != "", "ok"] = False
            frames.append(df)
        tmp = f.with_suffix(".tmp")
        pd.concat(frames).to_parquet(tmp, index=False)
        tmp.rename(f)
        done = min(s + args.part, len(keys))
        el = time.time() - t0
        print(f"part {pi} done {done}/{len(keys)}  {done / el:.0f} ECG/s (incl. skipped)", flush=True)
    print(f"finished in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
