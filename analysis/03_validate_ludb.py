"""Tune (odd records) and test (even records) the spatial-velocity QRS onset/offset thresholds
of ecgmeasure_gpu against LUDB manual annotations.

Reference per beat: global QRS onset = earliest '(' over the 12 leads, offset = latest ')',
for the annotated QRS whose peak lies within +-75 ms of our aligned fiducial.
Per record error = median over matched kept beats (ours - reference), in ms.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch
import wfdb

import ecgmeasure_gpu as G

LUDB = Path("/home/fatih/datasets/physionet/ludb/1.0.1/data")
OUT = Path(__file__).resolve().parent / "out"
EXT = ["i", "ii", "iii", "avr", "avl", "avf", "v1", "v2", "v3", "v4", "v5", "v6"]
GRID = [0.01, 0.015, 0.02, 0.03, 0.05]


def qrs_triplets(rec: str, ext: str) -> np.ndarray:
    a = wfdb.rdann(str(LUDB / rec), ext)
    out = []
    for k, s in enumerate(a.symbol):
        if s == "N" and 0 < k < len(a.symbol) - 1 and a.symbol[k - 1] == "(" and a.symbol[k + 1] == ")":
            out.append((a.sample[k - 1], a.sample[k], a.sample[k + 1]))
    return np.array(out).reshape(-1, 3)


recs = sorted({p.stem for p in LUDB.glob("*.hea")}, key=int)
X = np.stack([wfdb.rdrecord(str(LUDB / r)).p_signal.T for r in recs]).astype(np.float32)
refs = {r: [qrs_triplets(r, e) for e in EXT] for r in recs}
x = torch.from_numpy(X).cuda()


def errors(out) -> np.ndarray:
    pos, keep = out["pos"].cpu().numpy(), out["keep"].cpu().numpy()
    on, off = out["on_idx"].cpu().numpy() - G.PRE, out["off_idx"].cpu().numpy() - G.PRE
    res = np.full((len(recs), 2), np.nan)
    for i, r in enumerate(recs):
        if not out["ok"][i]:
            continue
        e_on, e_off = [], []
        for p in pos[i][keep[i]]:
            ons, offs = [], []
            for tr in refs[r]:
                if len(tr):
                    j = np.argmin(np.abs(tr[:, 1] - p))
                    if abs(tr[j, 1] - p) <= 37:
                        ons.append(tr[j, 0]); offs.append(tr[j, 2])
            if len(ons) >= 6:
                e_on.append(p + on[i] - min(ons)); e_off.append(p + off[i] - max(offs))
        if e_on:
            res[i] = [np.median(e_on) * 2, np.median(e_off) * 2]           # samples -> ms
    return res


train = np.array([int(r) % 2 == 1 for r in recs])
grid = {}
for a in GRID:
    out = G.measure(x, a_on=a, a_off=a, on_shift=0, off_shift=0)
    grid[a] = errors(out)


def summ(e):
    e = e[~np.isnan(e)]
    return {"n": int(len(e)), "mean": float(e.mean()), "sd": float(e.std()), "mae": float(np.abs(e).mean()),
            "within_10ms": float((np.abs(e) <= 10).mean())}


# criterion: smallest SD on train (a constant bias is removed separately, also on train only)
best_on = min(GRID, key=lambda a: summ(grid[a][train, 0])["sd"])
best_off = min(GRID, key=lambda a: summ(grid[a][train, 1])["sd"])
e0 = np.column_stack([grid[best_on][:, 0], grid[best_off][:, 1]])
on_shift = -int(round(np.nanmedian(e0[train, 0]) / 2))
off_shift = -int(round(np.nanmedian(e0[train, 1]) / 2))
out = G.measure(x, a_on=best_on, a_off=best_off, on_shift=on_shift, off_shift=off_shift)
e = errors(out)
res = {
    "grid_train_sd": {str(a): [summ(grid[a][train, 0])["sd"], summ(grid[a][train, 1])["sd"]] for a in GRID},
    "a_on": best_on, "a_off": best_off, "on_shift": on_shift, "off_shift": off_shift,
    "ok_rate": float(out["ok"].float().mean()),
    "test_onset": summ(e[~train, 0]), "test_offset": summ(e[~train, 1]),
    "train_onset": summ(e[train, 0]), "train_offset": summ(e[train, 1]),
}
OUT.mkdir(exist_ok=True)
(OUT / "ludb_validation.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
