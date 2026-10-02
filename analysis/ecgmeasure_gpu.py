"""Batched GPU ECG measurement: median beat + Kors VCG (protocol v3, concept.md §11.5).

Pipeline per record (all records of a batch processed together on the GPU):
  1. zero-phase FFT filter 0.5-100 Hz (Butterworth |H|^2, i.e. filtfilt-equivalent) + 60 Hz notch
  2. QRS detection on the 12-lead derivative energy, 250 ms refractory
  3. ectopic rule (RR < 0.8 x median RR -> drop beat and the next), cross-correlation alignment
     (+-20 ms) to a first template, template-correlation filter (r >= CORR_MIN)
  4. median beat (-300..+500 ms around the fiducial)
  5. QRS onset/offset from 12-lead spatial velocity of the median beat (thresholds tuned on
     LUDB, see 03_validate_ludb.py)
  6. features relative to the PR-segment baseline; Kors VCG early/late QRS magnitudes

Input: float32 tensor [B, 12, 5000] in mV, lead order LEADS, 500 Hz.
"""
from __future__ import annotations

import math

import torch
import torch.nn.functional as F

FS = 500
N = 5000
LEADS = ["I", "II", "III", "aVR", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]
PRE, POST = 150, 250            # median-beat window: -300 ms .. +500 ms
W = PRE + POST
K_MAX = 48                      # max beats per 10 s record (HR <= ~280)
CORR_MIN = 0.90
MIN_BEATS = 6
# spatial-velocity thresholds (fraction of the QRS peak velocity); set from LUDB tuning
A_ON = 0.05
A_OFF = 0.05
ON_SHIFT = -7    # samples; constant bias correction calibrated on LUDB odd records (03_validate_ludb.py)
OFF_SHIFT = 7
A_HIGH = 0.30
EARLY_FRAC, LATE_FRAC = 0.45, 0.55

# Kors regression matrix: rows X, Y, Z; columns I, II, V1..V6 (Kors et al., Eur Heart J 1990)
KORS = torch.tensor([
    [0.38, -0.07, -0.13, 0.05, -0.01, 0.14, 0.06, 0.54],
    [-0.07, 0.93, 0.06, -0.02, -0.05, 0.06, -0.17, 0.13],
    [0.11, -0.23, -0.43, -0.06, -0.14, -0.20, -0.11, 0.31],
])
KORS_IDX = [0, 1, 6, 7, 8, 9, 10, 11]


def _filter(x: torch.Tensor) -> torch.Tensor:
    pad = FS
    xp = F.pad(x, (pad, pad), mode="reflect")
    n = xp.shape[-1]
    f = torch.fft.rfftfreq(n, 1.0 / FS, device=x.device)
    hp = (f / 0.5) ** 4 / (1 + (f / 0.5) ** 4)          # 2nd-order Butterworth, squared
    lp = 1 / (1 + (f / 100.0) ** 4)
    d = (f ** 2 - 60.0 ** 2) ** 2
    notch = d / (d + (f * 2.0) ** 2)
    y = torch.fft.irfft(torch.fft.rfft(xp) * (hp * lp * notch), n=n)
    return y[..., pad:pad + x.shape[-1]]


def _gather(sig: torch.Tensor, pos: torch.Tensor, lo: int, hi: int) -> torch.Tensor:
    """sig [B,L,N], pos [B,K] -> windows [B,K,L,hi-lo] (indices clamped)."""
    B, L, n = sig.shape
    K = pos.shape[1]
    off = torch.arange(lo, hi, device=sig.device)
    idx = (pos[:, :, None] + off).clamp(0, n - 1)                       # [B,K,w]
    idx = idx[:, :, None, :].expand(B, K, L, hi - lo)
    return sig[:, None].expand(B, K, L, n).gather(3, idx)


def _corr(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    """Pearson correlation over the last two dims; a [B,K,L,w], b [B,1,L,w]."""
    a = a - a.mean(dim=(-1, -2), keepdim=True)
    b = b - b.mean(dim=(-1, -2), keepdim=True)
    num = (a * b).sum(dim=(-1, -2))
    den = torch.sqrt((a * a).sum(dim=(-1, -2)) * (b * b).sum(dim=(-1, -2))) + 1e-9
    return num / den


def _masked_median_beat(beats: torch.Tensor, keep: torch.Tensor) -> torch.Tensor:
    b = beats.masked_fill(~keep[:, :, None, None], float("nan"))
    return torch.nanmedian(b, dim=1).values                              # [B,L,W]


def _signed_extreme(x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    """x [B,L,W], mask [B,1,W] -> value with the largest |x| inside the mask (nan if empty)."""
    a = x.abs().masked_fill(~mask, -1.0)
    i = a.argmax(dim=-1, keepdim=True)
    v = x.gather(-1, i).squeeze(-1)
    return v.masked_fill(~mask.any(dim=-1), float("nan"))


@torch.no_grad()
def measure(x: torch.Tensor, a_on: float = A_ON, a_off: float = A_OFF,
            on_shift: int | None = None, off_shift: int | None = None) -> dict[str, torch.Tensor]:
    dev = x.device
    B = x.shape[0]
    nan = float("nan")
    sig = _filter(x)

    # --- 2. QRS detection -------------------------------------------------------------------
    d = torch.diff(sig, dim=-1, prepend=sig[..., :1])
    e = (d * d).sum(dim=1, keepdim=True)                                   # [B,1,N]
    e = F.avg_pool1d(e, 21, stride=1, padding=10)[:, 0]                    # 40 ms smoothing
    thr = 0.25 * torch.quantile(e, 0.99, dim=-1, keepdim=True)
    mp = F.max_pool1d(e[:, None], 251, stride=1, padding=125)[:, 0]
    t = torch.arange(N, device=dev)
    peak = (e == mp) & (e > thr) & (t >= PRE) & (t < N - POST)
    key = torch.where(peak, t.float(), torch.full_like(e, float("inf")))
    key, _ = key.sort(dim=-1)
    key = key[:, :K_MAX]
    valid = torch.isfinite(key)
    pos = torch.where(valid, key, torch.zeros_like(key)).long()           # [B,K]

    # RR and ectopic rule (all detected beats)
    rr = (pos[:, 1:] - pos[:, :-1]).float()
    rr_valid = valid[:, 1:] & valid[:, :-1]
    rr_n = rr.masked_fill(~rr_valid, nan)
    rr_med = torch.nanmedian(rr_n, dim=1).values                           # samples
    rr_cv = (torch.sqrt(torch.nanmean((rr_n - torch.nanmean(rr_n, 1, keepdim=True)) ** 2, 1))
             / torch.nanmean(rr_n, 1))
    prem = rr_valid & (rr < 0.8 * rr_med[:, None])                         # beat k+1 premature
    drop = torch.zeros_like(valid)
    drop[:, 1:] |= prem
    drop[:, 2:] |= prem[:, :-1]
    keep = valid & ~drop

    # --- 3. alignment to a first template, correlation filter --------------------------------
    qlo, qhi = -30, 50                                                     # -60..+100 ms
    tmpl0 = _masked_median_beat(_gather(sig, pos, qlo, qhi), keep)[:, None]
    best_c = torch.full(pos.shape, -2.0, device=dev)
    best_s = torch.zeros_like(pos)
    for s in range(-10, 11):
        c = _corr(_gather(sig, pos + s, qlo, qhi), torch.nan_to_num(tmpl0))
        upd = c > best_c
        best_c = torch.where(upd, c, best_c)
        best_s = torch.where(upd, torch.full_like(pos, s), best_s)
    pos = (pos + best_s).clamp(PRE, N - POST - 1)
    beats = _gather(sig, pos, -PRE, POST)                                  # [B,K,L,W]
    tmpl1 = _masked_median_beat(beats[..., PRE + qlo:PRE + qhi], keep)[:, None]
    c = _corr(beats[..., PRE + qlo:PRE + qhi], torch.nan_to_num(tmpl1))
    keep = keep & (c >= CORR_MIN)
    n_beats = keep.sum(1)
    ok = n_beats >= MIN_BEATS
    med = _masked_median_beat(beats, keep)                                 # [B,L,W]
    med = torch.nan_to_num(med)

    # --- 5. QRS onset / offset from spatial velocity ---------------------------------------
    sv = torch.sqrt((torch.diff(med, dim=-1) ** 2).sum(1))                 # [B,W-1]
    sv = F.avg_pool1d(sv[:, None], 5, stride=1, padding=2)[:, 0]
    tw = torch.arange(W - 1, device=dev)
    reg = (tw >= PRE - 60) & (tw <= PRE + 60)
    svmax = sv.masked_fill(~reg, 0).amax(-1, keepdim=True) + 1e-9
    high = (sv > A_HIGH * svmax) & (tw >= PRE - 60) & (tw <= PRE + 75)
    first_high = torch.where(high, tw, W).amin(-1)
    last_high = torch.where(high, tw, -1).amax(-1)
    # flattest 20 ms window before (PR segment) and after (ST segment) the high-velocity core;
    # thresholds are noise-adaptive and the search never crosses these windows
    rm = F.avg_pool1d(sv[:, None], 11, stride=1, padding=5)[:, 0]
    fh, lh = first_high[:, None], last_high[:, None]
    pre_reg = (tw >= PRE - 75) & (tw <= fh - 5)
    post_reg = (tw >= lh + 5) & (tw <= PRE + 110)
    flat_pre = rm.masked_fill(~pre_reg, float("inf")).argmin(-1)
    flat_post = rm.masked_fill(~post_reg, float("inf")).argmin(-1)
    n_pre = rm.gather(-1, flat_pre[:, None])
    n_post = rm.gather(-1, flat_post[:, None])
    thr_on = n_pre + a_on * (svmax - n_pre)
    thr_off = n_post + a_off * (svmax - n_post)
    cand_on = (sv < thr_on) & (tw >= flat_pre[:, None]) & (tw < fh)
    on = torch.where(cand_on, tw, -1).amax(-1) + 1
    on_flag = on <= 0
    on = torch.where(on_flag, flat_pre, on)
    cand_off = (sv < thr_off) & (tw > lh) & (tw <= flat_post[:, None])
    off = torch.where(cand_off, tw, W).amin(-1)
    off_flag = off >= W
    off = torch.where(off_flag, flat_post, off)
    on = (on + (ON_SHIFT if on_shift is None else on_shift)).clamp(1, PRE - 1)
    off = (off + (OFF_SHIFT if off_shift is None else off_shift)).clamp(PRE + 1, W - 41)

    # --- 6. features ----------------------------------------------------------------------
    tt = torch.arange(W, device=dev)[None, None, :]
    on3, off3 = on[:, None, None], off[:, None, None]
    bmask = (tt >= on3 - 20) & (tt <= on3 - 6)                             # -40..-12 ms
    base = torch.nanmedian(med.masked_fill(~bmask, nan), dim=-1).values    # [B,L]
    m = med - base[..., None]
    qmask = (tt >= on3) & (tt <= off3)
    r_amp = m.masked_fill(~qmask, -1e9).amax(-1).clamp(min=0)
    s_min = m.masked_fill(~qmask, 1e9).amin(-1).clamp(max=0)
    qrs_pp = r_amp - s_min
    j_amp = m.gather(-1, off[:, None, None].expand(-1, 12, 1)).squeeze(-1)
    st60 = m.gather(-1, (off + 30).clamp(max=W - 1)[:, None, None].expand(-1, 12, 1)).squeeze(-1)
    t_end = torch.minimum(torch.full_like(off, W - 1), (PRE + 0.7 * rr_med).nan_to_num(W - 1).long())
    tmask = (tt >= off3 + 40) & (tt <= t_end[:, None, None])
    t_amp = _signed_extreme(m, tmask)
    pmask = (tt >= (on3 - 125).clamp(min=0)) & (tt <= on3 - 25)
    p_amp = _signed_extreme(m, pmask)

    xyz = torch.einsum("pk,bkw->bpw", KORS.to(dev), m[:, KORS_IDX])
    mag = torch.sqrt((xyz ** 2).sum(1))                                    # [B,W]
    dur = (off - on).float()
    e_end = on + (EARLY_FRAC * dur).long()
    l_start = on + (LATE_FRAC * dur).long()
    t2 = tt[:, 0]
    m_early = mag.masked_fill(~((t2 >= on[:, None]) & (t2 <= e_end[:, None])), 0).amax(-1)
    m_late = mag.masked_fill(~((t2 >= l_start[:, None]) & (t2 <= off[:, None])), 0).amax(-1)
    m_qrs = mag.masked_fill(~((t2 >= on[:, None]) & (t2 <= off[:, None])), 0).amax(-1)

    out = {
        "ok": ok, "n_beats": n_beats, "n_detected": valid.sum(1), "tmpl_corr": c.masked_fill(~keep, nan).nanmean(1),
        "rr_ms": rr_med * 1000 / FS, "hr": 60 * FS / rr_med, "rr_cv": rr_cv,
        "qrs_on_ms": (on - PRE).float() * 1000 / FS, "qrs_off_ms": (off - PRE).float() * 1000 / FS,
        "qrs_dur_ms": dur * 1000 / FS, "on_flag": on_flag, "off_flag": off_flag,
        "m_early": m_early, "m_late": m_late, "m_qrs": m_qrs,
        "sum_r": r_amp.sum(1), "sum_qrs_pp": qrs_pp.sum(1), "sum_p_abs": p_amp.abs().nansum(1),
        "st60_min": st60.amin(1),
        "r_amp": r_amp, "s_min": s_min, "qrs_pp": qrs_pp, "j_amp": j_amp, "st60": st60,
        "t_amp": t_amp, "p_amp": p_amp,
        "pos": pos, "keep": keep, "on_idx": on, "off_idx": off,
    }
    return out


PER_LEAD = ["r_amp", "s_min", "qrs_pp", "j_amp", "st60", "t_amp", "p_amp"]
SCALARS = ["ok", "n_beats", "n_detected", "tmpl_corr", "rr_ms", "hr", "rr_cv", "qrs_on_ms", "qrs_off_ms",
           "qrs_dur_ms", "on_flag", "off_flag", "m_early", "m_late", "m_qrs", "sum_r", "sum_qrs_pp",
           "sum_p_abs", "st60_min"]


def to_rows(out: dict[str, torch.Tensor]) -> dict[str, list]:
    """Flatten to columns: scalars + per-lead `<feat>_<lead>`; invalid records get NaN features."""
    ok = out["ok"].cpu()
    cols = {}
    for k in SCALARS:
        if k in ("ok", "n_beats", "n_detected", "on_flag", "off_flag"):
            cols[k] = out[k].cpu().numpy()
            continue
        cols[k] = out[k].float().cpu().masked_fill(~ok, math.nan).numpy()
    for k in PER_LEAD:
        v = out[k].float().cpu().masked_fill(~ok[:, None], math.nan).numpy()
        for li, ln in enumerate(LEADS):
            cols[f"{k}_{ln}"] = v[:, li]
    return cols
