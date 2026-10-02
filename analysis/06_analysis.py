"""Pre-registered analysis, protocol v3.1 (concept.md §11.3-11.7). Identical code for both cohorts.

  python3 06_analysis.py --cohort mgh      (discovery)
  python3 06_analysis.py --cohort mimic    (external; only after this file's hash is logged in §11.9)

Input out/analysis_<cohort>.parquet with columns: pid, delta_hb (base-post, + = fall), base_*/post_* ECG
features and labs, gap_h, age, male, bleed_icd, base_paced/post_paced, post_hb.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer

HERE = Path(__file__).resolve().parent
MARGIN = 0.01          # mV per g/dL (H1 equivalence margin)
M_IMP = 5
ADJ = ["d_rr", "d_k", "d_ca", "gap_h", "age", "male"]
VOL = ["d_na", "d_bun", "d_alb"]


def prep(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    flow = {"pairs": len(df)}
    df = df[(df.base_ok == True) & (df.post_ok == True)]  # noqa: E712
    flow["both_measured"] = len(df)
    df = df[(df.post_qrs_dur_ms - df.base_qrs_dur_ms).abs() <= 20]
    flow["qrs_change_le20"] = len(df)
    irr_b, irr_p = df.base_rr_cv > 0.15, df.post_rr_cv > 0.15
    df = df[irr_b == irr_p]
    flow["no_rhythm_change"] = len(df)
    df = df[~(df.base_paced.fillna(False).astype(bool) | df.post_paced.fillna(False).astype(bool))]
    flow["not_paced"] = len(df)
    df = df.copy()
    for f in ["st60_min", "m_early", "m_late", "m_qrs", "sum_r", "sum_qrs_pp", "sum_p_abs", "s_min_V6",
              "rr_ms", "hr"] + [c[5:] for c in df.columns if c.startswith("base_st60_") or c.startswith("base_t_amp_")]:
        df[f"d_{f}"] = df[f"post_{f}"] - df[f"base_{f}"]
    for lab in ["k", "ca", "mg", "na", "bun", "alb", "cr"]:
        df[f"d_{lab}"] = df[f"post_{lab}"] - df[f"base_{lab}"]
    df["d_rr"] = df.d_rr_ms
    # one pair per patient, deterministic
    df["h"] = [hashlib.md5(f"v3:{p}:{k}".encode()).hexdigest() for p, k in zip(df.pid, df.post_key)]
    df["primary"] = df.h == df.groupby("pid").h.transform("min")
    flow["primary_one_per_patient"] = int(df.primary.sum())
    flow["patients"] = int(df.pid.nunique())
    return df, flow


def fit(df, y, x="delta_hb", adj=ADJ, extra=(), cluster=False, impute=True):
    cols = [x] + list(adj) + list(extra) + [f"base_{y[2:]}"] if f"base_{y[2:]}" in df else [x] + list(adj) + list(extra)
    d = df[[y] + cols + ["pid"]].replace([np.inf, -np.inf], np.nan)
    d = d[d[y].notna() & d[x].notna()]
    miss = d[cols].isna().any(axis=1).mean()
    ests, vars_ = [], []
    reps = M_IMP if (impute and miss > 0) else 1
    for m in range(reps):
        X = d[cols].to_numpy(float)
        if reps > 1:
            X = IterativeImputer(sample_posterior=True, random_state=m, max_iter=5).fit_transform(
                np.column_stack([d[y].to_numpy(float), X]))[:, 1:]
        else:
            keep = ~np.isnan(X).any(axis=1)
            X = X[keep]
        Xc = sm.add_constant(X)
        yy = d[y].to_numpy(float) if reps > 1 else d[y].to_numpy(float)[keep]
        mod = sm.OLS(yy, Xc)
        if cluster:
            g = d.pid.to_numpy() if reps > 1 else d.pid.to_numpy()[keep]
            r = mod.fit(cov_type="cluster", cov_kwds={"groups": g})
        else:
            r = mod.fit(cov_type="HC3")
        ests.append(r.params[1]); vars_.append(r.bse[1] ** 2); n = int(r.nobs)
    q = np.mean(ests)
    w = np.mean(vars_)
    b = np.var(ests, ddof=1) if reps > 1 else 0.0
    se = np.sqrt(w + (1 + 1 / reps) * b)
    z = q / se
    return {"beta": float(q), "se": float(se), "p": float(2 * stats.norm.sf(abs(z))),
            "ci95": [float(q - 1.96 * se), float(q + 1.96 * se)], "ci90": [float(q - 1.645 * se), float(q + 1.645 * se)],
            "n": n, "frac_imputed_rows": float(miss)}


def tost(r):
    p_lo = stats.norm.sf((r["beta"] + MARGIN) / r["se"])
    p_hi = stats.norm.cdf((r["beta"] - MARGIN) / r["se"])
    return float(max(p_lo, p_hi))


def holm(ps):
    order = np.argsort(ps)
    adj = np.empty(len(ps))
    run = 0
    for i, k in enumerate(order):
        run = max(run, (len(ps) - i) * ps[k])
        adj[k] = min(1.0, run)
    return adj.tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    args = ap.parse_args()
    raw = pd.read_parquet(HERE / f"out/analysis_{args.cohort}.parquet")
    if "post_key" not in raw:
        raw["post_key"] = raw.post_fname
    df, flow = prep(raw)
    P = df[df.primary]
    R = {"cohort": args.cohort, "flow": flow, "code_sha1": hashlib.sha1(Path(__file__).read_bytes()).hexdigest()}

    # --- primary -------------------------------------------------------------------------
    h1 = fit(P, "d_st60_min"); h1["p_tost"] = tost(h1); h1["equivalent"] = h1["p_tost"] < 0.05
    e = fit(P, "d_m_early"); l = fit(P, "d_m_late")
    adj = holm([h1["p_tost"], e["p"], l["p"]])
    R["H1_st60_min"] = h1
    R["H2a_m_early"] = e
    R["H2a_m_late"] = l
    R["holm_adjusted"] = {"H1_tost": adj[0], "H2a_early": adj[1], "H2a_late": adj[2]}
    R["H2a_supported"] = bool(e["beta"] > 0 and l["beta"] < 0 and adj[1] < 0.05 and adj[2] < 0.05)

    # --- secondary -----------------------------------------------------------------------
    sec = {"H1_no_hr_adjust": fit(P, "d_st60_min", adj=[a for a in ADJ if a != "d_rr"]),
           "H2a_sum_r": fit(P, "d_sum_r"), "H2a_s_min_V6": fit(P, "d_s_min_V6"),
           "H2b_sum_qrs_pp": fit(P, "d_sum_qrs_pp", extra=VOL), "H2b_sum_p_abs": fit(P, "d_sum_p_abs", extra=VOL),
           "H2b_m_qrs": fit(P, "d_m_qrs", extra=VOL),
           "H2a_m_early_vol_adj": fit(P, "d_m_early", extra=VOL), "H2a_m_late_vol_adj": fit(P, "d_m_late", extra=VOL)}
    for v in VOL:
        sec[f"H2b_sum_qrs_pp_on_{v}"] = fit(P, "d_sum_qrs_pp", x=v, extra=["delta_hb"] + [w for w in VOL if w != v])
    for c in [c for c in P.columns if c.startswith("d_st60_") and c != "d_st60_min"] + \
             [c for c in P.columns if c.startswith("d_t_amp_")]:
        sec[c] = fit(P, c, impute=False)
    sec["hr_on_dhb"] = fit(P, "d_hr", adj=["d_k", "d_ca", "gap_h", "age", "male"], impute=False)
    R["secondary"] = sec

    # --- negative controls ---------------------------------------------------------------
    st = P[P.delta_hb.abs() < 0.5]
    R["noise_floor_sd_stable_hb"] = {y: float(st[y].std()) for y in ["d_st60_min", "d_m_early", "d_m_late", "d_sum_r"]}
    R["standardized_effect_per_3gdl"] = {y: float(3 * R[k]["beta"] / R["noise_floor_sd_stable_hb"][y])
                                         for y, k in [("d_st60_min", "H1_st60_min"), ("d_m_early", "H2a_m_early"),
                                                      ("d_m_late", "H2a_m_late")]}
    rng = np.random.default_rng(20261002)
    W = P.copy()
    perm = rng.permutation(len(W))
    for f in ["st60_min", "m_early", "m_late"]:
        W[f"d_{f}"] = W[f"post_{f}"].to_numpy()[perm] - W[f"base_{f}"].to_numpy()
    R["wrong_pair"] = {f: fit(W, f"d_{f}", impute=False) for f in ["st60_min", "m_early", "m_late"]}
    R["creatinine_exposure"] = {f: fit(P, f"d_{f}", x="d_cr", extra=["delta_hb"], impute=False)
                                for f in ["m_early", "m_late", "st60_min"]}

    # --- subgroups & sensitivity -----------------------------------------------------------
    sub = {"bleed_icd": P[P.bleed_icd],
           "dilution_proxy": P[(~P.bleed_icd) & (P.d_alb <= -0.3)],
           "post_hb_le8": P[P.post_hb <= 8]}
    R["subgroups"] = {k: {"n": len(v), **{f: fit(v, f"d_{f}", impute=False) for f in ["m_early", "m_late", "st60_min",
                                                                                      "sum_qrs_pp", "sum_p_abs"]}}
                      for k, v in sub.items()}
    R["sensitivity"] = {"all_pairs_cluster": {f: fit(df, f"d_{f}", cluster=True, impute=False)
                                              for f in ["st60_min", "m_early", "m_late"]},
                        "complete_case": {f: fit(P, f"d_{f}", impute=False) for f in ["st60_min", "m_early", "m_late"]}}
    out = HERE / f"out/results_{args.cohort}.json"
    out.write_text(json.dumps(R, indent=1))
    print(json.dumps({k: R[k] for k in ["flow", "H1_st60_min", "H2a_m_early", "H2a_m_late", "holm_adjusted",
                                        "H2a_supported"]}, indent=1))


if __name__ == "__main__":
    main()
