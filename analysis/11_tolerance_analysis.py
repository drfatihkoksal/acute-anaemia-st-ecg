"""concept_tolerance.md v1.0 §7 — pre-registered analysis. Identical code for both cohorts.

  python3 11_tolerance_analysis.py --cohort mgh
  python3 11_tolerance_analysis.py --cohort mimic     (only after hashes are logged)

Input: out/tolerance_<cohort>.parquet from the cohort builder.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
import statsmodels.api as sm
from lifelines import CoxPHFitter
from scipy import stats
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
M = 10
NBOOT = 2000
SEED = 20261002
LABS = ["post_cr", "post_bun", "post_na", "post_k", "post_alb", "post_wbc", "post_plt", "post_bili", "post_inr",
        "post_lactate", "post_glucose", "post_hco3"]
# base_hb omitted: base_hb = delta_hb + post_hb exactly (collinear with the post_hb spline); logged §10
BASE = ["age", "male", "year", "delta_hb", "bleed_d62", "bleed_major"] + LABS + \
       ["lactate_measured", "inr_measured", "bnp_measured", "cad", "hf", "ckd", "dm", "elix_vw",
        "base_st60_min", "base_qrs_dur_ms"]
HR = ["post_hr"]
IMPUTE = [c for c in BASE + HR + ["post_hb"] if c.startswith("post_") or c in ("base_hb", "base_st60_min", "base_qrs_dur_ms")]


def rubin(ests, ses):
    q = np.mean(ests); w = np.mean(np.square(ses)); b = np.var(ests, ddof=1) if len(ests) > 1 else 0.0
    se = np.sqrt(w + (1 + 1 / len(ests)) * b)
    return float(q), float(se)


def or_out(q, se):
    return {"OR": float(np.exp(q)), "ci95": [float(np.exp(q - 1.96 * se)), float(np.exp(q + 1.96 * se))],
            "p": float(2 * stats.norm.sf(abs(q / se))), "log_or": q, "se": se}


def imputations(d, m=M):
    cols = [c for c in IMPUTE if c in d]
    out = []
    for i in range(m):
        X = d.copy()
        imp = IterativeImputer(sample_posterior=True, random_state=SEED + i, max_iter=10)
        aux = ["died30", "new_std", "age", "male", "anemic"] if "anemic" in d else ["died30", "new_std", "age", "male"]
        X[cols] = imp.fit_transform(d[cols + aux].astype(float))[:, :len(cols)]
        out.append(X)
    return out


def design(d, expo, covs, spline_df=3, extra_terms=""):
    """Design matrix; non-binary covariates standardised for numerical stability (exposure/interaction terms untouched)."""
    rhs = f"{expo} + bs(post_hb, df={spline_df}) + " + " + ".join(covs) + extra_terms
    X = patsy.dmatrix(rhs, d, return_type="dataframe")
    for c in covs:
        if c not in X:
            continue
        if X[c].nunique() > 2:
            X[c] = (X[c] - X[c].mean()) / (X[c].std() or 1.0)
        elif min((X[c] == 1).sum(), (X[c] != 1).sum()) < 10:   # outcome-blind sparsity rule (logged §10)
            X = X.drop(columns=c)
    return X


def logit_fit(y, X, w=None, robust=True):
    if w is None:
        r = sm.GLM(y, X, family=sm.families.Binomial()).fit()
    else:
        r = sm.GLM(y, X, family=sm.families.Binomial(), var_weights=w).fit(cov_type="HC0")
    return r


def term(r, name):
    return r.params[name], r.bse[name]


def primary(imps, expo="new_std", covs=BASE + HR, outcome="died30"):
    ests, ses = [], []
    for d in imps:
        X = design(d, expo, covs)
        r = logit_fit(d.loc[X.index, outcome].astype(float), X)
        e, s = term(r, expo); ests.append(e); ses.append(s)
    return or_out(*rubin(ests, ses))


def interaction(imps, var, covs=BASE + HR, expo="new_std"):
    """OR of exposure within strata of a binary modifier + interaction p (Rubin)."""
    res = {"e0": ([], []), "e1": ([], []), "int": ([], [])}
    for d in imps:
        X = design(d, expo, [c for c in covs if c != var], extra_terms=f" + {var} + {expo}:{var}")
        r = logit_fit(d.loc[X.index, "died30"].astype(float), X)
        V = r.cov_params()
        b0, b1 = r.params[expo], r.params[f"{expo}:{var}"]
        res["e0"][0].append(b0); res["e0"][1].append(np.sqrt(V.loc[expo, expo]))
        res["e1"][0].append(b0 + b1)
        res["e1"][1].append(np.sqrt(V.loc[expo, expo] + V.loc[f"{expo}:{var}", f"{expo}:{var}"] + 2 * V.loc[expo, f"{expo}:{var}"]))
        res["int"][0].append(b1); res["int"][1].append(np.sqrt(V.loc[f"{expo}:{var}", f"{expo}:{var}"]))
    return {"n_1": None, f"OR_{var}=0": or_out(*rubin(*res["e0"])), f"OR_{var}=1": or_out(*rubin(*res["e1"])),
            "ratio_of_OR": or_out(*rubin(*res["int"]))}


def specificity(T):
    """Anaemic first index vs contrast cohort, contrast reweighted to anaemic age/year distribution."""
    d = T[(T["first"]) & (T.fu30_complete == 1)].copy()
    Z = patsy.dmatrix("bs(age, df=4) + C(year)", d, return_type="dataframe")
    lr = LogisticRegression(max_iter=3000, C=1e6).fit(Z, d.anemic)
    p = lr.predict_proba(Z)[:, 1]
    w = np.where(d.anemic == 1, 1.0, p / (1 - p))
    w[d.anemic == 0] = w[d.anemic == 0] / w[d.anemic == 0].mean()
    d["w"] = w
    imps = imputations(d, m=M)
    covs = [c for c in BASE + HR if c not in ("bleed_d62", "bleed_major", "delta_hb")]  # undefined/degenerate in contrast
    ests, ses, e_an, e_ct = [], [], ([], []), ([], [])
    for x in imps:
        X = design(x, "new_std", covs, extra_terms=" + anemic + new_std:anemic")
        r = logit_fit(x.loc[X.index, "died30"].astype(float), X, w=x.loc[X.index, "w"].to_numpy())
        V = r.cov_params()
        b0, b1 = r.params["new_std"], r.params["new_std:anemic"]
        ests.append(b1); ses.append(np.sqrt(V.loc["new_std:anemic", "new_std:anemic"]))
        e_ct[0].append(b0); e_ct[1].append(np.sqrt(V.loc["new_std", "new_std"]))
        e_an[0].append(b0 + b1)
        e_an[1].append(np.sqrt(V.loc["new_std", "new_std"] + V.loc["new_std:anemic", "new_std:anemic"] + 2 * V.loc["new_std", "new_std:anemic"]))
    return {"n_anemic": int((d.anemic == 1).sum()), "n_contrast": int((d.anemic == 0).sum()),
            "ess_contrast": float(w[d.anemic == 0].sum() ** 2 / (w[d.anemic == 0] ** 2).sum()),
            "OR_contrast": or_out(*rubin(*e_ct)), "OR_anemic": or_out(*rubin(*e_an)),
            "ratio_of_OR_anemic_vs_contrast": or_out(*rubin(ests, ses))}


def incremental(d):
    """ΔC and continuous NRI, first imputation, bootstrap NBOOT."""
    Xb = design(d, "1", BASE + HR).drop(columns=["Intercept"], errors="ignore")
    Xf = Xb.copy(); Xf["new_std"] = d.loc[Xb.index, "new_std"].to_numpy()
    y = d.loc[Xb.index, "died30"].to_numpy()
    mu, sd = Xf.mean(), Xf.std().replace(0, 1)
    Zb, Zf = ((Xb - mu[Xb.columns]) / sd[Xb.columns]).to_numpy(), ((Xf - mu) / sd).to_numpy()

    def fitpred(idx):
        pb = LogisticRegression(C=1e6, max_iter=5000).fit(Zb[idx], y[idx]).predict_proba(Zb[idx])[:, 1]
        pf = LogisticRegression(C=1e6, max_iter=5000).fit(Zf[idx], y[idx]).predict_proba(Zf[idx])[:, 1]
        cb, cf = roc_auc_score(y[idx], pb), roc_auc_score(y[idx], pf)
        up = pf > pb
        ev, ne = y[idx] == 1, y[idx] == 0
        nri = (up[ev].mean() - (~up)[ev].mean()) + ((~up)[ne].mean() - up[ne].mean())
        return cb, cf, nri, pf
    cb, cf, nri, pf = fitpred(np.arange(len(y)))
    rng = np.random.default_rng(SEED)
    boot = np.array([fitpred(rng.integers(0, len(y), len(y)))[:3] for _ in range(NBOOT)])
    dC = boot[:, 1] - boot[:, 0]
    dec = pd.qcut(pf, 10, labels=False, duplicates="drop")
    calib = pd.DataFrame({"pred": pf, "obs": y, "dec": dec}).groupby("dec").mean().round(4).to_dict("list")
    return {"C_base": cb, "C_full": cf, "dC": cf - cb, "dC_ci95": [float(np.percentile(dC, 2.5)), float(np.percentile(dC, 97.5))],
            "NRI": nri, "NRI_ci95": [float(np.percentile(boot[:, 2], 2.5)), float(np.percentile(boot[:, 2], 97.5))],
            "calibration_deciles": calib}


def cox1y(imps):
    ests, ses = [], []
    for d in imps:
        X = design(d, "new_std", BASE + HR).drop(columns=["Intercept"])
        X = X.loc[:, X.std() > 0]
        X["T"], X["E"] = d.loc[X.index, "time365"].to_numpy(), d.loc[X.index, "event365"].to_numpy()
        cph = CoxPHFitter(penalizer=1e-4).fit(X, "T", "E")
        ests.append(cph.params_["new_std"]); ses.append(cph.standard_errors_["new_std"])
    q, se = rubin(ests, ses)
    return {"HR": float(np.exp(q)), "ci95": [float(np.exp(q - 1.96 * se)), float(np.exp(q + 1.96 * se))],
            "p": float(2 * stats.norm.sf(abs(q / se)))}


def injury_ipw(imps):
    ests, ses = [], []
    for d in imps:
        Xm = design(d, "new_std", BASE + HR)
        pm = logit_fit(d.loc[Xm.index, "trop_measured"].astype(float), Xm).predict(Xm)
        s = (d.loc[Xm.index, "trop_measured"] == 1) & d.loc[Xm.index, "mi_injury"].notna()
        X = Xm[s]
        r = logit_fit(d.loc[Xm.index, "mi_injury"][s].astype(float), X, w=(1 / np.clip(pm[s], 0.02, 1)).to_numpy())
        e, se = term(r, "new_std"); ests.append(e); ses.append(se)
    out = or_out(*rubin(ests, ses))
    out["n_measured"] = int(imps[0].trop_measured.sum())
    out["injury_rate_measured"] = float(imps[0].mi_injury.mean())
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--cohort", required=True); a = ap.parse_args()
    T = pd.read_parquet(HERE / f"out/tolerance_{a.cohort}.parquet")
    T = T[T.early_death_before_index == 0]
    T["first"] = T["first"].astype(bool); T["random"] = T["random"].astype(bool)
    an = T[(T.anemic == 1) & T["first"] & (T.fu30_complete == 1)].copy()
    R = {"cohort": a.cohort, "code_sha1": hashlib.sha1(Path(__file__).read_bytes()).hexdigest(),
         "n_index": len(an), "deaths30": int(an.died30.sum()), "new_std_n": int(an.new_std.sum()),
         "deaths30_new_std": int(an.died30[an.new_std == 1].sum()),
         "crude_mortality": {"new_std": float(an.died30[an.new_std == 1].mean()), "no_std": float(an.died30[an.new_std == 0].mean())}}
    imps = imputations(an)
    R["primary_OR_new_std"] = primary(imps)
    R["specificity"] = specificity(T)
    R["incremental"] = incremental(imps[0])
    R["subgroups"] = {}
    for v, f in [("hb_le7", lambda d: (d.post_hb <= 7).astype(int)), ("cad", None),
                 ("bleed_any", lambda d: ((d.bleed_d62 + d.bleed_major) > 0).astype(int)), ("hr_rise10", None)]:
        ii = [x.assign(**{v: f(x)}) if f else x for x in imps]
        R["subgroups"][v] = interaction(ii, v)
        R["subgroups"][v]["n_1"] = int(ii[0][v].sum())
        R["subgroups"][v]["cells_stratum_expo_deaths_n"] = (ii[0].groupby([v, "new_std"]).died30.agg(["sum", "size"])
                                                            .reset_index().astype(int).to_dict("records"))
    R["secondary"] = {"cox_1y": cox1y(imps), "injury_ipw": injury_ipw(imps),
                      "new_std_1mm": primary(imps, expo="new_std_1mm"), "std_burden_per_0.1mV": primary(imps, expo="std_burden")}
    R["sensitivity"] = {
        "no_hr_adjust": primary(imps, covs=BASE),
        "complete_case": primary([an.dropna(subset=BASE + HR + ["post_hb"])]),
        "only_500hz": primary([x[(x.base_fs == 500) & (x.post_fs == 500)] for x in imps]),
    }
    rnd = T[(T.anemic == 1) & T["random"] & (T.fu30_complete == 1)].copy()
    R["sensitivity"]["random_pair"] = primary(imputations(rnd))
    R["sensitivity"]["complete_case_n"] = int(len(an.dropna(subset=BASE + HR + ["post_hb"])))
    R["deep_learning_secondary"] = "not run in this script (separate pre-specified step)"
    out = HERE / f"out/tolerance_results_{a.cohort}.json"
    out.write_text(json.dumps(R, indent=1, default=float))
    print(json.dumps({k: R[k] for k in ["n_index", "deaths30", "new_std_n", "deaths30_new_std", "crude_mortality",
                                        "primary_OR_new_std"]}, indent=1, default=float))


if __name__ == "__main__":
    main()
