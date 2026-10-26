"""Addendum P (concept_addendum_posthoc.md): post hoc analyses anticipating peer review.

    python3 12_posthoc_reviewer.py

Reuses the frozen Analysis 1 code (06_analysis.py) and Analysis 2 code (11_tolerance_analysis.py) by
import, without modifying them. Writes aggregates only to out/posthoc_reviewer.json; counts below 10
are suppressed.
"""
from __future__ import annotations

import glob
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
LEADS = ["I", "II", "III", "aVL", "aVF", "V1", "V2", "V3", "V4", "V5", "V6"]     # non-aVR, as in the cohort builders
THR = 0.05
UPSLOPE = 0.01
SMALL = 10
CATS = [("rise", -np.inf, -0.5), ("stable", -0.5, 0.5), ("0.5-2", 0.5, 2.0), ("2-3", 2.0, 3.0), ("ge3", 3.0, np.inf)]


def load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


A = load("analysis06", "06_analysis.py")
T11 = load("tolerance11", "11_tolerance_analysis.py")


def cnt(n):
    """Small-cell rule: counts below 10 are not written."""
    return int(n) if n >= SMALL else None


def category(dhb):
    out = pd.Series(index=dhb.index, dtype=object)
    for name, lo, hi in CATS:
        # stable is the open interval (-0.5, 0.5); the others are [lo, hi)
        m = (dhb > lo) & (dhb < hi) if name == "stable" else ((dhb > lo) & (dhb <= hi) if name == "rise" else (dhb >= lo) & (dhb < hi))
        out[m] = name
    return out


def hits(d, kind):
    return sum(((d[f"post_{kind}_{l}"] <= -THR) & ((d[f"post_{kind}_{l}"] - d[f"base_{kind}_{l}"]) <= -THR)).astype(int)
               for l in LEADS)


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [float(c - h), float(c + h)]


def pool(items):
    """Fixed-effect inverse-variance pooling of log ORs; items = [(log_or, se), ...]."""
    b = np.array([x[0] for x in items]); s = np.array([x[1] for x in items])
    w = 1 / s ** 2
    q = float((w * b).sum() / w.sum()); se = float(np.sqrt(1 / w.sum()))
    Q = float((w * (b - q) ** 2).sum())
    return {**T11.or_out(q, se), "het_p": float(stats.chi2.sf(Q, len(b) - 1))}


def evalue(or_point, or_lo):
    def e(rr):
        rr = 1 / rr if rr < 1 else rr
        return float(rr + np.sqrt(rr * (rr - 1)))
    rr, rr_lo = np.sqrt(or_point), np.sqrt(or_lo)
    return {"rr_approx": float(rr), "E_point": e(rr), "E_ci": 1.0 if rr_lo <= 1 else e(rr_lo)}


# ------------------------------------------------------------------------------- Analysis 1
def analysis1(cohort):
    raw = pd.read_parquet(OUT / f"analysis_{cohort}.parquet")
    if "post_key" not in raw:
        raw["post_key"] = raw.post_fname
    df, _ = A.prep(raw)
    P = df[df.primary].copy()
    R = {}
    # P1 adults only
    ad = P[P.age >= 18]
    r = A.fit(ad, "d_st60_min"); r["p_tost"] = A.tost(r); r["equivalent"] = r["p_tost"] < 0.05
    R["P1_adults"] = r
    # P2 categories of Hb change
    P["cat"] = category(P.delta_hb)
    cols = list(A.ADJ) + ["base_st60_min"]
    d = P[["d_st60_min", "cat"] + cols].replace([np.inf, -np.inf], np.nan).dropna()
    X = pd.get_dummies(d.cat).astype(float).drop(columns="stable")
    X = sm.add_constant(pd.concat([X, d[cols].astype(float)], axis=1))
    fit = sm.OLS(d.d_st60_min.astype(float), X).fit(cov_type="HC3")
    cats = {}
    for name, *_ in CATS:
        n = int((d.cat == name).sum())
        if name == "stable":
            cats[name] = {"n": cnt(n)}
            continue
        b, se = float(fit.params[name]), float(fit.bse[name])
        cats[name] = {"n": cnt(n), "diff_vs_stable": b, "ci90": [b - 1.645 * se, b + 1.645 * se], "se": se}
    b, se = cats["ge3"]["diff_vs_stable"], cats["ge3"]["se"]
    m3 = 3 * A.MARGIN
    cats["ge3"]["p_tost_margin_0.03"] = float(max(stats.norm.sf((b + m3) / se), stats.norm.cdf((b - m3) / se)))
    cats["ge3"]["median_dhb"] = float(P.loc[d.index[d.cat == "ge3"], "delta_hb"].median())
    R["P2_categories"] = {"n_complete": int(len(d)), "categories": cats}
    # P3 very low nadir
    lo = P[P.post_hb <= 7]
    R["P3_post_hb_le7"] = {"n": cnt(len(lo)), "st60_min": A.fit(lo, "d_st60_min", impute=False)}
    # P4 frequency of new ST depression by category (adults)
    ad = P[P.age >= 18].copy()
    ad["new_std"] = (hits(ad, "st60") >= 2).astype(int)
    freq = {}
    for name, *_ in CATS:
        g = ad[ad.cat == name]
        k, n = int(g.new_std.sum()), len(g)
        ok = k >= SMALL and n - k >= SMALL
        freq[name] = {"n": cnt(n), "events": cnt(k), "pct": float(100 * k / n) if ok else None,
                      "ci95_pct": [100 * x for x in wilson(k, n)] if ok else None}
    R["P4_new_std_by_category_adults"] = freq
    # P5 noise floor, from the frozen results
    R["P5_noise_floor_sd_stable_hb_mV"] = json.loads((OUT / f"results_{cohort}.json").read_text())["noise_floor_sd_stable_hb"]["d_st60_min"]
    return R


# ------------------------------------------------------------------------------- Analysis 2
def j_amps_mgh(T):
    cols = [f"j_amp_{l}" for l in LEADS]
    meas = pd.concat([pd.read_parquet(f, columns=["key"] + cols) for f in sorted(glob.glob(str(OUT / "measure_mgh/part_*.parquet")))])
    meas = meas.drop_duplicates("key").set_index("key")
    T = T.drop(columns=[c for c in T.columns if c.startswith(("base_j_amp_", "post_j_amp_"))])
    for side in ("base", "post"):
        T = T.join(meas.add_prefix(f"{side}_"), on=f"{side}_fname")
    return T


def analysis2(cohort):
    T = pd.read_parquet(OUT / f"tolerance_{cohort}.parquet")
    T = T[T.early_death_before_index == 0]
    T["first"] = T["first"].astype(bool)
    if cohort == "mgh":
        T = j_amps_mgh(T)
    an = T[(T.anemic == 1) & T["first"] & (T.fu30_complete == 1)].copy()
    an["new_std_J"] = (hits(an, "j_amp") >= 2).astype(int)
    # P7 morphology among new ST depression (primary definition)
    n_hd = 0
    for l in LEADS:
        q = (an[f"post_st60_{l}"] <= -THR) & ((an[f"post_st60_{l}"] - an[f"base_st60_{l}"]) <= -THR)
        flat = (an[f"post_st60_{l}"] - an[f"post_j_amp_{l}"]) <= UPSLOPE
        n_hd = n_hd + (q & flat).astype(int)
    an["std_hd"] = ((an.new_std == 1) & (n_hd >= 2)).astype(int)
    an["std_up"] = ((an.new_std == 1) & (an.std_hd == 0)).astype(int)
    R = {"n_index": len(an), "j_amp_missing_frac": float(an[[f"post_j_amp_{l}" for l in LEADS]].isna().any(axis=1).mean()),
         "new_std_J_n": cnt(an.new_std_J.sum()), "new_std_n": cnt(an.new_std.sum()),
         "agreement_primary_vs_J": {"both": cnt(((an.new_std == 1) & (an.new_std_J == 1)).sum()),
                                    "primary_only": cnt(((an.new_std == 1) & (an.new_std_J == 0)).sum()),
                                    "J_only": cnt(((an.new_std == 0) & (an.new_std_J == 1)).sum())},
         "std_hd_n": cnt(an.std_hd.sum()), "std_up_n": cnt(an.std_up.sum()),
         "std_hd_deaths30": cnt(an.died30[an.std_hd == 1].sum()), "std_up_deaths30": cnt(an.died30[an.std_up == 1].sum())}
    imps = T11.imputations(an)
    R["P6_J_point_OR"] = T11.primary(imps, expo="new_std_J")
    est = {"std_hd": ([], []), "std_up": ([], [])}
    for d in imps:
        X = T11.design(d, "std_hd + std_up", T11.BASE + T11.HR)
        r = T11.logit_fit(d.loc[X.index, "died30"].astype(float), X)
        for k in est:
            est[k][0].append(r.params[k]); est[k][1].append(r.bse[k])
    R["P7_morphology_OR"] = {k: T11.or_out(*T11.rubin(*v)) for k, v in est.items()}
    return R


def main():
    res = {"addendum": "concept_addendum_posthoc.md", "label": "post hoc",
           "code_sha1": hashlib.sha1(Path(__file__).read_bytes()).hexdigest(),
           "addendum_sha256": hashlib.sha256((HERE.parent / "concept_addendum_posthoc.md").read_bytes()).hexdigest()}
    for c in ("mgh", "mimic"):
        res[c] = {"analysis1": analysis1(c), "analysis2": analysis2(c)}
        print(c, "done", flush=True)
    pooled = {}
    for key in ("P6_J_point_OR",):
        pooled[key] = pool([(res[c]["analysis2"][key]["log_or"], res[c]["analysis2"][key]["se"]) for c in ("mgh", "mimic")])
    for k in ("std_hd", "std_up"):
        pooled[f"P7_{k}"] = pool([(res[c]["analysis2"]["P7_morphology_OR"][k]["log_or"], res[c]["analysis2"]["P7_morphology_OR"][k]["se"])
                                  for c in ("mgh", "mimic")])
    tol = {c: json.loads((OUT / f"tolerance_results_{c}.json").read_text()) for c in ("mgh", "mimic")}
    pooled["P8_no_hr_adjust"] = pool([(tol[c]["sensitivity"]["no_hr_adjust"]["log_or"], tol[c]["sensitivity"]["no_hr_adjust"]["se"])
                                      for c in ("mgh", "mimic")])
    prim = json.loads((OUT / "tolerance_pooled.json").read_text())["primary"]
    res["pooled"] = pooled
    res["P9_evalues"] = {"primary": evalue(prim["OR"], prim["ci95"][0]),
                         "no_hr_adjust": evalue(pooled["P8_no_hr_adjust"]["OR"], pooled["P8_no_hr_adjust"]["ci95"][0])}
    (OUT / "posthoc_reviewer.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({"pooled": {k: {kk: v[kk] for kk in ("OR", "ci95", "het_p")} for k, v in pooled.items()},
                      "evalues": res["P9_evalues"]}, indent=1))


if __name__ == "__main__":
    main()
