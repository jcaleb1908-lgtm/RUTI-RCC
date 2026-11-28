"""analysis.py - estimation tools for the RUTI-RCC outcome cells (F05-F10). Study-agnostic numerics only.

Everything here runs on participant-level arrays inside the Workbench and returns aggregate numbers. Masking is done
by the cells through aou_core.show; nothing here prints or saves.

  cox()               weighted Cox model, Breslow ties, counting-process entry (delayed entry, time-varying rows),
                      strata, cluster-robust (sandwich) variance; optional Schoenfeld-type PH score test per term
  ps_weights()        stabilized inverse-probability-of-treatment weights from a binary or multinomial logistic model,
                      truncated at percentiles
  balance()           standardized differences before and after weighting, per design column
  aalen_johansen()    cumulative incidence with a competing event (weights, delayed entry)
  split_cumulative()  counting-process rows for time-varying cumulative counts (Q2b, post-T0 imaging)
  fine_gray_rows()    rows and weights that turn cox() into a Fine-Gray model (Geskus 2011 weighting)
  mice_pmm()          multiple imputation by chained equations with predictive mean matching
  reri()              relative excess risk due to interaction with a delta-method CI
  rubin_log()         Rubin's rules on the log scale

Statsmodels' PHReg has no case weights, so the Cox model is written here and tested against PHReg (unweighted)
and against row replication (integer weights) in test_synthetic.py.
"""
import warnings

import numpy as np
import pandas as pd
from scipy import sparse, stats

DIVERGENT_LOG = 10.0     # |log HR| above this on any term: not estimable (aou_core M12 rule)


# =============================================================================
# COX MODEL
# =============================================================================
def _codes(v, n):
    if v is None:
        return np.zeros(n, dtype=int)
    return pd.factorize(pd.Series(v).astype(str).values)[0]


def cox(time, event, X, w=None, entry=None, strata=None, cluster=None, names=None, ph=False, max_iter=60, tol=1e-9):
    """Weighted Cox partial likelihood (Breslow ties).
    time: exit time; entry: entry time (row at risk for event times u with entry < u <= time); event: 0/1.
    w: case weights (IPTW, Fine-Gray). strata: baseline hazard per stratum. cluster: rows of one person (robust SE).
    Returns dict: beta, se (robust), vcov (robust), vcov_model, names, loglik, converged, n_rows, n_events,
    divergent, and with ph=True 'ph' = DataFrame(term, chi2, p) for beta(t) = beta + theta*log(t)."""
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    t = np.asarray(time, dtype=float)
    n, p = X.shape
    d = np.asarray(event).astype(int) == 1
    w = np.ones(n) if w is None else np.asarray(w, dtype=float)
    en = np.zeros(n) if entry is None else np.asarray(entry, dtype=float)
    if np.any(en >= t):
        raise ValueError("cox: every row needs entry < time (drop rows that end before they start)")
    if not (np.isfinite(X).all() and np.isfinite(t).all() and np.isfinite(w).all()):
        raise ValueError("cox: non-finite values in X, time or weights")
    sc = _codes(strata, n)
    cl = _codes(cluster, n) if cluster is not None else np.arange(n)
    names = list(names) if names is not None else [f"x{i}" for i in range(p)]
    pre = []
    for k in np.unique(sc):
        idx = np.where(sc == k)[0]
        ti, di = t[idx], d[idx]
        tau = np.unique(ti[di])
        if len(tau) == 0:
            continue
        a = np.searchsorted(tau, en[idx], side="right")
        b = np.searchsorted(tau, ti, side="right")
        D = len(tau)
        ones = np.ones(len(idx))
        M = (sparse.csr_matrix((ones, (a, np.arange(len(idx)))), shape=(D + 1, len(idx)))
             - sparse.csr_matrix((ones, (b, np.arange(len(idx)))), shape=(D + 1, len(idx))))
        je = np.searchsorted(tau, ti[di])
        pre.append((idx, tau, a, b, di, je, M))

    def risk_sums(M, v, D):
        return np.cumsum(np.asarray(M @ v), axis=0)[:D]

    def evaluate(beta, resid=False, want_ph=False):
        ll, U, I = 0.0, np.zeros(p), np.zeros((p, p))
        R = np.zeros((n, p)) if resid else None
        phs = []
        for idx, tau, a, b, di, je, M in pre:
            Xs, ws = X[idx], w[idx]
            eta = Xs @ beta
            c = eta.max()
            r = ws * np.exp(eta - c)
            D = len(tau)
            S0 = risk_sums(M, r, D)
            S1 = risk_sums(M, r[:, None] * Xs, D)
            dj = np.bincount(je, weights=ws[di], minlength=D)
            xw = np.zeros((D, p))
            np.add.at(xw, je, ws[di, None] * Xs[di])
            xbar = S1 / S0[:, None]
            ll += float((ws[di] * eta[di]).sum() - (dj * (np.log(S0) + c)).sum())
            U += (xw - dj[:, None] * xbar).sum(axis=0)
            cj = dj / S0
            cc = np.concatenate([[0.0], np.cumsum(cj)])
            Ci = cc[b] - cc[a]
            I += (Xs * (r * Ci)[:, None]).T @ Xs - (xbar * dj[:, None]).T @ xbar
            if resid:
                G = np.vstack([np.zeros((1, p)), np.cumsum(cj[:, None] * xbar, axis=0)])
                Ri = -r[:, None] * (Xs * Ci[:, None] - (G[b] - G[a]))
                Ri[di] += ws[di, None] * (Xs[di] - xbar[je])
                R[idx] = Ri
            if want_ph:
                S2d = risk_sums(M, r[:, None] * Xs ** 2, D)
                phs.append((tau, dj, xw - dj[:, None] * xbar, S2d / S0[:, None] - xbar ** 2))
        return ll, U, I, R, phs

    beta = np.zeros(p)
    ll, U, I, _, _ = evaluate(beta)
    converged = False
    it = 0
    for it in range(max_iter):
        try:
            step = np.linalg.solve(I, U)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(I, U, rcond=None)[0]
        lam = 1.0
        while True:
            nb = beta + lam * step
            nll, nU, nI, _, _ = evaluate(nb)
            if np.isfinite(nll) and nll >= ll - 1e-10 or lam < 1e-4:
                break
            lam /= 2
        done = abs(nll - ll) < tol * (1 + abs(ll)) and np.max(np.abs(lam * step)) < 1e-6
        beta, ll, U, I = nb, nll, nU, nI
        if done:
            converged = True
            break
        if np.max(np.abs(beta)) > 50:
            break
    ll, U, I, R, phs = evaluate(beta, resid=True, want_ph=ph)
    Iinv = np.linalg.pinv(I)
    Rc = np.zeros((cl.max() + 1, p))
    np.add.at(Rc, cl, R)
    V = Iinv @ (Rc.T @ Rc) @ Iinv
    se = np.sqrt(np.clip(np.diag(V), 0, None))
    out = {"beta": beta, "se": se, "vcov": V, "vcov_model": Iinv, "names": names, "loglik": ll,
           "converged": converged, "iterations": it + 1, "n_rows": n, "n_events": int(d.sum()),
           "divergent": [names[i] for i in range(p) if not np.isfinite(beta[i]) or not np.isfinite(se[i])
                         or abs(beta[i]) > DIVERGENT_LOG]}
    if ph:
        tau = np.concatenate([x[0] for x in phs]) if phs else np.zeros(0)
        dj = np.concatenate([x[1] for x in phs]) if phs else np.zeros(0)
        res = np.vstack([x[2] for x in phs]) if phs else np.zeros((0, p))
        var = np.vstack([x[3] for x in phs]) if phs else np.zeros((0, p))
        g = np.log(np.clip(tau, 1e-6, None))
        gbar = (dj * g).sum() / dj.sum() if dj.sum() > 0 else 0.0
        num = ((g - gbar)[:, None] * res).sum(axis=0)
        den = ((dj * (g - gbar) ** 2)[:, None] * var).sum(axis=0)
        chi2 = np.where(den > 0, num ** 2 / np.where(den > 0, den, 1), np.nan)
        out["ph"] = pd.DataFrame({"term": names, "chi2": chi2, "p": stats.chi2.sf(chi2, 1)})
    return out


def wald_equal(fit, i, j):
    """Wald test that coefficients i and j are equal (robust vcov). Returns (chi2, p)."""
    V = fit["vcov"]
    diff = fit["beta"][i] - fit["beta"][j]
    v = V[i, i] + V[j, j] - 2 * V[i, j]
    if not v > 0:
        return np.nan, np.nan
    c2 = diff ** 2 / v
    return float(c2), float(stats.chi2.sf(c2, 1))


def wald_joint(fit, idx):
    """Joint Wald test that the coefficients in idx are all 0 (robust vcov). Returns (chi2, df, p)."""
    idx = list(idx)
    b = fit["beta"][idx]
    V = fit["vcov"][np.ix_(idx, idx)]
    try:
        c2 = float(b @ np.linalg.solve(V, b))
    except np.linalg.LinAlgError:
        return np.nan, len(idx), np.nan
    return c2, len(idx), float(stats.chi2.sf(c2, len(idx)))


# =============================================================================
# PROPENSITY SCORES AND BALANCE
# =============================================================================
def drop_constant(X):
    """Design columns that vary in these rows (a dummy for an absent level makes the model singular)."""
    keep = [c for c in X.columns if X[c].nunique(dropna=True) > 1]
    return X[keep], [c for c in X.columns if c not in keep]


def drop_collinear(X, tol=1e-8):
    """Drop design columns that are exact linear combinations of earlier ones (pivoted QR rank check)."""
    from scipy.linalg import qr
    if X.shape[1] == 0:
        return X, []
    Z = X.values.astype(float)
    Z = (Z - Z.mean(axis=0)) / np.where(Z.std(axis=0) > 0, Z.std(axis=0), 1)
    _, Rm, piv = qr(Z, mode="economic", pivoting=True)
    d = np.abs(np.diag(Rm))
    rank = int((d > tol * max(d.max(), 1)).sum())
    keep = sorted(piv[:rank])
    cols = [X.columns[i] for i in keep]
    return X[cols], [c for c in X.columns if c not in cols]


def ps_weights(A, X, levels, trunc=(1, 99)):
    """Stabilized ATE weights P(A = a) / P(A = a | X) from a logistic (2 levels) or multinomial logistic model.
    A: exposure labels; X: design DataFrame without constant; levels: [reference, other levels...].
    Weights are truncated at the trunc percentiles of the whole sample. Returns (weights, info dict)."""
    import statsmodels.api as sm
    A = pd.Series(A).reset_index(drop=True)
    y = A.map({lev: i for i, lev in enumerate(levels)}).values
    if np.isnan(y.astype(float)).any():
        raise ValueError("ps_weights: exposure has a level outside `levels`")
    y = y.astype(int)
    X, dropped = drop_constant(X.reset_index(drop=True))
    X, coll = drop_collinear(X)
    dropped += coll
    Xc = sm.add_constant(X, has_constant="add").astype(float)
    info = {"dropped_constant": dropped, "converged": True, "method": "logit" if len(levels) == 2 else "mnlogit"}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if len(levels) == 2:
            r = sm.GLM(y, Xc, family=sm.families.Binomial()).fit()
            p1 = np.asarray(r.predict(Xc))
            P = np.column_stack([1 - p1, p1])
            info["converged"] = bool(r.converged)
        else:
            try:
                r = sm.MNLogit(y, Xc).fit(method="newton", maxiter=100, disp=0)
            except Exception:
                r = sm.MNLogit(y, Xc).fit(method="bfgs", maxiter=1000, disp=0)
            P = np.asarray(r.predict(Xc))
            info["converged"] = bool(r.mle_retvals.get("converged", True))
    P = np.clip(P, 1e-8, 1)
    pa = P[np.arange(len(y)), y]
    marg = np.bincount(y, minlength=len(levels)) / len(y)
    sw = marg[y] / pa
    lo, hi = np.percentile(sw, trunc)
    info["share_truncated"] = float(np.mean((sw < lo) | (sw > hi)))
    sw = np.clip(sw, lo, hi)
    info["ess"] = {lev: float(sw[y == i].sum() ** 2 / (sw[y == i] ** 2).sum()) if (y == i).any() else np.nan
                   for i, lev in enumerate(levels)}
    info["pctl"] = {lev: np.percentile(sw[y == i], [1, 25, 50, 75, 99]) if (y == i).any() else None
                    for i, lev in enumerate(levels)}
    return sw, info


def _wstats(x, w):
    m = np.average(x, weights=w)
    return m, np.average((x - m) ** 2, weights=w)


def balance(X, A, level, ref, w):
    """Per design column: count of rows with a non-zero value in each arm (for masking binary columns),
    unweighted and weighted standardized difference (level vs ref)."""
    A = np.asarray(A)
    w = np.asarray(w, dtype=float)
    rows = []
    g1, g0 = A == level, A == ref
    for c in X.columns:
        x = X[c].values.astype(float)
        binary = set(np.unique(x[np.isfinite(x)])) <= {0.0, 1.0}
        out = {"term": c, "binary": binary,
               "nz_1": int((x[g1] != 0).sum()) if binary else np.nan, "nz_0": int((x[g0] != 0).sum()) if binary else np.nan}
        for lab, ww in (("unweighted", np.ones_like(w)), ("weighted", w)):
            m1, v1 = _wstats(x[g1], ww[g1])
            m0, v0 = _wstats(x[g0], ww[g0])
            sd = np.sqrt((v1 + v0) / 2)
            out[f"smd_{lab}"] = float((m1 - m0) / sd) if sd > 0 else 0.0
        rows.append(out)
    return pd.DataFrame(rows)


# =============================================================================
# CUMULATIVE INCIDENCE (AALEN-JOHANSEN)
# =============================================================================
def aalen_johansen(t, cause, w=None, entry=None, at=(5.0, 10.0)):
    """Cumulative incidence of cause 1 with cause 2 competing (0 = censored). Weighted, delayed entry allowed.
    Returns DataFrame(time, cif, n_at_risk, events_by) where n_at_risk and events_by are UNWEIGHTED counts
    (persons at risk just before `time`; cause-1 events by `time`), used only to decide masking."""
    t = np.asarray(t, dtype=float)
    cause = np.asarray(cause).astype(int)
    w = np.ones(len(t)) if w is None else np.asarray(w, dtype=float)
    en = np.zeros(len(t)) if entry is None else np.asarray(entry, dtype=float)
    ev = cause > 0
    tau = np.unique(t[ev])
    rows = []
    if len(tau):
        o = np.argsort(t)
        ts, wt = t[o], w[o]
        suf = np.concatenate([np.cumsum(wt[::-1])[::-1], [0.0]])
        oe = np.argsort(en)
        es, we = en[oe], w[oe]
        sufe = np.concatenate([np.cumsum(we[::-1])[::-1], [0.0]])
        Y = suf[np.searchsorted(ts, tau, side="left")] - sufe[np.searchsorted(es, tau, side="left")]
        j = np.searchsorted(tau, t[ev])
        d1 = np.bincount(j, weights=w[ev] * (cause[ev] == 1), minlength=len(tau))
        d2 = np.bincount(j, weights=w[ev] * (cause[ev] == 2), minlength=len(tau))
        Y = np.where(Y > 0, Y, np.nan)
        S = np.cumprod(1 - np.nan_to_num((d1 + d2) / Y))
        S_prev = np.concatenate([[1.0], S[:-1]])
        F1 = np.cumsum(S_prev * np.nan_to_num(d1 / Y))
    for x in at:
        k = np.searchsorted(tau, x, side="right") - 1 if len(tau) else -1
        rows.append({"time": x, "cif": float(F1[k]) if k >= 0 else 0.0,
                     "n_at_risk": int(((t >= x) & (en < x)).sum()),
                     "events_by": int(((cause == 1) & (t <= x)).sum())})
    return pd.DataFrame(rows)


def cif_bootstrap(t, cause, w, entry=None, at=(5.0, 10.0), B=200, seed=1):
    """Bootstrap SD of the CIF at each time (persons resampled; weights held fixed)."""
    rng = np.random.default_rng(seed)
    n = len(t)
    t, cause, w = np.asarray(t), np.asarray(cause), np.asarray(w)
    en = None if entry is None else np.asarray(entry)
    draws = np.zeros((B, len(at)))
    for b in range(B):
        i = rng.integers(0, n, n)
        draws[b] = aalen_johansen(t[i], cause[i], w[i], None if en is None else en[i], at)["cif"].values
    return draws


# =============================================================================
# COUNTING-PROCESS ROWS
# =============================================================================
def split_cumulative(F, changes, variables, id_col="person_id", time="t_years", event="event", entry=None, caps=None):
    """F: one row per person (time = exit in years from outcome-counting start; entry = optional start column).
    changes: rows (id_col, u, var): one increment of `var` at time u (same time scale; u <= start counts from the
    first row). Returns counting-process rows: start, stop, event (on the person's last row only) and, for each
    variable, the cumulative count valid on (start, stop], capped at caps[var]; plus F's other columns."""
    caps = caps or {}
    s0 = F[entry].astype(float).values if entry else np.zeros(len(F))
    P = F[[id_col, time]].assign(_s0=s0)
    ch = changes.loc[changes["var"].isin(variables), [id_col, "u", "var"]].merge(P, on=id_col)
    cuts = ch.loc[(ch["u"] > ch["_s0"]) & (ch["u"] < ch[time]), [id_col, "u"]].drop_duplicates()
    S = pd.concat([P[[id_col]].assign(start=P["_s0"].values), cuts.rename(columns={"u": "start"})], ignore_index=True)
    S["start"] = S["start"].astype(float)
    S = S.sort_values("start", kind="stable").reset_index(drop=True)
    for v in variables:
        cv = ch.loc[ch["var"].eq(v), [id_col, "u"]].sort_values([id_col, "u"], kind="stable")
        cv["_cum"] = cv.groupby(id_col).cumcount() + 1
        cv = cv.sort_values(["u", "_cum"], kind="stable").rename(columns={"u": "_u"})
        cv["_u"] = cv["_u"].astype(float)
        if len(cv):
            S = pd.merge_asof(S, cv, left_on="start", right_on="_u", by=id_col, direction="backward")
            S[v] = S.pop("_cum").fillna(0).astype(int)
            S = S.drop(columns="_u")
        else:
            S[v] = 0
        if v in caps:
            S[v] = S[v].clip(upper=caps[v])
    S = S.sort_values([id_col, "start"], kind="stable").reset_index(drop=True)
    S["stop"] = S.groupby(id_col)["start"].shift(-1)
    last = S["stop"].isna()
    rest = F.drop(columns=[c for c in variables if c in F.columns])
    S = S.merge(rest, on=id_col, how="left")
    S.loc[last, "stop"] = S.loc[last, time]
    S[event] = np.where(last, S[event], 0)
    return S.loc[S["stop"] > S["start"]].reset_index(drop=True)


def fine_gray_rows(F, time="t_years", cause="cause", id_col="person_id"):
    """Rows for a Fine-Gray model fitted by weighted Cox (Geskus 2011): a person with the competing event (cause 2)
    stays at risk after it, on pieces between later cause-1 event times, with weight G(u-) / G(t_i), G = Kaplan-Meier
    of censoring (cause 0). Returns rows with start, stop, event (cause-1 indicator), fg_w and the person columns."""
    t = F[time].values.astype(float)
    c = F[cause].values.astype(int)
    tau = np.unique(t[c == 1])
    # censoring KM: censoring at t (cause 0); events at the same time happen first
    ct = np.unique(t[c == 0])
    o = np.sort(t)
    at_risk = len(t) - np.searchsorted(o, ct, side="left")
    nc = pd.Series(t[c == 0]).value_counts().reindex(ct).values
    Gs = np.cumprod(1 - nc / np.maximum(at_risk, 1)) if len(ct) else np.zeros(0)

    def G(u, left=True):          # G(u-) if left else G(u)
        k = np.searchsorted(ct, u, side="left" if left else "right") - 1
        return np.where(k >= 0, Gs[np.clip(k, 0, None)] if len(Gs) else 1.0, 1.0)
    base = F.assign(start=0.0, stop=t, event=(c == 1).astype(int), fg_w=1.0)
    comp = F.loc[c == 2]
    pieces = []
    if len(comp) and len(tau):
        ti = comp[time].values.astype(float)
        Gt = G(ti, left=False)
        for i in range(len(comp)):
            later = tau[tau > ti[i]]
            if not len(later):
                continue
            starts = np.concatenate([[ti[i]], later[:-1]])
            wi = G(later) / max(Gt[i], 1e-12)
            pieces.append(pd.DataFrame({id_col: comp[id_col].values[i], "start": starts, "stop": later, "event": 0, "fg_w": wi}))
    if pieces:
        extra = pd.concat(pieces, ignore_index=True).merge(F.drop(columns=[time, "event"], errors="ignore"), on=id_col, how="left")
        extra[time] = extra["stop"]
        base = pd.concat([base, extra], ignore_index=True)
    return base.loc[base["stop"] > base["start"]].reset_index(drop=True)


# =============================================================================
# MULTIPLE IMPUTATION (chained equations, predictive mean matching)
# =============================================================================
def mice_pmm(T, P, m=20, iters=5, k=5, seed=1, fit_rows=None):
    """T: DataFrame of targets (numeric; NaN = missing; ordinal targets coded 0..L-1).
    P: complete numeric predictor matrix (same rows). Each target is imputed in turn from P and the other targets,
    by a Bayesian linear model draw and predictive mean matching (k donors), for `iters` cycles.
    Returns a list of m DataFrames (the targets, completed)."""
    rng = np.random.default_rng(seed)
    T = T.reset_index(drop=True).astype(float)
    P = np.asarray(P, dtype=float)
    miss = T.isna()
    out = []
    for _ in range(m):
        Z = T.copy()
        for c in Z.columns:
            obs = Z.loc[~miss[c], c].values
            if miss[c].any():
                Z.loc[miss[c], c] = rng.choice(obs, int(miss[c].sum())) if len(obs) else 0.0
        for _ in range(iters):
            for c in Z.columns:
                if not miss[c].any() or (~miss[c]).sum() < 2:
                    continue
                others = Z.drop(columns=c).values
                Xall = np.column_stack([np.ones(len(Z)), P, others])
                oi = np.where(~miss[c].values)[0]
                fi = oi if fit_rows is None or len(oi) <= fit_rows else rng.choice(oi, fit_rows, replace=False)
                Xo, yo = Xall[fi], Z[c].values[fi]
                XtX = Xo.T @ Xo + 1e-6 * np.eye(Xo.shape[1])
                bh = np.linalg.solve(XtX, Xo.T @ yo)
                resid = yo - Xo @ bh
                dfree = max(len(yo) - Xo.shape[1], 1)
                s2 = float(resid @ resid) / rng.chisquare(dfree)
                L = np.linalg.cholesky(np.linalg.inv(XtX) * s2 + 1e-12 * np.eye(Xo.shape[1]))
                bs = bh + L @ rng.standard_normal(Xo.shape[1])
                yhat_obs = Xall[oi] @ bh
                mi = np.where(miss[c].values)[0]
                yhat_mis = Xall[mi] @ bs
                order = np.argsort(yhat_obs)
                so = yhat_obs[order]
                pos = np.searchsorted(so, yhat_mis)
                lo = np.clip(pos - k, 0, len(so) - 1)
                cand = lo[:, None] + np.arange(2 * k)[None, :]
                cand = np.clip(cand, 0, len(so) - 1)
                dist = np.abs(so[cand] - yhat_mis[:, None])
                near = np.argsort(dist, axis=1)[:, :k]
                pick = cand[np.arange(len(mi)), near[np.arange(len(mi)), rng.integers(0, k, len(mi))]]
                Z.loc[mi, c] = Z[c].values[oi][order][pick]
        out.append(Z)
    return out


def nelson_aalen(t, event):
    """Nelson-Aalen cumulative hazard at each person's own time (White & Royston 2009 imputation predictor)."""
    t = np.asarray(t, float)
    e = np.asarray(event).astype(int)
    tau, cnt = np.unique(t[e == 1], return_counts=True)
    o = np.sort(t)
    at = len(t) - np.searchsorted(o, tau, side="left")
    H = np.cumsum(cnt / at)
    k = np.searchsorted(tau, t, side="right") - 1
    return np.where(k >= 0, H[np.clip(k, 0, None)] if len(H) else 0.0, 0.0)


# =============================================================================
# POOLING AND INTERACTION
# =============================================================================
def rubin_log(betas, ses):
    """Rubin's rules for log hazard ratios. Returns (beta, se, df)."""
    b, s = np.asarray(betas, float), np.asarray(ses, float)
    m = len(b)
    if m == 1:
        return float(b[0]), float(s[0]), np.inf
    qbar, ubar, bv = b.mean(), (s ** 2).mean(), b.var(ddof=1)
    tv = ubar + (1 + 1 / m) * bv
    r = (1 + 1 / m) * bv / ubar if ubar > 0 else np.inf
    dfree = (m - 1) * (1 + 1 / r) ** 2 if np.isfinite(r) and r > 0 else np.inf
    return float(qbar), float(np.sqrt(tv)), float(dfree)


def ratio_ci(beta, se, dfree=np.inf):
    """(HR, lo, hi, p) from a log HR and SE (t reference when the MI df is finite)."""
    q = stats.t.ppf(0.975, dfree) if np.isfinite(dfree) else 1.959963984540054
    z = beta / se if se > 0 else np.nan
    p = 2 * (stats.t.sf(abs(z), dfree) if np.isfinite(dfree) else stats.norm.sf(abs(z)))
    return float(np.exp(beta)), float(np.exp(beta - q * se)), float(np.exp(beta + q * se)), float(p)


def reri(beta, vcov):
    """beta, vcov for (b10, b01, b11) of a joint-category model with a common reference (00).
    RERI = HR11 - HR10 - HR01 + 1, delta-method SE. Returns (reri, lo, hi)."""
    b10, b01, b11 = beta
    g = np.array([-np.exp(b10), -np.exp(b01), np.exp(b11)])
    est = np.exp(b11) - np.exp(b10) - np.exp(b01) + 1
    se = float(np.sqrt(max(g @ np.asarray(vcov) @ g, 0)))
    return float(est), float(est - 1.96 * se), float(est + 1.96 * se)


def rubin_vec(betas, vcovs):
    """Multivariate Rubin's rules: pooled coefficient vector and total covariance (Ubar + (1 + 1/m) B)."""
    B_ = np.asarray(betas, float)
    m = len(B_)
    qbar = B_.mean(axis=0)
    ubar = np.mean(np.asarray(vcovs, float), axis=0)
    if m == 1:
        return qbar, ubar
    return qbar, ubar + (1 + 1 / m) * np.cov(B_.T, ddof=1).reshape(ubar.shape)
