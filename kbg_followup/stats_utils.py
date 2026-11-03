"""Small, dependency-free logistic-regression utilities (ML and Firth) used for speed in simulations."""
import numpy as np
from scipy import stats


def _sig(z):
    return 0.5 * (1 + np.tanh(0.5 * z))


def logit_ml(X, y, max_iter=100, tol=1e-9):
    """Maximum-likelihood logistic regression by Newton-Raphson with step halving.
    Returns (beta, loglik, converged, cov). Separation shows up as non-convergence / huge |beta|."""
    n, p = X.shape
    b = np.zeros(p)
    ll_old = -np.inf
    conv = False
    for _ in range(max_iter):
        eta = X @ b
        mu = _sig(eta)
        w = np.clip(mu * (1 - mu), 1e-12, None)
        g = X.T @ (y - mu)
        Hm = (X * w[:, None]).T @ X
        try:
            step = np.linalg.solve(Hm, g)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(Hm, g, rcond=None)[0]
        t = 1.0
        for _h in range(30):
            bn = b + t * step
            etan = X @ bn
            ll = np.sum(y * etan - np.logaddexp(0, etan))
            if ll >= ll_old - 1e-12:
                break
            t /= 2
        b = bn
        if abs(ll - ll_old) < tol:
            conv = True
            break
        ll_old = ll
    mu = _sig(X @ b)
    w = np.clip(mu * (1 - mu), 1e-12, None)
    try:
        cov = np.linalg.inv((X * w[:, None]).T @ X)
    except np.linalg.LinAlgError:
        cov = np.full((p, p), np.nan)
    ll = np.sum(y * (X @ b) - np.logaddexp(0, X @ b))
    return b, ll, conv and np.all(np.abs(b) < 15), cov


def logit_firth(X, y, max_iter=200, tol=1e-9, fixed=None):
    """Firth-penalised logistic regression (Jeffreys prior). `fixed` = dict {col_index: value} to
    constrain coefficients (used for the penalised likelihood-ratio test).
    Returns (beta, penalised loglik, cov)."""
    n, p = X.shape
    free = np.array([j not in (fixed or {}) for j in range(p)])
    b = np.zeros(p)
    for j, v in (fixed or {}).items():
        b[j] = v

    def pll(beta):
        eta = X @ beta
        mu = _sig(eta)
        w = mu * (1 - mu)
        sign, logdet = np.linalg.slogdet((X * w[:, None]).T @ X)
        return np.sum(y * eta - np.logaddexp(0, eta)) + 0.5 * logdet

    ll_old = pll(b)
    for _ in range(max_iter):
        eta = X @ b
        mu = _sig(eta)
        w = np.clip(mu * (1 - mu), 1e-12, None)
        Xw = X * np.sqrt(w)[:, None]
        I = Xw.T @ Xw
        Iinv = np.linalg.pinv(I)
        h = np.einsum("ij,jk,ik->i", Xw, Iinv, Xw)            # hat diagonal
        U = X.T @ (y - mu + h * (0.5 - mu))                  # modified score
        If = I[np.ix_(free, free)]
        step = np.zeros(p)
        step[free] = np.linalg.solve(If, U[free])
        mx = np.max(np.abs(step))
        if mx > 5:
            step *= 5 / mx
        t = 1.0
        for _h in range(30):
            bn = b + t * step
            ll = pll(bn)
            if ll >= ll_old - 1e-12:
                break
            t /= 2
        b = bn
        if abs(ll - ll_old) < tol:
            ll_old = ll
            break
        ll_old = ll
    mu = _sig(X @ b)
    w = mu * (1 - mu)
    cov = np.linalg.pinv((X * w[:, None]).T @ X)
    return b, ll_old, cov


def lrt_ml(X0, X1, y):
    """Nested ML likelihood-ratio test. Returns (stat, df, p, fit1, converged both)."""
    b0, ll0, c0, _ = logit_ml(X0, y)
    b1, ll1, c1, cov1 = logit_ml(X1, y)
    df = X1.shape[1] - X0.shape[1]
    stat = max(0.0, 2 * (ll1 - ll0))
    return stat, df, stats.chi2.sf(stat, df), (b1, cov1), (c0 and c1)


def lrt_firth(X1, y, test_cols):
    """Penalised likelihood-ratio test for H0: beta[test_cols]=0 (Heinze & Schemper 2002)."""
    b1, pl1, cov1 = logit_firth(X1, y)
    _, pl0, _ = logit_firth(X1, y, fixed={j: 0.0 for j in test_cols})
    stat = max(0.0, 2 * (pl1 - pl0))
    return stat, stats.chi2.sf(stat, len(test_cols)), b1, cov1


def profile_ci_firth(X, y, j, level=0.95, span=8.0):
    """Profile penalised-likelihood CI for coefficient j (root finding)."""
    from scipy.optimize import brentq
    b, pl, cov = logit_firth(X, y)
    crit = stats.chi2.ppf(level, 1) / 2

    def f(v):
        _, plv, _ = logit_firth(X, y, fixed={j: v})
        return pl - plv - crit

    se = np.sqrt(cov[j, j]) if np.isfinite(cov[j, j]) else 1.0
    lo = hi = np.nan
    try:
        lo = brentq(f, b[j] - span * max(se, 0.5), b[j])
    except ValueError:
        pass
    try:
        hi = brentq(f, b[j], b[j] + span * max(se, 0.5))
    except ValueError:
        pass
    return b[j], lo, hi


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


if __name__ == "__main__":
    # validation 1: Firth on a 2x2 table == Haldane-corrected log OR (exact property of Firth's method)
    a, b_, c, d = 9, 1, 3, 12
    y = np.r_[np.ones(a), np.zeros(b_), np.ones(c), np.zeros(d)]
    x = np.r_[np.ones(a + b_), np.zeros(c + d)]
    X = np.c_[np.ones_like(x), x]
    bf, _, _ = logit_firth(X, y)
    print("Firth logOR", bf[1], "Haldane", np.log((a + .5) * (d + .5) / ((b_ + .5) * (c + .5))))
    # validation 2: ML vs statsmodels
    import statsmodels.api as sm
    rng = np.random.default_rng(1)
    Xr = np.c_[np.ones(200), rng.normal(size=(200, 2))]
    yr = (rng.random(200) < _sig(Xr @ [0.3, 0.8, -0.5])).astype(float)
    print("ML", logit_ml(Xr, yr)[0], "statsmodels", sm.Logit(yr, Xr).fit(disp=0).params)
    # validation 3: separation -> Firth finite
    ys = np.r_[np.ones(10), np.zeros(10)]
    xs = np.c_[np.ones(20), np.r_[np.ones(10), np.zeros(10)]]
    print("separated: ML conv", logit_ml(xs, ys)[2], "Firth", logit_firth(xs, ys)[0], "Haldane", np.log(10.5 * 10.5 / .25))
