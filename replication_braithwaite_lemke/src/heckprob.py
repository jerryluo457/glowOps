"""Probit with sample selection (censored / Heckman probit), equivalent to Stata's `heckprob`.

Selection:  s* = Z g + u,   s = 1[s* > 0]
Outcome:    y* = X b + e,   y = 1[y* > 0], observed only when s = 1
corr(u, e) = rho,  rho = tanh(athrho).

Log-likelihood
  s = 0:  log Phi(-Zg)
  s = 1:  log Phi2(q Xb, Zg, q rho),  q = 2y - 1

Standard errors: cluster-robust sandwich with Stata's G/(G-1) small-sample factor.
"""
from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats
from scipy.special import log_ndtr, ndtr

# 20-point Gauss-Legendre nodes/weights on [-1, 1].
_GL_X, _GL_W = np.polynomial.legendre.leggauss(20)


def bvn_cdf(a, b, r):
    """P(X < a, Y < b) for standard bivariate normal with correlation r (scalar).

    Uses Phi2 = Phi(a)Phi(b) + (1/2pi) * int_0^{asin r} exp(-(a^2+b^2-2ab sin t)/(2cos^2 t)) dt
    (Genz 2004, low-|r| branch) with 20-point Gauss-Legendre; accurate to ~1e-14 for |r| < 0.925.
    Falls back to scipy for larger |r|.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if abs(r) < 0.925:
        if r == 0:
            return ndtr(a) * ndtr(b)
        asr = np.arcsin(r)
        t = 0.5 * asr * (_GL_X + 1.0)                        # nodes on [0, asr]
        sn = np.sin(t)
        hs = (a * a + b * b)[:, None] / 2.0
        hk = (a * b)[:, None]
        integrand = np.exp((sn * hk - hs) / (1.0 - sn * sn))
        val = integrand @ (_GL_W * 0.5 * asr) / (2 * np.pi)
        return np.clip(ndtr(a) * ndtr(b) + val, 1e-300, 1.0)
    cov = np.array([[1.0, r], [r, 1.0]])
    pts = np.column_stack([a, b])
    return np.clip(stats.multivariate_normal(cov=cov).cdf(pts), 1e-300, 1.0)


def _phi(x):
    return np.exp(-0.5 * x * x) / np.sqrt(2 * np.pi)


@dataclass
class HeckprobResult:
    params: np.ndarray
    vcov: np.ndarray
    names: list
    k_out: int
    k_sel: int
    loglik: float
    n: int
    n_selected: int
    n_clusters: int
    converged: bool

    @property
    def se(self):
        return np.sqrt(np.diag(self.vcov))

    @property
    def rho(self):
        return np.tanh(self.params[-1])

    def wald_outcome(self, exclude_const=True):
        """Wald chi2 that all outcome-equation slopes are zero (Stata's model chi2 for heckprob)."""
        idx = [i for i in range(self.k_out) if not (exclude_const and self.names[i].endswith(":_cons"))]
        b = self.params[idx]
        V = self.vcov[np.ix_(idx, idx)]
        w = float(b @ np.linalg.solve(V, b))
        return w, len(idx), stats.chi2.sf(w, len(idx))

    def table(self):
        z = self.params / self.se
        p = 2 * stats.norm.sf(np.abs(z))
        import pandas as pd
        return pd.DataFrame({"coef": self.params, "se": self.se, "z": z, "p": p}, index=self.names)


class Heckprob:
    def __init__(self, y, X, s, Z, cluster=None, x_names=None, z_names=None):
        self.s = np.asarray(s, int)
        sel = self.s == 1
        self.X1 = np.asarray(X, float)[sel]
        self.y1 = np.asarray(y, float)[sel]
        self.Z = np.asarray(Z, float)
        self.Z1, self.Z0 = self.Z[sel], self.Z[~sel]
        self.sel = sel
        self.kx, self.kz = self.X1.shape[1], self.Z.shape[1]
        self.cluster = None if cluster is None else np.asarray(cluster)
        self.x_names = x_names or [f"x{i}" for i in range(self.kx)]
        self.z_names = z_names or [f"z{i}" for i in range(self.kz)]

    def _unpack(self, theta):
        return theta[: self.kx], theta[self.kx: self.kx + self.kz], theta[-1]

    def obs_loglik_and_scores(self, theta):
        """Per-observation log-likelihood and score matrix (n x k)."""
        b, g, ath = self._unpack(theta)
        rho = np.tanh(ath)
        n = len(self.s)
        ll = np.empty(n)
        S = np.zeros((n, len(theta)))

        # Non-selected: log Phi(-Zg)
        zg0 = self.Z0 @ g
        ll[~self.sel] = log_ndtr(-zg0)
        mills = _phi(zg0) / np.clip(ndtr(-zg0), 1e-300, None)
        S[np.ix_(~self.sel, range(self.kx, self.kx + self.kz))] = -mills[:, None] * self.Z0

        # Selected: log Phi2(q xb, zg, q rho), evaluated separately for y=1 and y=0 (r = +/- rho)
        q = 2 * self.y1 - 1
        xb, zg = self.X1 @ b, self.Z1 @ g
        sel_idx = np.flatnonzero(self.sel)
        ll1 = np.empty(len(q))
        da, db, dr = np.empty(len(q)), np.empty(len(q)), np.empty(len(q))
        for sign in (1, -1):
            m = q == sign
            if not m.any():
                continue
            r = sign * rho
            a, c = sign * xb[m], zg[m]
            P = bvn_cdf(a, c, r)
            s1r = np.sqrt(1 - r * r)
            ll1[m] = np.log(P)
            da[m] = _phi(a) * ndtr((c - r * a) / s1r) / P
            db[m] = _phi(c) * ndtr((a - r * c) / s1r) / P
            dens = np.exp(-(a * a - 2 * r * a * c + c * c) / (2 * (1 - r * r))) / (2 * np.pi * s1r)
            dr[m] = dens / P
        ll[sel_idx] = ll1
        S[sel_idx, : self.kx] = (q * da)[:, None] * self.X1
        S[sel_idx, self.kx: self.kx + self.kz] = db[:, None] * self.Z1
        S[sel_idx, -1] = q * dr * (1 - rho * rho)
        return ll, S

    def negll(self, theta):
        ll, S = self.obs_loglik_and_scores(theta)
        return -ll.sum(), -S.sum(axis=0)

    def start_values(self):
        import statsmodels.api as sm
        g = sm.Probit(self.s, self.Z).fit(disp=0).params
        b = sm.Probit(self.y1, self.X1).fit(disp=0).params
        return np.concatenate([b, g, [0.0]])

    def fit(self, start=None, maxiter=2000):
        theta0 = self.start_values() if start is None else start
        res = optimize.minimize(self.negll, theta0, jac=True, method="BFGS",
                                options={"maxiter": maxiter, "gtol": 1e-7})
        # Polish with Newton steps using a numerical Hessian of the analytic gradient.
        theta = res.x
        for _ in range(20):
            H = self.hessian(theta)
            g = self.negll(theta)[1]
            step = np.linalg.solve(H, g)
            theta = theta - step
            if np.max(np.abs(step)) < 1e-9:
                break
        H = self.hessian(theta)
        ll, S = self.obs_loglik_and_scores(theta)
        Hinv = np.linalg.inv(H)
        if self.cluster is not None:
            import pandas as pd
            G = pd.DataFrame(S).groupby(self.cluster).sum().to_numpy()
            nc = G.shape[0]
            meat = G.T @ G * nc / (nc - 1)
        else:
            nc = len(S)
            meat = S.T @ S * nc / (nc - 1)
        V = Hinv @ meat @ Hinv
        names = [f"outcome:{n}" for n in self.x_names] + [f"select:{n}" for n in self.z_names] + ["athrho"]
        # Convergence criterion in the spirit of Stata's scaled gradient: g' H^-1 g.
        grad = S.sum(axis=0)
        scaled = float(grad @ Hinv @ grad)
        return HeckprobResult(theta, V, names, self.kx, self.kz, ll.sum(), len(self.s),
                              int(self.sel.sum()), nc, bool(scaled < 1e-8))

    def hessian(self, theta, eps=1e-6):
        k = len(theta)
        H = np.empty((k, k))
        for j in range(k):
            e = np.zeros(k)
            e[j] = eps * max(1.0, abs(theta[j]))
            H[:, j] = (self.negll(theta + e)[1] - self.negll(theta - e)[1]) / (2 * e[j])
        return (H + H.T) / 2
