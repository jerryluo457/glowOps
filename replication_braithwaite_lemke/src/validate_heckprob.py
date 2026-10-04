"""Checks for the heckprob implementation: CDF accuracy, gradient, and parameter recovery."""
import sys
from pathlib import Path

import numpy as np
from scipy import optimize, stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from heckprob import Heckprob, bvn_cdf  # noqa: E402

rng = np.random.default_rng(620)

# 1. Bivariate normal CDF against scipy.
worst = 0.0
for r in (-0.9, -0.6, -0.3, 0.05, 0.4, 0.8, 0.95):
    a, b = rng.normal(0, 1.5, 300), rng.normal(0, 1.5, 300)
    ref = stats.multivariate_normal(cov=[[1, r], [r, 1]]).cdf(np.column_stack([a, b]))
    worst = max(worst, np.max(np.abs(bvn_cdf(a, b, r) - ref)))
print(f"bvn_cdf max abs error vs scipy: {worst:.2e}")
assert worst < 1e-5

# 2. Simulated selection model; check gradient and recovery.
n = 60000
Z = np.column_stack([np.ones(n), rng.normal(size=n), rng.binomial(1, 0.3, n), rng.normal(size=n)])
X = Z[:, :3]
g_true = np.array([-1.2, 0.5, 0.8, 0.7])          # last column is the exclusion restriction
b_true = np.array([0.3, -0.4, 0.6])
rho_true = -0.4
u, e = rng.multivariate_normal([0, 0], [[1, rho_true], [rho_true, 1]], n).T
s = (Z @ g_true + u > 0).astype(int)
y = (X @ b_true + e > 0).astype(int)
cl = rng.integers(0, 3000, n)

m = Heckprob(y, X, s, Z, cluster=cl)
theta = np.concatenate([b_true, g_true, [np.arctanh(rho_true)]]) + 0.05
num = optimize.approx_fprime(theta, lambda t: m.negll(t)[0], 1e-6)
print(f"max |analytic - numeric gradient|: {np.max(np.abs(num - m.negll(theta)[1])):.2e}")

res = m.fit()
print(f"converged={res.converged}  selected={res.n_selected}  rho_hat={res.rho:.3f} (true {rho_true})")
truth = np.concatenate([b_true, g_true, [np.arctanh(rho_true)]])
z = (res.params - truth) / res.se
for nm, est, se_, t, zz in zip(res.names, res.params, res.se, truth, z):
    print(f"  {nm:12s} est {est:7.3f}  se {se_:.3f}  true {t:6.3f}  z {zz:5.2f}")
assert np.all(np.abs(z) < 4), "parameter recovery failed"
print("OK")
