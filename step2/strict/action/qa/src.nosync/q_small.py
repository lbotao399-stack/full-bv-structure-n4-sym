"""QA: numerical check of Prop at-small (ii): on the simplex of n_t Schwinger fractions (Lebesgue measure on the first
n_t-1 coordinates, volume 1/(n_t-1)!), region alpha_s smallest:
  int log(min_{f!=s} a_f / a_s) = log(n/(n-1)) / (n-1)!,   int (min_{f!=s} a_f - a_s) = 1/(n^2 (n-1) (n-1)!)."""
import numpy as np, math
rng = np.random.default_rng(1)
out = {}
for n in (2, 3, 4, 5):
    N = 4_000_000
    a = rng.dirichlet(np.ones(n), size=N)          # uniform on the simplex
    vol = 1 / math.factorial(n - 1)
    s = a[:, 0]; m = a[:, 1:].min(axis=1)
    reg = s <= m
    I0 = vol * np.mean(np.where(reg, np.log(np.where(reg, m / s, 1.0)), 0.0))
    I1 = vol * np.mean(np.where(reg, m - s, 0.0))
    out[n] = (I0, math.log(n / (n - 1)) / math.factorial(n - 1), I1, 1 / (n * n * (n - 1) * math.factorial(n - 1)))
    print(n, 'log-int MC %.5f  claimed %.5f | pole-int MC %.6f claimed %.6f' % out[n])
