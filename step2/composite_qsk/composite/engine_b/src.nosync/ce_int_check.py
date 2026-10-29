"""Independent check of ce_int.J for n_t = 4, 5: quasi-Monte Carlo on the alpha simplex of
int alpha^n W_m(alpha), W_0 = log(min_f(alpha_f/c_f)/alpha_sigma)_+, W_m = (1/m)[min_f(alpha_f/c_f)^m - alpha_sigma^m]_+
(the tau-integral done first), compared with the Feynman-parameter evaluation of ce_int.J."""
import numpy as np, math, sys
from scipy.stats import qmc
import ce_int as CI

def J_qmc(m, a, nvec, cut, npow=22, seed=0):
    nt = len(nvec) + 1
    s = qmc.Sobol(d=nt - 1, scramble=True, seed=seed)
    x = s.random_base2(m=npow)
    # uniform on simplex via sorted uniforms (density (nt-1)! on the simplex)
    xs = np.sort(x, axis=1)
    al = np.diff(np.concatenate([np.zeros((len(xs), 1)), xs, np.ones((len(xs), 1))], axis=1))
    asg = al[:, 0]; af = al[:, 1:]
    mn = np.min(af / np.array(cut, dtype=float), axis=1)
    val = asg ** a
    for j, n in enumerate(nvec):
        val = val * af[:, j] ** n
    with np.errstate(divide='ignore', invalid='ignore'):
        if m == 0:
            W = np.where(mn > asg, np.log(mn / asg), 0.0)
        else:
            W = np.where(mn > asg, (mn ** m - asg ** m) / m, 0.0)
    W = np.nan_to_num(W)
    return float(np.mean(val * W)) / math.factorial(nt - 1)

if __name__ == '__main__':
    out = []
    for (m, a, nvec, cut) in [(0, 0, (0, 0), (1, 1)), (0, 0, (0, 0, 0), (1, 1, 1)), (0, 2, (1, 0, 1), (1, 1, 1)),
                              (1, 1, (0, 2, 0), (1, 1, 1)), (0, 1, (1, 0), (2, 1)), (1, 0, (0, 1, 1), (2, 1, 1)),
                              (0, 1, (2, 1, 0, 1), (1, 1, 1, 1)), (2, 0, (1, 1, 0), (1, 1, 1)), (0, 3, (0, 0, 2), (1, 2, 1))]:
        ex = CI.J(m, a, nvec, cut)
        q = [J_qmc(m, a, nvec, cut, seed=s) for s in range(3)]
        mq = sum(q) / 3; sd = (max(q) - min(q))
        print(m, a, nvec, cut, ex, '%.10g' % float(ex), '%.10g' % mq, 'rel diff %.2e (qmc spread %.1e)' % (abs(float(ex) - mq) / abs(mq), sd / abs(mq)), flush=True)
