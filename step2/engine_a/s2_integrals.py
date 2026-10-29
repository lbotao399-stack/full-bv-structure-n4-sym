"""Step 2, engine A: exact Schwinger integrals of the special-edge graphs.

J(n; b; N) = int_0^1 dtau_1 int_{[1,infty)^{n-1}} dtau_2..dtau_n  tau_2^{b_2} ... tau_n^{b_n} (tau_1+...+tau_n)^{-N},
n = number of edges (special edge tau_1 in [0,1], the others in [1,infty)), b = (b_2,...,b_n).
Convergence: the integrand is at most T^{|b|-N}; for N - |b| = m the region tau_1 <= 1, tau_i >= 1 (i >= 2) gives
int (tau_2+...+tau_n)^{-m} over [1,infty)^{n-1}, finite iff m > n-1; T >= n-1 >= 1 on the region (no small-T singularity).
Exact evaluation: iterated integration of rational functions (sympy), innermost tau_n first.  Second, independent
evaluation: mpmath at 50 digits (iterated tanh-sinh quadrature), and an integer relation (PSLQ) on the basis
{1, log 2, log 3}."""
import json, os, sys
import sympy as sp
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'logs', 's2_integrals.json')
t1, t2, t3 = sp.symbols('t1 t2 t3', positive=True)
L2, L3 = sp.log(2), sp.log(3)


def split_logs(r):
    def rep(x):
        x = sp.nsimplify(x)
        if x.is_Rational and x > 0:
            out = 0
            for pr, e in sp.factorint(x.p).items():
                out += e * sp.log(pr)
            for pr, e in sp.factorint(x.q).items():
                out -= e * sp.log(pr)
            return out
        return sp.log(x)
    return sp.expand(r.replace(sp.log, rep))


def exact(n, b, N):
    assert N - sum(b) > n - 1, ('divergent', n, b, N)
    if n == 2:
        f = t2 ** b[0] * (t1 + t2) ** (-N)
        g = sp.integrate(f, (t2, 1, sp.oo))
        r = sp.integrate(sp.apart(sp.simplify(g), t1), (t1, 0, 1))
    elif n == 3:
        f = t2 ** b[0] * t3 ** b[1] * (t1 + t2 + t3) ** (-N)
        g = sp.integrate(f, (t3, 1, sp.oo))
        g = sp.apart(sp.simplify(g), t2)
        h = sp.integrate(g, (t2, 1, sp.oo))
        h = sp.apart(sp.simplify(h), t1)
        r = sp.integrate(h, (t1, 0, 1))
    else:
        raise ValueError(n)
    r = split_logs(sp.expand(sp.expand_log(sp.simplify(r), force=True)))
    c0 = r.subs({L2: 0, L3: 0}); c2 = r.coeff(L2); c3 = r.coeff(L3)
    chk = sp.simplify(r - (c0 + c2 * L2 + c3 * L3))
    assert chk == 0 and c0.is_rational and c2.is_rational and c3.is_rational, (n, b, N, r)
    return (sp.Rational(c0), sp.Rational(c2), sp.Rational(c3))


def numeric(n, b, N, dps=50):
    """independent evaluation in the variables (T, u): tau_i = T u_i, d^n tau = T^{n-1} dT d^{n-1}u; the region
    tau_1 <= 1, tau_{i>=2} >= 1 is T in [1/min_{i>=2} u_i, 1/u_1] with u_1 <= min_{i>=2} u_i; the T-integral of a power is
    elementary, the u-integral is done by mpmath (tanh-sinh, split along the kinks) at dps digits"""
    mp.mp.dps = dps + 10
    def Tpow(e, lo, hi):          # int_lo^hi T^e dT
        if e == -1:
            return mp.log(hi / lo)
        return (hi ** (e + 1) - lo ** (e + 1)) / (e + 1)
    if n == 2:
        e = 1 + b[0] - N
        f = lambda u1: (1 - u1) ** b[0] * Tpow(e, 1 / (1 - u1), 1 / u1)
        val = mp.quad(f, [0, mp.mpf(1) / 4, mp.mpf(1) / 2])
    else:
        e = 2 + b[0] + b[1] - N
        def f(u2, u1):
            u3 = 1 - u2 - u1
            return u2 ** b[0] * u3 ** b[1] * Tpow(e, 1 / min(u2, u3), 1 / u1)
        # inner variable u1 (u3 = 1 - u2 - u1); piece u2 <= u3: u2 in [0,1/2], u1 in [0, min(u2, 1-2u2)];
        # piece u3 < u2: u2 in [1/3,1], u1 in [max(1-2u2,0), (1-u2)/2]
        p1 = mp.quad(lambda u2: mp.quad(lambda u1: f(u2, u1), [0, min(u2, 1 - 2 * u2)]), [0, mp.mpf(1) / 3, mp.mpf(1) / 2])
        p2 = mp.quad(lambda u2: mp.quad(lambda u1: f(u2, u1), [max(1 - 2 * u2, 0), (1 - u2) / 2]),
                     [mp.mpf(1) / 3, mp.mpf(1) / 2, 1])
        val = p1 + p2
    mp.mp.dps = dps
    return +val


_cache = None
def J(n, b, N):
    global _cache
    if _cache is None:
        _cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    key = '%d|%s|%d' % (n, ','.join(map(str, b)), N)
    if key not in _cache:
        c = exact(n, tuple(b), N)
        _cache[key] = {'exact': [str(x) for x in c]}
        json.dump(_cache, open(CACHE, 'w'), indent=1)
    return tuple(sp.Rational(x) for x in _cache[key]['exact'])


def verify_all(dps=50):
    """second independent evaluation of every cached integral, and PSLQ identification"""
    global _cache
    J(2, (0,), 2)
    out = {}
    for key, v in _cache.items():
        n, b, N = key.split('|'); n = int(n); N = int(N); b = tuple(int(x) for x in b.split(','))
        val = numeric(n, b, N, dps)
        c0, c2, c3 = (sp.Rational(x) for x in v['exact'])
        ex = mp.mpf(c0.p) / c0.q + mp.mpf(c2.p) / c2.q * mp.log(2) + mp.mpf(c3.p) / c3.q * mp.log(3)
        mp.mp.dps = dps
        rel = mp.pslq([val, 1, mp.log(2), mp.log(3)], maxcoeff=10 ** 8, maxsteps=10 ** 6)
        v['numeric_%ddigits' % dps] = mp.nstr(val, dps - 5)
        v['abs_diff'] = mp.nstr(abs(val - ex), 3)
        v['pslq_[val,1,log2,log3]'] = None if rel is None else [int(x) for x in rel]
        out[key] = v
    json.dump(_cache, open(CACHE, 'w'), indent=1)
    return out

if __name__ == '__main__':
    for (n, b, N) in [(2, (0,), 2), (2, (0,), 3), (3, (0, 0), 3), (3, (0, 0), 4)]:
        print(n, b, N, J(n, b, N), flush=True)
    if len(sys.argv) > 1 and sys.argv[1] == 'verify':
        for k, v in verify_all().items():
            print(k, v, flush=True)
