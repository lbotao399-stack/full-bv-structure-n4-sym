"""Step 2, Part D4, engine A, task (1): internal-insertion bubbles Pi_int in all sectors and the one-loop self-energies,
fitted exactly on the invariant structures (several random polarisation samples, exact linear solve, residual checked).
Canonical form: alpha delta(Sigma) + beta^mu (d_mu delta)(Sigma), field leg first (momentum p, held fixed)."""
import os, sys, json, time, itertools
from fractions import Fraction
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import s2_core as C
from r2_core import Q2
import d4_two_point as TP
import d4_integrals as DI
from cliff import s as SIG, sb as SIGB

SB = np.array(SIGB); SS = np.array(SIG)
with np.errstate(all='ignore'):
    CHI = np.rint((SB[0] @ SS[1] @ SB[2] @ SS[3]).real)     # 4d chirality on Sigma^* (Sigma^* -> Sigma^*)
assert np.allclose(CHI @ CHI, np.eye(16))
EPS4 = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])
GR = C.G.real


def cdot(x, y):
    return np.einsum('...a,...b,ab->...', x, y, GR)


def structures(sector, a, b, p):
    """invariant structures for alpha (list) and beta^mu (list of 4-vectors) for legs field a (momentum p), antifield b"""
    p = np.asarray(p, dtype=float); p2 = float(p @ p)
    if sector in ('AsP,AP',):
        s = sum(cdot(a[M], b[M]) for M in range(4, 10))
        return [s], ['(a,b)'], [p * s], ['p^mu (a,b)']
    if sector == 'AsV,AV':
        ab = sum(cdot(a[m], b[m]) for m in range(4))
        pa = sum(p[m] * a[m] for m in range(4)); pb = sum(p[m] * b[m] for m in range(4))
        papb = cdot(pa, pb)
        al = [ab, papb / p2]; aln = ['(a.b)', '(p.a)(p.b)/p^2']
        be = [p * ab, np.array([cdot(a[m], pb) for m in range(4)]), np.array([cdot(b[m], pa) for m in range(4)]),
              p * papb / p2,
              np.array([sum(EPS4[mu, m, n, s] * cdot(a[m], b[n]) * p[s] for m in range(4) for n in range(4) for s in range(4))
                        for mu in range(4)])]
        ben = ['p^mu (a.b)', 'a^mu (p.b)', 'b^mu (p.a)', 'p^mu (p.a)(p.b)/p^2', 'eps^{mu m n s} (a_m,b_n) p_s']
        return al, aln, be, ben
    if sector == 'ls,l':
        uv = np.einsum('Aa,Ab,ab->', a, b, GR)
        uCv = np.einsum('Aa,AB,Bb,ab->', a, CHI, b, GR)
        al = [uv, uCv]; aln = ['(u,v)', '(u,chi v)']
        rp = np.einsum('m,mAB->AB', p, SS[:4])          # rho(p): Sigma^* -> Sigma
        be = [p * uv, p * uCv]
        ben = ['p^mu (u,v)', 'p^mu (u,chi v)']
        v1 = np.array([np.einsum('Aa,AB,Bb,ab->', a, SB[mu] @ rp, b, GR) for mu in range(4)])
        v2 = np.array([np.einsum('Aa,AB,Bb,ab->', a, SB[mu] @ rp @ CHI, b, GR) for mu in range(4)])
        be += [v1, v2]; ben += ['(u, rho(e^mu) rho(p) v)', '(u, rho(e^mu) rho(p) chi v)']
        return al, aln, be, ben
    if sector == 'cs,c':
        gw = cdot(a, b)
        return [gw], ['(gamma,w)'], [p * gw], ['p^mu (gamma,w)']
    raise ValueError(sector)


def tosym(z):
    return sp.Rational(int(z.re.numerator), int(z.re.denominator)) + sp.I * sp.Rational(int(z.im.numerator), int(z.im.denominator))


def fit(rows, rhs):
    M = sp.Matrix(rows); b = sp.Matrix(rhs)
    if M.shape[1] == 0:
        return [], all(x == 0 for x in rhs)
    rk = M.rank()
    if all(x == 0 for x in rhs):
        return [0] * M.shape[1], True
    sol, params = M.gauss_jordan_solve(b)
    if params.shape[0] > 0:
        sol = sol.subs({q: 0 for q in params})
    resid = (M * sol - b).applyfunc(sp.nsimplify)
    return [sp.nsimplify(x) for x in sol], all(sp.simplify(r) == 0 for r in resid), rk


def pol_pair(rng, sector):
    tf, wf, ta, wa = TP.SECT[sector]
    return C.rpol(rng, tf, wf), C.rpol(rng, ta, wa)


def run_sector(sector, nsamp, seed=1, p=(1, 2, 0, 1)):
    rng = np.random.default_rng(seed)
    samples = []
    for it in range(nsamp):
        a, b = pol_pair(rng, sector)
        t0 = time.time()
        alpha, beta, n0 = TP.pi_int(sector, p, a, b)
        al, aln, be, ben = structures(sector, a, b, p)
        samples.append((alpha, beta, al, be))
        print(sector, 'sample', it, 'nseq', n0, 'time %.1f' % (time.time() - t0), flush=True)
    keys = sorted(set(k for s_ in samples for k in s_[0]) | set(k for s_ in samples for bm in s_[1] for k in bm))
    out = {'alpha': {}, 'beta': {}, 'structures_alpha': aln, 'structures_beta': ben, 'momentum': list(p)}
    for key in keys:
        rows = []; rhs = []
        for (alpha, beta, al, be) in samples:
            rows.append([sp.nsimplify(complex(x).real) + sp.I * sp.nsimplify(complex(x).imag) for x in al])
            rhs.append(tosym(alpha.get(key, Q2(0))))
        sol = fit(rows, rhs)
        out['alpha'][key] = {'coeffs': [str(x) for x in sol[0]], 'exact_fit': bool(sol[1])}
        rows = []; rhs = []
        for (alpha, beta, al, be) in samples:
            for mu in range(4):
                rows.append([sp.nsimplify(complex(v[mu]).real) + sp.I * sp.nsimplify(complex(v[mu]).imag) for v in be])
                rhs.append(tosym(beta[mu].get(key, Q2(0))))
        sol = fit(rows, rhs)
        out['beta'][key] = {'coeffs': [str(x) for x in sol[0]], 'exact_fit': bool(sol[1])}
    return out


if __name__ == '__main__':
    sector = sys.argv[1]
    nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    pstr = sys.argv[3] if len(sys.argv) > 3 else '1,2,0,1'
    p = tuple(int(x) for x in pstr.split(','))
    res = run_sector(sector, nsamp, p=p)
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'logs', 'd4_pi_%s_%s_p%s.json' % (sector.replace(',', '_'), C.ALG,
                                                                             pstr.replace(',', ''))), 'w'), indent=1)
