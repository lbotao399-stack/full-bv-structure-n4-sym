"""QA: analysis of the third evaluation of Ob^(2)(Delta,xi)|_2 (q_dxi_raw.json).
(1) poles; (2) rational part: intrinsic fit on T0..T4, position part; (3) log part: fit of (intrinsic, position) data
onto Q of the degree -2, weight 9/2, position-degree-1 functionals of the census (Table at-census2, F_1^(2)(Delta,xi_+)):
  psi1 = int x^mu (Gamma(xi,lam)^mu, c*)            (lam_-, c*)_0
  psi2 = int x^mu (A*_mu, (xi,lam*))                (A*_V, lam*_+)_0
  psi3 = int x^mu (A*_nu, (rho_mu rho_nu xi,lam*)), nu != mu   (A*_V, lam*_+)_0, self-dual part
  psi4 = int x^mu (A*_I, (rho_I rho_mu xi, lam*))   (A*_perp, lam*_-)_0
using engine B's Q on densities (left-derivation rule, fixed by its G-sign and confirmed by engine A's mode model) and its
canonical-density rule.  theory = -(raw) (compare.json)."""
import sys, os, json
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
BSRC = os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b', 'src.nosync'))
sys.path.insert(0, BSRC)
import eb_qb as Q
QA = os.path.abspath(os.path.join(HERE, '..'))
CHIR = sys.argv[1] if len(sys.argv) > 1 else '+'
raw = json.load(open(os.path.join(QA, 'out.nosync', 'q_dxi_raw.json' if CHIR == '+' else 'q_dxi_raw_minus.json')))
MR, MP, R, RB = Q.MR, Q.MP, Q.R, Q.RB

def gi(z):  # CLog json -> dict key -> sympy Gaussian rational, negated (theory = -raw)
    out = {}
    for key in ('1', 'log2', 'log3'):
        re = sp.Rational(z['re'].get(key, '0')); im = sp.Rational(z['im'].get(key, '0'))
        out[key] = -(re + sp.I * im)
    return out

def rebuild(rec):
    legs0 = rec['legs']
    p = [Fr(x) for x in legs0[0]['p']]
    u0 = np.array(legs0[0]['u'], complex); u1 = np.array(legs0[1]['u'], complex)
    s0, s1 = legs0[0]['sp'], legs0[1]['sp']
    def legs_at(r):
        return [dict(sp=s0, u=u0, p=p), dict(sp=s1, u=u1, p=[-p[m] - r[m] for m in range(4)])]
    return legs_at, np.array(rec['xi'], complex)

def psis(xi):
    def p1(mu):
        return Q.dens(1, ('lam', 'cs'), lambda s: (xi @ MR[mu] @ s[0]['u']) * s[1]['u'][0])
    def p2(mu):
        return Q.dens(1, ('As', 'lams'), lambda s: s[0]['u'][mu] * (xi @ MP @ s[1]['u']))
    def p3(mu):
        return Q.dens(1, ('As', 'lams'), lambda s: sum(s[0]['u'][nu] * ((RB[mu] @ R[nu] @ xi) @ MP @ s[1]['u'])
                                                       for nu in range(4) if nu != mu))
    def p4(mu):
        return Q.dens(1, ('As', 'lams'), lambda s: sum(s[0]['u'][I] * ((RB[I] @ R[mu] @ xi) @ MP @ s[1]['u'])
                                                       for I in range(4, 10)))
    return [(None, [f(mu) for mu in range(4)]) for f in (p1, p2, p3, p4)]

def gz(t):
    return sp.Rational(t[0]) + sp.I * sp.Rational(t[1])

report = {'poles_zero': True, 'top_zero': True, 'rational': {}, 'rational_beta_zero': True}
rows_log, rhs_log, rows_log3 = [], [], []
for sec, recs in raw.items():
    rat_rows, rat_rhs = [], []
    labels = list(recs[0]['basis'].keys())
    for rec in recs:
        report['top_zero'] &= rec['top_zero']
        a0 = gi(rec['alpha']['k']['0']); a1 = gi(rec['alpha']['k']['1'])
        report['poles_zero'] &= all(v == 0 for v in a1.values())
        betas = {}
        for mu in range(4):
            b = rec['beta'][str(mu)]['k']
            report['poles_zero'] &= all(v == 0 for v in gi(b['1']).values())
            betas[mu] = gi(b['0'])
            if betas[mu]['1'] != 0:
                report['rational_beta_zero'] = False
        rat_rows.append([gz(rec['basis'][l]) for l in labels]); rat_rhs.append(a0['1'])
        legs_at, xi = rebuild(rec)
        cands = [Q.canon(Q.Qfn(fn), legs_at) for fn in psis(xi)]
        # intrinsic equation
        rows_log.append([gz(c[0]) for c in cands]); rhs_log.append(a0['log2']); rows_log3.append(a0['log3'])
        for mu in range(4):
            rows_log.append([gz(c[1][mu]) for c in cands]); rhs_log.append(betas[mu]['log2']); rows_log3.append(betas[mu]['log3'])
    M = sp.Matrix(rat_rows); bvec = sp.Matrix(rat_rhs)
    sol, params = M.gauss_jordan_solve(bvec)
    sol = sol.subs({p: 0 for p in params})
    resid = (M * sol - bvec).applyfunc(sp.simplify)
    report['rational'][sec] = {l: str(sp.nsimplify(x)) for l, x in zip(labels, sol)}
    report['rational'][sec]['consistent'] = all(r == 0 for r in resid)
# log part: every log entry proportional to log(4/3) = 2 log2 - log3  <=>  log3 coefficient = -(1/2) log2 coefficient
report['log_entries_prop_log4/3'] = all(sp.simplify(l3 + l2 / 2) == 0 for l2, l3 in zip(rhs_log, rows_log3))
M = sp.Matrix(rows_log); bvec = sp.Matrix(rhs_log)
sol, params = M.gauss_jordan_solve(bvec)
sol = sol.subs({p: 0 for p in params})
resid = (M * sol - bvec).applyfunc(sp.simplify)
report['log_part_exact'] = {'coefficients of log2 on Q psi1..psi4 (free params set 0)': [str(sp.nsimplify(x)) for x in sol],
                            'free parameters': len(params), 'consistent (residual zero)': all(r == 0 for r in resid),
                            'equations': len(rhs_log), 'rank': M.rank()}
print(json.dumps(report, indent=1))
json.dump(report, open(os.path.join(QA, 'q_dxi_fit.json' if CHIR == '+' else 'q_dxi_fit_minus.json'), 'w'), indent=1)
