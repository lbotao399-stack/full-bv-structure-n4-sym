"""QA: fit of the eta_+ intrinsic two-leg obstruction (q_eta_raw.json) on the bases of engine A's task T6, realised in
engine B's conventions: b1 = int(eta,lam) d*A_V, b2 = int(eta,rho(dA_V)lam), b3 = int(eta,rho(dA_perp)lam),
b4 = int(rho(dc)eta,lam*); (eta,lam) = lam^T MP eta.  Theory = raw (odd generator)."""
import sys, os, json
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b', 'src.nosync')))
import eb_qb as Q
QA = os.path.abspath(os.path.join(HERE, '..'))
raw = json.load(open(os.path.join(QA, 'out.nosync', 'q_eta_raw.json')))
MP, RB = Q.MP, Q.RB
def bases(eta):
    b1 = Q.dens(1, ('A', 'lam'), lambda s: (-np.dot(1j * Q.pv(s[0]), s[0]['u'][:4])) * (s[1]['u'] @ MP @ eta))
    b2 = Q.dens(1, ('A', 'lam'), lambda s: Q.rhoF(Q.dA(s[0], range(4)), s[1]['u']) @ MP @ eta)
    b3 = Q.dens(1, ('A', 'lam'), lambda s: Q.rhoF(Q.dA(s[0], range(4, 10)), s[1]['u']) @ MP @ eta)
    b4 = Q.dens(1, ('c', 'lams'), lambda s: sum(1j * Q.pv(s[0])[m] * s[0]['u'][0] * ((RB[m] @ eta) @ MP @ s[1]['u']) for m in range(4)))
    return dict(b1=b1, b2=b2, b3=b3, b4=b4)
LAB = {'AV_lam+': ['b1', 'b2'], 'Ap_lam-': ['b3'], 'c_lams-': ['b4']}
def gi(z, key):
    return sp.Rational(z['re'].get(key, '0')) + sp.I * sp.Rational(z['im'].get(key, '0'))
out = {'sigma(s0)': raw.get('sigma'), 'closedness residuals': raw.get('closedness')}
for sec, recs in raw['sectors'].items():
    labs = LAB[sec]
    rows, rhs = [], {k: [] for k in ('1', 'log2', 'log3')}
    poles, betanz = [], False
    for rec in recs:
        p = [Fr(x) for x in rec['p']]
        sp0, sp1 = {'AV_lam+': ('A', 'lam'), 'Ap_lam-': ('A', 'lam'), 'c_lams-': ('c', 'lams')}[sec]
        legs = [dict(sp=sp0, u=np.array(rec['u0'], complex), p=p), dict(sp=sp1, u=np.array(rec['u1'], complex), p=[-x for x in p])]
        B = bases(np.array(rec['eta'], complex))
        rows.append([sp.Rational(str(Q.to_fr(B[l](legs))[0])) + sp.I * sp.Rational(str(Q.to_fr(B[l](legs))[1])) for l in labs])
        a = rec['alpha']['k']
        for k in rhs:
            rhs[k].append(gi(a['0'], k))
        poles.append({k: str(gi(a['1'], k)) for k in ('1', 'log2', 'log3')})
        betanz |= any(gi(rec['beta'][mu]['k']['0'], k) != 0 for mu in rec['beta'] for k in ('1', 'log2', 'log3'))
    M = sp.Matrix(rows)
    res = {}
    for k in rhs:
        try:
            sol, params = M.gauss_jordan_solve(sp.Matrix(rhs[k]))
            sol = sol.subs({q: 0 for q in params})
            res[k] = [str(sp.nsimplify(x)) for x in sol]
        except ValueError:
            res[k] = 'inconsistent'
    out[sec] = dict(bases=labs, coefficients=res, poles_first_config=poles[0], position_part_nonzero=betanz)
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(QA, 'q_eta_fit.json'), 'w'), indent=1)
