"""Exact fits of the B1 results (flowed arity-one defect at two legs) to local structures, per sector.
Canonical form: alpha delta + beta^mu d/dp2 delta, leg 1 at p (kept), leg 2 eliminated."""
import json, sys, glob
import sympy as sp
import delta_engine as D
from realization import R, Rb

V = D.V
def S(x):
    return sp.sympify(x)

def num(d):
    """{'1':..,'L2':..,'L3':..} -> (rational part, dict of log parts)"""
    return S(d.get('1', 0)), {k: S(v) for k, v in d.items() if k != '1' and S(v) != 0}

def load(pattern):
    recs = []
    for f in glob.glob(pattern):
        d = json.load(open(f))
        for sec, lst in d.items():
            for r in lst:
                recs.append((sec, r))
    return recs

def eps_tensor(a, b, c, mu):
    return sum(sp.LeviCivita(mu, i, j, k) * a[i] * b[j] * c[k] for i in range(4) for j in range(4) for k in range(4))

def fit_sector(sec, recs):
    out = {}
    if sec in ('ApAp', 'AVAV'):
        cs = list(range(4, 10)) if sec == 'ApAp' else list(range(4))
        names_a0 = ['int(dA,dA)', 'int(d*A,d*A)'] if sec == 'AVAV' else ['int(dA,dA)']
        names_am1 = ['int(A,A)']
        names_b0 = ['p^mu p^2 (u1,u2)', 'p^mu (p,u1)(p,u2)', 'u1^mu p^2 (p,u2)', 'u2^mu p^2 (p,u1)', 'eps^mu(u1,u2,p) p^2']
        names_bm1 = ['p^mu (u1,u2)', 'u1^mu (p,u2)', 'u2^mu (p,u1)', 'eps^mu(u1,u2,p)']
    elif sec == 'cAsV':
        names_a0 = ['int(dc,A*)']; names_am1 = []
        names_b0 = ['p^mu (p,a) c', 'a^mu p^2 c']; names_bm1 = ['a^mu c']
    else:
        names_a0 = ['int(lambda,dslash lambda)']; names_am1 = []
        names_b0 = []; names_bm1 = []
    unk = {n: sp.symbols(f'x0_0:{len(names_a0)}') for n in ['a0']}
    xa0 = sp.symbols(f'xa0_0:{max(1,len(names_a0))}'); xam1 = sp.symbols(f'xam1_0:{max(1,len(names_am1))}')
    xb0 = sp.symbols(f'xb0_0:{max(1,len(names_b0))}'); xbm1 = sp.symbols(f'xbm1_0:{max(1,len(names_bm1))}')
    eqs = {'a0': [], 'am1': [], 'b0': [], 'bm1': []}
    logs = []
    for s, r in recs:
        if s != sec:
            continue
        u1 = [S(x) for x in r['u1']]; u2 = [S(x) for x in r['u2']]; p = [S(x) for x in r['p1']]
        p2 = sum(x * x for x in p)
        res = r['res']
        a0, l0 = num(res['0']['alpha']); am1, lm1 = num(res['-1']['alpha'])
        b0 = [num(b) for b in res['0']['beta']]; bm1 = [num(b) for b in res['-1']['beta']]
        logs.append({'a0': {k: str(v) for k, v in l0.items()}, 'am1': {k: str(v) for k, v in lm1.items()},
                     'b0': [{k: str(v) for k, v in b[1].items()} for b in b0], 'bm1': [{k: str(v) for k, v in b[1].items()} for b in bm1]})
        legs = [dict(sp=('A' if sec in ('ApAp', 'AVAV') else ('c' if sec == 'cAsV' else 'lam')), u=u1, p=p),
                dict(sp=('A' if sec in ('ApAp', 'AVAV') else ('As' if sec == 'cAsV' else 'lam')), u=u2, p=[-x for x in p])]
        if sec in ('ApAp', 'AVAV'):
            uu = sum(a * b for a, b in zip(u1, u2)); pu1 = sum(p[m] * u1[m] for m in range(4)); pu2 = sum(p[m] * u2[m] for m in range(4))
            ref_a0 = [2 * (p2 * uu - pu1 * pu2), 2 * pu1 * pu2] if sec == 'AVAV' else [2 * p2 * uu]
            ref_am1 = [2 * uu]
            ref_b0 = lambda m: [p[m] * p2 * uu, p[m] * pu1 * pu2, (u1[m] if m < 4 else 0) * p2 * pu2, (u2[m] if m < 4 else 0) * p2 * pu1,
                                eps_tensor(u1, u2, p, m) * p2]
            ref_bm1 = lambda m: [p[m] * uu, (u1[m] if m < 4 else 0) * pu2, (u2[m] if m < 4 else 0) * pu1, eps_tensor(u1, u2, p, m)]
        elif sec == 'cAsV':
            pa = sum(p[m] * u2[m] for m in range(4)); c = u1[0]
            ref_a0 = [V.vertex('kingh', legs)['k']]; ref_am1 = []
            ref_b0 = lambda m: [p[m] * pa * c, u2[m] * p2 * c]
            ref_bm1 = lambda m: [u2[m] * c]
        else:
            u1p = [u1[i] if i < 8 else 0 for i in range(16)]; u1m = [u1[i] - u1p[i] for i in range(16)]
            u2p = [u2[i] if i < 8 else 0 for i in range(16)]; u2m = [u2[i] - u2p[i] for i in range(16)]
            def kl(a, b):
                return 2 * V.vertex('kinlam', [dict(sp='lam', u=a, p=p), dict(sp='lam', u=b, p=[-x for x in p])]).get('k', 0)
            ref_a0 = [kl(u1, u2)]
            ref_am1 = []
            ref_b0 = lambda m: []; ref_bm1 = lambda m: []
        eqs['a0'].append(a0 - sum(x * y for x, y in zip(xa0, ref_a0)))
        if names_am1:
            eqs['am1'].append(am1 - sum(x * y for x, y in zip(xam1, ref_am1)))
        else:
            eqs['am1'].append(am1)
        for m in range(4):
            if names_b0:
                eqs['b0'].append(b0[m][0] - sum(x * y for x, y in zip(xb0, ref_b0(m))))
            else:
                eqs['b0'].append(None)
            if names_bm1:
                eqs['bm1'].append(bm1[m][0] - sum(x * y for x, y in zip(xbm1, ref_bm1(m))))
            else:
                eqs['bm1'].append(bm1[m][0])
    def solve(eq, xs, names):
        eq = [e for e in eq if e is not None]
        if not names:
            return {'all_zero': all(sp.simplify(e) == 0 for e in eq)}
        sol = sp.solve(eq, list(xs[:len(names)]), dict=True)
        if not sol:
            return {'no_solution': True}
        s0 = sol[0]
        return {names[i]: str(s0.get(xs[i], 'free')) for i in range(len(names))}
    out['alpha eps^0'] = solve(eqs['a0'], xa0, names_a0)
    out['alpha eps^-1'] = solve(eqs['am1'], xam1, names_am1)
    out['beta eps^0'] = solve(eqs['b0'], xb0, names_b0) if names_b0 else {'raw': 'not fitted'}
    out['beta eps^-1'] = solve(eqs['bm1'], xbm1, names_bm1)
    out['n_samples'] = len([1 for s, r in recs if s == sec])
    out['log_constant_parts'] = logs
    return out

if __name__ == '__main__':
    recs = load('delta_B1s_*_s*.json') + load('delta_B1_ApAp_s23.json')
    res = {}
    for sec in ('ApAp', 'AVAV', 'cAsV', 'lamlam'):
        if any(s == sec for s, r in recs):
            res[sec] = fit_sector(sec, recs)
    json.dump(res, open('delta_fit_B1.json', 'w'), indent=1)
    print(json.dumps(res, indent=1))
