"""Fit of K5 (A|_Wperp, A|_Wperp) results: alpha = c1 T1 + c2 T2, beta^mu = (d1 p1 + d2 k)^mu T1 + (d3 p1 + d4 k)^mu T2,
T1 = (w1,u1)(w2,u2) + (w1,u2)(w2,u1) (vertex value of O), T2 = (w1,w2)(u1,u2); eps^-1 parts must vanish."""
import json, sys, glob
import sympy as sp

def load(files):
    recs = []
    for f in files:
        d = json.load(open(f))
        for sec, lst in d.items():
            recs += lst
    return recs

def fit(recs):
    c1, c2, d1, d2, d3, d4 = sp.symbols('c1 c2 d1 d2 d3 d4')
    eqs = []; checks = []
    for r in recs:
        w1 = [sp.sympify(x) for x in r['w1']]; w2 = [sp.sympify(x) for x in r['w2']]
        u1 = [sp.sympify(x) for x in r['u1']]; u2 = [sp.sympify(x) for x in r['u2']]
        p1 = [sp.sympify(x) for x in r['p1']]; k = [sp.sympify(x) for x in r['k']]
        dot = lambda a, b: sum(x * y for x, y in zip(a, b))
        T1 = dot(w1, u1) * dot(w2, u2) + dot(w1, u2) * dot(w2, u1); T2 = dot(w1, w2) * dot(u1, u2)
        res = r['res']
        for n in ('-1',):
            checks.append(all(sp.sympify(v) == 0 for v in res[n]['alpha'].values()) and
                          all(sp.sympify(v) == 0 for b in res[n]['beta'] for v in b.values()))
        al = res['0']['alpha']
        logs = {kk: sp.sympify(v) for kk, v in al.items() if kk != '1'}
        eqs.append(sp.sympify(al.get('1', 0)) - c1 * T1 - c2 * T2)
        for mu in range(4):
            bm = sp.sympify(res['0']['beta'][mu].get('1', 0))
            eqs.append(bm - (d1 * p1[mu] + d2 * k[mu]) * T1 - (d3 * p1[mu] + d4 * k[mu]) * T2)
        checks.append(('log constants', {kk: str(v) for kk, v in logs.items()}))
    sol = sp.solve(eqs, [c1, c2, d1, d2, d3, d4], dict=True)
    return sol, checks, len(eqs)

if __name__ == '__main__':
    files = sys.argv[1:]
    sol, checks, ne = fit(load(files))
    print(sol, checks, ne)
