"""F0: general SO(6) scalar pair (A_perp, A_perp) -> A_perp A_perp at W0 with symbolic letters and outputs;
decomposition of the planar density into identity, exchange and trace (Prop prop:cc-scalar)."""
import os, sys, time, json, pickle
import sympy as sp
from sympy.polys.domains import QQ_I
import ce_job as CJ
HERE = os.path.dirname(os.path.abspath(__file__))
I = sp.I

def run():
    t0 = time.time()
    job = CJ.CJob('AW', 'AW', ['AW', 'AW'], pol='sym', za=None, zb=None)
    res = job.kernels(mlist=(0, 1), verbose=True)
    R = job.R; g = job.g
    za = [g[f'za_{c}'] for c in range(4, 10)]; zb = [g[f'zb_{c}'] for c in range(4, 10)]
    Y0 = [g[f'o0_{c}'] for c in range(4, 10)]; Y1 = [g[f'o1_{c}'] for c in range(4, 10)]
    dot = lambda x, y: sum((a * b for a, b in zip(x, y)), R.zero)
    p0 = [g[f'r0_{m}'] for m in range(4)]; p1 = [g[f'r1_{m}'] for m in range(4)]
    a = [g[f'a{m}'] for m in range(4)]
    S = {'(a.Y0)(b.Y1)': dot(za, Y0) * dot(zb, Y1), '(a.Y1)(b.Y0)': dot(za, Y1) * dot(zb, Y0), '(a.b)(Y0.Y1)': dot(za, zb) * dot(Y0, Y1)}
    out = {}
    for tab in ('Pb', 'D0'):
        for w, dm in res[tab].items():
            for m, d in dm.items():
                for b, P in d.items():
                    if tab == 'D0':
                        basis = dict(S)
                    else:
                        basis = {}
                        for kname, sv in S.items():
                            for j, pj in ((0, p0), (1, p1)):
                                basis[f'i(b.p{j}) {kname}'] = sv * dot(a, pj) * R(QQ_I(0, 1))
                    # least squares over exact coefficients: solve by matching monomials
                    names = list(basis)
                    syms = sp.symbols(f'c0:{len(names)}')
                    expr = (P - sum((R(0) for _ in [0]), R.zero)).as_expr() - sum(s * basis[nm].as_expr() for s, nm in zip(syms, names))
                    eqs = sp.Poly(sp.expand(expr), *[sp.Symbol(x) for x in job.names if x.startswith(('za', 'zb', 'o0', 'o1', 'r0', 'r1', 'a'))]).coeffs()
                    sol = sp.solve(eqs, syms, dict=True)
                    out[f'{tab} word {w} eps^-{m} {b}'] = {nm: str(sol[0].get(s, 'free')) for s, nm in zip(syms, names)} if sol else 'NOT IN SPAN'
    out['seconds'] = time.time() - t0
    return out

if __name__ == '__main__':
    out = run()
    for k, v in out.items():
        print(k, v)
    json.dump(out, open(os.path.join(HERE, '..', 'out.nosync', 'so6_AWAW.json'), 'w'), indent=1)
