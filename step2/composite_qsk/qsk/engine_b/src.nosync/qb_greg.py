"""Gate G-reg along P_b (route B): the defect cycles with S^(1)(P_b) anywhere on the cycle and one special edge, for the
scalar pair (A_perp, A_perp) -> A_perp A_perp at the pair weight, symbolic letters, outputs and b; decomposition on
the SO(6) structures and comparison with the composite kernel C^(2)_{P_b} of composite engine B (out.nosync/so6_AWAW.json:
word (0,1): -1/3 i(b.p0)(a.Y1)(b.Y0) - 1/3 i(b.p1)(a.Y1)(b.Y0); word (1,0): the same with Y0 <-> Y1), Cor qs-DP.
Also the eps^0 part of the pair (A_V, A_V) is reported (gauge-variant letters, no comparison claimed)."""
import os, sys, json, time
import sympy as sp
from sympy.polys.domains import QQ_I
import qb_core as B
import qb_job as QJ

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')


def decompose(job, P, basis):
    names = list(basis)
    syms = sp.symbols(f'c0:{len(names)}')
    expr = P.as_expr() - sum(s * basis[nm].as_expr() for s, nm in zip(syms, names))
    gens = [sp.Symbol(x) for x in job.names if x.startswith(('za', 'zb', 'o0', 'o1', 'r0', 'r1', 'b'))]
    eqs = sp.Poly(sp.expand(expr), *gens).coeffs()
    sol = sp.solve(eqs, syms, dict=True)
    if not sol:
        return 'NOT IN SPAN'
    return {nm: str(sol[0].get(s, 'free')) for s, nm in zip(syms, names)}


def run(slot='AW'):
    t0 = time.time()
    bsym = list(sp.symbols('b0:4'))
    job = QJ.DJob(slot, slot, [slot, slot], 'Pb', bsym, tag='greg')
    res = job.kernel(mlist=(0, 1, 2), verbose=True)
    R = job.R; g = job.g
    comps = range(4, 10) if slot == 'AW' else range(4)
    za = [g[f'za_{c}'] for c in comps]; zb = [g[f'zb_{c}'] for c in comps]
    Y0 = [g[f'o0_{c}'] for c in comps]; Y1 = [g[f'o1_{c}'] for c in comps]
    dot = lambda x, y: sum((a * b for a, b in zip(x, y)), R.zero)
    p0 = [g[f'r0_{m}'] for m in range(4)]; p1 = [g[f'r1_{m}'] for m in range(4)]
    bb = [g[f'b{m}'] for m in range(4)]
    S = {'(a.Y0)(b.Y1)': dot(za, Y0) * dot(zb, Y1), '(a.Y1)(b.Y0)': dot(za, Y1) * dot(zb, Y0), '(a.b)(Y0.Y1)': dot(za, zb) * dot(Y0, Y1)}
    basis = {}
    for kname, sv in S.items():
        for j, pj in ((0, p0), (1, p1)):
            basis[f'i(b.p{j}) {kname}'] = sv * dot(bb, pj) * R(QQ_I(0, 1))
    out = {'raw': QJ.table_str(job, res['total']), 'decomposition': {}}
    for w, dm in res['total'].items():
        for m, d in dm.items():
            for bk, P in d.items():
                out['decomposition'][f'word {list(w)} eps^-{m} {bk}'] = decompose(job, P, basis) if slot == 'AW' else 'not decomposed'
    out['seconds'] = round(time.time() - t0, 1)
    out['configs'] = res['configs']
    out['peak_MB'] = round(QJ.peak_mb(), 1)
    return out


if __name__ == '__main__':
    slot = sys.argv[1] if len(sys.argv) > 1 else 'AW'
    out = run(slot)
    json.dump(out, open(os.path.join(OUT, f'greg_Pb_{slot}.json'), 'w'), indent=1)
    for kk, v in out['decomposition'].items():
        print(kk, v)
    print('seconds', out['seconds'], 'peak MB', out['peak_MB'])
