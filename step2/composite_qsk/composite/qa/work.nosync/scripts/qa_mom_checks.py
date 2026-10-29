"""QA third method, step 1: fix F0 on the k=1 chain (A_V, A_perp) -> A_perp at one engine point, then predict every
other number: the second k=1 point, the k=1 C^0 values, the bilinear k=2 chain (GC3) and the AW|AW off-shell table."""
import json, os, sys
import qa_mom as M, qa_mom_chains as CH
from qa_mom import C, I, cc
from fractions import Fraction as Fr
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
TAB = json.load(open(os.path.join(HERE, '..', '..', '..', 'engine_a', 'kernels_C2.json')))
def cnum(s):
    z = complex(s.replace('(', '').replace(')', ''))
    return C(Fr(z.real).limit_denominator(1000), Fr(z.imag).limit_denominator(1000))
def tosym(v):
    tot = 0
    for k, x in v.items():
        f = 1 if k == '1' else sp.log(int(k[3:]))
        tot += (sp.Rational(x.r.numerator, x.r.denominator) + sp.I * sp.Rational(x.i.numerator, x.i.denominator)) * f
    return sp.nsimplify(sp.expand(tot))
def asym(s):
    return sp.nsimplify(sp.expand(sp.sympify(s.replace('i', 'I').replace('+-', '-').replace('log2', 'log(2)').replace('log3', 'log(3)'))))
def dot(a, b):
    t = C()
    for x, y in zip(a, b):
        t = t + cc(x) * cc(y)
    return t
beta = [1, 2, -1, 3]
F0 = M.F0
out = {'F0': str(F0)}
# k=1
rows = TAB['AV|AW']['leg_sets']['AW']['rows']
res = []
for r in rows:
    za = [cnum(x) for x in r['zeta_a']]; zb = [cnum(x) for x in r['zeta_b']]; e = [cnum(x) for x in r['legs'][0]['pol']]
    q = r['legs'][0]['mom']
    fib = dot(e, zb)
    vP = CH.az(0, 0, q, ('P', beta), za); vD = CH.az(0, 0, q, 'D0', za)
    mineP = tosym(M.vscale(vP, F0 * fib)); mineD = tosym(M.vscale(vD, F0 * fib))
    res.append({'A_Pbeta': r['C2_Pbeta']['letter'], 'QA_Pbeta': str(mineP), 'agree_P': sp.simplify(mineP - asym(r['C2_Pbeta']['letter'])) == 0,
                'A_D0': r['C2_0_Delta']['letter'], 'QA_D0': str(mineD), 'agree_D0': sp.simplify(mineD - asym(r['C2_0_Delta']['letter'])) == 0})
out['k1_AV_AW_to_AW'] = res
print(json.dumps(res, indent=1))
# bilinear k=2 (GC3): per word (1,0) with generic momenta: P_b coefficient of i b.(p0+p1); D0 constant
q0 = [1, 2, -1, 1]; q1 = [2, -1, 0, 3]
def word10(kind, a=0, b=0, qa=q0, qb=q1):
    v = M.vadd(CH.zz(a, b, qa, qb, kind), CH.zz(b, a, qb, qa, kind))
    return M.vscale(v, F0)
vP = tosym(word10(('P', beta))); vD = tosym(word10('D0'))
ib = sp.I * sum(bb * (x + y) for bb, x, y in zip(beta, q0, q1))
out['bilinear_word10'] = {'C_Pb / (i b.(p0+p1))': str(sp.nsimplify(vP / ib)), 'C0': str(vD)}
print(out['bilinear_word10'])
# AW|AW table (word ['AW','AW'] = B word (1,0)): fibre (zeta_a.Y0)(zeta_b.Y1)
res = []
for r in TAB['AW|AW']['leg_sets']['AW AW']['rows']:
    za = [cnum(x) for x in r['zeta_a']]; zb = [cnum(x) for x in r['zeta_b']]
    Y0 = [cnum(x) for x in r['legs'][0]['pol']]; Y1 = [cnum(x) for x in r['legs'][1]['pol']]
    p0 = r['legs'][0]['mom']; p1 = r['legs'][1]['mom']
    fib = dot(za, Y0) * dot(zb, Y1)
    vP = tosym(M.vscale(word10(('P', beta), qa=p0, qb=p1), fib)); vD = tosym(M.vscale(word10('D0', qa=p0, qb=p1), fib))
    res.append({'A_Pbeta': r['C2_Pbeta']['letter'], 'QA_Pbeta': str(vP), 'agree_P': sp.simplify(vP - asym(r['C2_Pbeta']['letter'])) == 0,
                'A_D0': r['C2_0_Delta']['letter'], 'QA_D0': str(vD), 'agree_D0': sp.simplify(vD - asym(r['C2_0_Delta']['letter'])) == 0})
out['k2_AW_AW_table'] = res
print(json.dumps(res, indent=1))
json.dump(out, open(os.path.join(HERE, '..', 'qa_mom_checks.json'), 'w'), indent=1)
