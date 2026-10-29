"""T2 cross-checks (engine A, round 2):
(1) Q on the relevant local functionals in the two-mode model (Q = {S_free,-}, [C07]/[ESW] bracket, classical_spec S3):
    Q int(A|_V,A|_V), Q int(d^*A|_V,d^*A|_V), Q int(dA|_V,dA|_V), Q of a (c,A^*) functional proportional to p,
    Q int(dc,A^*), Q of a general A|_V-A|_V kernel; read off the (c at +p, A|_V at -p) coefficient.
(2) Route (b): exact eps^{-1}, log eps and finite parts of every A|_V-A|_V self-energy class (cross coefficient
    W_1(a1 e^{ipx} + a2 e^{-ipx})), J = e^{-1} x/eps + y log(eps p^2 e^{gamma_E}) + z + o(1); longitudinal part
    gives Q W_1 in the (c, A|_V) sector."""
import json
import numpy as np
import sympy as sp
import free as Fm
from grass import GA
import r2_core as C
from r2_T2 import selfenergy_finite
from fractions import Fraction

p = Fm.p; I = sp.I; y = Fm.y
Sf = Fm.S_free()
Q = lambda X: Fm.bracket(Sf, X, 0).map(sp.expand)

def coeff_cA(Fq):
    """coefficient vector v_N with (QF)(r_c gamma e^{ipx} + a e^{-ipx}) = r_c gamma sum_N v_N a_N"""
    v = []
    for N in range(10):
        key = tuple(sorted((Fm.VAR[('c', 0, 1)], Fm.VAR[('A', N, -1)])))
        c_ = Fq.t.get(key, 0)
        # monomial stored with sorted keys; c is odd, A even: order irrelevant
        v.append(sp.expand(c_))
    return v

out = {}
F_AA = sum((y('A', M, 1) * y('A', M, -1) * 2 for M in range(4)), GA())
F_dd = GA()
for M in range(4):
    for N in range(4):
        F_dd = F_dd + y('A', M, 1) * y('A', N, -1) * (2 * p[M] * p[N])
F_FF = GA()
for M in range(4):
    for N in range(M + 1, 4):
        def Fmn(sg):
            k = [sg * x for x in p]
            return y('A', N, sg) * (I * k[M]) - y('A', M, sg) * (I * k[N])
        F_FF = F_FF + Fmn(1) * Fmn(-1) * 2
F_gh = GA()
for mu in range(4):
    F_gh = F_gh + y('c', 0, 1) * y('As', mu, -1) * p[mu] + y('c', 0, -1) * y('As', mu, 1) * (-p[mu])
F_dcAs = GA()
for (s1, s2) in ((1, -1), (-1, 1)):
    k1 = [s1 * x for x in p]
    for mu in range(4):
        F_dcAs = F_dcAs + y('c', 0, s1) * y('As', mu, s2) * (I * k1[mu])
x0, x1 = sp.symbols('x0 x1')
F_Pi = GA()
for M in range(4):
    for N in range(4):
        F_Pi = F_Pi + y('A', M, 1) * y('A', N, -1) * (x0 * (p[0]**2 + p[1]**2 + p[2]**2 + p[3]**2) * (1 if M == N else 0) + x1 * p[M] * p[N])
for nm, Fx in (('int(A|_V,A|_V)', F_AA), ('int(d^*A|_V,d^*A|_V)', F_dd), ('int(dA|_V,dA|_V)', F_FF),
               ('(c,A^*) functional propto p', F_gh), ('int(dc,A^*)', F_dcAs), ('A|_V kernel x0 p^2 delta + x1 p p', F_Pi)):
    v = coeff_cA(Q(Fx))
    out[nm] = [str(sp.factor(c)) for c in v[:4]]
    print('Q', nm, ' (c,+p ; A_N,-p) coefficients N=0..3:', out[nm], flush=True)

# ---- route (b): self-energy finite parts (exact)
rng = np.random.default_rng(77)
rows = []; data = []
for it in range(4):
    pp = rng.integers(-1, 2, size=4) * 2
    if not pp.any(): pp[0] = 2
    a1 = np.zeros((10, 3), complex); a1[:4] = rng.integers(-2, 3, size=(4, 3))
    a2 = np.zeros((10, 3), complex); a2[:4] = rng.integers(-2, 3, size=(4, 3))
    S0 = int(pp @ pp) * int(round(np.einsum('ma,ma->', a1[:4].real, a2[:4].real)))
    S1 = int(round(np.einsum('m,ma,n,na->', pp, a1[:4].real, pp, a2[:4].real)))
    S3 = int(round(np.einsum('ma,ma->', a1[:4].real, a2[:4].real)))
    se = selfenergy_finite(a1, a2, pp)
    rows.append((S0, S1, S3)); data.append(se)
    print('sample', it, flush=True)
keys = sorted(set().union(*[d.keys() for d in data]))
res = {}
tot = {'epsinv': [0, 0, 0], 'log': [0, 0, 0], 'fin': [0, 0, 0]}
M = sp.Matrix([[r[0], r[1], r[2]] for r in rows])
for k in keys:
    ent = {}
    for part in ('epsinv', 'log', 'fin'):
        b = sp.Matrix([d.get(k, {}).get(part, 0) for d in data])
        sol, params = M.gauss_jordan_solve(b)
        sol = sol.subs({q: 0 for q in params})
        assert (M * sol - b).applyfunc(sp.simplify) == sp.zeros(len(rows), 1)
        sol = [sp.nsimplify(sp.simplify(x / C.C_ADJ), [sp.log(2)]) for x in sol]
        ent[part] = {'p^2(a1.a2)': str(sol[0]), '(p.a1)(p.a2)': str(sol[1]), '(a1.a2)': str(sol[2])}
        tot[part] = [sp.simplify(t + x) for t, x in zip(tot[part], sol)]
    res[k] = ent
    print(k, ent, flush=True)
totals = {part: {'p^2(a1.a2)': str(v[0]), '(p.a1)(p.a2)': str(v[1]), '(a1.a2)': str(v[2]),
                 'longitudinal x0+x1': str(sp.simplify(v[0] + v[1]))} for part, v in tot.items()}
print('TOTAL (units u):', totals, flush=True)
json.dump({'Q_two_mode': out, 'selfenergy_classes': res, 'selfenergy_totals': totals},
          open('logs/r2_T2b.json', 'w'), indent=1)
