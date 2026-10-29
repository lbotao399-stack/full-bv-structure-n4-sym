"""Quadratic-level check that the log-eps counterterm L_2 = -(B1+B4) + 4 B6 - B7 (units hbar c_adj/16pi^2) is BV-exact:
L_2 = Q Psi_2,  Psi_2 = zA int(A,A*) + zl int(lambda*,lambda) + zc int(c*,c)  ([C07] pairing order of the spec).
Also the brackets {O, Psi} for single letters and the scalar bilinear (two-mode model, Vol stripped)."""
import json
import sympy as sp
from grass import GA, PAR
import free as F
p, I = F.p, sp.I
y = F.y
Sf = F.S_free()
def mode_pair(n1, i1, n2, i2):
    return y(n1, i1, 1) * y(n2, i2, -1) + y(n1, i1, -1) * y(n2, i2, 1)
PsiA = sum((mode_pair('A', M, 'As', M) for M in range(10)), GA())
Psil = sum((mode_pair('ls', a, 'l', a) for a in range(16)), GA())
Psic = mode_pair('cs', 0, 'c', 0)
Q = lambda X: F.bracket(Sf, X, 0).map(sp.expand)
def mom(sg): return [sg * x for x in p]
def B(name):
    t = GA()
    if name in ('B1', 'B4'):
        for M in range(10):
            for N in range(M + 1, 10):
                if name == 'B1' and N >= 4: continue
                if name == 'B4' and not (M < 4 <= N): continue
                def Fm(sg):
                    k = mom(sg); r = GA()
                    if M < 4: r = r + y('A', N, sg) * (I * k[M])
                    if N < 4: r = r - y('A', M, sg) * (I * k[N])
                    return r
                t = t + Fm(1) * Fm(-1) * 2
    if name == 'B6':
        for (s1, s2) in ((1, -1), (-1, 1)):
            k2 = mom(s2)
            for mu in range(4):
                for a in range(16):
                    for b in range(16):
                        c = F.SB[mu][a, b]
                        if c != 0:
                            t = t + y('l', a, s1) * y('l', b, s2) * (c * I * k2[mu])
    if name == 'B7':
        for (s1, s2) in ((1, -1), (-1, 1)):
            k1 = mom(s1)
            for mu in range(4):
                t = t + y('c', 0, s1) * y('As', mu, s2) * (I * k1[mu])
    return t.map(sp.expand)
zA, zl, zc = sp.symbols('zA zl zc')
L2 = (B('B1') + B('B4')) * (-1) + B('B6') * 4 - B('B7')
QPsi = Q(PsiA) * zA + Q(Psil) * zl + Q(Psic) * zc
diff = (QPsi - L2).map(sp.expand)
eqs = []
for c in diff.t.values():
    eqs += sp.Poly(sp.expand(c), *p).coeffs()
sol = sp.solve(eqs, [zA, zl, zc], dict=True)
print('Psi coefficients solving Q Psi_2 = L_2:', sol)
s0 = sol[0]
Psi = PsiA * s0[zA] + Psil * s0[zl] + Psic * s0[zc]
out = {'Psi': {str(k): str(v) for k, v in s0.items()}}
# brackets {O, Psi}
tests = {'letter_Aperp': (y('A', 4, 1), 0), 'letter_AV': (y('A', 1, 1), 0), 'letter_lambda': (y('l', 3, 1), 1),
         'letter_c': (y('c', 0, 1), 1), 'bilinear_Aperp': (y('A', 4, 1) * y('A', 5, -1), 0)}
for nm, (O, par) in tests.items():
    br = F.bracket(O, Psi, par).map(sp.expand)
    ratio = None
    for k, c in br.t.items():
        ratio = sp.simplify(c / O.t[k]) if k in O.t else ('other', k, c)
    print(nm, '{O,Psi} = ', ratio, '* O')
    out[nm] = str(ratio)
json.dump(out, open('logs/bv_exact.json', 'w'), indent=1)
