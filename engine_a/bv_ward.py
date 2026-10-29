"""Two-mode check of the tree-level identity used in route (i):
   cross_P(Q Psi, O) = Q(cross_P(Psi, O)) - {O, Psi}   (P = P(0,infty), QO = 0 on the tested O)
and the explicit leg-correction cross_P(L_2, O) for single letters and the scalar bilinear."""
import json
import sympy as sp
from grass import GA, PAR, dP
import free as F
from bv_exact import Q, Psi, L2, y
p = F.p; I = sp.I; p2 = F.p2
inv = {v: k for k, v in F.VAR.items()}
def Phat(a, b, sg):
    (n1, i1), (n2, i2) = a, b
    q = [sg * x for x in p]
    if n1 == 'A' and n2 == 'A': return 2 if i1 == i2 else 0
    if n1 == 'c' and n2 == 'As': return -I * q[i2] if i2 < 4 else 0
    if n1 == 'As' and n2 == 'c': return -I * q[i1] if i1 < 4 else 0
    if n1 == 'l' and n2 == 'l': return I * sum(F.S_[mu][i1, i2] * q[mu] for mu in range(4))
    return 0
Pm = {}
for va, ka in inv.items():
    for vb, kb in inv.items():
        if ka[2] == -kb[2]:
            v = Phat((ka[0], ka[1]), (kb[0], kb[1]), ka[2])
            if v != 0: Pm[(va, vb)] = v / p2
VARS = list(inv.keys())
def cross(Fa, Fb):
    return (dP(Fa * Fb, Pm, VARS) - dP(Fa, Pm, VARS) * Fb - Fa * dP(Fb, Pm, VARS)).map(lambda c: sp.simplify(c))
tests = {'letter_Aperp': y('A', 4, 1), 'letter_lambda': y('l', 3, 1), 'letter_c': y('c', 0, 1),
         'bilinear_Aperp': y('A', 4, 1) * y('A', 5, -1), 'letter_AV': y('A', 1, 1)}
out = {}
for nm, O in tests.items():
    par = sum(PAR[v] for v in list(O.t.keys())[0]) % 2
    lhs = cross(Q(Psi), O)
    cP = cross(Psi, O)
    QcP = Q(cP).map(sp.simplify)
    brk = F.bracket(O, Psi, par).map(sp.simplify)
    chk = (lhs - QcP + brk).map(sp.simplify)
    leg = cross(L2, O)
    def ratio(X):
        r = {}
        for k, c in X.t.items():
            r[str([inv[v] for v in k])] = str(sp.simplify(c / O.t[k])) if k in O.t else str(sp.factor(c))
        return r
    out[nm] = {'identity_residual_terms': len(chk.t), 'cross(L2,O)/O': ratio(leg), 'Q(cross(Psi,O))/O': ratio(QcP), '{O,Psi}/O': ratio(brk)}
    print(nm, json.dumps(out[nm]))
json.dump(out, open('logs/bv_ward.json', 'w'), indent=1)
