"""Q of the candidate quadratic arity-one finite counterterms (canonical form, x on the antifield leg):
   int(A,A*) on each sector, int(E A, A*), and the field-counting functional pieces."""
import json, random, sympy as sp
import delta_local as DL
import fa
V = DL.V
P = list(sp.symbols('P0:4')); KAP = list(sp.symbols('K0:4'))
rnd = random.Random(9)
def rp(sp_, comps):
    u = [sp.Integer(0)] * V.DIM[sp_]
    for c in comps: u[c] = sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3]))
    return u
cands = {
 'int(A_perp,A*_perp)': DL.canon('A', 'As', lambda u, w, p: DL.dot(u, w, DL.APc), None, 'F_Ap'),
 'int(E A_perp,A*_perp)': DL.canon('A', 'As', None, lambda u, w, p, m: p[m] * DL.dot(u, w, DL.APc), 'F_EAp'),
 'int(A_V,A*_V)': DL.canon('A', 'As', lambda u, w, p: DL.dot(u, w, DL.AVc), None, 'F_AV'),
 'int(E A_V,A*_V)': DL.canon('A', 'As', None, lambda u, w, p, m: p[m] * DL.dot(u, w, DL.AVc), 'F_EAV'),
 'int(lambda*,lambda)': DL.canon('lam', 'lams', lambda u, w, p: DL.lam_pair(u, w), None, 'F_l'),
 'int(c*,c)': DL.canon('c', 'cs', lambda u, w, p: u[0] * w[0], None, 'F_c'),
}
out = {}
for nm, F in cands.items():
    QF = fa.Q(F); res = {}
    for (s1, s2, c1, c2) in [('A', 'A', DL.APc, DL.APc), ('A', 'A', DL.AVc, DL.AVc), ('c', 'As', [0], DL.AVc), ('lam', 'lam', range(16), range(16))]:
        u1 = rp(s1, c1); u2 = rp(s2, c2)
        a, b = DL.to_canon(QF, s1, s2, u1, u2, P, KAP)
        if a == 0 and all(x == 0 for x in b):
            continue
        # express in terms of plane-wave values of standard functionals
        PP = sum(x * x for x in P)
        if s1 == 'A' and list(c1) == DL.APc:
            ref = 2 * PP * DL.dot(u1, u2, DL.APc)        # int(dA_perp, dA_perp)
            res['(A_perp,A_perp) in units of int(dA_perp,dA_perp)'] = [str(sp.simplify(a / ref)), [str(x) for x in b]]
        elif s1 == 'A':
            pu1 = DL.dot(P, u1, DL.AVc); pu2 = DL.dot(P, u2, DL.AVc)
            ref1 = 2 * (PP * DL.dot(u1, u2, DL.AVc) - pu1 * pu2)    # int(dA_V,dA_V)
            ref2 = 2 * pu1 * pu2                                     # int(d*A_V, d*A_V)
            x1, x2 = sp.symbols('x1 x2')
            sol = sp.solve(sp.Poly(sp.expand(a - x1 * ref1 - x2 * ref2), *P).coeffs(), [x1, x2], dict=True)
            res['(A_V,A_V) in units of [int(dA_V,dA_V), int(d*A_V,d*A_V)]'] = [str(sol), [str(x) for x in b]]
        elif s1 == 'c':
            ref = sp.I * DL.dot(P, u2, DL.AVc) * u1[0]           # plane-wave value of int(dc, A*) on (c at P, A* at -P)
            res['(c,A*_V) in units of int(dc,A*)'] = [str(sp.simplify(a / ref)), [str(x) for x in b]]
        else:
            ref = 2 * V.vertex('kinlam', [dict(sp='lam', u=u1, p=list(P)), dict(sp='lam', u=u2, p=[-x for x in P])])['k']
            res['(lambda,lambda) in units of int(lambda,dslash lambda)'] = [str(sp.simplify(a / ref)), [str(x) for x in b]]
    out[nm] = res
json.dump(out, open('delta_Qbasis.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
