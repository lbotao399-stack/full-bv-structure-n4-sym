"""classical part: {S_M, O} for O = int omega (iota_w1 A, iota_w2 A), canonical form (leg 2 eliminated, omega momentum k).
Expected 2(E + w_O) O = -2(2 + k.d/dk) O_k:  alpha = -4 v, beta^mu = -2 k^mu v  (v = vertex value of O)."""
import json, random, sympy as sp
import delta_local as DL
import fa
D = DL.D; V = DL.V; I = sp.I
P = list(sp.symbols('P0:4')); KAP = list(sp.symbols('K0:4'))
rnd = random.Random(7)
def rp(c):
    u = [sp.Integer(0)] * 10
    for i in c: u[i] = sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3]))
    return u
out = []
for s in range(2):
    w1 = rp(range(4, 10)); w2 = rp(range(4, 10)); u1 = rp(range(4, 10)); u2 = rp(range(4, 10))
    kv = [sp.Integer(rnd.choice([-2, -1, 1, 2])) for _ in range(4)]
    D.set_O(w1, w2)
    O = fa.local([(1, 'O')], 'O', 0)
    SM = DL.S_M_functional()
    br = fa.bracket(SM, O)
    leg1 = dict(sp='A', u=u1, p=list(P)); om = dict(sp='om', u=[sp.Integer(1)], p=list(kv))
    leg2 = dict(sp='A', u=u2, p=[-P[m] - kv[m] for m in range(4)])
    alpha = sp.expand(br([leg1, leg2, om]).get('k', 0)); beta = []
    for mu in range(4):
        l2 = dict(sp='A', u=u2, p=[-P[m] - kv[m] - KAP[m] for m in range(4)])
        lk = dict(sp=f'k{mu}', u=[sp.Integer(1)], p=list(KAP))
        v = sp.expand(br([leg1, l2, om, lk]).get('k', 0))
        alpha += -I * sp.diff(v, KAP[mu]).xreplace({x: 0 for x in KAP})
        beta.append(sp.expand(-I * v.xreplace({x: 0 for x in KAP})))
    vO = O([leg1, leg2, om]).get('k', 0)
    out.append(dict(alpha_over_v=str(sp.simplify(alpha / vO)), beta_over_v=[str(sp.simplify(b / vO)) for b in beta], k=[str(x) for x in kv]))
print(json.dumps(out, indent=1))
json.dump(dict(statement='{S_M, O} in canonical form (x on leg 2): expected alpha = -4 v, beta^mu = -2 k^mu v (= 2(E+w_O)O, w_O = 2)', samples=out),
          open('delta_classical.json', 'w'), indent=1)
