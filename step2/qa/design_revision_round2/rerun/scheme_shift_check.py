# Audit (odd_rev_audit_graphs), script 3: what a nonzero P_1 forces on K_1, Q_1, S_1 (engine finding of Part D4,
# step2/engine_b/delta_results.json part2.K5 'beta_part' and step2/qa/compare_delta.json: translation component of the
# scalar bilinears = (1 - (1/3 + f_P) unit) D_P O, i.e. P_1 = c D_P with c = -(1/3 + f_P) unit, letters c = -f_P).
#
# (A) Even part, explicit.  The conformal Killing fields of R^4 act on scalar densities of weight w by
#     rho0(v) f = v.grad f + (w/4) div(v) f, [v,u] = v(u) - u(v) (ch02 def:c2-even).  One-loop ansatz on a family of
#     scalar operators: rho1(P_b) = cP rho0(P_b), rho1(M) = cM rho0(M), rho1(Delta) = cD rho0(Delta) + 2 gam,
#     rho1(K_b) = cK rho0(K_b) + kap (1/4) div(v_K) (the K_1 structures of ch07 cor:c7-K for scalar operators).
#     Impose the order-hbar homomorphism  [rho0(x),rho1(y)] + [rho1(x),rho0(y)] = rho1([x,y])  on all basis pairs,
#     applied to a generic cubic polynomial f, and solve for (cM, cD, cK, kap) given (cP, gam).
# (B) Odd part, block level.  For multiplier ansaetze rho1(X) = c_X rho0(X) the condition is c_X + c_Y = c_Z for every
#     block bracket [X,Y] with a nonzero component along Z (ch02 prop:c2-special, def:c2-bracket).  Solve the
#     homogeneous system; compare with the ad_Delta grading (SW (17)) and the outer derivation B (ch08 def:c8-B).
import json, os, itertools
import sympy as sp

x = sp.symbols('x0:4'); w = sp.Symbol('w')
def vf_P(b): return [sp.Integer(bi) for bi in b]
def vf_M(i, j):
    v = [0] * 4; v[i] = x[j]; v[j] = -x[i]; return [sp.sympify(t) for t in v]
def vf_D(): return [2 * xi for xi in x]
def vf_K(b):
    xx = sum(xi * xi for xi in x); bx = sum(bi * xi for bi, xi in zip(b, x))
    return [sp.expand(xx * bi - 2 * bx * xi) for bi, xi in zip(b, x)]
def div(v): return sp.expand(sum(sp.diff(v[i], x[i]) for i in range(4)))
def bracket(v, u):
    return [sp.expand(sum(v[j] * sp.diff(u[i], x[j]) - u[j] * sp.diff(v[i], x[j]) for j in range(4))) for i in range(4)]
def rho0(v):
    return lambda f: sp.expand(sum(v[i] * sp.diff(f, x[i]) for i in range(4)) + w / 4 * div(v) * f)

# generic cubic test function
mons = [m for d in range(4) for m in itertools.combinations_with_replacement(x, d)]
coef = sp.symbols('f0:%d' % len(mons))
f = sum(cf * sp.Mul(*m) for cf, m in zip(coef, mons))

cP, cM, cD, cK, kap, gam = sp.symbols('cP cM cD cK kap gam')
basis = {}
for k in range(4):
    e = [0] * 4; e[k] = 1
    basis['P%d' % k] = ('P', vf_P(e)); basis['K%d' % k] = ('K', vf_K(e))
for i, j in itertools.combinations(range(4), 2):
    basis['M%d%d' % (i, j)] = ('M', vf_M(i, j))
basis['D'] = ('D', vf_D())

def decompose(v):
    """write a conformal Killing field as a combination of the basis fields (exact linear solve)."""
    names = list(basis); cs = sp.symbols('a0:%d' % len(names))
    expr = [sp.expand(v[i] - sum(c * basis[n][1][i] for c, n in zip(cs, names))) for i in range(4)]
    eqs = []
    for e in expr:
        eqs += sp.Poly(e, *x).coeffs() if e != 0 else []
    sol = sp.solve(eqs, cs, dict=True)
    assert sol, 'not in the span'
    return {n: sol[0].get(c, 0) for c, n in zip(cs, names) if sol[0].get(c, 0) != 0}

def rho1(name):
    kind, v = basis[name]
    r0 = rho0(v)
    if kind == 'P': return lambda g: cP * r0(g)
    if kind == 'M': return lambda g: cM * r0(g)
    if kind == 'D': return lambda g: cD * r0(g) + 2 * gam * g
    if kind == 'K': return lambda g: cK * r0(g) + kap * sp.Rational(1, 4) * div(v) * g

# check rho0 is a homomorphism (sanity) and collect the order-hbar equations
eqs = []; hom_ok = True
names = list(basis)
for a, b in itertools.combinations(names, 2):
    va, vb = basis[a][1], basis[b][1]
    br = bracket(va, vb)
    lhs0 = sp.expand(rho0(va)(rho0(vb)(f)) - rho0(vb)(rho0(va)(f)))
    rhs0 = rho0(br)(f)
    if sp.expand(lhs0 - rhs0) != 0:
        hom_ok = False
    dec = decompose(br) if any(t != 0 for t in br) else {}
    lhs1 = sp.expand(rho0(va)(rho1(b)(f)) - rho1(b)(rho0(va)(f)) + rho1(a)(rho0(vb)(f)) - rho0(vb)(rho1(a)(f)))
    rhs1 = sp.expand(sum(cc * rho1(n)(f) for n, cc in dec.items()))
    diff = sp.expand(lhs1 - rhs1)
    if diff != 0:
        eqs += sp.Poly(diff, *(list(x) + list(coef))).coeffs()
eqs = list({sp.factor(e) for e in eqs if e != 0})
solA = sp.solve(eqs, [cM, cD, cK, kap], dict=True)

# (B) block level
blocks = ['P', 'K', 'D', 'M', 'R', 'Q', 'Qb', 'S', 'Sb']
lam = {'P': -2, 'K': 2, 'D': 0, 'M': 0, 'R': 0, 'Q': -1, 'Qb': -1, 'S': 1, 'Sb': 1}      # ad_Delta, SW (17), ch02
bB = {'P': 0, 'K': 0, 'D': 0, 'M': 0, 'R': 0, 'Q': 1, 'Qb': -1, 'S': -1, 'Sb': 1}       # ch08 def:c8-B
# nonzero block brackets [X,Y] -> components (ch02 prop:c2-special; Q: S+ (x) W, Qb: S- (x) W*, S: S+ (x) W*, Sb: S- (x) W)
BR = [('D', 'P', ['P']), ('D', 'K', ['K']), ('D', 'Q', ['Q']), ('D', 'Qb', ['Qb']), ('D', 'S', ['S']), ('D', 'Sb', ['Sb']),
      ('M', 'P', ['P']), ('M', 'K', ['K']), ('M', 'Q', ['Q']), ('M', 'Qb', ['Qb']), ('M', 'S', ['S']), ('M', 'Sb', ['Sb']),
      ('M', 'M', ['M']), ('R', 'R', ['R']), ('R', 'Q', ['Q']), ('R', 'Qb', ['Qb']), ('R', 'S', ['S']), ('R', 'Sb', ['Sb']),
      ('Q', 'Qb', ['P']), ('S', 'Sb', ['K']), ('Q', 'S', ['D', 'M', 'R']), ('Qb', 'Sb', ['D', 'M', 'R']),
      ('K', 'Q', ['Sb']), ('K', 'Qb', ['S']), ('P', 'S', ['Qb']), ('P', 'Sb', ['Q']), ('K', 'P', ['D', 'M'])]
c = {X: sp.Symbol('c_' + X) for X in blocks}
eqB = []
for X, Y, Zs in BR:
    for Z in Zs:
        eqB.append(c[X] + c[Y] - c[Z])
solB = sp.linsolve(eqB, [c[X] for X in blocks])
s, t = sp.symbols('s t')
gen = {X: s * lam[X] + t * bB[X] for X in blocks}
check_gen = all(sp.simplify(e.subs({c[X]: gen[X] for X in blocks})) == 0 for e in eqB)
sol_list = list(solB)[0]
free = sorted({str(sy) for e in sol_list for sy in sp.sympify(e).free_symbols})
cc = sp.Symbol('c')       # P_1 = c D_P
shift = {X: sp.simplify(gen[X].subs(s, -cc / 2)) for X in blocks}
out = {
    'A_rho0_homomorphism': hom_ok,
    'A_solution_given_cP_gam': [{str(k): str(v) for k, v in sol.items()} for sol in solA],
    'B_general_solution': {str(X): str(e) for X, e in zip(blocks, sol_list)},
    'B_free_parameters': free,
    'B_solution_space_is_span(ad_Delta grading, B)': bool(check_gen) and len(free) == 2,
    'B_multipliers_when_P1_equals_c_DP (t free)': {X: str(shift[X]) for X in blocks},
    'consequences': {
        'K_1 multiplier': str(shift['K']),
        'S_1(eta) - Q_1(rho(x)eta), eta in S+ (x) W*': str(sp.simplify(shift['S'] - shift['Qb'])),
        'S_1(etabar) - Q_1(rho(x)etabar), etabar in S- (x) W': str(sp.simplify(shift['Sb'] - shift['Q'])),
    },
    'statement': ('If P_1 = c {S^(1)(P_b),-} on a family (engine B: c = -(1/3 + f_P) unit on the scalar bilinears, -f_P on '
                  'letters), the order-hbar realisation forces K_1 = gamma (1/4) div v - c {S^(1)(K_b),-} (kap = gam, cK = -cP) '
                  'and S_1(eta) = Q_1(rho(x)eta) - c {S^(1)(eta),-} (+ U_eta), with M_1 and the Delta rescaling zero; '
                  'the outer direction t cancels in these differences.  Design thm:o-KD and eq:o-SQ hold only for c = 0.'),
}
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scheme_shift_check.json'), 'w'), indent=1)
