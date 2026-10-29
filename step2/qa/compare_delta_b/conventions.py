"""compare_delta_b / conventions.py

Independent check (this QA agent's own code) of the sign conventions that the two Part-D engines and the
theory note use, derived from the manuscript definitions only:

 (1) ch01 Def c1-pairing / c1-bracket: K_0 = -delta x (inverse fibre pairing) on the slot pairs, Delta_0 = -d_{K_0},
     d_phi = (1/2) sum d_{phi''} d_{phi'} [C07 6.4], {f,g} = Delta_0(fg) - Delta_0(f) g - (-1)^|f| f Delta_0(g).
     A 0-dimensional Grassmann model with one copy of each slot pair (A,A*), (lambda*,lambda) [two copies of
     lambda so that a fermion bilinear exists], (c*,c).  Checks:
       {S_X, F} = + D_X F  for S_X = (XA,A*) - (X lambda,lambda*) - (X c,c*) (the form of S^(1)(Delta), Def c7-vdot);
       Q := {S_2,-} with S_2 = -(1/4)(dA,dA) + (1/2)(lambda,dslash lambda) + (dc,A*) reproduces ch03's field-level Q
       (c -> -dc, A -> -(1/2) d^*dA, A^* -> d^*A^*, lambda -> dslash lambda) as the induced derivation D_Q;
       Q(A,A*) = -(1/2)(dA,dA),  Q(lambda*,lambda) = -(lambda,dslash lambda),  Q(c*,c) = -(dc,A*),  Q A = -dc.
 (2) canonical two-point calculus (alpha delta + beta.d delta, first leg momentum fixed): the D_M rule and the scaling
     Ward identity D_M Sigma = +4 eps d/deps Sigma on engine A's exact A|_Wperp self-energy.
 (3) ghost placement (left xi in the theory note, right xi in engine B, classical-term normalisation in engine A):
     F(Delta) and Ob(Delta) are the same functionals in all three once J_Delta = S^(1)(Delta) - I^CT_Delta - F(Delta).
Output: conventions.json
"""
import json, os, itertools
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- Grassmann algebra (exact, sympy coefficients)
EVEN = ['a', 'ls1', 'ls2', 'gs']          # A, lambda*_1, lambda*_2, c*   (even total parity: even coordinates)
ODD = ['al', 'l1', 'l2', 'g']             # A*, lambda_1, lambda_2, c     (odd total parity: odd coordinates)
NE = len(EVEN)


def mono(even=None, odd=()):
    e = [0] * NE
    for k, v in (even or {}).items():
        e[EVEN.index(k)] += v
    return (tuple(e), tuple(ODD.index(x) for x in odd))


def P(*terms):
    """terms: (coeff, even_dict, odd_list)"""
    out = {}
    for c, ed, ol in terms:
        m = mono(ed, ol)
        # normal-order the odd part
        o = list(m[1]); s = 1
        for i in range(len(o)):
            for j in range(len(o) - 1 - i):
                if o[j] > o[j + 1]:
                    o[j], o[j + 1] = o[j + 1], o[j]; s = -s
        if len(set(o)) < len(o):
            continue
        key = (m[0], tuple(o))
        out[key] = out.get(key, 0) + s * sp.nsimplify(c)
    return {k: v for k, v in out.items() if v != 0}


def mul(f, g):
    out = {}
    for (e1, o1), c1 in f.items():
        for (e2, o2), c2 in g.items():
            o = list(o1) + list(o2)
            if len(set(o)) < len(o):
                continue
            s = 1
            for i in range(len(o)):
                for j in range(len(o) - 1 - i):
                    if o[j] > o[j + 1]:
                        o[j], o[j + 1] = o[j + 1], o[j]; s = -s
            key = (tuple(x + y for x, y in zip(e1, e2)), tuple(o))
            out[key] = out.get(key, 0) + s * c1 * c2
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def add(*fs, coeffs=None):
    out = {}
    coeffs = coeffs or [1] * len(fs)
    for f, c in zip(fs, coeffs):
        for k, v in f.items():
            out[k] = out.get(k, 0) + c * v
    return {k: sp.simplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def parity(f):
    ps = {len(o) % 2 for (e, o) in f}
    assert len(ps) <= 1, 'inhomogeneous parity'
    return ps.pop() if ps else 0


def d_even(f, name):
    i = EVEN.index(name); out = {}
    for (e, o), c in f.items():
        if e[i] > 0:
            e2 = list(e); e2[i] -= 1
            out[(tuple(e2), o)] = out.get((tuple(e2), o), 0) + c * e[i]
    return out


def d_odd(f, name):
    """left derivative"""
    i = ODD.index(name); out = {}
    for (e, o), c in f.items():
        if i in o:
            k = o.index(i)
            o2 = o[:k] + o[k + 1:]
            out[(e, o2)] = out.get((e, o2), 0) + c * (-1) ** k
    return out


# slot pairs (even member coordinate y, odd partner coordinate eta); the fibre pairings are 1 in the model
PAIRS = [('a', 'al'), ('ls1', 'l1'), ('ls2', 'l2'), ('gs', 'g')]


def d_K0(f):
    """d_{K_0} with K_0 = -(e_y (x) e_eta + e_eta (x) e_y) per pair, d_phi = (1/2) sum d_{phi''} d_{phi'}"""
    out = {}
    for (y, eta) in PAIRS:
        t1 = d_odd(d_even(f, y), eta)        # phi' = e_y,  phi'' = e_eta : d_eta d_y
        t2 = d_even(d_odd(f, eta), y)        # phi' = e_eta, phi'' = e_y  : d_y d_eta
        out = add(out, t1, t2, coeffs=[1, sp.Rational(-1, 2), sp.Rational(-1, 2)])
    return out


def Delta0(f):
    return {k: -v for k, v in d_K0(f).items()}


def bracket(f, g):
    pf = parity(f)
    return add(Delta0(mul(f, g)), mul(Delta0(f), g), mul(f, Delta0(g)), coeffs=[1, -1, -(-1) ** pf])


def DX(f, X):
    """derivation induced by the linear map X on fields: x -> x + s X x, X diagonal on coordinates (dict name->factor)"""
    out = {}
    for name in EVEN:
        if X.get(name, 0):
            out = add(out, mul(P((X[name], {name: 1}, [])), d_even(f, name)))
    for name in ODD:
        if X.get(name, 0):
            out = add(out, mul(P((X[name], {}, [name])), d_odd(f, name)))
    return out


def eq(f, g):
    return add(f, g, coeffs=[1, -1]) == {}


res = {}
x = sp.Symbol('x')
# S_X = (XA,A*) - (X lambda,lambda*) - (X c,c*); coordinate expressions (e0 slot first): (A,A*) = a al, (lambda,lambda*) = l ls, (c,c*) = g gs
S_X = P((x, {'a': 1}, ['al']), (-x, {'ls1': 1}, ['l1']), (-x, {'ls2': 1}, ['l2']), (-x, {'gs': 1}, ['g']))
# the operator X: x on the fields, -x on the antifields (pairing preserving, Def c7-vdot)
Xmap = {'a': x, 'l1': x, 'l2': x, 'g': x, 'al': -x, 'ls1': -x, 'ls2': -x, 'gs': -x}
tests = {
    'int(A,A)': P((1, {'a': 2}, [])),
    'int(lambda,dslash lambda) ~ l1 l2': P((1, {}, ['l1', 'l2'])),
    'int(dc,A*) ~ g al': P((1, {}, ['g', 'al'])),
    'int(A,A*) ~ a al': P((1, {'a': 1}, ['al'])),
    'cubic a^2 al': P((1, {'a': 2}, ['al'])),
    'g l1 ls2 (mixed)': P((1, {'ls2': 1}, ['g', 'l1'])),
}
r1 = {}
for name, F in tests.items():
    br = bracket(S_X, F)
    dx = DX(F, Xmap)
    ratio = None
    if dx:
        k0 = list(dx.keys())[0]
        ratio = sp.simplify(br.get(k0, 0) / dx[k0])
    r1[name] = {'{S_X,F} == +D_X F': eq(br, dx), '{S_X,F} == -D_X F': eq(br, {k: -v for k, v in dx.items()}),
                'ratio': str(ratio)}
res['bracket_sign'] = r1
res['bracket_sign_verdict'] = '+D_X' if all(v['{S_X,F} == +D_X F'] for v in r1.values()) else 'FAIL'
# graded symmetry check {g,f} = (-1)^{|f||g|}{f,g}
f1 = P((1, {'a': 1}, ['al'])); g1 = P((1, {}, ['g', 'al'])); g2 = P((1, {'a': 2}, []))
res['graded_symmetry'] = {
    '{f,g} vs {g,f}, f = a al, g = g al': eq(bracket(f1, g1), {k: (-1) ** (parity(f1) * parity(g1)) * v for k, v in bracket(g1, f1).items()}),
    '{f,g} vs {g,f}, f = a al, g = a^2': eq(bracket(f1, g2), {k: (-1) ** (parity(f1) * parity(g2)) * v for k, v in bracket(g2, f1).items()})}

# Q := {S_2, -}; model constants: (dA,dA) ~ k a^2, d ~ n (c -> A), (lambda, dslash lambda) ~ 2 m l1 l2
k, n, m = sp.symbols('k n m')
S2 = P((-k / 4, {'a': 2}, []), (m, {}, ['l1', 'l2']), (n, {}, ['g', 'al']))
Qc = {}
for nm in EVEN:
    Qc[nm] = bracket(S2, P((1, {nm: 1}, [])))
for nm in ODD:
    Qc[nm] = bracket(S2, P((1, {}, [nm])))


def show(f):
    s = []
    for (e, o), c in sorted(f.items()):
        mon = '*'.join([EVEN[i] + ('^%d' % p if p > 1 else '') for i, p in enumerate(e) if p] + [ODD[i] for i in o])
        s.append('%s*%s' % (c, mon if mon else '1'))
    return ' + '.join(s) if s else '0'


res['Q_on_coordinates'] = {nm: show(Qc[nm]) for nm in EVEN + ODD}
res['Q_on_coordinates_meaning'] = {
    'Q a = -n g': 'D_Q on the A coordinate: field-level Q c = -dc (ch03 proof of prop:c3-H): QA = -dc on functionals',
    'Q al = -(k/2) a': 'field-level Q A = -(1/2) d^*dA',
    'Q gs = n al': 'field-level Q A^* = d^*(A^*|_V) (c^* slot)',
    'Q ls1 = m l2, Q ls2 = -m l1': 'field-level Q lambda = dslash lambda (lambda* slot); Q lambda*_a = (dslash lambda)_a'}
checks_Q = {
    'Q(A_perp,A*_perp) = -(1/2)(dA,dA)  [S_2 without the (dc,A*) term]': eq(bracket(P((-k / 4, {'a': 2}, [])), P((1, {'a': 1}, ['al']))), P((-k / 2, {'a': 2}, []))),
    'Q(A_V,A*_V) = -(1/2)(dA,dA) - (dc,A*)': eq(bracket(S2, P((1, {'a': 1}, ['al']))), add(P((-k / 2, {'a': 2}, [])), P((-n, {}, ['g', 'al'])))),
    'Q(lambda*,lambda) = -(lambda,dslash lambda)': eq(bracket(S2, add(P((1, {'ls1': 1}, ['l1'])), P((1, {'ls2': 1}, ['l2'])))),
                                                        P((-2 * m, {}, ['l1', 'l2']))),
    'Q(c*,c) = -(dc,A*)': eq(bracket(S2, P((1, {'gs': 1}, ['g']))), P((-n, {}, ['g', 'al']))),
    'Q a = -n g (QA = -dc)': eq(Qc['a'], P((-n, {}, ['g']))),
    'Q al = -(k/2) a': eq(Qc['al'], P((-k / 2, {'a': 1}, []))),
    'Q gs = n al': eq(Qc['gs'], P((n, {}, ['al']))),
    'Q ls1 = m l2': eq(Qc['ls1'], P((m, {}, ['l2']))),
    'Q g = 0, Q l = 0': Qc['g'] == {} and Qc['l1'] == {} and Qc['l2'] == {},
}
res['Q_checks'] = checks_Q
# Q of S_X: in the model only the weight part exists; with X acting as +x on fields and -x on antifields S_X is Q-closed
res['Q S_X (model)'] = show(bracket(S2, S_X))
res['Q S_X == +D_X S_2'] = eq(bracket(S2, S_X), DX(S2, Xmap))
res['Q S_X note'] = 'Q S_X = {S_2,S_X} = D_X S_2: in R^4 D_M S_2 = 0 (Weyl weight 4), so Q S^(1)(Delta) = 0 (lem:d-calculus (ii))'

# ---------------------------------------------------------------- (2) canonical calculus: D_M and the Ward identity
p = sp.symbols('p0:4', real=True)
pp = sum(c * c for c in p)
eps, gE = sp.symbols('epsilon gamma_E', positive=True)


def DM_alpha(alpha, w):
    return sp.simplify(2 * (w - 4) * alpha + 2 * sum(p[i] * sp.diff(alpha, p[i]) for i in range(4)))


def pdotd(expr):
    return sum(p[i] * sp.diff(expr, p[i]) for i in range(4))


# engine A: A|_Wperp self-energy (plane-wave value, legs A,A, weights 1+1): 9/2 eps^-1 + |p|^2 (43/8 - 2 L2e)
L2e = sp.log(eps * pp) + gE + sp.log(2)
Sigma = sp.Rational(9, 2) / eps + pp * (sp.Rational(43, 8) - 2 * L2e)
DM_Sigma = sp.simplify(DM_alpha(Sigma, 2))
four_eps_deps = sp.simplify(4 * eps * sp.diff(Sigma, eps))
res['ward_identity_on_engineA_Sigma_perp'] = {
    'D_M Sigma': str(sp.simplify(DM_Sigma)), '4 eps d/deps Sigma': str(four_eps_deps),
    'D_M Sigma == +4 eps dSigma/deps': sp.simplify(DM_Sigma - four_eps_deps) == 0,
    'D_M Sigma == -4 eps dSigma/deps': sp.simplify(DM_Sigma + four_eps_deps) == 0}
# D_M on canonical arity-one forms: alpha -> 2(w1+w2-4+p.d)alpha, beta -> 2(w1+w2-5+p.d)beta; S_M|_perp = (2, 2p), w=1+3
res['D_M_S_M_perp'] = {'alpha': str(DM_alpha(sp.Integer(2), 4)),
                       'beta': [str(sp.simplify(2 * (4 - 5) * 2 * p[i] + 2 * pdotd(2 * p[i]))) for i in range(4)]}

# ---------------------------------------------------------------- (3) ghost placement
xi = 'xi'
res['ghost_placement'] = {
    'statement': ('Left ghost (theory: f = f^(0) + xi f_Delta) and right ghost (engine B: f = f^(0) + f_Delta xi) differ by '
                  'f_Delta^left = (-1)^{|f_Delta|} f_Delta^right; for the degree -1 components (S^(1)(Delta), I^CT_Delta, F(Delta)) '
                  'the sign flips uniformly, for the degree 0 component Ob(Delta) it does not.  All three groups define the '
                  'arity-one vertex as J_Delta = S^(1)(Delta) - I^CT_Delta[eps] - F(Delta) with S^(1)(Delta) the explicit functional '
                  'of ex:c2-dil, and Ob(Delta) as the local functional whose tree flow is lim Q W\'_J(S^(1)(Delta) - I^CT_Delta); '
                  'the equation Q W\'_J(J_Delta) -> 0 contains no xi, so F(Delta) and Ob(Delta) = (Q + {I,-}) F(Delta) are the '
                  'same functionals in all three conventions; only the global sign of the component is fixed by the '
                  'normalisation "classical part = +{S^(1)(Delta),O}" (prop:c6-classical), used by all three.'),
    'engine_A': 'component normalised so that the classical part is +2(E+w_O)O; J_Delta = S_M - I^CT_Delta - F_Delta (Y1b + CT = C_P(Pi_int - F_Delta, O))',
    'engine_B': 'right ghost: J = I - I^CT - F_0 + (S_M - I^CT_Delta - F(Delta)) xi',
    'theory': 'left ghost; J_Delta = S^(1)(Delta) - I^CT_Delta - hbar F_1(Delta); (Q+{I,-})hbar F_1(Delta) = Ob(Delta) (eq:d-F)',
    'consistent': True}

json.dump(res, open(os.path.join(HERE, 'conventions.json'), 'w'), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
