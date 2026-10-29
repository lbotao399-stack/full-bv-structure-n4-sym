"""Engine B, Step 2, Part D4: one-loop graphs with one vertex S_M = S^(1)(Delta) (and optionally one O vertex), with edges
P(eps,infty) and at most one special edge Khat = K_eps - K_0.

S_M = int (v.A, A*) - int (v.lambda, lambda*) - int (v.c, c*),  v = 2E (Example ex:c2-dil), v.A = 2E(A)+2A, v.lambda = 2E(lambda)
+3 lambda, v.c = 2E(c).  Explicit position dependence: S_M = S_w + S_E with
   S_w = 2 int(A,A*) - 3 int(lambda,lambda*)            (translation invariant),
   S_E = sum_mu (-i d/dkappa_mu)|_{kappa=0} T_mu[kappa],  T_mu[kappa] = 2 int e^{i kappa.x} [(d_mu A,A*) - (d_mu lambda,lambda*) - (d_mu c, c*)]
(int x^mu f = -i d/dkappa_mu int e^{i kappa.x} f at kappa = 0).  T_mu[kappa] is a translation-invariant vertex injecting kappa.
Ghost of Delta: an odd parameter xi carried by an extra odd 'leg' of species 'xi' (dimension 1, momentum = injection),
placed LAST in the leg order: the term is S_M xi (right-ghost convention = cochain evaluation of the manuscript, Def c2-dce),
so every coefficient returned below is directly the component along Delta.
O = int omega (iota_w1 A, iota_w2 A): even 'leg' of species 'om' (dimension 1) carrying the momentum k of omega.

Word formula (toy-verified, engine_b/toy_r2.py, step2 ia_engine.py):
  P only      : coefficient of (odd r's in leg order) of W_1 = sum_words (1/2n) rsign colour STr[P Pi h_1 ... P Pi h_n]
  one Khat    : coefficient in hbar Delta_Khat W + 1/2{W,W}_Khat (1PI part) =
                sum_words sum_pos (1/2n) rsign (-1)^{odd legs before pos} colour STr[..Khat at pos..].
Explicit component bases live only in code (engine_b/realization.py)."""
import itertools, sys, os
from fractions import Fraction as Fr
from math import comb
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
EB = os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b'))
for pth in (EB, HERE):
    if pth not in sys.path:
        sys.path.insert(0, pth)
import vertices as V
import loops as L
import r2_engine as R
import ia_engine as IA

I = sp.I
q = L.q
V.DIM['xi'] = 1; V.PAR['xi'] = 1
V.DIM['om'] = 1; V.PAR['om'] = 0
INT6 = ['A', 'lam', 'c', 'As', 'lams', 'cs']
PARD = {s: V.PAR[s] for s in INT6}
Mp = V.Mp

def _lampair(l, ls):
    return (sp.Matrix(l['u']).T * Mp * sp.Matrix(ls['u']))[0]

# ------------------------------------------------------------------ new terms (slot order: fields, antifields, then xi)
def _SwA(s):
    A, As, x = s
    return sum(A['u'][M] * As['u'][M] for M in range(10)) * x['u'][0]
def _Swlam(s):
    l, ls, x = s
    return _lampair(l, ls) * x['u'][0]
V.TERMS['SwA'] = (sp.Integer(2), ('A', 'As', 'xi'), _SwA, ('k', (0, 1)))
V.TERMS['Swlam'] = (sp.Integer(-3), ('lam', 'lams', 'xi'), _Swlam, ('k', (0, 1)))
for mu in range(4):
    V.TERMS[f'TA{mu}'] = (sp.Integer(2), ('A', 'As', 'xi'),
                         (lambda m: (lambda s: I * s[0]['p'][m] * sum(s[0]['u'][M] * s[1]['u'][M] for M in range(10)) * s[2]['u'][0]))(mu),
                         ('k', (0, 1)))
    V.TERMS[f'Tlam{mu}'] = (sp.Integer(-2), ('lam', 'lams', 'xi'),
                           (lambda m: (lambda s: I * s[0]['p'][m] * _lampair(s[0], s[1]) * s[2]['u'][0]))(mu), ('k', (0, 1)))
    V.TERMS[f'Tc{mu}'] = (sp.Integer(-2), ('c', 'cs', 'xi'),
                         (lambda m: (lambda s: I * s[0]['p'][m] * s[0]['u'][0] * s[1]['u'][0] * s[2]['u'][0]))(mu), ('k', (0, 1)))
AD = list(sp.symbols('AD0:4'))           # symbolic direction a of the translation letter T_a = sum_mu a_mu T_mu
def _adot(p):
    return sum(AD[m] * p[m] for m in range(4))
V.TERMS['TAa'] = (sp.Integer(2), ('A', 'As', 'xi'),
                  lambda s: I * _adot(s[0]['p']) * sum(s[0]['u'][M] * s[1]['u'][M] for M in range(10)) * s[2]['u'][0], ('k', (0, 1)))
V.TERMS['Tlama'] = (sp.Integer(-2), ('lam', 'lams', 'xi'), lambda s: I * _adot(s[0]['p']) * _lampair(s[0], s[1]) * s[2]['u'][0], ('k', (0, 1)))
V.TERMS['Tca'] = (sp.Integer(-2), ('c', 'cs', 'xi'), lambda s: I * _adot(s[0]['p']) * s[0]['u'][0] * s[1]['u'][0] * s[2]['u'][0], ('k', (0, 1)))
SW_TERMS = ['SwA', 'Swlam']
TA_TERMS = ['TAa', 'Tlama', 'Tca']
def T_TERMS(mu):
    return [f'TA{mu}', f'Tlam{mu}', f'Tc{mu}']

W1 = None; W2 = None
def set_O(w1, w2):
    """O = int omega (iota_w1 A, iota_w2 A), w1, w2 in W_R^perp (10-vectors with V-part 0)."""
    global W1, W2
    W1 = [sp.Integer(x) for x in w1]; W2 = [sp.Integer(x) for x in w2]
    def _O(s):
        a, b, o = s
        return sum(W1[M] * a['u'][M] for M in range(10)) * sum(W2[M] * b['u'][M] for M in range(10)) * o['u'][0]
    V.TERMS['O'] = (sp.Integer(1), ('A', 'A', 'om'), _O, ('k', (0, 1)))

# ------------------------------------------------------------------ letters
KS = sp.symbols('hk0:4'); KPS = sp.symbols('hkp0:4')
class Letter:
    """Hessian of the given terms with the given external legs, colour stripped (vertices.hessian_blocks)."""
    def __init__(self, kind, legs, leg_ids, terms):
        self.kind = kind; self.leg_ids = tuple(leg_ids)
        pairs = [(a, b) for a in INT6 for b in INT6]
        bl = V.hessian_blocks(terms, [legs[i] for i in leg_ids], list(KS), list(KPS), pairs)
        self.blocks = {}
        for (a, b), d in bl.items():
            for key, Mx in d.items():
                if not Mx.is_zero_matrix:
                    self.blocks.setdefault((a, b), {})[key] = Mx
    def block(self, a, rr, b, cc, key, kv, kpv):
        Mx = self.blocks[(a, b)][key]
        sub = Mx.extract(list(rr), list(cc))
        rep = {KS[m]: kv[m] for m in range(4)}; rep.update({KPS[m]: kpv[m] for m in range(4)})
        sub = sub.xreplace(rep).applyfunc(sp.expand)
        return sub * (-1) ** PARD[a]

def colour_M(kind, ids, key):
    if kind in ('SM', 'O'):
        assert key == 'k'
        return {(): Fr(1)}
    return R.letter_M(kind, ids, key)

# ------------------------------------------------------------------ word enumeration
def enumerate_terms(legs, partition, mode, colour_legs, special_filter=None):
    """legs: list of dicts (sp, u, p); partition: list of (kind, ids, terms) letters; mode 'P' (all edges P) or 'K'
    (exactly one Khat edge).  Returns list of records dict(fac, Num, curs, s, word, edges)."""
    odd = [V.PAR[l['sp']] for l in legs]
    lets = {}
    for (kind, ids, terms) in partition:
        lets[(kind, ids)] = Letter(kind, legs, ids, terms)
    n = len(partition)
    out = []
    for perm in itertools.permutations(range(n)):
        word = [partition[i] for i in perm]
        curs = [[sp.Integer(0)] * 4]
        for (kind, ids, terms) in word:
            inj = [sum(legs[i]['p'][m] for i in ids) for m in range(4)]
            curs.append([sp.expand(curs[-1][m] + inj[m]) for m in range(4)])
        assert all(sp.expand(x) == 0 for x in curs[-1]), curs[-1]
        curs = curs[:n]
        odd_seq = [i for (kind, ids, terms) in word for i in ids if odd[i]]
        rsign = R.perm_parity(odd_seq)
        L_ = [lets[(kind, ids)] for (kind, ids, terms) in word]
        positions = [None] if mode == 'P' else list(range(n))
        for pos in positions:
            if pos is None:
                w = sp.Rational(1, 2 * n) * rsign
            else:
                odd_before = sum(odd[i] for (kind, ids, terms) in word[:pos] for i in ids)
                w = sp.Rational(1, 2 * n) * rsign * (-1) ** odd_before
            choices = [list(IA.KEDGE) if e == pos else list(IA.PEDGE) for e in range(n)]
            for edges in itertools.product(*choices):
                spec = [IA.KEDGE[edges[e]] if e == pos else IA.PEDGE[edges[e]] for e in range(n)]
                ok = all((spec[kk][2], spec[(kk + 1) % n][0]) in L_[kk].blocks for kk in range(n))
                if not ok:
                    continue
                if special_filter is not None and not special_filter(word, edges):
                    continue
                keysets = [list(L_[kk].blocks[(spec[kk][2], spec[(kk + 1) % n][0])].keys()) for kk in range(n)]
                for keys in itertools.product(*keysets):
                    Ms = [colour_M(word[kk][0], word[kk][1], keys[kk]) for kk in range(n)]
                    col = R.reduce_traces(R.colour_of_word(Ms), colour_legs)
                    ckey = {2: 'cadj k', 3: 'cadj f'}[len(colour_legs)]
                    assert set(col) <= {ckey}, col
                    cf = col.get(ckey, 0)
                    if cf == 0:
                        continue
                    mats = []
                    for kk in range(n):
                        line_mom = [q[m] - curs[kk][m] for m in range(4)]
                        if kk == pos:
                            mats.append(IA.K_fibre(edges[kk]))
                        else:
                            mats.append(IA.P_numerator(edges[kk], line_mom))
                        kv = [-(q[m] - curs[kk][m]) for m in range(4)]
                        nxt = curs[kk + 1] if kk + 1 < n else [0] * 4
                        kpv = [q[m] - nxt[m] for m in range(4)]
                        a, rr = spec[kk][2], spec[kk][3]
                        b, cc = spec[(kk + 1) % n][0], spec[(kk + 1) % n][1]
                        mats.append(L_[kk].block(a, rr, b, cc, keys[kk], kv, kpv))
                    Mx = mats[0]
                    for X in mats[1:]:
                        Mx = (Mx * X).applyfunc(sp.expand)
                        if Mx.is_zero_matrix:
                            break
                    if Mx.is_zero_matrix:
                        continue
                    Num = sp.expand((-1) ** PARD[spec[0][0]] * Mx.trace())
                    if Num == 0:
                        continue
                    out.append(dict(fac=w * sp.Rational(cf.numerator, cf.denominator), Num=Num, curs=[list(c) for c in curs],
                                    s=pos, word=[(k_, i_) for (k_, i_, t_) in word], edges=edges))
    return out

# ------------------------------------------------------------------ exact sphere averages and the log eps coefficient of P-only loops
def _dfact(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r

def sphere_avg(poly_expr):
    """exact average over the unit sphere S^3 of a polynomial in q0..q3."""
    poly = sp.Poly(sp.expand(poly_expr), *q)
    tot = 0
    for mon, cf in poly.terms():
        if any(e % 2 for e in mon):
            continue
        num = 1
        for e in mon:
            num *= _dfact(e - 1) if e > 0 else 1
        n2 = sum(mon) // 2
        den = 1
        for j in range(n2):
            den *= (4 + 2 * j)
        tot += cf * sp.Rational(num, den)
    return sp.expand(tot)

LAM = sp.Symbol('lam_scale')
def log_coeff(Num, curs, omega, ext_syms):
    """coefficient of log eps (units (16 pi^2)^-1) of int d^4q/(2pi)^4 Num prod_e e^{-eps|q-curs_e|^2}/|q-curs_e|^2:
    -< [lambda^omega] Num(q, lambda ext) prod_e |q - lambda curs_e|^-2 >_{|q|=1} (Lemma c5-log, with all external momenta
    -- every symbol in ext_syms and every number in curs -- scaled by lambda)."""
    sub = {s: LAM * s for s in ext_syms}
    N = sp.expand(Num.xreplace(sub))
    q2 = sum(x * x for x in q)
    E = len(curs)
    Dtot = E + omega                       # multiply everything by q2^Dtot: 1/|q - lam c|^2 -> sum_k X^k q2^(-(k+1))
    # series per edge as list of (power of q2 removed, polynomial)
    ser_list = []
    for c in curs:
        cl = [sp.expand(sp.sympify(x).xreplace(sub)) for x in c]
        X = sp.expand(2 * sum(q[m] * cl[m] for m in range(4)) - sum(x * x for x in cl))
        ser_list.append([(k + 1, sp.expand(X ** k)) for k in range(omega + 1)])
    tot = 0
    for combo in itertools.product(*ser_list):
        pw = sum(c_[0] for c_ in combo)
        if pw > Dtot:
            continue
        prod = 1
        for c_ in combo:
            prod = sp.expand(prod * c_[1])
        tot += prod * q2 ** (Dtot - pw)
    tot = sp.expand(N * sp.expand(tot))
    htil = tot.coeff(LAM, omega)
    htil = sp.expand(htil)
    if htil == 0:
        return sp.Integer(0)
    degs = set(sum(m_) for m_ in sp.Poly(htil, *q).monoms())
    assert degs == {2 * Dtot - 4}, (degs, Dtot)
    return -sphere_avg(htil)

def sphere_avg_ratio(num, den):
    q2 = sum(x * x for x in q)
    # den = c * q2^N
    N = sp.Poly(den, *q).total_degree() // 2
    c = sp.expand(den / q2 ** N)
    c = sp.simplify(c)
    assert not any(c.has(x) for x in q), den
    numh = sp.Poly(sp.expand(num), *q)
    degs = set(sum(m) for m in numh.monoms())
    assert all(d == 2 * N - 4 for d in degs), (degs, N)
    return sp.expand(sphere_avg(num) / c)

# ------------------------------------------------------------------ exact masters for one special edge
def master_special(a_s, a_rest, m):
    """int_0^1 du_s u_s^a_s int_{[1,oo)^r} prod u_i^a_i (u_s + sum u_i)^-m ;  returns {'1': Fr, 'L': Fr} with L =
    log(c+1/c), c = r (number of [1,oo) variables)."""
    terms = {(a_s, tuple(a_rest), 0, m): Fr(1)}     # (a_s, rest exps, c, e): coeff * u_s^a_s prod u^a (c + u_s + sum u)^-e
    rest = list(a_rest)
    for idx in range(len(rest) - 1, -1, -1):
        new = {}
        for (as_, ar, c, e), cf in terms.items():
            a = ar[idx]
            ar2 = ar[:idx]
            for i in range(a + 1):
                assert e - i - 1 >= 1, ('divergent master', a_s, a_rest, m)
                ci = Fr(comb(a, i) * 1, e - i - 1)
                for j in range(a - i + 1):
                    cj = Fr(comb(a - i, j) * (-1) ** j)
                    key = (as_, ar2, c + 1, e - i - j - 1)
                    new[key] = new.get(key, Fr(0)) + cf * ci * cj
        terms = {k: v for k, v in new.items() if v != 0}
    out = {'1': Fr(0), 'L': Fr(0)}
    for (as_, ar, c, e), cf in terms.items():
        # int_0^1 u^as (c+u)^-e du, v = c+u : sum_k C(as,k) (-c)^{as-k} int_c^{c+1} v^{k-e} dv
        for k in range(as_ + 1):
            ck = Fr(comb(as_, k)) * Fr(-c) ** (as_ - k)
            pw = k - e
            if pw == -1:
                out['L'] += cf * ck
            else:
                out['1'] += cf * ck * (Fr(c + 1) ** (pw + 1) - Fr(c) ** (pw + 1)) / (pw + 1)
    return out, len(rest)

BET = sp.symbols('bb0:6')
def special_coeffs(Num, curs, s, kap_op=None):
    """one special edge at position s (Schwinger time in [0,eps]), others P(eps,infty).  Returns {n: {'1','L'}: value}
    for n = -1, 0 (coefficients of eps^n, units (16 pi^2)^-1), L = log((E)/(E-1)).  kap_op: function applied to the
    polynomial integrand (in the injection symbols) before the Schwinger integration (e.g. d/dkappa at 0)."""
    E = len(curs)
    ks = [q[m] - curs[s][m] for m in range(4)]
    Nt = sp.expand(-Num * sum(x * x for x in ks))
    shift = [sum(BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
    C = L.gaussian_moments(Nt, shift)
    P2 = lambda v: sum(x * x for x in v)
    rho = sp.expand(sum(BET[e] * P2(curs[e]) for e in range(E)) - P2(shift))
    order = [s] + [e for e in range(E) if e != s]
    out = {}
    for n in (-1, 0):
        acc = {'1': 0, 'L': 0}
        for k, Ck in C.items():
            j = k + 2 - E + n
            if j < 0:
                continue
            expr = sp.expand(Ck * rho ** j) * sp.Rational((-1) ** j, 2 ** k * sp.factorial(j))
            if kap_op is not None:
                expr = sp.expand(kap_op(expr))
            if expr == 0:
                continue
            poly = sp.Poly(expr, *BET[:E])
            for mon, cf in poly.terms():
                al = [mon[e] for e in order]
                mm = E - n + sum(al)
                J, _ = master_special(al[0], al[1:], mm)
                acc['1'] += cf * sp.Rational(J['1'].numerator, J['1'].denominator)
                acc['L'] += cf * sp.Rational(J['L'].numerator, J['L'].denominator)
        out[n] = {kk: sp.expand(v) for kk, v in acc.items()}
    return out

def special_coeffs_multi(Num, curs, s, ops):
    """as special_coeffs for several kappa-operations at once (one Gaussian integration): returns {opname: {n: {...}}}."""
    E = len(curs)
    ks = [q[m] - curs[s][m] for m in range(4)]
    Nt = sp.expand(-Num * sum(x * x for x in ks))
    shift = [sum(BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
    C = L.gaussian_moments(Nt, shift)
    P2 = lambda v: sum(x * x for x in v)
    rho = sp.expand(sum(BET[e] * P2(curs[e]) for e in range(E)) - P2(shift))
    order = [s] + [e for e in range(E) if e != s]
    res = {nm: {} for nm in ops}
    cache = {}
    for n in (-1, 0):
        base = 0
        for k, Ck in C.items():
            j = k + 2 - E + n
            if j < 0:
                continue
            base += sp.expand(Ck * rho ** j) * sp.Rational((-1) ** j, 2 ** k * sp.factorial(j))
        base = sp.expand(base)
        for nm, op in ops.items():
            expr = sp.expand(op(base))
            acc = {'1': 0, 'L': 0}
            if expr != 0:
                poly = sp.Poly(expr, *BET[:E])
                for mon, cf in poly.terms():
                    al = tuple(mon[e] for e in order)
                    mm = E - n + sum(al)
                    key = (al, mm)
                    if key not in cache:
                        cache[key] = master_special(al[0], list(al[1:]), mm)[0]
                    J = cache[key]
                    acc['1'] += cf * sp.Rational(J['1'].numerator, J['1'].denominator)
                    acc['L'] += cf * sp.Rational(J['L'].numerator, J['L'].denominator)
            res[nm][n] = {kk: sp.expand(v) for kk, v in acc.items()}
    return res

def kap_derivative_op(kap, mu):
    """-i d/dkappa_mu at kappa = 0."""
    def op(expr):
        return sp.expand(-I * sp.diff(expr, kap[mu]).xreplace({k: 0 for k in kap}))
    return op

def kap_zero_op(kap):
    def op(expr):
        return sp.expand(expr.xreplace({k: 0 for k in kap}))
    return op
