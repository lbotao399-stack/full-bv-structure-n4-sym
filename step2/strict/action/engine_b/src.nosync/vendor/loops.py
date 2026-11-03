"""Engine B loop machinery.
One-loop formula (0-dim verified, toy_verify.py, with the agreed contraction rule d_P(f_X f_Y) = P^{XY}):
   W_1 = (hbar/2) sum_n (1/n) STr[(P Pi H)^n],  H_jk = D_j D_k I,  Pi = parity operator,  STr = sum (-1)^{|i|}.
Momentum space: STr[P Pi H_1 P Pi H_2 ...] -> int d^4q/(2pi)^4 STr[ P(q) Pi h_1(-q, q-p_1) P(q-p_1) Pi h_2(...) ... ].
Propagator: P = -(Pi H_free)^{-1} (Feynman-gauge A sector, L-restricted ghost sector) times e^{-eps q^2}:
   P(q) = N(q) * S_eps(q),  S_eps(q) = e^{-eps q^2}/q^2 = int_eps^infty e^{-t q^2} dt.
UV boundary (method B): eps d/deps of a graph = - sum_e (graph with t_e pinned at eps); pinned graphs are expanded in eps
(Taylor in external momenta) with t_f = eps*u_f, u_f in [1,infty), u-integrals exact."""
import itertools
import sympy as sp
import vertices as V
from realization import Rb

I = sp.I
q = sp.symbols('q0:4'); l = sp.symbols('l0:4')
eps = sp.Symbol('epsilon', positive=True)
INT = ['A', 'lam', 'c', 'As']          # internal species with nonzero propagator
OFF = {}; _o = 0
for s in INT:
    OFF[s] = _o; _o += V.DIM[s]
NTOT = _o

# ---------- propagator numerators from inversion of the free Hessian ----------
def _free_hessian(spcs, terms, qv):
    return V.hessian_blocks(terms, [], [-x for x in qv], list(qv), spcs)

def derive_propagators():
    qs = sp.symbols('Q0:4'); q2 = sum(x * x for x in qs)
    out = {}
    # lambda
    h = _free_hessian([('lam', 'lam')], ['kinlam'], qs)[('lam', 'lam')]['k']
    Plam = sp.simplify(-((-1) * h).inv())
    expect = I * sum((qs[m] * Rb[m] for m in range(4)), sp.zeros(16, 16)) * V.Mpinv / q2
    out['lam_matches_iRbar_Minv'] = sp.simplify(Plam - expect) == sp.zeros(16, 16)
    out['N_lam'] = sp.simplify(Plam * q2)
    # A (Feynman gauge: free -1/4(dA,dA) plus -1/4 (d^*A_V)^2, i.e. -1/4 |dA|^2 on all ten components)
    V.TERMS['gfA'] = (sp.Rational(-1, 4), ('A', 'A'),
                      lambda s: (I * V.dot4(s[0]['p'], s[0]['u'])) * (I * V.dot4(s[1]['p'], s[1]['u'])), ('k', (0, 1)))
    h = _free_hessian([('A', 'A')], ['kinA', 'gfA'], qs)[('A', 'A')]['k']
    PA = sp.simplify(-(h).inv())
    out['A_is_2_over_q2'] = sp.simplify(PA - 2 * sp.eye(10) / q2) == sp.zeros(10, 10)
    del V.TERMS['gfA']
    # ghosts: on L = {c} + {A*_V = d cbar}, S_2 = (dc, d cbar)
    kin = lambda s: sum((I * s[0]['p'][m]) * (I * s[1]['p'][m]) for m in range(4)) * s[0]['u'][0] * s[1]['u'][0]
    V.DIM['cb'] = 1; V.PAR['cb'] = 1
    V.TERMS['ghL'] = (sp.Integer(1), ('c', 'cb'), kin, ('k', (0, 1)))
    hb = _free_hessian([('c', 'c'), ('c', 'cb'), ('cb', 'c'), ('cb', 'cb')], ['ghL'], qs)
    Hm = sp.Matrix([[hb[('c', 'c')].get('k', sp.zeros(1, 1))[0], hb[('c', 'cb')]['k'][0]],
                    [hb[('cb', 'c')]['k'][0], hb[('cb', 'cb')].get('k', sp.zeros(1, 1))[0]]])
    Pg = sp.simplify(-((-1) * Hm).inv())
    P_c_cb = Pg[0, 1]; P_cb_c = Pg[1, 0]
    P_c_As = [sp.simplify(P_c_cb * (-I * qs[m])) for m in range(4)]     # A*_mu(y) = d_mu cbar(y)
    P_As_c = [sp.simplify(I * qs[m] * P_cb_c) for m in range(4)]        # A*_mu(x) = d_mu cbar(x)
    out['ghost_c_As'] = [sp.simplify(x * q2) for x in P_c_As]
    out['ghost_As_c'] = [sp.simplify(x * q2) for x in P_As_c]
    del V.TERMS['ghL']; del V.DIM['cb']; del V.PAR['cb']
    out['qs'] = qs
    return out

_PROP = None
def prop_numerator(qv, mask=None):
    """Full NTOT x NTOT numerator N(q) with P(q) = N(q) e^{-eps q^2}/q^2.  mask: set of allowed blocks
    ('AV','AV'), ('Ap','Ap'), ('lam','lam'), ('c','As'), ('As','c')."""
    global _PROP
    if _PROP is None:
        _PROP = derive_propagators()
    d = dict(zip(_PROP['qs'], qv))
    N = sp.zeros(NTOT, NTOT)
    def allowed(b):
        return mask is None or b in mask
    for M in range(10):
        b = ('AV', 'AV') if M < 4 else ('Ap', 'Ap')
        if allowed(b):
            N[OFF['A'] + M, OFF['A'] + M] = 2
    if allowed(('lam', 'lam')):
        Nl = _PROP['N_lam'].subs(d)
        for i in range(16):
            for j in range(16):
                N[OFF['lam'] + i, OFF['lam'] + j] = sp.expand(Nl[i, j])
    for m in range(4):
        if allowed(('c', 'As')):
            N[OFF['c'], OFF['As'] + m] = sp.expand(_PROP['ghost_c_As'][m].subs(d))
        if allowed(('As', 'c')):
            N[OFF['As'] + m, OFF['c']] = sp.expand(_PROP['ghost_As_c'][m].subs(d))
    return N

PIdiag = []
for s in INT:
    PIdiag += [(-1) ** V.PAR[s]] * V.DIM[s]

def full_hessian(terms, ext_legs, k, kp, key_weights):
    """NTOT x NTOT colour-stripped Hessian, colour keys combined with the given weights (colour factor / (c_adj (c1,c2)))."""
    pairs = [(a, b) for a in INT for b in INT]
    blocks = V.hessian_blocks(terms, ext_legs, k, kp, pairs)
    H = sp.zeros(NTOT, NTOT)
    for (a, b), bl in blocks.items():
        for key, M in bl.items():
            w = key_weights.get(key, None)
            if w is None:
                raise ValueError('unexpected colour key ' + key)
            if w == 0:
                continue
            for i in range(M.rows):
                for j in range(M.cols):
                    if M[i, j] != 0:
                        H[OFF[a] + i, OFF[b] + j] += w * M[i, j]
    return H.applyfunc(sp.expand)

def PiH(H):
    return sp.Matrix(H.rows, H.cols, lambda i, j: PIdiag[i] * H[i, j])

def str_product(mats):
    """STr of the product of the given NTOT x NTOT matrices."""
    M = mats[0]
    for X in mats[1:-1]:
        M = (M * X).applyfunc(sp.expand)
    last = mats[-1]
    return sp.expand(sum(PIdiag[i] * M[i, j] * last[j, i] for i in range(NTOT) for j in range(NTOT)
                         if M[i, j] != 0 and last[j, i] != 0))

# ---------- Gaussian integration and UV-boundary expansion ----------
def gaussian_moments(Num, shift):
    """Num(q) with q = l + shift; returns {k: coeff} with int d^4l e^{-T l^2} Num / int e^{-T l^2} = sum_k coeff_k (2T)^{-k}."""
    sub = {q[m]: l[m] + shift[m] for m in range(4)}
    expr = sp.expand(sp.sympify(Num).subs(sub, simultaneous=True))
    poly = sp.Poly(expr, *l)
    out = {}
    for mon, cf in poly.terms():
        if any(e % 2 for e in mon):
            continue
        fac = 1
        for e in mon:
            fac *= sp.factorial2(e - 1) if e > 0 else 1
        k = sum(mon) // 2
        out[k] = sp.expand(out.get(k, 0) + fac * cf)
    return out

def boundary_singular(Num, lines, keep_check=True):
    """Graph X = int d^4q/(2pi)^4 Num(q) prod_e S_eps(q - P_e).
    Returns (a, b) with X = a/eps + b log eps + O(1), computed from the pinned-edge (UV boundary) graphs:
    eps dX/deps = -sum_e eps^E int_{[1,inf)^{E-1}} f(t_e = eps, t_f = eps u_f) = -(A_{-1}/eps + A_0) + o(1),  a = A_{-1}, b = -A_0."""
    E = len(lines)
    ts = sp.symbols(f't0:{E}', positive=True)
    betas = sp.symbols(f'beta0:{E}')
    T = sp.Symbol('T', positive=True)
    shift = [sum(betas[e] * lines[e][m] for e in range(E)) for m in range(4)]
    C = gaussian_moments(Num, shift)               # f = (16 pi^2)^{-1} T^{-2} sum_k C_k (2T)^{-k} e^{-R}
    P2 = lambda v: sum(x * x for x in v)
    # R/T = sum_e beta_e P_e^2 - (sum beta_e P_e)^2
    RoverT = sp.expand(sum(betas[e] * P2(lines[e]) for e in range(E)) - P2(shift))
    us = sp.symbols(f'u0:{E}', positive=True)
    Am1 = 0; A0 = 0
    for e in range(E):
        that = [sp.Integer(1) if f == e else us[f] for f in range(E)]
        Th = sum(that)
        bsub = {betas[f]: that[f] / Th for f in range(E)}
        free_u = [us[f] for f in range(E) if f != e]
        # eps^E f(eps t) = (16pi^2)^{-1} sum_k C_k 2^{-k} eps^{E-2-k} Th^{-2-k} (1 - eps Th RoverT + ...)
        for k, Ck in C.items():
            base = sp.expand(Ck.subs(bsub)) / 2**k * Th**(-2 - k)
            pw = E - 2 - k
            # order eps^{pw}: base ; order eps^{pw+1}: -base*Th*RoverT
            terms = [(pw, base), (pw + 1, -base * Th * RoverT.subs(bsub))]
            for (pp, expr) in terms:
                if pp == -1:
                    Am1 += _uint(expr, free_u)
                elif pp == 0:
                    A0 += _uint(expr, free_u)
                elif pp < -1:
                    raise ValueError('unexpected eps power %d' % pp)
            if pw + 2 <= 0:
                raise ValueError('need higher eps expansion (k=%d, E=%d)' % (k, E))
    fac = 1 / (16 * sp.pi**2)
    return sp.expand(sp.simplify(Am1 * fac)), sp.expand(sp.simplify(-A0 * fac))

def _uint(expr, free_u):
    expr = sp.together(sp.expand(expr))
    if not free_u:
        return expr
    # integrate a rational function of u over [1, oo)^n, after splitting into coefficient monomials of the other symbols
    num, den = sp.fraction(expr)
    num = sp.expand(num)
    others = sorted(num.free_symbols - set(free_u), key=lambda s: s.name)
    if others:
        poly = sp.Poly(num, *others)
        tot = 0
        for mon, cf in poly.terms():
            val = cf / den
            for u in free_u:
                val = sp.integrate(sp.apart(sp.together(val), u) if len(free_u) == 1 else val, (u, 1, sp.oo))
                if val.has(sp.oo) or val.has(sp.zoo) or val.has(sp.nan) or isinstance(val, sp.Integral):
                    raise ValueError('divergent/unevaluated u-integral: ' + str(cf / den))
            tot += sp.Mul(*[s**e for s, e in zip(others, mon)]) * sp.simplify(val)
        return sp.expand(tot)
    val = num / den
    for u in free_u:
        val = sp.integrate(val, (u, 1, sp.oo))
        if val.has(sp.oo) or val.has(sp.zoo) or val.has(sp.nan) or isinstance(val, sp.Integral):
            raise ValueError('divergent/unevaluated u-integral: ' + str(num / den))
    return sp.simplify(val)
