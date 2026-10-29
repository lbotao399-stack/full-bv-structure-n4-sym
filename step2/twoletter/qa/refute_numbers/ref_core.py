"""Referee engine R (adversarial numbers check of R8 two-letter kernels).

Independent of engine A / engine B code.  Built only from the manuscript data:
  * vertices: Ch1 eq c1-I densities / Def c1-vertex table (AAA eq c1-AAA, AAAA eq c1-AAAA, lambda lambda A,
    A c A*, c lambda lambda*, c c c*), read here as ordered densities in coordinate fields;
  * propagator: Rem tl-P / Prop c3-P, p incoming at the FIRST factor:
        A x A: 2 delta,  c x A*_V and A*_V x c: -i p_mu,  lambda x lambda: i rho(p)|_{Sigma* -> Sigma}
    (each times (-,-)^{-1}, times int_eps^T e^{-t p^2} dt);
  * kernel: K(X_a, X_b) = hbar Sing <X_a X_b>_e  (Wick contraction with the dressed propagator, Def def:tl-dressed),
    Grassmann signs by explicit Koszul reordering of the ordered product X_a X_b V_1 ... V_k (each vertex in its
    density order); colour reduced to C_w = (1 (x) w)(-,-)^{-1}, w = ad_{Y_k} ... ad_{Y_1} (Y_1 nearest X_a).
Integration (method different from both engines): after the Gaussian loop integral, every term is a homogeneous
rational function s^a tau^{-N} of the Schwinger times; t = eps s;
  [eps^{-kappa}] (kappa >= 1) = int_{[1,oo)^E} s^a tau^{-N} ds, evaluated EXACTLY by the Laplace representation
      tau^{-N} = Gamma(N)^{-1} int_0^oo x^{N-1} e^{-x tau} dx and int_1^oo s^a e^{-xs} ds = a! e^{-x} sum_k x^{k-a-1}/k!,
  [log eps]  = - (log-Lambda coefficient of int_{[1,Lambda]^E}) = - lim_{delta->0} delta int_{[1,oo)^E} s^a tau^{-N-delta},
      read off the same Laplace representation (the Gamma(delta) pole).
No dimensional regularisation, T = infinity at generic momenta (justified at finite T, Prop prop:tl-map (iii)).
Units of every returned coefficient: hbar/(16 pi^2) times C_w times the ordered product of the output letters.
"""
import itertools
from fractions import Fraction
from math import factorial
import sympy as sp

L = sp.symbols('l0:4')
I = sp.I


def dot4(u, v):
    return sum(u[i] * v[i] for i in range(4))


def emb10(p4):
    """momentum (4-vector along V) as a 10-vector"""
    return list(p4) + [0] * 6


# ----------------------------------------------------------------------------------------------- Schwinger integrals
_JC = {}


def J_orthant(a, N):
    """int_{[1,oo)^E} prod s_e^{a_e} (sum s)^{-N} ds, exact (requires N > |a| + E)."""
    key = (tuple(a), N)
    if key in _JC:
        return _JC[key]
    E = len(a)
    assert N > sum(a) + E
    tot = Fraction(0)
    for ks in itertools.product(*[range(ae + 1) for ae in a]):
        coef = Fraction(1)
        for ae, ke in zip(a, ks):
            coef *= Fraction(factorial(ae), factorial(ke))
        p = sum(ke - ae - 1 for ae, ke in zip(a, ks))
        q = N + p                      # >= 1 by convergence
        assert q >= 1
        tot += coef * Fraction(factorial(q - 1), E ** q)
    tot /= factorial(N - 1)
    _JC[key] = tot
    return tot


def log_coef(a, N):
    """log-Lambda coefficient of int_{[1,Lambda]^E} s^a tau^{-N} for N = |a| + E: the Gamma(delta) pole of the
    Laplace representation, i.e. prod a_e!/(N-1)!."""
    E = len(a)
    assert N == sum(a) + E
    c = Fraction(1)
    for ae in a:
        c *= factorial(ae)
    return c / factorial(N - 1)


def laplacian(P):
    return sum(sp.diff(P, L[m], 2) for m in range(4))


def singular_parts(N, Rs, maxdeg=None):
    """N: sympy polynomial in L (loop momentum l) and other symbols; Rs: list of E shifts R_0=0, R_1, ..., R_k
    (edge j carries l + R_j).  Returns {'log': expr, 1: expr, 2: expr, ...} = coefficients of log eps, eps^{-m}
    times 16 pi^2 (per-edge cutoff t_e >= eps, T = infinity)."""
    E = len(Rs)
    s = sp.symbols('s0:%d' % E)
    w = sp.Symbol('w')                # w = 1/tau
    N = sp.expand(N)
    degl = sp.Poly(N, *L).total_degree() if N != 0 else 0
    if N == 0:
        return {}
    c = [-w * sum(s[j] * Rs[j][mu] for j in range(E)) for mu in range(4)]
    X = sum(s[j] * dot4(Rs[j], Rs[j]) for j in range(E)) - w * sum(
        (sum(s[j] * Rs[j][mu] for j in range(E))) ** 2 for mu in range(4))
    out = {}
    lapN = N
    for j in range(0, degl // 2 + 1):
        if j > 0:
            lapN = sp.expand(laplacian(lapN))
        if lapN == 0:
            break
        nmax = 2 + j - E
        if nmax < 0:
            continue
        numj = sp.expand(lapN.subs({L[m]: c[m] for m in range(4)}, simultaneous=True)) / (factorial(j) * 4 ** j)
        base = sp.expand(numj * w ** (2 + j))
        for n in range(0, nmax + 1):
            term = sp.expand(base * (-X) ** n / factorial(n))
            if term == 0:
                continue
            kappa = 2 + j - E - n
            P = sp.Poly(term, *s, w)
            for monom, coef in P.terms():
                a = monom[:E]
                Nw = monom[E]
                # t-degree |a| - Nw must equal -(E + kappa)
                assert sum(a) - Nw == -(E + kappa), (a, Nw, kappa)
                if kappa >= 1:
                    val = J_orthant(a, Nw)
                    out[kappa] = out.get(kappa, 0) + coef * sp.Rational(val.numerator, val.denominator)
                else:
                    val = log_coef(a, Nw)
                    out['log'] = out.get('log', 0) - coef * sp.Rational(val.numerator, val.denominator)
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


# ----------------------------------------------------------------------------------------------- fields and lines
SIZE = {'A': 10, 'c': 1, 'As': 10, 'l': 16, 'ls': 16, 'cs': 1}
PARITY = {'A': 0, 'c': 1, 'As': 1, 'l': 1, 'ls': 0, 'cs': 0}


class Cliff:
    """sig[M]: rho(e_M)|_Sigma : Sigma -> Sigma^* (bilinear (s_a, rho s_b) = sig[M][a,b]);
    sigb[M]: rho(e_M)|_Sigma^* : Sigma^* -> Sigma (upper indices)."""

    def __init__(self, sig, sigb):
        def ex(m):
            return sp.Matrix(16, 16, lambda i, j: sp.Integer(int(round(complex(m[i][j]).real)))
                             + sp.I * sp.Integer(int(round(complex(m[i][j]).imag))))
        self.sig = [ex(m) for m in sig]
        self.sigb = [ex(m) for m in sigb]

    def sig_v(self, v):
        return sum((v[M] * self.sig[M] for M in range(10)), sp.zeros(16, 16))

    def sigb_v(self, v):
        return sum((v[M] * self.sigb[M] for M in range(10)), sp.zeros(16, 16))


def prop_hat(t1, t2, p, cl=None):
    """P-hat^{t1 t2}(p) as a SIZE[t1] x SIZE[t2] matrix, p = incoming momentum at the first factor (4-vector)."""
    if t1 == 'A' and t2 == 'A':
        return 2 * sp.eye(10)
    if t1 == 'c' and t2 == 'As':
        return sp.Matrix([[(-I * p[m]) if m < 4 else 0 for m in range(10)]])
    if t1 == 'As' and t2 == 'c':
        return sp.Matrix([[(-I * p[m]) if m < 4 else 0] for m in range(10)])
    if t1 == 'l' and t2 == 'l':
        return I * cl.sigb_v(emb10(p))
    return None


PARTNER = {'A': ['A'], 'c': ['As'], 'As': ['c'], 'l': ['l']}


# ----------------------------------------------------------------------------------------------- vertices
def perm_sign(seq, ref):
    """sign of the permutation taking ref -> seq (lists of distinct labels)"""
    idx = [ref.index(x) for x in seq]
    sgn = 1
    for i in range(len(idx)):
        for j in range(i + 1, len(idx)):
            if idx[i] > idx[j]:
                sgn = -sgn
    return sgn


# ordered densities of the Grassmann vertices: (slot types in density order, coefficient tensor builder)
# lambda lambda A: (1/2) lambda^a sig_M[a,b] A_M lambda^b f(lam1, A, lam2)
# A c A*:         A_M c A*_M f(A, c, A*)
# c lambda lambda*: c lambda^a lambda*_a f(c, lam, lam*)
# c c c*:         (1/2) c c c* f(c, c, c*)
DENS = {
    'lAl': ['l', 'A', 'l'],
    'AcAs': ['A', 'c', 'As'],
    'clls': ['c', 'l', 'ls'],
    'cccs': ['c', 'c', 'cs'],
}


def dens_tensor(name, cl):
    """returns function T(i0, i1, i2) -> coefficient for slot components (density order)"""
    if name == 'lAl':
        return lambda a, M, b: sp.Rational(1, 2) * cl.sig[M][a, b]
    if name == 'AcAs':
        return lambda M, z, N: 1 if M == N else 0
    if name == 'clls':
        return lambda z, a, b: 1 if a == b else 0
    if name == 'cccs':
        return lambda x, y, z: sp.Rational(1, 2)
    raise KeyError(name)


def grass_vertex_terms(name, tin, tout, tbg, ybg, cl):
    """All role assignments of a Grassmann-relevant cubic vertex with in-type tin, out-type tout, background type tbg
    (background selector vector ybg of length SIZE[tbg]).  Returns list of (roles_in_density_order, Matrix[in,out],
    colour_sign) with colour_sign relative to (ad_Y Z_in, Z_out)."""
    slots = DENS[name]
    T = dens_tensor(name, cl)
    out = []
    for perm in itertools.permutations(range(3)):
        roles = [None] * 3
        # perm[0] = slot of in, perm[1] = slot of out, perm[2] = slot of bg
        si, so, sb = perm
        if slots[si] != tin or slots[so] != tout or slots[sb] != tbg:
            continue
        roles[si], roles[so], roles[sb] = 'in', 'out', 'bg'
        Mx = sp.zeros(SIZE[tin], SIZE[tout])
        for i in range(SIZE[tin]):
            for o in range(SIZE[tout]):
                acc = 0
                for b in range(SIZE[tbg]):
                    if ybg[b] == 0:
                        continue
                    comp = [None] * 3
                    comp[si], comp[so], comp[sb] = i, o, b
                    acc += ybg[b] * T(*comp)
                Mx[i, o] = acc
        if Mx == sp.zeros(SIZE[tin], SIZE[tout]):
            continue
        csign = perm_sign(roles, ['bg', 'in', 'out'])
        out.append((roles, Mx, csign))
    return out


def AAA_matrix(y, p_in, p_Y, p_out):
    """vertex function eq c1-AAA with legs (1,2,3) = (Y, in, out); f(Y,in,out) = (ad_Y Z_in, Z_out).
    returns 10x10 matrix [in, out] (without the colour element)."""
    P_in, P_Y, P_out = emb10(p_in), emb10(p_Y), emb10(p_out)
    M = sp.zeros(10, 10)
    yd = lambda v: sum(y[k] * v[k] for k in range(10))
    for i in range(10):
        for o in range(10):
            t = y[i] * (P_in[o] - P_Y[o]) + (1 if i == o else 0) * yd([P_out[k] - P_in[k] for k in range(10)]) \
                + (P_Y[i] - P_out[i]) * y[o]
            M[i, o] = -I / 2 * t
    return M


def AAAA_matrix(y1, y2):
    """coefficient of the word ad_{y2} ad_{y1} (y1 nearer letter a) in eq c1-AAAA, legs (in, out, y1, y2):
    (1/2)[y1 y2^T - 2 y2 y1^T + (y1.y2) Id]"""
    y1m, y2m = sp.Matrix(y1), sp.Matrix(y2)
    yy = sum(y1[k] * y2[k] for k in range(10))
    return sp.Rational(1, 2) * (y1m * y2m.T - 2 * y2m * y1m.T + yy * sp.eye(10))


# ----------------------------------------------------------------------------------------------- chains
def wick_sign(parities, pairs):
    """sign of bringing an ordered product to (pair)(pair)...(rest, original order); pairs: list of (i, j) i<j"""
    order = list(range(len(parities)))
    sgn = 1
    partner = {}
    for i, j in pairs:
        partner[i] = j; partner[j] = i
    work = list(order)
    out = []
    while True:
        k = next((x for x in work if x in partner), None)
        if k is None:
            break
        pk = partner[k]
        a, b = work.index(k), work.index(pk)
        assert a < b
        between = work[a + 1:b]
        if parities[pk] % 2 == 1:
            nodd = sum(parities[x] % 2 for x in between)
            if nodd % 2:
                sgn = -sgn
        work.remove(pk); work.remove(k)
    return sgn


def chain_numerator(ta, za, tb, zb, verts, rs, cl=None, fa=None, fb=None):
    """Numerator N(l) of one chain from letter a (type ta, selector za) to letter b (type tb, selector zb).
    verts: list of vertex specs, in chain order from a:
        ('AAA', y)                    cubic gauge, all A lines, output A with selector y (10)
        ('AAAA', y1, y2)              quartic, word ad_{y2} ad_{y1}
        (name, tin, tout, tbg, ybg)   Grassmann density vertex (name in DENS)
    rs: list of output momenta r_j (4-vectors) per vertex (quartic: sum of its two outputs).
    fa, fb: letter derivative factors as functions of the incoming letter momentum (4-vector) -> sympy expr.
    Returns (N(l), Rs) with Rs the edge shifts (R_0 = 0)."""
    k = len(verts)
    Rs = [[0, 0, 0, 0]]
    for j in range(k):
        Rs.append([Rs[-1][m] + rs[j][m] for m in range(4)])
    lvec = list(L)
    edge_p = [[lvec[m] + Rs[j][m] for m in range(4)] for j in range(k + 1)]   # momentum along edge j (a -> b)
    # expand each vertex into its role-assignment terms (Grassmann ones) -> product over vertices
    vterms = []
    types_in, types_out = [], []
    for j, v in enumerate(verts):
        p_in = edge_p[j]
        p_out = [-x for x in edge_p[j + 1]]
        if v[0] == 'AAA':
            vterms.append([(['in', 'out', 'bg'], AAA_matrix(v[1], p_in, rs[j], p_out), 1, ['A', 'A', 'A'])])
            types_in.append('A'); types_out.append('A')
        elif v[0] == 'AAAA':
            vterms.append([(['in', 'out', 'bg', 'bg2'], AAAA_matrix(v[1], v[2]), 1, ['A', 'A', 'A', 'A'])])
            types_in.append('A'); types_out.append('A')
        else:
            name, tin, tout, tbg, ybg = v
            ts = grass_vertex_terms(name, tin, tout, tbg, ybg, cl)
            lst = []
            for roles, Mx, cs in ts:
                slot_types = [None] * 3
                for r, t in zip(roles, DENS[name]):
                    slot_types[roles.index(r)] = t
                lst.append((roles, Mx, cs, DENS[name]))
            vterms.append(lst)
            types_in.append(tin); types_out.append(tout)
    # edge propagators
    total = 0
    for combo in itertools.product(*vterms):
        # build ordered factor list: X_a, X_b, then vertex factors in density order
        parities = [PARITY[ta], PARITY[tb]]
        pos = {}
        for j, (roles, Mx, cs, stypes) in enumerate(combo):
            for r, t in zip(roles, stypes):
                pos[(j, r)] = len(parities)
                parities.append(PARITY[t])
        pairs = [(0, pos[(0, 'in')])]
        for j in range(k - 1):
            pairs.append((pos[(j, 'out')], pos[(j + 1, 'in')]))
        pairs.append((1, pos[(k - 1, 'out')]))
        sgn = wick_sign(parities, pairs)
        # matrices
        P0 = prop_hat(ta, types_in[0], [-x for x in lvec], cl)                     # p_a = -l at the letter
        mats = [sp.Matrix([za]) * P0]
        csign = 1
        for j, (roles, Mx, cs, stypes) in enumerate(combo):
            mats.append(Mx)
            csign *= cs
            if j < k - 1:
                pj = [-x for x in edge_p[j + 1]]                                     # incoming at V_j.out
                mats.append(prop_hat(types_out[j], types_in[j + 1], pj, cl))
        Pb = prop_hat(tb, types_out[-1], edge_p[k], cl)                              # p_b = l + R at the letter
        mats.append(Pb.T * sp.Matrix(zb))
        prod = mats[0]
        for m in mats[1:]:
            prod = prod * m
        total += sgn * csign * prod[0, 0]
    fa_v = fa([-x for x in lvec]) if fa else 1
    fb_v = fb(edge_p[k]) if fb else 1
    return sp.expand(total * fa_v * fb_v), Rs


def kernel_chain(ta, za, tb, zb, verts, rs, cl=None, fa=None, fb=None):
    N, Rs = chain_numerator(ta, za, tb, zb, verts, rs, cl, fa, fb)
    return singular_parts(N, Rs)


def kernel_k0(ta, za, tb, zb, fa=None, fb=None, cl=None):
    """self-contraction (k = 0): <X_a X_b> with one edge"""
    lvec = list(L)
    P0 = prop_hat(ta, tb, [-x for x in lvec], cl)
    if P0 is None:
        return {}
    N = (sp.Matrix([za]) * P0 * sp.Matrix(zb))[0, 0]
    fa_v = fa([-x for x in lvec]) if fa else 1
    fb_v = fb(lvec) if fb else 1
    return singular_parts(sp.expand(N * fa_v * fb_v), [[0, 0, 0, 0]])


def deriv_factor(alpha):
    """letter derivative d^alpha (alpha: tuple of 4 nonneg ints) -> factor (i p)^alpha"""
    def f(p):
        out = 1
        for m in range(4):
            out *= (I * p[m]) ** alpha[m]
        return out
    return f


def unit(n, i):
    v = [0] * n
    v[i] = 1
    return v
