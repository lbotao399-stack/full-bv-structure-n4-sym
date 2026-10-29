"""Engine B, Step 2: tree-level algebra of functionals on plane waves: the contraction B_K(f,g) along a kernel K
(K_0 = local bracket, K_eps - K_0 = special edge, P(eps,infty) = propagator), Q = {S_2, -}, and the checks.

Conventions (identical to engine_b/vertices.py and the toy-verified loop formula):
 * a functional f is given by its vertex functions: vf(legs) = {colour key: value}, value = coefficient of the product
   of the odd r's in leg order, f evaluated on sum_j r_j u_j e^{i p_j.x}, (2pi)^4 delta stripped.
   Colour keys: 'k' = (c_0,c_1), 'f' = ([c_0,c_1],c_2), 'K12','K13','K14' = ([c_0,c_1],[c_2,c_3]), ([c_0,c_2],[c_1,c_3]),
   ([c_0,c_3],[c_1,c_2]) (legs in leg order).
 * contraction operator d_K = 1/2 sum_{ij} K_{ij} d_j d_i (left derivatives), graded-symmetric kernel, d_K(x_a x_b) = K_{ab}.
   B_K(f,g) := cross part of d_K(f g) (one derivative on f, one on g).
   For local brackets: {f,g}_K = -B_K(f,g)  (Delta_K = -d_K; Delta_eps = -d_{K_eps}, DESIGN.md).
 * kernel in momentum space: K^{XY}(q) with the slot X (on the f side) carrying the plane wave e^{iq.x}, the slot Y (g side)
   e^{-iq.x}; this is the convention of the propagator blocks of the loop formula (engine_b/loops.py, r2_T2.py).
 * Grassmann sign of one contraction (f = R_f z_f V_f, g = z_g R_g V_g; z = internal coordinates, R = odd r's in leg order):
   d_{z_g} d_{z_f} (f g) = (-1)^{|R_f| (|X| + |Y|)} R_f R_g V_f V_g, then R_f R_g is reordered to leg order."""
import itertools, os, sys
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
EB = os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b'))
if EB not in sys.path:
    sys.path.insert(0, EB)
import vertices as V
import loops as L

I = sp.I
PAR = dict(V.PAR); DIM = dict(V.DIM)

def perm_sign(seq):
    s = 1; seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s

# ------------------------------------------------------------------ colour
def colour_canon(struct, legs):
    """struct: ('k', (x, y)) or ('f', (x, y, z)) or ('bb', (x, y, z, w)) with global leg labels; returns (key, sign) in the
    basis of the module docstring (legs relabelled to their rank)."""
    kind, ls = struct
    rank = {l: i for i, l in enumerate(sorted(ls))}
    sig = [rank[l] for l in ls]
    if kind == 'k':
        return 'k', 1
    if kind == 'f':
        return 'f', perm_sign(sig)
    key, s = V.colour_key(('bb', (0, 1, 2, 3)), sig)
    return key, s

def contract_colour(kf, legs_f, kg, legs_g):
    """f has colour key kf on its legs (ext legs in global labels, internal last, label 'i'); g has kg on legs ('i' first,
    then ext).  Contraction of the internal index with (-,-)^{-1}.  Returns (struct over global labels, sign) or None."""
    # express f's key as a structure on its own legs
    def struct_of(key, ls):
        if key == 'k':
            return ('k', tuple(ls))
        if key == 'f':
            return ('f', tuple(ls))
        a, b, c, d = ls
        if key == 'K12':
            return ('bb', (a, b, c, d))
        if key == 'K13':
            return ('bb', (a, c, b, d))
        if key == 'K14':
            return ('bb', (a, d, b, c))
        raise ValueError(key)
    sf = struct_of(kf, legs_f); sg = struct_of(kg, legs_g)
    # internal label 'i' is last in f, first in g
    if sf[0] == 'k':                      # (x, i) . g  -> g with i -> x
        x = sf[1][0]
        gl = tuple(x if l == 'i' else l for l in sg[1])
        return (sg[0], gl), 1
    if sg[0] == 'k':                      # f . (i, y) -> f with i -> y
        y = sg[1][1]
        fl = tuple(y if l == 'i' else l for l in sf[1])
        return (sf[0], fl), 1
    if sf[0] == 'f' and sg[0] == 'f':     # ([x,y],t) (t',[z,w]) ... g = ([t,z],w) = (t,[z,w])
        x, y, _ = sf[1]; _, z, w = sg[1]
        return ('bb', (x, y, z, w)), 1
    raise NotImplementedError((kf, kg))

# ------------------------------------------------------------------ functionals
class Functional:
    """vf(legs) -> {colour key: value}; nlegs: allowed leg numbers (None = any)."""
    def __init__(self, vf, name='', parity=0, nlegs=None):
        self.vf = vf; self.name = name; self.parity = parity; self.nlegs = nlegs
    def __call__(self, legs):
        if self.nlegs is not None and len(legs) not in self.nlegs:
            return {}
        return self.vf(legs)
    def __add__(self, o):
        return Functional(lambda legs: addd(self(legs), o(legs)), f'({self.name}+{o.name})', self.parity)
    def scale(self, c, name=None):
        return Functional(lambda legs: {k: sp.expand(c * v) for k, v in self(legs).items()}, name or f'{c}*{self.name}', self.parity, self.nlegs)

def addd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = sp.expand(out.get(k, 0) + s * v)
    return {k: v for k, v in out.items() if v != 0}

def local(terms, name='', parity=None):
    """local functional: a list of (coefficient, term name in V.TERMS)."""
    nl = set(len(V.TERMS[t][1]) for c, t in terms)
    def vf(legs):
        out = {}
        for c, t in terms:
            if len(V.TERMS[t][1]) != len(legs):
                continue
            for k, v in V.vertex(t, legs).items():
                out[k] = sp.expand(out.get(k, 0) + c * v)
        return {k: v for k, v in out.items() if v != 0}
    if parity is None:
        slots = V.TERMS[terms[0][1]][1]
        parity = sum(PAR[s] for s in slots) % 2
    return Functional(vf, name, parity, nl)

def add_term(name, coef, slots, Lf, colour):
    V.TERMS[name] = (sp.sympify(coef), tuple(slots), Lf, colour)

# ------------------------------------------------------------------ kernels: list of (X, Y, matrix-function of q)
MPINV = V.Mpinv
def _fibre(X, Y):
    if (X, Y) in (('A', 'As'), ('As', 'A')):
        return lambda q: -sp.eye(10)
    if (X, Y) == ('lam', 'lams'):
        return lambda q: sp.Matrix(16, 16, lambda i, j: -MPINV[j, i])
    if (X, Y) == ('lams', 'lam'):
        return lambda q: sp.Matrix(16, 16, lambda j, i: -MPINV[j, i])
    if (X, Y) in (('c', 'cs'), ('cs', 'c')):
        return lambda q: -sp.eye(1)
    raise ValueError
FPAIRS = [('A', 'As'), ('As', 'A'), ('lam', 'lams'), ('lams', 'lam'), ('c', 'cs'), ('cs', 'c')]

def kernel_K0():
    return [(X, Y, _fibre(X, Y)) for X, Y in FPAIRS]

def kernel_scaled(base, scal):
    """base kernel times the scalar function scal(q2)."""
    return [(X, Y, (lambda f: (lambda q: f(q) * scal(sum(x * x for x in q))))(f)) for X, Y, f in base]

def kernel_P(eps):
    """P(eps,infty): N(q) e^{-eps q^2}/q^2 (eps may be 0 or a sympy symbol)."""
    out = []
    o4 = L.OFF
    def blk(X, Y):
        def f(q):
            N = L.prop_numerator(list(q))
            q2 = sum(x * x for x in q)
            M = sp.Matrix(DIM[X], DIM[Y], lambda i, j: N[o4[X] + i, o4[Y] + j] if (i < DIM[X] and j < DIM[Y]) else 0)
            return M * (sp.exp(-eps * q2) / q2)
        return f
    for X, Y in (('A', 'A'), ('lam', 'lam'), ('c', 'As'), ('As', 'c')):
        out.append((X, Y, blk(X, Y)))
    return out

# ------------------------------------------------------------------ the contraction B_K(f, g)
def sym_leg(sp_, mom, tag):
    return dict(sp=sp_, u=list(sp.symbols(f'{tag}0:{DIM[sp_]}')), p=list(mom))

def grad(vals, leg):
    """vals: {key: expr linear in leg['u']} -> {key: [coefficient of u_a]}"""
    out = {}
    for k, v in vals.items():
        v = sp.expand(v)
        out[k] = [v.coeff(s) for s in leg['u']]
    return out

_DEPTH = [0]
def B(K, f, g, legs, ext_sizes=None):
    """cross part of d_K(f g) evaluated on legs (global order 0..n-1).  Returns {colour key: value}.
    (internal-leg symbols carry the recursion depth, so nested contractions never share symbols)"""
    _DEPTH[0] += 1
    try:
        return _B(K, f, g, legs, _DEPTH[0])
    finally:
        _DEPTH[0] -= 1

def _B(K, f, g, legs, depth):
    n = len(legs)
    out = {}
    for r in range(1, n):
        for Lf in itertools.combinations(range(n), r):
            Lg = tuple(i for i in range(n) if i not in Lf)
            if f.nlegs is not None and len(Lf) + 1 not in f.nlegs:
                continue
            if g.nlegs is not None and len(Lg) + 1 not in g.nlegs:
                continue
            qf = [-sum(legs[i]['p'][m] for i in Lf) for m in range(4)]
            Rf = sum(PAR[legs[i]['sp']] for i in Lf) % 2
            oddseq = [i for i in Lf if PAR[legs[i]['sp']]] + [i for i in Lg if PAR[legs[i]['sp']]]
            ps = perm_sign(oddseq)
            for X, Y, Kf in K:
                lf = sym_leg(X, qf, f'zf{depth}_'); lg = sym_leg(Y, [-x for x in qf], f'zg{depth}_')
                vf = f([legs[i] for i in Lf] + [lf])
                if not vf:
                    continue
                vg = g([lg] + [legs[i] for i in Lg])
                if not vg:
                    continue
                Km = Kf(qf)
                gf = grad(vf, lf); gg = grad(vg, lg)
                sgn = (-1) ** (Rf * (PAR[X] + PAR[Y])) * ps
                for kf, af in gf.items():
                    for kg, ag in gg.items():
                        val = sum(af[a] * Km[a, b] * ag[b] for a in range(len(af)) for b in range(len(ag)) if af[a] != 0 and ag[b] != 0 and Km[a, b] != 0)
                        if val == 0:
                            continue
                        (st, cs) = contract_colour(kf, list(Lf) + ['i'], kg, ['i'] + list(Lg))
                        key, s2 = colour_canon(st, None)
                        out[key] = sp.expand(out.get(key, 0) + sgn * cs * s2 * val)
    return {k: v for k, v in out.items() if v != 0}

def Bfun(K, f, g, name=''):
    par = (f.parity + g.parity + (1 if K and PAR[K[0][0]] != PAR[K[0][1]] else 0)) % 2
    return Functional(lambda legs: B(K, f, g, legs), name or f'B({f.name},{g.name})', par)

def bracket(f, g, K=None, name=''):
    """{f,g}_K = -B_K(f,g)  (K = K_0 by default)."""
    K = kernel_K0() if K is None else K
    par = (f.parity + g.parity + 1) % 2
    return Functional(lambda legs: {k: -v for k, v in B(K, f, g, legs).items()}, name or '{%s,%s}' % (f.name, g.name), par)

def reduce_K13(d):
    """use the Jacobi identity K13 = K12 + K14 to write 4-leg colour in the basis K12, K14."""
    d = dict(d)
    if 'K13' in d:
        v = d.pop('K13')
        d['K12'] = sp.expand(d.get('K12', 0) + v); d['K14'] = sp.expand(d.get('K14', 0) + v)
    return {k: v for k, v in d.items() if v != 0}

# ------------------------------------------------------------------ the action and standard functionals
S2 = local([(1, 'kinA'), (1, 'kinlam'), (1, 'kingh')], 'S_2', 0)
I3 = local([(1, 'AAA'), (1, 'lamlamA'), (1, 'AcAs'), (1, 'clamlams'), (1, 'cccs')], 'I_3', 0)
I4 = local([(1, 'AAAA')], 'I_4', 0)
I3_TERMS = {'AAA': local([(1, 'AAA')], 'AAA', 0), 'lamlamA': local([(1, 'lamlamA')], 'lamlamA', 0),
            'AcAs': local([(1, 'AcAs')], 'AcAs', 0), 'clamlams': local([(1, 'clamlams')], 'clamlams', 0),
            'cccs': local([(1, 'cccs')], 'cccs', 0)}

def Q(F):
    """Q F = {S_2, F} (Prop. c1-Q: the linear part of {S,-} is the derivation induced by Q)."""
    return bracket(S2, F, name=f'Q({F.name})')

# Z pieces: int(A,A^*) (all ten components), int(lambda^*,lambda), int(c^*,c)
add_term('ZA', 1, ('A', 'As'), lambda s: sum(s[0]['u'][M] * s[1]['u'][M] for M in range(10)), ('k', (0, 1)))
add_term('Zlam', 1, ('lam', 'lams'), lambda s: (sp.Matrix(s[0]['u']).T * V.Mp * sp.Matrix(s[1]['u']))[0], ('k', (0, 1)))
add_term('Zc', 1, ('c', 'cs'), lambda s: s[0]['u'][0] * s[1]['u'][0], ('k', (0, 1)))
ZA = local([(1, 'ZA')], 'int(A,A*)', 1)
Zlam_raw = local([(1, 'Zlam')], 'lambda^T Mp lambda*', 1)
Zc_raw = local([(1, 'Zc')], 'c c*', 1)
# mass term int(A,A) (all ten) and the log-part functionals of I^CT_2 (Theorem c4-quadratic)
add_term('mass', 1, ('A', 'A'), lambda s: sum(s[0]['u'][M] * s[1]['u'][M] for M in range(10)), ('k', (0, 1)))
MASS = local([(1, 'mass')], 'int(A,A)', 0)
