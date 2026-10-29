"""Engine A, part C (route A): special-edge composite cycles through two letters, exact.

Objects (composite_theory.tex, Def def:cc-sp, Def def:cc-cycle, Thm thm:cc-simplex).  A chain from letter a to letter b
with k vertices of I, one composite line (two segments P with the internal vertex of S^(1)(x) between them, Lemma
lem:cc-kernels (iv)) and one special edge K_eps - K_0 (field-antifield components of K_0, numerator -|l+q|^2, Schwinger
time in [0, eps]); every other edge P(eps, infty).  The cycle has E = k + 2 Schwinger times.

Conventions are those of the two-letter engine A (step2/twoletter/engine_a/src.nosync/ea_core.py, imported read-only):
colour-stripped blocks (cubic: (ad_Y t_in, t_out), quartic words), plane-wave legs with all momenta incoming, Hessian
H_jl = d_j d_l (left derivatives, l first), chain T = L_1 S h_1 L_2 S h_2 ... L_n S h_n S L_0 S with the rows in the slot
of letter a and the columns in the slot of letter b, line i entering block i, line 0 at letter b (incoming -k), letter a
incoming k - P_0, P_0 = -(sum of the block momenta).  Here the space is the full 54-component stripped space (the special
edge reaches the antifield slots A*|_perp, lambda*, c*).

Internal vertex blocks (no colour factor, (t_in, t_out)):
  'P'  (direction b): int (d_b A|_V, A*|_V)           H[A_M, A*_M] = i (kj)_b, H[A*_M, A_M] = i (kl)_b  (M < 4)
  's0'               : int 2 (A|_V, A*|_V)             H = 2 on both
  's'  (direction mu): int 2 (d_mu A|_V, A*|_V) at the injected momentum r (the position factor x^mu = -i d/dr_mu)
Koszul reading: the internal vertex (odd, degree -1) and the special edge (odd kernel) are read as odd symbols at their
place in the chain; the value is brought to the order (vertex, special edge, odd outputs in chain order) -- see
odd_sign().  The overall normalisation (one global constant) is fixed by gate GC3 and recorded in gates.json.

Exact arithmetic: matrices are Gaussian integers in complex128 with the 2^52 bound asserted on every product (r2_core.mm);
simplex weights are Fractions; the finite part of a cycle with E Schwinger times lies in Q + Q log(E/(E-1))."""
import os, sys, itertools, math, hashlib
from fractions import Fraction
import numpy as np

SRC = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(SRC)                                            # composite/engine_a
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))     # project root
TLSRC = os.path.join(ROOT, 'step2', 'twoletter', 'engine_a', 'src.nosync')
sys.path.insert(0, TLSRC)
import ea_core as E                                                    # two-letter engine A core (read-only)
import ea_kern as EK

NS = E.NS
OFS = E.OFS
SD = E.SDIAG.copy()
SLOT = E.SLOT
WEYL = E.WEYL
DEG = E.DEG
PARITY_T = E.PARITY_T
BASIS = E.BASIS
I4 = np.eye(4, dtype=int)
Z4 = np.zeros(4, dtype=int)


def provenance():
    out = {}
    for nm in ('ea_core.py', 'ea_kern.py', 'vendor/r2_core.py', 'vendor/cliff.py'):
        p = os.path.join(TLSRC, nm)
        out[nm] = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    return out

# ------------------------------------------------------------------------------------------------ rows of a letter
ROWS = {
    'AV': [OFS['A'] + M for M in range(4)], 'AW': [OFS['A'] + M for M in range(4, 10)],
    'AsV': [OFS['As'] + M for M in range(4)], 'AsW': [OFS['As'] + M for M in range(4, 10)],
    'lamP': [OFS['l'] + a for a in range(16)], 'lamM': [OFS['l'] + a for a in range(16)],
    'lamsP': [OFS['ls'] + a for a in range(16)], 'lamsM': [OFS['ls'] + a for a in range(16)],
    'c': [OFS['c']], 'cs': [OFS['cs']],
}

# ------------------------------------------------------------------------------------------------ lines (full space)
def _lines():
    P0 = np.zeros((NS, NS), complex)
    for M in range(10):
        P0[OFS['A'] + M, OFS['A'] + M] = 2
    P1 = []
    for mu in range(4):
        X = np.zeros((NS, NS), complex)
        X[OFS['c'], OFS['As'] + mu] = -1j
        X[OFS['As'] + mu, OFS['c']] = -1j
        X[OFS['l']:OFS['l'] + 16, OFS['l']:OFS['l'] + 16] = 1j * np.asarray(E.SG[mu])
        P1.append(X)
    Kn = np.zeros((NS, NS), complex)          # -Khat: the special edge numerator is -|l+q|^2 Khat = +|l+q|^2 Kn
    for M in range(10):
        Kn[OFS['A'] + M, OFS['As'] + M] = 1; Kn[OFS['As'] + M, OFS['A'] + M] = 1
    for a in range(16):
        Kn[OFS['l'] + a, OFS['ls'] + a] = 1; Kn[OFS['ls'] + a, OFS['l'] + a] = 1
    Kn[OFS['c'], OFS['cs']] = 1; Kn[OFS['cs'], OFS['c']] = 1
    return P0, P1, Kn


LP0, LP1, KN = _lines()
# consistency with the restricted line matrices of ea_core
_A = E.ACT
assert np.array_equal(LP0[np.ix_(_A, _A)], E.LINES['AV'][0] + E.LINES['AP'][0])
assert all(np.array_equal(LP1[m][np.ix_(_A, _A)], E.LINES['gh'][1][m] + E.LINES['l'][1][m]) for m in range(4))
# consistency with r2_core's Khat (colour index 0): Kn = -Khat on the stripped components
_KH = sum(E.C.Khat(L) for L in E.C.SLINES)
_idx = np.array([E.full_index(t, a, 0) for t in E.TYPES for a in range(E.SZ[t])])
assert np.array_equal(-_KH[np.ix_(_idx, _idx)], KN)


RULE = {'koszul': 'R2', 'global': 1, 'spec': 'K'}     # sign conventions (gates GC3, GS)


def line_poly(kind, q):
    q = np.asarray(q)
    if kind == 'P':
        const = LP0 + sum(q[mu] * LP1[mu] for mu in range(4))
        return E.affine_poly(const, LP1)
    if kind == 'S':
        KM = KN if RULE['spec'] == 'K' else (SD[:, None] * KN if RULE['spec'] == 'SK' else KN * SD[None, :])
        out = {}
        qq = int(q @ q)
        if qq:
            out[(0, 0, 0, 0)] = qq * KM
        for mu in range(4):
            if q[mu]:
                out[tuple(I4[mu])] = (2 * int(q[mu])) * KM
            e = [0, 0, 0, 0]; e[mu] = 2
            out[tuple(e)] = KM.copy()
        return out
    raise ValueError(kind)

# ------------------------------------------------------------------------------------------------ blocks (full space)
class Aff:
    __slots__ = ['h0', 'Hj', 'Hl', 'kind']

    def __init__(self, h0, Hj, Hl, kind):
        self.h0, self.Hj, self.Hl, self.kind = h0, Hj, Hl, kind


def cubic_full(ltype, pol, mom):
    B = E.blocks()['cub']
    slot = SLOT[ltype]
    h0 = np.zeros((NS, NS), complex); Hj = [np.zeros((NS, NS), complex) for _ in range(4)]
    Hl = [np.zeros((NS, NS), complex) for _ in range(4)]
    for comp in range(E.SZ[slot]):
        c = pol[comp]
        if c == 0 or (slot, comp) not in B:
            continue
        b0, Hb, bj, bl = B[(slot, comp)]
        h0 = h0 + c * (b0 + sum(mom[m] * Hb[m] for m in range(4)))
        for m in range(4):
            Hj[m] = Hj[m] + c * bj[m]; Hl[m] = Hl[m] + c * bl[m]
    return Aff(h0, Hj, Hl, 'cub')


def quartic_full(polY, polY2):
    """word ad_Y ad_Y' (Y' nearer to letter a)"""
    Q = E.blocks()['qua']
    h0 = np.zeros((NS, NS), complex)
    for M in range(10):
        if polY[M] == 0:
            continue
        for N in range(10):
            if polY2[N] == 0:
                continue
            h0 = h0 + polY[M] * polY2[N] * Q[(M, N)][0]
    z = [np.zeros((NS, NS), complex) for _ in range(4)]
    return Aff(h0, z, [x.copy() for x in z], 'qua')


def v2_block(kind, d=None):
    """internal vertex of S^(1)(x) between the two segments of the composite line"""
    h0 = np.zeros((NS, NS), complex); Hj = [np.zeros((NS, NS), complex) for _ in range(4)]
    Hl = [np.zeros((NS, NS), complex) for _ in range(4)]
    if kind == 's0':
        for M in range(4):
            h0[OFS['A'] + M, OFS['As'] + M] = 2; h0[OFS['As'] + M, OFS['A'] + M] = 2
    elif kind in ('P', 's'):
        sc = 1 if kind == 'P' else 2
        for M in range(4):
            Hj[d][OFS['A'] + M, OFS['As'] + M] = 1j * sc      # j-leg A: (i kj)_d
            Hl[d][OFS['As'] + M, OFS['A'] + M] = 1j * sc      # l-leg A: (i kl)_d
    elif kind == 'Pv':                                        # int (beta.d A|_V, A*|_V), beta = d (a 4-vector)
        for mu in range(4):
            if d[mu] == 0:
                continue
            for M in range(4):
                Hj[mu][OFS['A'] + M, OFS['As'] + M] = 1j * d[mu]
                Hl[mu][OFS['As'] + M, OFS['A'] + M] = 1j * d[mu]
    else:
        raise ValueError(kind)
    return Aff(h0, Hj, Hl, 'v2')


def block_poly(blk, qin, qout):
    """S h(kj = -(k+qin), kl = k+qout) as a polynomial matrix in k"""
    const = blk.h0 - sum(qin[mu] * blk.Hj[mu] for mu in range(4)) + sum(qout[mu] * blk.Hl[mu] for mu in range(4))
    lin = [blk.Hl[mu] - blk.Hj[mu] for mu in range(4)]
    S = SD[:, None]
    return E.affine_poly(S * const, [S * x for x in lin])


def chain(blks, Ps, lines, ta, tb):
    """polynomial matrix K(k)[rows of ta, cols of tb] of T = L_1 S h_1 ... L_n S h_n S L_0 S; lines[i] in {'P','S'}"""
    n = len(blks)
    P0 = -sum((np.asarray(p, dtype=int) for p in Ps), np.zeros(4, dtype=int))
    allP = [P0] + [np.asarray(p, dtype=int) for p in Ps]
    qs = [np.zeros(4, dtype=int)]
    for i in range(1, n + 1):
        qs.append(qs[-1] - allP[i - 1])
    ra = ROWS[ta]; rb = ROWS[tb]
    T = {(0, 0, 0, 0): np.eye(NS, dtype=complex)[ra]}
    for i in range(1, n + 1):
        T = E.pmul(T, line_poly(lines[i], qs[i]))
        qn = qs[i + 1] if i < n else qs[0]
        T = E.pmul(T, block_poly(blks[i - 1], qs[i], qn))
        if not T:
            return {}, qs, P0
    T = {e: M * SD[None, :] for e, M in T.items()}
    T = E.pmul(T, line_poly(lines[0], qs[0]))
    T = {e: (M * SD[None, :])[:, rb] for e, M in T.items()}
    T = {e: M for e, M in T.items() if np.any(M)}
    return T, qs, P0

# ------------------------------------------------------------------------------------------------ finite part weights
_FIN = {}


def _H(n):
    return sum((Fraction(1, j) for j in range(1, n + 1)), Fraction(0))


def _beta_log(a, b, E_):
    """int_0^1 x^a (1-x)^b log((x + c(1-x))/x) dx / B(a+1,b+1), c = E/(E-1): returned as (rational, coeff of log c)"""
    c = Fraction(E_, E_ - 1); d = 1 - c
    # R(m) = int_0^1 x^m/(c + d x) dx = (1/d)(1/m - c R(m-1)),  R(0) = -(log c)/d
    R = [(Fraction(0), -1 / d)]
    for m in range(1, a + b + 3):
        r0, r1 = R[-1]
        R.append(((Fraction(1, m) - c * r0) / d, (-c * r1) / d))
    def I(nn):              # int_0^1 x^nn log(c + d x) dx = -(d/(nn+1)) R(nn+1)
        r0, r1 = R[nn + 1]
        return (-d / (nn + 1) * r0, -d / (nn + 1) * r1)
    J0, J1 = Fraction(0), Fraction(0)
    for i in range(b + 1):
        cf = math.comb(b, i) * (-1) ** i
        x0, x1 = I(a + i)
        J0 += cf * x0; J1 += cf * x1
    Bv = Fraction(math.factorial(a) * math.factorial(b), math.factorial(a + b + 1))
    elogx = _H(a) - _H(a + b + 1)
    return (J0 / Bv - elogx, J1 / Bv)


def fin_monomial(exps, s):
    """int_{R_s} alpha^exps log(min_{f != s} alpha_f / alpha_s) d alpha = (rational, coeff of log(E/(E-1)))"""
    key = (exps, s)
    if key in _FIN:
        return _FIN[key]
    E_ = len(exps)
    others = [g for g in range(E_) if g != s]
    t0, t1 = Fraction(0), Fraction(0)
    for order in itertools.permutations(others):
        seq = [s] + list(order)
        verts = []
        for j in range(E_):
            v = [Fraction(0)] * E_
            for g in seq[j:]:
                v[g] = Fraction(1, E_ - j)
            verts.append(v)
        poly = {tuple([0] * E_): Fraction(1)}
        for g in range(E_):
            lin = {}
            for j in range(E_):
                if verts[j][g] != 0:
                    e = [0] * E_; e[j] = 1
                    lin[tuple(e)] = verts[j][g]
            poly = E.sp_mul(poly, E.sp_pow(lin, exps[g], E_))
        Mx = [[verts[j][g] for g in range(E_ - 1)] + [Fraction(1)] for j in range(E_)]
        det = abs(E._det(Mx))
        for b, cf in poly.items():
            dv = E.dirichlet(b)
            r0, r1 = _beta_log(b[0], b[1], E_)
            t0 += det * cf * dv * r0; t1 += det * cf * dv * r1
    _FIN[key] = (t0, t1)
    return (t0, t1)


def log_basis(E_):
    """log(E/(E-1)) in the basis (log2, log3, log5, log7)"""
    def fac(n):
        out = {2: 0, 3: 0, 5: 0, 7: 0}
        for p in (2, 3, 5, 7):
            while n % p == 0:
                out[p] += 1; n //= p
        assert n == 1
        return out
    a, b = fac(E_), fac(E_ - 1)
    return [a[p] - b[p] for p in (2, 3, 5, 7)]


BASIS_NAMES = ['1', 'log2', 'log3', 'log5', 'log7']


class CSimplex(E.Simplex):
    """ea_core.Simplex plus the finite part of a cycle with one special edge (Thm thm:cc-simplex)"""

    def fin(self, N, s):
        degk = max((sum(e) for e in N), default=0)
        jmax = degk // 2
        laps = [N]
        for j in range(1, jmax + 1):
            laps.append(E.lap_poly(laps[-1]))
        E_ = self.E
        acc0 = None; acc1 = None
        for j, Lj in enumerate(laps):
            mj = j + 2 - E_
            if mj < 0:
                continue
            pref = Fraction(1, math.factorial(mj) * math.factorial(j) * 4 ** j)
            for kexp, coef in Lj.items():
                key = ('fin', kexp, j, s)
                if key not in self._w:
                    poly = self.phipow(mj)
                    for mu in range(4):
                        if kexp[mu]:
                            poly = E.sp_mul(poly, self.cpow(mu, kexp[mu]))
                    w0, w1 = Fraction(0), Fraction(0)
                    for e_, v in poly.items():
                        a0, a1 = fin_monomial(e_, s)
                        w0 += v * a0; w1 += v * a1
                    self._w[key] = (w0 * pref, w1 * pref)
                w0, w1 = self._w[key]
                if w0:
                    acc0 = E._acc(acc0, coef, w0)
                if w1:
                    acc1 = E._acc(acc1, coef, w1)
        return E._finish(acc0), E._finish(acc1)


# ------------------------------------------------------------------------------------------------ configurations
def odd_sign(symbols, rule='front'):
    """symbols: chain-ordered list of 'v' (internal vertex), 's' (special edge), 't' (odd output).
    rule 'front': sign of the permutation to the order (v, s, t...);  'none': +1"""
    if rule == 'none':
        return 1
    if rule == 'R2':                 # internal vertex and special edge odd relative to the odd outputs only
        sgn = 1
        for sym in ('v', 's'):
            if sym in symbols:
                i = symbols.index(sym)
                if sum(1 for x in symbols[:i] if x == 't') % 2:
                    sgn = -sgn
        return sgn
    seq = list(symbols)
    sgn = 1
    for target in ('v', 's'):
        if target not in seq:
            continue
        i = seq.index(target)
        n_before = sum(1 for x in seq[:i])
        if n_before % 2:
            sgn = -sgn
        seq.pop(i)
    return sgn




def make_config(outs, v2pos, spos, v2kind, v2dir=None, r=None):
    """outs: blocks of I as in ea_kern.make_blocks (('cub', ltype, pol, mom) or ('qua', (p1,m1), (p2,m2)), ordered from
    letter a); v2pos in 0..len(outs): the internal vertex is inserted before outs[v2pos]; spos: index of the special
    line in the final chain numbering (lines 1..n, 0), must not be adjacent to the internal vertex.
    Returns (blks, Ps, lines, scale, symbols) or None."""
    blks, Ps, scale, syms = [], [], 1, []
    seq = list(outs[:v2pos]) + [('v2',)] + list(outs[v2pos:])
    n = len(seq)
    lines = ['P'] * (n + 1)        # lines[i]: line i (1..n), lines[0] the last one
    # line numbering in chain order: line 1, block 1, line 2, ..., block n, line 0
    order_lines = list(range(1, n + 1)) + [0]
    vb = v2pos + 1                 # block index of the internal vertex (1-based)
    adj = {vb, vb + 1 if vb < n else 0}     # lines adjacent to it: line vb (before) and the line after
    if spos in adj:
        return None
    lines[spos] = 'S'
    for b in seq:
        if b[0] == 'cub':
            _, lt, pol, mom = b
            blks.append(cubic_full(lt, pol, mom)); Ps.append(np.asarray(mom, dtype=int)); scale *= 2
        elif b[0] == 'qua':
            _, (p1, m1), (p2, m2) = b
            blks.append(quartic_full(p2, p1)); Ps.append(np.asarray(m1, dtype=int) + np.asarray(m2, dtype=int)); scale *= 8
        else:
            blks.append(v2_block(v2kind, v2dir)); Ps.append(np.asarray(r if r is not None else Z4, dtype=int))
    # odd symbols in chain order
    for pos, li in enumerate(order_lines):
        if lines[li] == 'S':
            syms.append('s')
        if pos < n:
            b = seq[pos]
            if b[0] == 'v2':
                syms.append('v')
            elif b[0] == 'cub' and PARITY_T[b[1]] == 1:
                syms.append('t')
    return blks, Ps, lines, scale, syms


def config_value(outs, v2pos, spos, v2kind, ta, tb, fa=None, fb=None, v2dir=None, r=None, want=('fin', 'poles')):
    """value of one configuration: {'fin': ((re,im) rational part, (re,im) log(E/(E-1)) part, E), m: (re,im) eps^-m}
    as letter-component matrices (rows: slot of ta, cols: slot of tb), units hbar/(16 pi^2), normalised."""
    cfg = make_config(outs, v2pos, spos, v2kind, v2dir, r)
    if cfg is None:
        return None
    blks, Ps, lines, scale, syms = cfg
    K, qs, P0 = chain(blks, Ps, lines, ta, tb)
    if not K:
        return None
    if fa is not None or fb is not None:
        K = EK.apply_letters(K, fa or {(0, 0, 0, 0): 1}, fb or {(0, 0, 0, 0): 1}, P0)
        K = {e: M for e, M in K.items() if np.any(M)}
        if not K:
            return None
    sgn = odd_sign(syms, RULE['koszul']) * RULE['global']
    S = CSimplex(qs)
    s_idx = lines.index('S')            # Schwinger index = line index (qs[i] is line i)
    n = len(blks)
    vv = s_idx not in (1, 0)            # special edge between two vertices (not at a letter)
    if vv:                              # direction sign of the vertex-vertex special edges (convention VV_a)
        if syms.index('s') < syms.index('v'):
            sgn = -sgn
    out = {'vv': vv}
    if 'fin' in want:
        f0, f1 = S.fin(K, s_idx)
        out['fin'] = (_scale(f0, sgn, scale), _scale(f1, sgn, scale), len(qs))
    if 'poles' in want:
        for m, v in S.special(K, s_idx).items():
            if v is not None:
                out[m] = _scale(v, sgn, scale)
    return out


def _scale(v, sgn, scale):
    if v is None:
        return None
    re, im = v
    return (re * Fraction(sgn, scale), im * Fraction(sgn, scale))


def plain_config_value(outs, v2pos, v2kind, ta, tb, fa=None, fb=None, v2dir=None, r=None):
    """plain composite cycle (no special edge): log eps and eps^-m coefficients (for ct^(1)_x, lc^(1)_x)"""
    blks, Ps, scale, syms = [], [], 1, []
    seq = list(outs[:v2pos]) + [('v2',)] + list(outs[v2pos:])
    for b in seq:
        if b[0] == 'cub':
            _, lt, pol, mom = b
            blks.append(cubic_full(lt, pol, mom)); Ps.append(np.asarray(mom, dtype=int)); scale *= 2
            if PARITY_T[lt] == 1:
                syms.append('t')
        elif b[0] == 'qua':
            _, (p1, m1), (p2, m2) = b
            blks.append(quartic_full(p2, p1)); Ps.append(np.asarray(m1, dtype=int) + np.asarray(m2, dtype=int)); scale *= 8
        else:
            blks.append(v2_block(v2kind, v2dir)); Ps.append(np.asarray(r if r is not None else Z4, dtype=int))
            syms.append('v')
    lines = ['P'] * (len(seq) + 1)
    K, qs, P0 = chain(blks, Ps, lines, ta, tb)
    if not K:
        return None
    if fa is not None or fb is not None:
        K = EK.apply_letters(K, fa or {(0, 0, 0, 0): 1}, fb or {(0, 0, 0, 0): 1}, P0)
        K = {e: M for e, M in K.items() if np.any(M)}
        if not K:
            return None
    sgn = odd_sign(syms, RULE['koszul']) * RULE['global']
    S = E.Simplex(qs)
    res = S.singular(K, kinds=('log', 'eps'))
    return {k: _scale(v, sgn, scale) for k, v in res.items() if v is not None}
