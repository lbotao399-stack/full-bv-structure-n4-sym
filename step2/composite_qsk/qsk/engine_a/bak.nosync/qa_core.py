"""Engine A, part Q (route A): defect cycles through two letters along a generator x, exact.

A defect cycle (qsk_theory.tex Def def:qs-defect): chain X_a - line 1 - block 1 - ... - block n - line 0 - X_b with
blocks of I (cubic: one output, quartic: two A outputs), EXACTLY ONE block of S^(1)(x) (a quadratic two-leg block, or the
cubic term T5 with one output) and EXACTLY ONE special line K_eps - K_0 (numerator -|l+q|^2 x the field-antifield
components of K_0, Schwinger time in [0, eps]); all other lines P(eps, T=infinity).  Every placement is generated: the
S-block at every position, the special line at every line (on a line / at the special edge / cubic), k_I >= 0.

Chain algebra, momenta, Gaussian loop integral and the exact finite part (Thm thm:qs-simplex):
step2/composite_qsk/composite/engine_a/src.nosync/cq_gcore.py (chain, GSimplex, apply_letters; read-only import):
  16 pi^2 [eps^0] = int_{alpha_s smallest} F_0 log(min_{f!=s} alpha_f/alpha_s),
  16 pi^2 [eps^-m] = (1/m) int_{alpha_s smallest} [(min_{f!=s} alpha_f)^m - alpha_s^m] F_m ;  no log eps.
Signs.  Raw chain T = L_1 S h_1 ... L_n S h_n S L_0 S (two-letter engine).  On top of it every configuration carries
   (-1)^{n(<sigma) + |S1(x)| (n(<s) + [s<sigma])} (-1)^{[sigma = line 0]} * GLOBAL,
n(<j) the number of odd outputs before position j along the chain from X_a, |S1(x)| = 0 for odd x, 1 for even x;
(-1)^{[sigma = line 0]} removes the factor that S L_0 S produces on a special line at X_b (opposite-parity slots), so
that the factor is that of Def def:qs-defect (ii) relative to Lemma c10-sign; for even x on a line this is literally
the rule of the Part-C engine (R2 x vertex-vertex sign -1), checked in gate E0-partC.  The intrinsic pair sign
(-1)^{|X_a|+|X_b|} of ea_kern.kernel is applied (INTR).  GLOBAL = +1 (Part-C normalisation; fixed by G-reg).
Values: list of five Gaussian rationals on the basis (1, log2, log3, log5, log7), units hbar/(16 pi^2), colour word
ad_{Y_r} ... ad_{Y_1} of the outputs in chain order (the output of T5 included)."""
import itertools
from fractions import Fraction
import numpy as np
import qa_s1 as S1
import ea_core as E
import ea_kern as EK
import ea_ops as EO
import cq_core as Q
import cq_gcore as G
import cq_eval as CV

NB = 5
vzero, vadd, vscale, vstr, viszero = CV.vzero, CV.vadd, CV.vscale, CV.vstr, CV.viszero
RULE = {'global': 1, 'L0fix': True, 'intr': True}
STATS = {'configs': 0, 'nonzero': 0}


class Gen:
    """a generator's vertex: quad = Aff block (or None); cub = callable(ltype, pol, mom) -> Aff (or None) for T5;
    par = |S^(1)(x)| (0 odd x, 1 even x); name for caching"""

    def __init__(self, name, quad=None, cub=None, par=0, cubtypes=('AV', 'AW', 'lamsP', 'lamsM'), zero=False):
        self.name, self.quad, self.cub, self.par, self.cubtypes = name, quad, cub, par, cubtypes
        self.zero = zero          # arity zero: no S block (special-edge cycles of L_1^(0))


ZERO = None


def zero_gen():
    global ZERO
    if ZERO is None:
        ZERO = Gen('arity0', None, None, 0, zero=True)
    return ZERO


def odd_gen(xi, name='xi', which=('T1', 'T2', 'T3', 'T4', 'T5')):
    xi = np.asarray(xi)
    quad_terms = S1.odd_terms(xi, tuple(w for w in which if w != 'T5'))
    quad = S1.quad_block(quad_terms) if quad_terms else None
    cub = None
    if 'T5' in which:
        t5 = S1.odd_terms(xi, ('T5',))
        key = ('T5', tuple(np.asarray(xi).tolist()))
        cub = lambda lt, pol, mom: S1.cubic_block(key, t5, lt, pol, mom)
    return Gen(name, quad, cub, 0)


def s0_gen(eta, name='s0'):
    return Gen(name, S1.quad_block(S1.s0_term(eta)), None, 0)


def even_gen(name, **kw):
    return Gen(name, S1.quad_block(S1.even_terms(**kw)), None, 1)


def _legblock(leg):
    lt, pol, q = leg
    return EO.slotvec_to_stripped(lt, pol), np.asarray(q)


def structures(legs, gen):
    """partitions of the ordered legs into blocks: ('I1', i), ('I2', i, i+1), ('S3', i) (T5 output), plus the quadratic
    S block ('S2',) inserted at a gap; exactly one S block.  Yields lists of block descriptors in chain order."""
    n = len(legs)
    def rec(i, used):
        if i == n:
            yield [], used; return
        lt = legs[i][0]
        for tail, u in rec(i + 1, used):
            yield [('I1', i)] + tail, u
        if gen.cub is not None and not used and lt in gen.cubtypes:
            for tail, u in rec(i + 1, True):
                yield [('S3', i)] + tail, u
        if i + 1 < n and E.SLOT[lt] == 'A' and E.SLOT[legs[i + 1][0]] == 'A':
            for tail, u in rec(i + 2, used):
                yield [('I2', i, i + 1)] + tail, u
    for st, used in rec(0, False):
        if gen.zero:
            if not used:
                yield st
            continue
        if used:
            yield st
        elif gen.quad is not None:
            for g in range(len(st) + 1):
                yield st[:g] + [('S2',)] + st[g:]


def build(st, legs, gen, r0=None, r1=None):
    blks, P0s, P1s, scale, odd, spos = [], [], [], 1, [], None
    for bi, b in enumerate(st):
        if b[0] == 'I1':
            lt, pol, q = legs[b[1]]
            blks.append(Q.cubic_full(lt, EO.slotvec_to_stripped(lt, pol), np.asarray(q)))
            P0s.append(np.asarray(q, dtype=complex)); P1s.append(np.zeros(4, dtype=int)); scale *= 2
            odd.append(Q.PARITY_T[lt])
        elif b[0] == 'I2':
            (l1, p1, q1), (l2, p2, q2) = legs[b[1]], legs[b[2]]
            blks.append(Q.quartic_full(EO.slotvec_to_stripped(l2, p2), EO.slotvec_to_stripped(l1, p1)))
            P0s.append(np.asarray(q1, dtype=complex) + np.asarray(q2, dtype=complex)); P1s.append(np.zeros(4, dtype=int))
            scale *= 8; odd.append(0)
        elif b[0] == 'S3':
            lt, pol, q = legs[b[1]]
            q = np.asarray(q, dtype=complex)
            blks.append(gen.cub(lt, EO.slotvec_to_stripped(lt, pol), q))
            extra = np.zeros(4, dtype=complex) if r0 is None else np.asarray(r0, dtype=complex)
            P0s.append(q + extra); P1s.append(np.zeros(4, dtype=int) if r1 is None else np.asarray(r1, dtype=int))
            scale *= 2; odd.append(Q.PARITY_T[lt]); spos = bi + 1
        else:
            blks.append(gen.quad)
            P0s.append(np.zeros(4, dtype=complex) if r0 is None else np.asarray(r0, dtype=complex))
            P1s.append(np.zeros(4, dtype=int) if r1 is None else np.asarray(r1, dtype=int))
            scale *= 2; odd.append(0); spos = bi + 1
    return blks, P0s, P1s, scale, odd, spos


def config_sign(odd, spos, sig, n, par):
    """sigma: line index (1..n, or 0 = last line).  Blocks 1..n; line i (i >= 1) precedes block i."""
    before_sigma = n if sig == 0 else sig - 1          # number of blocks before line sigma
    n_sig = sum(odd[:before_sigma])
    if spos is None:
        par, spos, n_s = 0, 1, 0
    n_s = sum(odd[:spos - 1])
    s_before = 1 if (sig == 0 or spos < sig) else 0
    e = n_sig + par * (n_s + s_before)
    if RULE['L0fix'] and sig == 0:
        e += 1
    return (-1) ** e * RULE['global']


def configs(legs, gen, r0=None, r1=None, only_sigma=None):
    for st in structures(legs, gen):
        blks, P0s, P1s, scale, odd, spos = build(st, legs, gen, r0, r1)
        n = len(blks)
        for sig in range(n + 1):
            if only_sigma is not None and not only_sigma(sig, spos, n):
                continue
            lines = ['P'] * (n + 1)
            lines[sig] = 'S'
            yield st, blks, P0s, P1s, scale, odd, spos, sig, lines


def _from_rl(rl, E_, za, zb):
    r0, r1 = rl
    out = vzero()
    if r0 is not None:
        out[0] = EK.contract(r0, za, zb)
    if r1 is not None:
        c = EK.contract(r1, za, zb)
        for bi, nn in enumerate(Q.log_basis(E_)):
            if nn:
                out[1 + bi] = (out[1 + bi][0] + nn * c[0], out[1 + bi][1] + nn * c[1])
    return out


def pair_value(Xa, Xb, legs, gen, what=('fin',), r0=None, r1=None, tord=0, split=False, only_sigma=None):
    """sum over all defect configurations of the pair (Xa, Xb) with the ordered outputs legs.
    Xa = (ltype, zeta, deriv-poly or None); legs = [(ltype, pol (letter rows), momentum (4-vector))].
    what: 'fin' (eps^0) and/or ('pole', m).  Returns {key: value} (value = 5 Gaussian rationals);
    split=True: {key: {(placement): value}} keyed by ('k', kI, 'S2'/'S3', sigma-adjacency)."""
    ta, za, fa = Xa; tb, zb, fb = Xb
    tot = {w: vzero() for w in what}
    parts = {}
    intr = -1 if (RULE['intr'] and (Q.PARITY_T[ta] + Q.PARITY_T[tb]) % 2 == 1) else 1
    for st, blks, P0s, P1s, scale, odd, spos, sig, lines in configs(legs, gen, r0, r1, only_sigma):
        STATS['configs'] += 1
        K, qs, P0 = G.chain(blks, P0s, P1s, lines, ta, tb)
        if not K:
            continue
        if fa is not None or fb is not None:
            K = G.apply_letters(K, fa or {(0, 0, 0, 0): 1}, fb or {(0, 0, 0, 0): 1}, P0)
            if not K:
                continue
        STATS['nonzero'] += 1
        n = len(blks)
        sgn = config_sign(odd, spos, sig, n, gen.par) * intr
        Sx = G.GSimplex(qs)
        E_ = len(qs[0])
        for w in what:
            if w == 'fin':
                r = Sx.evaluate(K, 'fin', s=sig)
            else:
                r = Sx.evaluate(K, 'spec', s=sig, m=w[1])
            rl = r[tord]
            rl = [None if v is None else (v[0] * Fraction(sgn, scale), v[1] * Fraction(sgn, scale)) for v in rl]
            v = _from_rl(rl, E_, za, zb)
            tot[w] = vadd(tot[w], v)
            if split:
                kI = sum(1 for b in st if b[0] in ('I1', 'I2'))
                adj = 'line'
                if spos is None:
                    adj = 'none'
                elif sig == spos or sig == (spos + 1 if spos < n else 0):
                    adj = 'special'
                key = (w, kI, 'S3' if any(b[0] == 'S3' for b in st) else 'S2', adj)
                parts[key] = vadd(parts.get(key, vzero()), v)
    if split:
        return tot, parts
    return tot


def plain_value(Xa, Xb, legs, gen, m, r0=None):
    """counterterm cycles along x (Def def:qs-ct1): one S block on the cycle, no special line: eps^-m and log eps"""
    ta, za, fa = Xa; tb, zb, fb = Xb
    out = {'log': vzero(), m: vzero()}
    intr = -1 if (RULE['intr'] and (Q.PARITY_T[ta] + Q.PARITY_T[tb]) % 2 == 1) else 1
    for st in structures(legs, gen):
        blks, P0s, P1s, scale, odd, spos = build(st, legs, gen, r0, None)
        n = len(blks)
        lines = ['P'] * (n + 1)
        K, qs, P0 = G.chain(blks, P0s, P1s, lines, ta, tb)
        if not K:
            continue
        if fa is not None or fb is not None:
            K = G.apply_letters(K, fa or {(0, 0, 0, 0): 1}, fb or {(0, 0, 0, 0): 1}, P0)
            if not K:
                continue
        # Koszul: the S vertex relative to the odd outputs before it (even x), as R2 for the internal vertex
        e = gen.par * sum(odd[:spos - 1])
        sgn = (-1) ** e * intr
        Sx = G.GSimplex(qs)
        for kd, key in (('log', 'log'), ('eps', m)):
            r = Sx.evaluate(K, kd, m=(m if kd == 'eps' else 0))
            rl = [None if v is None else (v[0] * Fraction(sgn, scale), v[1] * Fraction(sgn, scale)) for v in r[0]]
            out[key] = vadd(out[key], _from_rl(rl, len(qs[0]), za, zb))
    return out


def twoletter_value(Xa, Xb, legs, kind):
    """arity-zero pair kernel (two-letter engine), kind 'log' or m (eps^-m); stripped, units hbar/16pi^2"""
    z = EO.pair_kernel_on_legs(Xa, Xb, legs, kind=kind)
    out = vzero(); out[0] = (Fraction(z[0]), Fraction(z[1]))
    return out


# ------------------------------------------------------------------------------------------------ helpers
def gq(z):
    z = complex(z)
    r, i = round(z.real), round(z.imag)
    assert abs(z.real - r) < 1e-9 and abs(z.imag - i) < 1e-9, z
    return (Fraction(r), Fraction(i))


def deriv_r(f, mu, dmax, order=1):
    """exact derivative d/dr_mu (order 1) or d^2/dr_mu^2 at r = 0 of a polynomial value-function f(r) of degree <= dmax
    (f: callable taking an integer 4-vector), by exact Lagrange interpolation on the points r = j e_mu, |j| <= dmax"""
    h = (dmax + 2) // 2
    pts = list(range(-h, h + 1))
    vals = {}
    for j in pts:
        r = np.zeros(4, dtype=int); r[mu] = j
        vals[j] = f(r)
    # Lagrange basis derivative at 0
    out = vzero()
    for j in pts:
        num = [Fraction(1)]       # polynomial coefficients of prod_{i != j} (x - i)/(j - i)
        den = Fraction(1)
        for i in pts:
            if i == j:
                continue
            num = [Fraction(0)] + num
            for k in range(len(num) - 1):
                num[k] = num[k] - i * num[k + 1]
            den *= (j - i)
        c = num[order] * (1 if order == 1 else 2) / den
        if c:
            out = vadd(out, vscale(vals[j], (c, Fraction(0))))
    return out


def deriv_r2(f, mu, nu, dmax):
    """mixed second derivative d^2/dr_mu dr_nu at 0 (mu != nu) by polarisation of second derivatives"""
    if mu == nu:
        return deriv_r(f, mu, dmax, order=2)
    def g(rv, mu=mu, nu=nu):
        return None
    # d^2/dmu dnu = 1/2 (d^2_{e_mu+e_nu} - d^2_mu - d^2_nu) along the diagonal direction
    def fdiag(r):
        j = r[0]
        rr = np.zeros(4, dtype=int); rr[mu] += j; rr[nu] += j
        return f(rr)
    dd = deriv_r(lambda r: fdiag(r), 0, dmax, order=2)
    d1 = deriv_r(f, mu, dmax, order=2)
    d2 = deriv_r(f, nu, dmax, order=2)
    half = (Fraction(1, 2), Fraction(0))
    return vscale(vadd(dd, vscale(vadd(d1, d2), (Fraction(-1), Fraction(0)))), half)
