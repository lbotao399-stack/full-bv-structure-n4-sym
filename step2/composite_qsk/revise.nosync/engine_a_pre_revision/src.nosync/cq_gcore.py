"""Engine A, part C: general kinematics core.  Same chains as cq_core, with
  * Gaussian-integer momenta (complex null momenta for on-shell evaluations), the Gaussian loop integral continued
    analytically (the exponent sum_e t_e (l+q_e).(l+q_e) uses the bilinear dot product);
  * a first-order parameter t: block momenta P = P0 + t P1, so that the injected momentum r = t e_mu of the internal
    vertex s_mu is differentiated exactly (coefficient of t^1), giving the translation-invariant part along Delta
    without finite differences (Prop prop:cc-position: x^mu = -i d/dr_mu).
Polynomials in (k, t) are dicts {(k0,k1,k2,k3,t): Gaussian-integer complex matrix}, truncated at t^2.
Scalar polynomials in (alpha_0..alpha_{E-1}, t) have Gaussian-rational coefficients (r2_core.Q2)."""
import math, itertools
from fractions import Fraction
import numpy as np
import cq_core as Q
import ea_core as E

Q2 = E.C.Q2
NS = Q.NS
SD = Q.SD


def _key(e4, t):
    return (e4[0], e4[1], e4[2], e4[3], t)


def pmul(A, B):
    out = {}
    for e1, M1 in A.items():
        for e2, M2 in B.items():
            t = e1[4] + e2[4]
            if t > 1:
                continue
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2], e1[3] + e2[3], t)
            P = E.mm(M1, M2)
            out[e] = out[e] + P if e in out else P
    return {e: M for e, M in out.items() if np.any(M)}


def affine(const0, const1, lin):
    """const0 + t const1 + sum_mu k_mu lin[mu]"""
    P = {}
    if np.any(const0):
        P[(0, 0, 0, 0, 0)] = const0
    if const1 is not None and np.any(const1):
        P[(0, 0, 0, 0, 1)] = const1
    for mu in range(4):
        if np.any(lin[mu]):
            e = [0, 0, 0, 0, 0]; e[mu] = 1
            P[tuple(e)] = lin[mu]
    return P


def line_poly(kind, q0, q1):
    """q = q0 + t q1 (q0 Gaussian integer vector, q1 integer vector)"""
    if kind == 'P':
        c0 = Q.LP0 + sum(q0[mu] * Q.LP1[mu] for mu in range(4))
        c1 = sum(q1[mu] * Q.LP1[mu] for mu in range(4)) if np.any(q1) else None
        return affine(c0, c1, Q.LP1)
    if kind == 'S':
        KM = Q.KN
        out = {}
        # |k + q|^2 = k.k + 2 q.k + q.q, q.q = q0.q0 + 2t q0.q1 (+ t^2 dropped)
        qq0 = complex(np.dot(q0, q0)); qq1 = complex(2 * np.dot(q0, q1))
        if qq0:
            out[(0, 0, 0, 0, 0)] = qq0 * KM
        if qq1:
            out[(0, 0, 0, 0, 1)] = qq1 * KM
        for mu in range(4):
            e = [0, 0, 0, 0, 0]; e[mu] = 1
            if q0[mu]:
                out[tuple(e)] = (2 * q0[mu]) * KM
            if q1[mu]:
                e2 = list(e); e2[4] = 1
                out[tuple(e2)] = (2 * q1[mu]) * KM
            e = [0, 0, 0, 0, 0]; e[mu] = 2
            out[tuple(e)] = KM.copy()
        return out
    raise ValueError(kind)


def block_poly(blk, qin0, qin1, qout0, qout1):
    c0 = blk.h0 - sum(qin0[mu] * blk.Hj[mu] for mu in range(4)) + sum(qout0[mu] * blk.Hl[mu] for mu in range(4))
    c1 = None
    if np.any(qin1) or np.any(qout1):
        c1 = -sum(qin1[mu] * blk.Hj[mu] for mu in range(4)) + sum(qout1[mu] * blk.Hl[mu] for mu in range(4))
    lin = [blk.Hl[mu] - blk.Hj[mu] for mu in range(4)]
    S = SD[:, None]
    return affine(S * c0, None if c1 is None else S * c1, [S * x for x in lin])


def chain(blks, Ps0, Ps1, lines, ta, tb):
    n = len(blks)
    z = np.zeros(4, dtype=complex)
    P00 = -sum((np.asarray(p, dtype=complex) for p in Ps0), z.copy())
    P01 = -sum((np.asarray(p, dtype=int) for p in Ps1), np.zeros(4, dtype=int))
    a0 = [P00] + [np.asarray(p, dtype=complex) for p in Ps0]
    a1 = [P01] + [np.asarray(p, dtype=int) for p in Ps1]
    q0 = [z.copy()]; q1 = [np.zeros(4, dtype=int)]
    for i in range(1, n + 1):
        q0.append(q0[-1] - a0[i - 1]); q1.append(q1[-1] - a1[i - 1])
    ra = Q.ROWS[ta]; rb = Q.ROWS[tb]
    T = {(0, 0, 0, 0, 0): np.eye(NS, dtype=complex)[ra]}
    for i in range(1, n + 1):
        T = pmul(T, line_poly(lines[i], q0[i], q1[i]))
        j = i + 1 if i < n else 0
        T = pmul(T, block_poly(blks[i - 1], q0[i], q1[i], q0[j], q1[j]))
        if not T:
            return {}, (q0, q1), (P00, P01)
    T = {e: M * SD[None, :] for e, M in T.items()}
    T = pmul(T, line_poly(lines[0], q0[0], q1[0]))
    T = {e: (M * SD[None, :])[:, rb] for e, M in T.items()}
    T = {e: M for e, M in T.items() if np.any(M)}
    return T, (q0, q1), (P00, P01)


def apply_letters(K, fa, fb, P0):
    """letter a: incoming k - P0 (P0 = P00 + t P01), letter b: -k; fa, fb: {exp4: complex}"""
    P00, P01 = P0
    def shift(f, sgn, c0, c1):
        res = {}
        for e, c in f.items():
            term = {(0, 0, 0, 0, 0): complex(c)}
            for mu in range(4):
                for _ in range(e[mu]):
                    lin = {}
                    ee = [0, 0, 0, 0, 0]; ee[mu] = 1
                    lin[tuple(ee)] = sgn
                    if c0[mu]:
                        lin[(0, 0, 0, 0, 0)] = complex(c0[mu])
                    if c1[mu]:
                        lin[(0, 0, 0, 0, 1)] = complex(c1[mu])
                    new = {}
                    for e1, v1 in term.items():
                        for e2, v2 in lin.items():
                            if e1[4] + e2[4] > 1:
                                continue
                            k = tuple(x + y for x, y in zip(e1, e2))
                            new[k] = new.get(k, 0) + v1 * v2
                    term = new
            for k, v in term.items():
                res[k] = res.get(k, 0) + v
        return {k: v for k, v in res.items() if v != 0}
    ga = shift(fa, 1, -P00, -P01)
    gb = shift(fb, -1, np.zeros(4), np.zeros(4, dtype=int))
    out = {}
    for e1, M in K.items():
        for e2, c2 in ga.items():
            for e3, c3 in gb.items():
                if e1[4] + e2[4] + e3[4] > 1:
                    continue
                ee = tuple(x + y + z for x, y, z in zip(e1, e2, e3))
                c = c2 * c3
                cc = complex(round(c.real), round(c.imag))
                assert abs(c - cc) < 1e-9
                out[ee] = out[ee] + M * cc if ee in out else M * cc
    return {e: M for e, M in out.items() if np.any(M)}


def lap(P):
    out = {}
    for e, v in P.items():
        for mu in range(4):
            if e[mu] >= 2:
                e2 = list(e); e2[mu] -= 2; e2 = tuple(e2)
                w = v * (e[mu] * (e[mu] - 1))
                out[e2] = out[e2] + w if e2 in out else w
    return out


# ------------------------------------------------------------------------- scalar polynomials over Q2 in (alpha, t)
def g_of(z):
    z = complex(z)
    r, i = round(z.real), round(z.imag)
    assert abs(z.real - r) < 1e-9 and abs(z.imag - i) < 1e-9
    return Q2(r, i)


def smul(A, B, E_):
    out = {}
    for e1, a in A.items():
        for e2, b in B.items():
            if e1[E_] + e2[E_] > 1:
                continue
            e = tuple(x + y for x, y in zip(e1, e2))
            v = a * b
            out[e] = out[e] + v if e in out else v
    return {e: v for e, v in out.items() if not v.iszero()}


def spow(A, n, E_):
    out = {tuple([0] * (E_ + 1)): Q2(1)}
    for _ in range(n):
        out = smul(out, A, E_)
    return out


class GSimplex:
    """E Schwinger times, shifts q_e = q0_e + t q1_e; weights for the finite part with special edge s, the special-edge
    poles, and (plain cycles) the log eps / eps^-m coefficients, each to first order in t."""

    def __init__(self, qs):
        q0, q1 = qs
        self.E = E_ = len(q0)
        self.cp = []
        for mu in range(4):
            P = {}
            for e in range(E_):
                for tt, qq in ((0, q0[e][mu]), (1, q1[e][mu])):
                    if qq != 0:
                        ex = [0] * (E_ + 1); ex[e] = 1; ex[E_] = tt
                        P[tuple(ex)] = P.get(tuple(ex), Q2(0)) + g_of(-qq)
            self.cp.append({k: v for k, v in P.items() if not v.iszero()})
        Phi = {}
        for e in range(E_):
            v0 = complex(np.dot(q0[e], q0[e])); v1 = complex(2 * np.dot(q0[e], q1[e]))
            for tt, v in ((0, v0), (1, v1)):
                if v != 0:
                    ex = [0] * (E_ + 1); ex[e] = 1; ex[E_] = tt
                    Phi[tuple(ex)] = Phi.get(tuple(ex), Q2(0)) + g_of(v)
        for mu in range(4):
            for k_, v in smul(self.cp[mu], self.cp[mu], E_).items():
                Phi[k_] = Phi.get(k_, Q2(0)) - v
        self.mPhi = {k: -v for k, v in Phi.items() if not v.iszero()}
        self._cpow = {}; self._ppow = {}; self._w = {}

    def cpow(self, mu, n):
        if (mu, n) not in self._cpow:
            self._cpow[(mu, n)] = spow(self.cp[mu], n, self.E)
        return self._cpow[(mu, n)]

    def ppow(self, n):
        if n not in self._ppow:
            self._ppow[n] = spow(self.mPhi, n, self.E)
        return self._ppow[n]

    def _poly(self, kexp, p, tdeg):
        """c^kexp (-Phi)^p as a polynomial in (alpha, t), times t^tdeg"""
        poly = self.ppow(p)
        for mu in range(4):
            if kexp[mu]:
                poly = smul(poly, self.cpow(mu, kexp[mu]), self.E)
        if tdeg:
            poly = {e[:self.E] + (e[self.E] + 1,): v for e, v in poly.items() if e[self.E] == 0}
        return poly

    def _weights(self, kind, kexp, j, m, tdeg, s):
        key = (kind, kexp, j, m, tdeg, s)
        if key in self._w:
            return self._w[key]
        E_ = self.E
        mj = j + 2 - E_
        res = {0: [Q2(0), Q2(0)], 1: [Q2(0), Q2(0)]}   # t-order -> [rational, log(E/(E-1)) or unused]
        p = mj - m
        if p < 0:
            self._w[key] = res; return res
        pref = Fraction(1, math.factorial(p) * math.factorial(j) * 4 ** j)
        poly = self._poly(kexp, p, tdeg)
        for e_, v in poly.items():
            a = e_[:E_]; tt = e_[E_]
            if kind == 'fin':
                w0, w1 = Q.fin_monomial(a, s)
                res[tt][0] = res[tt][0] + v * w0; res[tt][1] = res[tt][1] + v * w1
            elif kind == 'spec':
                res[tt][0] = res[tt][0] + v * (E.specialregion_monomial(a, s, m) / m)
            elif kind == 'log':
                res[tt][0] = res[tt][0] - v * E.dirichlet(a)
            elif kind == 'eps':
                acc = Fraction(0)
                for e in range(E_):
                    ex = list(a); ex[e] += m
                    acc += E.region_monomial(tuple(ex), e)
                res[tt][0] = res[tt][0] + v * (acc / m)
        for tt in (0, 1):
            res[tt] = [x * pref for x in res[tt]]
        self._w[key] = res
        return res

    def evaluate(self, N, kind, s=None, m=0):
        """N: {(k0..k3, t): matrix}.  Returns {t-order: [acc_rational, acc_log]} (each (re, im) Fraction arrays or None)"""
        degk = max((sum(e[:4]) for e in N), default=0)
        laps = [N]
        for j in range(1, degk // 2 + 1):
            laps.append(lap(laps[-1]))
        acc = {0: [None, None], 1: [None, None]}
        for j, Lj in enumerate(laps):
            for e, coef in Lj.items():
                kexp = e[:4]; tdeg = e[4]
                w = self._weights(kind, kexp, j, m, tdeg, s)
                for tt in (0, 1):
                    for b in (0, 1):
                        z = w[tt][b]
                        if z.iszero():
                            continue
                        if z.re:
                            acc[tt][b] = E._acc(acc[tt][b], coef, z.re)
                        if z.im:
                            acc[tt][b] = E._acc(acc[tt][b], 1j * coef, z.im)
        return {tt: [E._finish(acc[tt][0]), E._finish(acc[tt][1])] for tt in (0, 1)}


# ------------------------------------------------------------------------- configurations
def gconfig(outs, v2pos, spos, v2kind, ta, tb, fa=None, fb=None, v2dir=None, r1=None, kinds=('fin',), plain=False):
    """one configuration (as cq_core.make_config); returns {'vv': bool, 'E': E, kind: {t-order: [rat, log]}}.
    plain=True: no special edge (spos ignored), kinds among ('log', ('eps', m))."""
    seq = list(outs[:v2pos]) + [('v2',)] + list(outs[v2pos:]) if v2pos is not None else list(outs)
    n = len(seq)
    lines = ['P'] * (n + 1)
    vb = v2pos + 1 if v2pos is not None else -1
    if not plain:
        adj = {vb, vb + 1 if vb < n else 0}
        if spos in adj:
            return None
        lines[spos] = 'S'
    blks, P0s, P1s, scale, syms = [], [], [], 1, []
    for b in seq:
        if b[0] == 'cub':
            _, lt, pol, mom = b
            blks.append(Q.cubic_full(lt, pol, mom)); P0s.append(np.asarray(mom, dtype=complex)); P1s.append(np.zeros(4, dtype=int))
            scale *= 2
        elif b[0] == 'qua':
            _, (p1, m1), (p2, m2) = b
            blks.append(Q.quartic_full(p2, p1)); P0s.append(np.asarray(m1, dtype=complex) + np.asarray(m2, dtype=complex))
            P1s.append(np.zeros(4, dtype=int)); scale *= 8
        else:
            blks.append(Q.v2_block(v2kind, v2dir)); P0s.append(np.zeros(4, dtype=complex))
            P1s.append(np.zeros(4, dtype=int) if r1 is None else np.asarray(r1, dtype=int))
    order_lines = list(range(1, n + 1)) + [0]
    for pos, li in enumerate(order_lines):
        if lines[li] == 'S':
            syms.append('s')
        if pos < n:
            b = seq[pos]
            if b[0] == 'v2':
                syms.append('v')
            elif b[0] == 'cub' and Q.PARITY_T[b[1]] == 1:
                syms.append('t')
    K, qs, P0 = chain(blks, P0s, P1s, lines, ta, tb)
    if not K:
        return None
    if fa is not None or fb is not None:
        K = apply_letters(K, fa or {(0, 0, 0, 0): 1}, fb or {(0, 0, 0, 0): 1}, P0)
        if not K:
            return None
    sgn = Q.RULE['global']
    if Q.RULE['koszul'] == 'R2':            # internal vertex and special edge odd relative to the odd outputs
        for sym in ('v', 's'):
            if sym in syms:
                i = syms.index(sym)
                if sum(1 for x in syms[:i] if x == 't') % 2:
                    sgn = -sgn
    vv = False
    if not plain:
        s_idx = lines.index('S')
        vv = s_idx not in (1, 0)
        if vv and syms.index('s') < syms.index('v'):
            sgn = -sgn
    if v2pos is None and Q.PARITY_T[ta] + Q.PARITY_T[tb] == 1:
        sgn = -sgn                        # pure two-letter chain: intrinsic sign of ea_kern.kernel
    S = GSimplex(qs)
    out = {'vv': vv, 'E': len(qs[0])}
    for kd in kinds:
        if kd == 'fin':
            r = S.evaluate(K, 'fin', s=s_idx)
        elif isinstance(kd, tuple) and kd[0] == 'spec':
            r = S.evaluate(K, 'spec', s=s_idx, m=kd[1])
        elif kd == 'log':
            r = S.evaluate(K, 'log')
        elif isinstance(kd, tuple) and kd[0] == 'eps':
            r = S.evaluate(K, 'eps', m=kd[1])
        else:
            raise ValueError(kd)
        out[kd] = {tt: [None if v is None else (v[0] * Fraction(sgn, scale), v[1] * Fraction(sgn, scale)) for v in r[tt]]
                   for tt in (0, 1)}
    return out
