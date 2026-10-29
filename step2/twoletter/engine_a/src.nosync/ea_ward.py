"""Engine A, gate G3 (Ward W1, Prop prop:tl-W1) in the planar reading with a constant traceless spectator zeta:
  M = (zeta, [X_a, X_b]) = kappa [tr(zeta X_a X_b) -+ tr(zeta X_b X_a)] (graded commutator), pairs adjacent once.
  A_perp: lc((zeta,[A_perp,c])) = u (zeta,[A_perp,c])
  lambda: lc((zeta,[c,lambda])) = u (zeta,[c,lambda])
  c:      lc((zeta,[c,c]))      = u (zeta,[c,c])
  A_V:    lc((zeta,[A_V,c]))    = u (zeta,[A_V,c]) - u (zeta, dc)      (s A = -d_A c, Q A = -dc, Prop c1-Q)
Evaluation: coefficient of tr(zeta T^{c_1}..T^{c_n}) with the legs in trace order; units u (planar factor (-1)^r/2)."""
import itertools
from fractions import Fraction
import numpy as np
import ea_core as E, ea_ops as O, ea_kern as K

ONE = O.ONE


def lc_const_trace(terms, legs):
    """terms: list of (coef, [letters]) meaning coef * tr(zeta X_1 ... X_L) (zeta fixed, not cyclic)"""
    n = len(legs)
    tot = E.Q0
    for coef, letters in terms:
        L = len(letters)
        nout = n - (L - 2)
        if nout < 1:
            continue
        for p in range(L - 1):
            Xa, Xb = letters[p], letters[p + 1]
            pre, post = letters[:p], letters[p + 2:]
            sv = 1
            for j, X in enumerate(pre):
                sv *= O.letter_on_leg(X, legs[j])
            for j, X in enumerate(post):
                sv *= O.letter_on_leg(X, legs[p + nout + j])
            if sv == 0:
                continue
            outs = legs[p:p + nout]
            kv = O.pair_kernel_on_legs(Xa, Xb, outs)
            f = Fraction((-1) ** nout, 2)
            c = complex(coef) * complex(sv)
            kv = (kv[0] * f, kv[1] * f)
            tot = E.q_add(tot, E.q_mul(kv, (Fraction(round(c.real)), Fraction(round(c.imag)))))
    return tot


def tree_const_trace(terms, legs):
    tot = 0j
    for coef, letters in terms:
        if len(letters) != len(legs):
            continue
        v = complex(coef)
        for X, l in zip(letters, legs):
            v *= O.letter_on_leg(X, l)
        tot += v
    return tot


def rv(rng, t):
    if t == 'AV': return rng.integers(-2, 3, 4).astype(complex)
    if t == 'AW': return rng.integers(-2, 3, 6).astype(complex)
    if t == 'c': return np.array([1 + 0j])
    b = E.BASIS[t]
    return sum(int(x) * v for x, v in zip(rng.integers(-1, 2, len(b)), b)).astype(complex)


def gate(seed=1, ntrials=3):
    rng = np.random.default_rng(seed)
    out = {}
    D = {(0, 0, 0, 0): 1}
    for name, ta, tb, sign in (('A_perp', 'AW', 'c', -1), ('lambda+', 'c', 'lamP', +1), ('lambda-', 'c', 'lamM', +1),
                               ('c', 'c', 'c', +1), ('A_V', 'AV', 'c', -1)):
        rows = []
        for it in range(ntrials):
            za, zb = rv(rng, ta), rv(rng, tb)
            Xa, Xb = (ta, za, D), (tb, zb, D)
            # M = tr(zeta Xa Xb) + sign * tr(zeta Xb Xa)  (sign -1: [A,c] = Ac - cA; +1: [c,lam] = c lam + lam c)
            if ta == tb:
                M = [(1, [Xa, Xb]), (sign, [Xb, Xa])]
            else:
                M = [(1, [Xa, Xb]), (sign, [Xb, Xa])]
            q1, q2 = rng.integers(-2, 3, 4), rng.integers(-2, 3, 4)
            ea, eb = rv(rng, ta), rv(rng, tb)
            res = {}
            # two-field outputs, both leg orders
            for legs in ([(ta, ea, q1), (tb, eb, q2)], [(tb, eb, q1), (ta, ea, q2)]):
                lhs = lc_const_trace(M, legs)
                rhs = tree_const_trace(M, legs)
                res[str([l[0] for l in legs])] = (E.qstr(lhs), str(rhs))
                rows.append(('two-field', lhs, rhs))
            if name == 'A_V':
                # one-field output c: expected -u (zeta, d c) evaluated on the leg c(q): -(i q . za) * (zeta.e)
                legs = [('c', np.array([1 + 0j]), q1)]
                lhs = lc_const_trace(M, legs)
                rhs = -(1j * complex(np.dot(za, q1)))
                rows.append(('one-field', lhs, rhs))
                res['one-field c'] = (E.qstr(lhs), str(rhs))
            out.setdefault(name, []).append(res)
        ok = all(abs(complex(float(l[0]), float(l[1])) - r) < 1e-12 for _, l, r in rows)
        out[name + '_pass'] = ok
    return out


if __name__ == '__main__':
    import json
    r = gate()
    print(json.dumps(r, indent=1))
