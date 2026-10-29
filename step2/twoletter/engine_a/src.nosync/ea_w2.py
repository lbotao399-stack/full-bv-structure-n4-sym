"""Engine A, gate G-W2: the two-body Ward identity W2 (Cor cor:tl-ward), tested as the full BRST covariance
lc(sM) = s lc(M) + {L,M} (Thm thm:tl-brst) on two-letter operators M = tr(zeta X_a X_b), zeta a constant traceless
matrix, in the planar colour-ordered reading of ea_ops (units u).  Given W1 (gate G3), this identity on all two-letter
M is equivalent to W2 (proof of Cor cor:tl-ward).

Letter-level BRST images (matrix letters, c odd), fixed by Prop c1-Q, by the W1 gate of ea_ward.py
({I,N} = -(zeta,[A,c]) and QN = -(zeta,dc) for N = (zeta,A_V)) and by s^2 = 0:
    s1 A_V = -dc - (A c - c A),   s1 A_perp = -(A c - c A),   s1 c = c c,   s1 lambda = c lambda + lambda c.
One-body bracket (Prop prop:tl-W1): {L,X} = u (gauge part of s1 X) for A_perp, lambda, c, and
{L,A_V} = u (gauge part of s1 A_V) - u (Q part of s1 A_V).
s acts on functionals from the left with Koszul signs; on a functional F evaluated on ordered legs,
  (Q F)[..., c(q), ...]  = sign * F[..., A_V(pol = -i q), ...]
  ({I,F})[..., Y(q1), c(q2), ...] = -sign F[..., Y(q1+q2), ...] for Y = A,   + for Y = lambda,
  ({I,F})[..., c(q1), Y(q2), ...] = +sign F[..., Y(q1+q2), ...] for Y = A or lambda,  ({I,F})[c(q1), c(q2)] = +sign F[c(q1+q2)],
sign = (-1)^(parity of the legs to the left).  Targets: output leg sets of field type (A_V, A_perp, lambda, c) on which
no antifield output of lc(M) can land under s (checked against census.json)."""
import itertools, json, math, os, sys, time
from fractions import Fraction
import numpy as np
import ea_core as E, ea_ops as O

SRC = os.path.dirname(os.path.abspath(__file__))
CENSUS = os.path.join(SRC, '..', '..', 'theory', 'census.json')
PAR = {'AV': 0, 'AW': 0, 'c': 1, 'lamP': 1, 'lamM': 1}
NROW = {'AV': 4, 'AW': 6, 'c': 1, 'lamP': 16, 'lamM': 16}
ONE = {(0, 0, 0, 0): 1}
NEG = {'Q': 1, 'A': 1, 'c': 1, 'lam': 1, 'oddsign': 1}   # negative controls flip one of these to -1

# ----------------------------------------------------------------------------------------------- polynomials
def pmul_lin(f, vec):
    """f(p) * sum_mu vec_mu p_mu"""
    out = {}
    for e, c in f.items():
        for mu in range(4):
            if vec[mu] == 0:
                continue
            ee = list(e); ee[mu] += 1; ee = tuple(ee)
            out[ee] = out.get(ee, 0) + c * vec[mu]
    return {e: c for e, c in out.items() if c != 0}


def psplit(f):
    """f(p1 + p2) = sum_k g_k(p1) h_k(p2), as a list of (g, h) monomial pairs"""
    out = []
    for e, c in f.items():
        for ks in itertools.product(*[range(x + 1) for x in e]):
            b = 1
            for x, k in zip(e, ks):
                b *= math.comb(x, k)
            out.append(({tuple(ks): c * b}, {tuple(x - k for x, k in zip(e, ks)): 1}))
    return out

# ----------------------------------------------------------------------------------------------- BRST on letters
def s1(X, part='all'):
    """image of a letter: list of (coef, [letters]); part 'Q', 'gauge' or 'all'"""
    t, z, f = X
    out = []
    if t == 'AV' and part in ('Q', 'all'):
        out.append((-NEG['Q'], [('c', np.array([1 + 0j]), pmul_lin(f, [1j * complex(x) for x in z]))]))
    if part in ('gauge', 'all'):
        for g, h in psplit(f):
            cl = ('c', np.array([1 + 0j]), h)
            if t in ('AV', 'AW'):
                out.append((-NEG['A'], [(t, z, g), cl])); out.append((NEG['A'], [cl, (t, z, g)]))
            elif t == 'c':
                out.append((NEG['c'], [('c', z, g), cl]))
            elif t in ('lamP', 'lamM'):
                out.append((NEG['lam'], [cl, (t, z, g)])); out.append((NEG['lam'], [(t, z, g), cl]))
    return out


def s_op(op, onebody_L=False):
    """s (or {L,-} in units u when onebody_L) applied to an operator (list of (coef, letters)) from the left"""
    out = []
    for coef, letters in op:
        sg = 1
        for i, X in enumerate(letters):
            if onebody_L:
                frag = [(c, ls) for c, ls in s1(X, 'gauge')]
                if X[0] == 'AV':
                    frag += [(-c, ls) for c, ls in s1(X, 'Q')]
            else:
                frag = s1(X, 'all')
            for c, ls in frag:
                out.append((coef * sg * c, letters[:i] + ls + letters[i + 1:]))
            if PAR[X[0]] and NEG['oddsign'] == 1:
                sg = -sg
    return out

# ----------------------------------------------------------------------------------------------- evaluations
def qc(z):
    z = complex(z)
    return (Fraction(z.real).limit_denominator(10 ** 9), Fraction(z.imag).limit_denominator(10 ** 9))


def lc_on_legs(op, legs):
    """planar lc of op = sum coef tr(zeta X_1..X_L) on ordered legs (pairs adjacent inside the word), units u"""
    n = len(legs)
    tot = E.Q0
    for coef, letters in op:
        L = len(letters)
        nout = n - (L - 2)
        if nout < 1 or L < 2:
            continue
        for p in range(L - 1):
            pre, post = letters[:p], letters[p + 2:]
            sv = complex(coef)
            for j, X in enumerate(pre):
                sv *= O.letter_on_leg(X, legs[j])
            for j, X in enumerate(post):
                sv *= O.letter_on_leg(X, legs[p + nout + j])
            if sv == 0:
                continue
            kv = O.pair_kernel_on_legs(letters[p], letters[p + 1], legs[p:p + nout])
            f = Fraction((-1) ** nout, 2)
            tot = E.q_add(tot, E.q_mul((kv[0] * f, kv[1] * f), qc(sv)))
    return tot


def tree_on_legs(op, legs):
    tot = 0j
    for coef, letters in op:
        if len(letters) != len(legs):
            continue
        v = complex(coef)
        for X, l in zip(letters, legs):
            v *= O.letter_on_leg(X, l)
        tot += v
    return qc(tot)


def s_functional_on_legs(F, legs):
    """(s F)[legs] for a functional F given as a callable on ordered field legs (no antifield outputs relevant)"""
    tot = E.Q0
    sg = 1
    for j, (t, pol, q) in enumerate(legs):
        if t == 'c':
            nl = legs[:j] + [('AV', -1j * np.asarray(q, dtype=complex), q)] + legs[j + 1:]
            tot = E.q_add(tot, E.q_scale(F(nl), Fraction(sg)))
        if j + 1 < len(legs):
            t2, pol2, q2 = legs[j + 1]
            qs = np.asarray(q) + np.asarray(q2)
            merged, c = None, 0
            if t in ('AV', 'AW') and t2 == 'c':
                merged, c = (t, pol, qs), -1
            elif t in ('lamP', 'lamM') and t2 == 'c':
                merged, c = (t, pol, qs), +1
            elif t == 'c' and t2 in ('AV', 'AW', 'lamP', 'lamM'):
                merged, c = (t2, pol2, qs), +1
            elif t == 'c' and t2 == 'c':
                merged, c = ('c', pol, qs), +1
            if merged is not None:
                nl = legs[:j] + [merged] + legs[j + 2:]
                tot = E.q_add(tot, E.q_scale(F(nl), Fraction(sg * c)))
        if PAR[t]:
            sg = -sg
    return tot

# ----------------------------------------------------------------------------------------------- targets
_CEN = None


def census():
    global _CEN
    if _CEN is None:
        _CEN = json.load(open(CENSUS))['output_leg_sets_by_pair_weight']
    return _CEN


ANTI = ('AsV', 'AsW', 'lamsP', 'lamsM', 'cs')


def blocked(pair, W, target):
    """True if an antifield-containing log leg set S of the pair at weight W could reach the target under s"""
    key = pair if pair in census() else '|'.join(reversed(pair.split('|')))
    ent = census()[key].get(str(Fraction(W)))
    if ent is None:
        raise KeyError('census has no entry %s at W=%s' % (key, W))
    tl = list(target)
    for S in ent['log']:
        legs = S.split()
        if not any(x in ANTI for x in legs):
            continue
        fields = [x for x in legs if x not in ANTI]
        rest = list(tl)
        ok = True
        for x in fields:
            if x in rest:
                rest.remove(x)
            else:
                ok = False; break
        if ok and len(tl) >= len(legs):
            return True
    return False


def rand_pol(rng, t):
    if t == 'AV': return rng.integers(-2, 3, 4).astype(complex)
    if t == 'AW': return rng.integers(-2, 3, 6).astype(complex)
    if t == 'c': return np.array([1 + 0j])
    b = E.BASIS[t]
    return sum(int(x) * v for x, v in zip(rng.integers(-1, 2, len(b)), b)).astype(complex)


def rand_mom(rng):
    while True:
        q = rng.integers(-2, 3, 4)
        if q.any():
            return q


def targets(ta, tb, nmax=3):
    """ordered field leg-type tuples with ghost number = number of c letters + 1 and the fermion parity of the pair"""
    types = ['AV', 'AW', 'lamP', 'lamM', 'c']
    nc = [ta, tb].count('c') + 1
    nl = sum(1 for x in (ta, tb) if x in ('lamP', 'lamM'))
    out = []
    for n in range(1, nmax + 1):
        for tup in itertools.product(types, repeat=n):
            if tup.count('c') != nc:
                continue
            if sum(1 for x in tup if x in ('lamP', 'lamM')) % 2 != nl % 2:
                continue
            out.append(tup)
    return out


def check_pair(ta, tb, fa, fb, W, rng, nmax=3, ntrial=1, log=None):
    rows = []
    pair = '%s|%s' % (ta, tb)
    for tup in targets(ta, tb, nmax):
        if blocked(pair, W, tup):
            rows.append({'target': ' '.join(tup), 'skipped': 'antifield output of lc(M) reaches it'})
            continue
        for it in range(ntrial):
            za, zb = rand_pol(rng, ta), rand_pol(rng, tb)
            M = [(1, [(ta, za, fa), (tb, zb, fb)])]
            legs = [(t, rand_pol(rng, t), rand_mom(rng)) for t in tup]
            lhs = lc_on_legs(s_op(M), legs)
            F = lambda lg: lc_on_legs(M, lg)
            rhs = E.q_add(s_functional_on_legs(F, legs), tree_on_legs(s_op(M, onebody_L=True), legs))
            ok = lhs == rhs
            rows.append({'target': ' '.join(tup), 'lhs': E.qstr(lhs), 'rhs': E.qstr(rhs), 'pass': ok,
                         'nontrivial': lhs != E.Q0 or rhs != E.Q0})
            if log:
                log('  %s -> %s : lhs %s rhs %s %s' % (pair, ' '.join(tup), E.qstr(lhs), E.qstr(rhs), 'ok' if ok else 'FAIL'))
    return rows


def deriv_poly(mu):
    e = [0, 0, 0, 0]; e[mu] = 1
    return {tuple(e): 1j}


def deriv2_poly(mu, nu):
    e = [0, 0, 0, 0]; e[mu] += 1; e[nu] += 1
    return {tuple(e): -1}


PAIRS = [('AV', 'AV'), ('AV', 'AW'), ('AW', 'AV'), ('AW', 'AW'), ('AV', 'c'), ('c', 'AV'), ('AW', 'c'), ('c', 'AW'),
         ('c', 'c'), ('lamP', 'lamM'), ('lamM', 'lamP'), ('lamP', 'lamP'), ('lamM', 'lamM'), ('lamP', 'c'), ('c', 'lamM'),
         ('AV', 'lamP'), ('lamM', 'AV'), ('AW', 'lamP'), ('lamM', 'AW')]


def run(out_path, seed=11, nmax=3, deriv_pairs=None, log=print):
    rng = np.random.default_rng(seed)
    res = {'description': __doc__, 'cases': []}
    t0 = time.time()
    W0 = {'AV': 1, 'AW': 1, 'c': 0, 'lamP': Fraction(3, 2), 'lamM': Fraction(3, 2)}
    for ta, tb in PAIRS:
        for dname, fa, fb, dW in (('none', ONE, ONE, 0), ('d_a', deriv_poly(int(rng.integers(0, 4))), ONE, 1),
                                  ('d_b', ONE, deriv_poly(int(rng.integers(0, 4))), 1)):
            W = W0[ta] + W0[tb] + dW
            log('pair (%s,%s) derivatives %s W=%s  [%.0fs]' % (ta, tb, dname, W, time.time() - t0))
            rows = check_pair(ta, tb, fa, fb, Fraction(W), rng, nmax=nmax, log=log)
            res['cases'].append({'pair': [ta, tb], 'derivatives': dname, 'W': str(W), 'rows': rows})
            json.dump(res, open(out_path, 'w'), indent=1, default=str)
    for ta, tb, mu, nu in (deriv_pairs or []):
        fa, fb = deriv_poly(mu), deriv_poly(nu)
        W = W0[ta] + W0[tb] + 2
        log('pair (%s,%s) derivatives d_%d, d_%d W=%s' % (ta, tb, mu, nu, W))
        rows = check_pair(ta, tb, fa, fb, Fraction(W), rng, nmax=nmax, log=log)
        res['cases'].append({'pair': [ta, tb], 'derivatives': 'd_a=%d d_b=%d' % (mu, nu), 'W': str(W), 'rows': rows})
        json.dump(res, open(out_path, 'w'), indent=1, default=str)
    allrows = [r for c in res['cases'] for r in c['rows'] if 'pass' in r]
    res['summary'] = {'tested': len(allrows), 'pass': sum(r['pass'] for r in allrows),
                      'nontrivial': sum(r['nontrivial'] for r in allrows),
                      'skipped': sum(1 for c in res['cases'] for r in c['rows'] if 'skipped' in r),
                      'seconds': round(time.time() - t0)}
    json.dump(res, open(out_path, 'w'), indent=1, default=str)
    return res


def negative_controls(seed=13):
    """flip one convention at a time; the gate must fail somewhere"""
    out = {}
    for key in ('Q', 'A', 'c', 'lam', 'oddsign'):
        NEG[key] = -1
        rng = np.random.default_rng(seed)
        nfail, ntest = 0, 0
        for ta, tb, Wv in (('AV', 'AW', 2), ('AW', 'c', 1), ('lamP', 'lamM', 3), ('c', 'c', 0), ('lamP', 'c', Fraction(3, 2))):
            for r in check_pair(ta, tb, ONE, ONE, Wv, rng, nmax=3):
                if 'pass' in r:
                    ntest += 1; nfail += (not r['pass'])
        out['flip ' + key] = {'tested': ntest, 'failed': nfail, 'detected': nfail > 0}
        NEG[key] = 1
    return out


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(SRC, '..', 'out.nosync', 'w2_gate.json')
    r = run(out, deriv_pairs=[('AV', 'AV', 0, 1), ('AV', 'AV', 2, 2), ('AV', 'AV', 1, 1), ('lamP', 'lamM', 1, 3),
                              ('lamP', 'lamM', 2, 2), ('AW', 'AW', 0, 0), ('AW', 'AW', 0, 3), ('AV', 'AW', 1, 2)])
    r['negative_controls'] = negative_controls()
    json.dump(r, open(out, 'w'), indent=1, default=str)
    print(json.dumps(r['summary']), json.dumps(r['negative_controls']))
