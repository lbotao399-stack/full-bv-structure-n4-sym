"""Holomorphic side: component entries of the user's one-loop operator, All-loop article (7.4) at n = 1
(block eq-Qnzero), with the superfield (B.6) (eq-A-13) and the fibre derivative (B.9) (eq-A-16).
  Q_1 = kappa int d^3theta C_{ADR} C_{BRE} (d_{+a}C^D)(d_+^a C^E) delta/delta C^A(theta) delta/delta C^B(theta),
  kappa = (-1)^1 i^2 / (2 (32 pi^2) (1!)^3 2) = 1/(128 pi^2).
Output normal form: coefficient of C_{A'PR} C_{B'RQ} W1^P . W2^Q, W1 attached to the first input X^{A'},
W2 to the second Y^{B'}, '.' = d_{+a} W1 d_+^a W2 in this order.  Grassmann algebra: ordered odd generators."""
import sympy as sp, itertools
g = sp.Symbol('g'); I = sp.I
# odd generators: 'th1','th2','th3' (integration theta), input letters, output letters with slot tags
PAR = {'c': 1, 'beta': 1, 'gamma': 0, 'b': 0}
class Gr:
    def __init__(self, terms=None): self.t = dict(terms or {})   # key: tuple of odd gens ; value: sympy coef (even symbols inside)
    def __add__(s, o):
        r = dict(s.t)
        for k, v in o.t.items(): r[k] = r.get(k, 0) + v
        return Gr({k: v for k, v in r.items() if sp.simplify(v) != 0})
    def scale(s, c): return Gr({k: c * v for k, v in s.t.items()})
    def __mul__(s, o):
        r = {}
        for k1, v1 in s.t.items():
            for k2, v2 in o.t.items():
                if set(k1) & set(k2): continue
                seq = list(k1) + list(k2)
                # sort with sign (bubble)
                sign = 1; seq2 = seq[:]
                for i in range(len(seq2)):
                    for j in range(len(seq2) - 1 - i):
                        if seq2[j] > seq2[j + 1]:
                            seq2[j], seq2[j + 1] = seq2[j + 1], seq2[j]; sign = -sign
                key = tuple(seq2)
                r[key] = r.get(key, 0) + sign * v1 * v2
        return Gr({k: v for k, v in r.items() if v != 0})
def odd(name): return Gr({(name,): 1})
def even(expr): return Gr({(): expr})
def dL(el, gen):
    """left derivative w.r.t. odd generator gen"""
    r = {}
    for k, v in el.t.items():
        if gen in k:
            pos = k.index(gen); sign = (-1) ** pos
            key = k[:pos] + k[pos + 1:]
            r[key] = r.get(key, 0) + sign * v
    return Gr(r)
def dE(el, sym):
    return Gr({k: sp.diff(v, sym) for k, v in el.t.items() if sp.diff(v, sym) != 0})
TH = ['th1', 'th2', 'th3']
def theta(S):
    out = even(1)
    for i in S: out = out * odd(TH[i - 1])
    return out
EPS = lambda a, b, c: sp.LeviCivita(a, b, c)
def superfield(slot):
    """d C^D(theta) for the output slot 'D' or 'E': components as generators/symbols tagged by slot"""
    out = odd('c_' + slot).scale(I * g)
    for Iq in (1, 2, 3):
        out = out + theta([Iq]) * even(sp.Symbol('gam%d_%s' % (Iq, slot)))
    for Iq in (1, 2, 3):
        for J in (1, 2, 3):
            for K in (1, 2, 3):
                e = EPS(Iq, J, K)
                if e != 0:
                    out = out + (theta([Iq, J]) * odd('beta%d_%s' % (K, slot))).scale(-sp.Rational(1, 2) * e)
    out = out + (theta([1, 2, 3]) * even(sp.Symbol('b_' + slot))).scale(1 / (I * g))
    return out
def delta_apply(el, inputs_tag):
    """delta/delta C(theta) acting on el, derivative w.r.t. the input letters with tag ('X' or 'Y'): returns Gr.
    delta = th1th2th3/(ig) dL/dc - 1/2 eps^{IJK} th_J th_K d/dgamma^I - th_I dL/dbeta_I - ig d/db"""
    tg = inputs_tag
    out = Gr()
    out = out + (theta([1, 2, 3]) * dL(el, 'c_' + tg)).scale(1 / (I * g))
    for Iq in (1, 2, 3):
        for J in (1, 2, 3):
            for K in (1, 2, 3):
                e = EPS(Iq, J, K)
                if e != 0:
                    out = out + (theta([J, K]) * dE(el, sp.Symbol('gam%d_%s' % (Iq, tg)))).scale(-sp.Rational(1, 2) * e)
    for Iq in (1, 2, 3):
        out = out + (theta([Iq]) * dL(el, 'beta%d_%s' % (Iq, tg))).scale(-1)
    out = out + dE(el, sp.Symbol('b_' + tg)).scale(-I * g)
    return out
def letter(name, tag):
    """input letter as an element: name in c, gamma1.., beta1.., b"""
    if name == 'c': return odd('c_' + tag)
    if name.startswith('beta'): return odd(name + '_' + tag)
    if name.startswith('gamma'): return even(sp.Symbol('gam' + name[5:] + '_' + tag))
    return even(sp.Symbol(name + '_' + tag))
def berezin(el):
    """int d^3theta with int d^3theta th1 th2 th3 = 1 (theta's moved to the left)"""
    r = {}
    for k, v in el.t.items():
        ths = [x for x in k if x in TH]
        if sorted(ths) != TH: continue
        rest = [x for x in k if x not in TH]
        # move the thetas to the left: count transpositions in k (thetas sorted first)
        seq = list(k); sign = 1
        target = TH + rest
        # compute permutation sign from seq to target
        perm = [target.index(x) for x in seq]
        inv = sum(1 for i in range(len(perm)) for j in range(i + 1, len(perm)) if perm[i] > perm[j])
        sign = (-1) ** inv
        key = tuple(rest)
        r[key] = r.get(key, 0) + sign * v
    return Gr(r)
def Q1_pair(Xn, Yn):
    """user's Q_1 on X Y (X with colour A', Y with colour B'); returns dict {(W1, W2): coefficient} in the normal form"""
    XY = letter(Xn, 'X') * letter(Yn, 'Y')
    KD, KE = superfield('D'), superfield('E')
    K = KD * KE
    res = {}
    for (ta, tb) in (('X', 'Y'), ('Y', 'X')):
        # delta_A(theta) delta_B(theta) with delta_A hitting tag ta, delta_B hitting tb (colour A <- ta)
        el = delta_apply(delta_apply(XY, tb), ta)
        full = berezin(K * el)
        for k, v in full.t.items():
            # identify output letters: one D-slot and one E-slot generator/symbol
            syms = list(k) + [str(s) for s in sp.sympify(v).free_symbols if str(s) not in ('g',)]
            if not full.t: continue
            # split coefficient into even symbols
            vv = sp.expand(v)
            for term in sp.Add.make_args(vv):
                evs = [s for s in term.free_symbols if s != g]
                letters = list(k) + [str(s) for s in evs]
                Dl = [x for x in letters if x.endswith('_D')]; El = [x for x in letters if x.endswith('_E')]
                assert len(Dl) == 1 and len(El) == 1, (letters, term)
                coef = term
                for s in evs: coef = coef / s
                coef = sp.simplify(coef)
                # odd generators order in key: need D-slot before E-slot for the '.' order; key is sorted lexicographically
                oddl = list(k)
                sgn = 1
                if len(oddl) == 2 and oddl[0].endswith('_E'):
                    sgn = -1      # reorder E D -> D E
                U, V = Dl[0][:-2], El[0][:-2]
                if ta == 'X':      # D-slot attached to X:  C_{A'DR} C_{B'RE} U^D . V^E
                    key2 = (U, V); c2 = sgn * coef
                else:              # D-slot attached to Y:  -> -(-1)^{|U||V|} C_{A'PR}C_{B'RQ} V^P . U^Q
                    pu = 1 if (U.startswith('c') or U.startswith('beta')) else 0
                    pv = 1 if (V.startswith('c') or V.startswith('beta')) else 0
                    key2 = (V, U); c2 = -((-1) ** (pu * pv)) * sgn * coef
                res[key2] = sp.simplify(res.get(key2, 0) + c2)
    return {k: v for k, v in res.items() if v != 0}
if __name__ == '__main__':
    kap = 1 / (128 * sp.pi ** 2)
    pairs = [('gamma1', 'beta1'), ('beta1', 'gamma1'), ('beta1', 'beta2'), ('b', 'c'), ('b', 'gamma1'), ('beta1', 'b'), ('b', 'b'),
             ('gamma1', 'gamma2'), ('c', 'c'), ('beta1', 'c'), ('gamma1', 'c'), ('beta1', 'gamma2')]
    import json
    out = {}
    for X, Y in pairs:
        r = Q1_pair(X, Y)
        out['%s %s' % (X, Y)] = {'%s.%s' % k: str(sp.simplify(kap * v)) for k, v in r.items()}
        print(X, Y, out['%s %s' % (X, Y)])
    json.dump(out, open('holo_entries.json', 'w'), indent=1)
