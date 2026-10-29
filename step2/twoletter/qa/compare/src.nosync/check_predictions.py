"""Task 2: predictions.json and the Chapter 5/6/9 values against the agreed W_0 kernels (engine A tables, equal to
engine B's entry by entry, cmp_W0.py).  Independent assembly (comparator's own code):
  P1  Ch5 graph table from the AW|AW block tables, contracted with letters iota_w1 A, iota_w2 A;
  P2  lc and lc^natural on the scalar bilinears (20', Konishi);
  P6  su(2) sector: planar two-site density from the agreed AW|AW -> AW AW tables, cyclic spectra L <= 8 vs P6;
  P5  twist two: (i) low-order entries of engine B's exported sl(2) data (two_site.json) re-derived from the agreed
      derivative tables; (ii) the comparator's own assembly of lc^natural on tr(D^k Z D^{N-k} Z) (N <= 8) from those
      data, spectrum vs -8 h(S) u; (iii) eigenvalues of engine A's exported twist-two matrices.
Colour reading (Lemma lem:tl-colour, checked by hand in compare_report): adjacent pair of a single trace, r outputs:
(-1)^r (value/2) u tr(... Y_1 ... Y_r ...), Y_1 nearest to letter a."""
import os, json, itertools
from fractions import Fraction
import sympy as sp
import numpy as np
import cmp_lib as L

A = json.load(open(os.path.join(L.EA, 'out.nosync', 'kernels_W0_raw.json')))
P = json.load(open(os.path.join(L.TL, 'theory', 'predictions.json')))
out = {}


def tensor(pair, ls, key, kind='log'):
    """exact (sympy) constant tensor T[ia,ib,io...] of an A table entry (no momenta)"""
    e = A[pair]['leg_sets'][ls].get(key, {}).get(kind, {})
    T = {}
    for idx, s in e.items():
        parts = [int(x) for x in idx.split(',') if x != '']
        T[tuple(parts)] = sp.nsimplify(s.replace('^', '**').replace('i', 'I'))
    return T


def tadd(*Ts):
    R = {}
    for T in Ts:
        for k, v in T.items():
            R[k] = R.get(k, 0) + v
    return R


# ---------------------------------------------------------------- P1: Ch5 graph table on O = int(iota_w1 A, iota_w2 A)
def bilinear_fit(T, n):
    """T[I,J,o1,o2] (letters I,J; outputs o1 (nearest a), o2): fit the symmetric output functional
    F = sum_{o1,o2} T[I,J,o1,o2] (A_o1, A_o2) as  c1 (w1-w2 swap-symmetric) + c2 delta_IJ delta_o1o2.
    Returns (c1, c2, residual)."""
    c1, c2 = sp.symbols('c1 c2')
    eqs = []
    for I, J, o1, o2 in itertools.product(range(6), range(6), range(n), range(n)):
        val = (T.get((I, J, o1, o2), 0) + T.get((I, J, o2, o1), 0)) / 2
        pred = 0
        if n == 6:
            pred += c1 * (sp.KroneckerDelta(I, o1) * sp.KroneckerDelta(J, o2) + sp.KroneckerDelta(I, o2) * sp.KroneckerDelta(J, o1)) / 2
        pred += c2 * sp.KroneckerDelta(I, J) * sp.KroneckerDelta(o1, o2)
        eqs.append(sp.expand(val - pred))
    sol = sp.solve([e for e in eqs if e != 0][:40], [c1, c2], dict=True)
    s = sol[0] if sol else {}
    res = max(abs(sp.N(e.subs(s).subs({c1: 0, c2: 0}))) for e in eqs)
    return s.get(c1, 0), s.get(c2, 0), res


e_Q = tensor('AW|AW', 'AW AW', 'AW AW | Q')
f_CC = tensor('AW|AW', 'AW AW', 'AW AW | CC')
c_Q = tensor('AW|AW', 'AV AV', 'AV AV | Q')
d_CC = tensor('AW|AW', 'AV AV', 'AV AV | CC')
p1 = {}
for name, T, n in (('(e) quartic bubble -> A_perp A_perp', e_Q, 6), ('(f) triangle -> A_perp A_perp', f_CC, 6),
                   ('(c) quartic bubble -> A_V A_V', c_Q, 4), ('(d) triangle -> A_V A_V', d_CC, 4)):
    c1, c2, res = bilinear_fit(T, n)
    p1[name] = dict(O=str(c1), w1w2_int=str(c2), residual=float(res))
b_C = A['AW|AW']['leg_sets'].get('AV', {})
p1['(b) cubic -> A_V'] = 'zero table' if not any(v for d in b_C.values() for v in d.values()) else 'nonzero'
tot6 = tadd(e_Q, f_CC); tot4 = tadd(c_Q, d_CC)
c1, c2, _ = bilinear_fit(tot6, 6); d1, d2, _ = bilinear_fit(tot4, 4)
# traceless: (w1,w2) = 0 -> lc = c1 O ; trace O = sum_I (A_I,A_I): lc = c1 O + 6 c2 O (c2 multiplies (w1,w2) int(A_perp,A_perp))
lc_traceless = c1; lc_trace = c1 + 6 * c2
p1['totals'] = dict(lc_traceless=str(lc_traceless), lc_trace=str(lc_trace), AV_AV_total=str((d1, d2)))
p1['predicted'] = P['P1_scalar_bilinears_Ch5']
p1['pass'] = bool(p1['(e) quartic bubble -> A_perp A_perp']['O'] == '2' and p1['(e) quartic bubble -> A_perp A_perp']['w1w2_int'] == '-2'
                  and p1['(f) triangle -> A_perp A_perp']['O'] == '2' and p1['(f) triangle -> A_perp A_perp']['w1w2_int'] == '0'
                  and p1['(c) quartic bubble -> A_V A_V']['w1w2_int'] == '-2' and p1['(d) triangle -> A_V A_V']['w1w2_int'] == '2'
                  and lc_traceless == 4 and lc_trace == -8)
# P2: lc^natural = lc - u{F,O}, {F,-} = 2 per A letter -> -4
p2 = dict(lc_natural_traceless=str(lc_traceless - 4), lc_natural_trace=str(lc_trace - 4),
          gamma_Konishi_units=str(2 * (lc_trace - 4)), predicted=P['P2_lc_natural']['values_units'])
p2['pass'] = bool(lc_traceless - 4 == 0 and lc_trace - 4 == -12)
out['P1_Ch5_graph_table'] = p1
out['P2_lc_natural'] = p2
print('P1', p1['pass'], {k: v for k, v in p1.items() if k not in ('predicted',)})
print('P2', p2)

# ---------------------------------------------------------------- P6: su(2) sector from the agreed tables
T6 = tadd(e_Q, f_CC)          # word AW AW (Y1 nearest a), all block structures
# density h: X_a X_b -> sum_{o1,o2} (T/2) Y_o1 Y_o2   (r = 2, sign +)
h = {}
for (I, J, o1, o2), v in T6.items():
    if v != 0:
        h[(I, J, o1, o2)] = v / 2
dens_ok = all(sp.simplify(h.get((I, J, o1, o2), 0) - (2 * (o1 == J) * (o2 == I) - (I == J) * (o1 == o2))) == 0
              for I, J, o1, o2 in itertools.product(range(6), repeat=4))
# complex null letters Z = e_0 + i e_1, X = e_2 + i e_3 (SO(6) components of engine A, any orthonormal pair)
Zv = [1, sp.I, 0, 0, 0, 0]; Xv = [0, 0, 1, sp.I, 0, 0]
# output projections: coefficient of Z (resp X) in A_o: dual vectors zbar/2
proj = {'Z': [sp.Rational(1, 2), -sp.I / 2, 0, 0, 0, 0], 'X': [0, 0, sp.Rational(1, 2), -sp.I / 2, 0, 0]}
vec = {'Z': Zv, 'X': Xv}
hs = {}
for a, b in itertools.product('ZX', repeat=2):
    for c, d in itertools.product('ZX', repeat=2):
        s = 0
        for (I, J, o1, o2), v in h.items():
            s += vec[a][I] * vec[b][J] * v * proj[c][o1] * proj[d][o2]
        hs[(a, b, c, d)] = sp.nsimplify(sp.expand(s))
sl_ok = all(hs[(a, b, c, d)] == 2 * (c == b) * (d == a) for a, b, c, d in hs)


def cyclic_charpoly(Lc, M):
    """E = lc^natural/(-2u) on cyclic (trace) states; trace words = necklaces; we act on the full word space and
    project with the cyclic symmetriser."""
    words = [''.join(p) for p in itertools.product('ZX', repeat=Lc) if p.count('X') == M]
    wi = {w: i for i, w in enumerate(words)}
    H = sp.zeros(len(words), len(words))
    for w in words:
        for l in range(Lc):
            a, b = w[l], w[(l + 1) % Lc]
            for c, d in itertools.product('ZX', repeat=2):
                co = hs[(a, b, c, d)]
                if co == 0:
                    continue
                nw = list(w); nw[l] = c; nw[(l + 1) % Lc] = d
                H[wi[''.join(nw)], wi[w]] += co
        H[wi[w], wi[w]] -= 2 * Lc            # - u n_F, n_F = 2 per scalar letter
    H = H / (-2)                             # E = lc^natural / (-2u)
    # cyclic invariant subspace
    def rot(w):
        return w[1:] + w[:1]
    reps, basis = [], []
    seen = set()
    for w in words:
        if w in seen:
            continue
        orb = []
        x = w
        for _ in range(Lc):
            orb.append(x); x = rot(x)
        seen |= set(orb)
        v = sp.zeros(len(words), 1)
        for x in orb:
            v[wi[x]] += 1
        basis.append(v)
    B = sp.Matrix.hstack(*basis)
    HB = H * B
    # solve B C = H B
    C = (B.T * B).inv() * B.T * HB
    assert sp.simplify(B * C - HB) == sp.zeros(*HB.shape), 'cyclic subspace not invariant'
    E = sp.symbols('E')
    return sp.factor(C.charpoly(E).as_expr())


E = sp.symbols('E')
p6 = dict(density_equals_2P_minus_K=bool(dens_ok), su2_two_site_equals_2P=bool(sl_ok), sectors={})
allok = True
for key, val in P['P6_su2_Heisenberg']['cyclic_spectra_E'].items():
    Lc = int(key.split(',')[0][2:]); M = int(key.split(',')[1][2:])
    cp = cyclic_charpoly(Lc, M)
    pred = 1
    for f in val['charpoly_factors']:
        pred *= sp.sympify(f['factor'].replace('E', '(E)'), locals={'E': E}) ** f['multiplicity']
    ok = sp.expand(cp - pred) == 0
    allok &= ok
    p6['sectors'][key] = dict(charpoly=str(cp), match=bool(ok))
p6['pass'] = bool(allok and dens_ok and sl_ok)
out['P6_su2_from_agreed_tables'] = p6
print('P6', p6['pass'], dens_ok, sl_ok)
json.dump(out, open(os.path.join(L.QA, 'out.nosync', 'predictions_check_part1.json'), 'w'), indent=1)
