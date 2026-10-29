"""Independent implementation of Beisert's harmonic action [B03c, (4.27)-(4.29)] on two-site graded cyclic states,
to check Prop c11-targets (iv), (v) (and (i)-(iii), the 20'/Konishi values) of Chapter 11.
Oscillators: a1,a2,b1,b2 bosonic, c1..c4 fermionic; central charge n_a - n_b + n_c - 2 = 0 per site."""
import itertools, json, sys
from fractions import Fraction as Fr
from math import gamma
import sympy as sp
from pathlib import Path
TYPES = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'c3', 'c4']
FERM = {t: t[0] == 'c' for t in TYPES}
CH = {'a': 1, 'b': -1, 'c': 1}
def hnum(n):
    return sum(Fr(1, j) for j in range(1, n+1))
def cfun(n, n12, n21):
    if n12 == 0 and n21 == 0:
        assert n % 2 == 0
        return sp.Rational(hnum(n//2))
    return (-1)**(1+n12*n21)*sp.gamma(sp.Rational(n12+n21, 2))*sp.gamma(1+sp.Rational(n-n12-n21, 2))/sp.gamma(1+sp.Rational(n, 2))
# a site state: tuple of counts per type (fermions 0/1); two-site basis state: (site1, site2)
def site_charge(cnt):
    return sum(CH[t[0]]*k for t, k in zip(TYPES, cnt)) - 2
def expand(state):
    """list of oscillators (site, typeindex) in canonical order: site 1 then site 2, types in TYPES order"""
    ops = []
    for s, cnt in ((1, state[0]), (2, state[1])):
        for ti, k in enumerate(cnt):
            ops += [(s, ti)]*k
    return ops
def canon(ops):
    """stable bubble sort by (site, type); sign from transpositions of two fermions"""
    ops = list(ops); sign = 1
    for i in range(len(ops)):
        for j in range(len(ops)-1-i):
            if ops[j] > ops[j+1]:
                if FERM[TYPES[ops[j][1]]] and FERM[TYPES[ops[j+1][1]]]:
                    sign = -sign
                ops[j], ops[j+1] = ops[j+1], ops[j]
    c1 = [0]*len(TYPES); c2 = [0]*len(TYPES)
    for s, ti in ops:
        (c1 if s == 1 else c2)[ti] += 1
    if any(c1[i] > 1 or c2[i] > 1 for i, t in enumerate(TYPES) if FERM[t]):
        return None, 0
    return (tuple(c1), tuple(c2)), sign
def H12(state):
    ops = expand(state); n = len(ops)
    out = {}
    for sp_ in itertools.product((1, 2), repeat=n):
        n12 = sum(1 for (s, _), t in zip(ops, sp_) if s == 1 and t == 2)
        n21 = sum(1 for (s, _), t in zip(ops, sp_) if s == 2 and t == 1)
        new = [(t, ti) for (s, ti), t in zip(ops, sp_)]
        st, sg = canon(new)
        if st is None:
            continue
        if site_charge(st[0]) != 0 or site_charge(st[1]) != 0:
            continue
        out[st] = out.get(st, 0) + sg*cfun(n, n12, n21)
    return out
def parity(cnt):
    return sum(k for t, k in zip(TYPES, cnt) if FERM[t]) % 2
def spectrum(states):
    """states: list of two-site basis states spanning a space closed under H12 and the graded swap"""
    idx = {s: i for i, s in enumerate(states)}
    N = len(states)
    H = sp.zeros(N, N)
    for s in states:
        for t, v in H12(s).items():
            if t not in idx:
                raise SystemExit(f'leaves family: {t}')
            H[idx[t], idx[s]] += v
    Pm = sp.zeros(N, N)
    for s in states:
        sw = (s[1], s[0]); sg = -1 if parity(s[0]) and parity(s[1]) else 1
        Pm[idx[sw], idx[s]] += sg
    Ht = H + Pm*H*Pm
    # graded cyclic states: P v = v
    Proj = (sp.eye(N) + Pm)/2
    cols = Proj.columnspace()
    if not cols:
        return []
    B = sp.Matrix.hstack(*cols)
    M = (B.T*B).inv()*B.T*Ht*B
    assert sp.simplify(Ht*B - B*M) == sp.zeros(N, B.shape[1])
    ev = []
    for e, m in M.eigenvals().items():
        ev += [sp.nsimplify(sp.simplify(e))]*m
    return sorted(ev, key=lambda e: float(e))
def cnt(**kw):
    return tuple(kw.get(t, 0) for t in TYPES)
def add(x, y):
    return tuple(i+j for i, j in zip(x, y))
D = cnt(a1=1, b1=1)
def Dk(k, X):
    return add(X, tuple(k*i for i in D))
Zl = cnt(c3=1, c4=1)
F = cnt(a1=2)
psi = cnt(a1=1, c4=1)
psi3 = cnt(a1=1, c3=1); psi4 = cnt(a1=1, c4=1)
def fam(pairs, Sp):
    sts = set()
    for X, Y in pairs:
        for k in range(Sp+1):
            sts.add((Dk(k, X), Dk(Sp-k, Y))); sts.add((Dk(Sp-k, Y), Dk(k, X)))
    return sorted(sts)
res = {}
# scalars without derivatives: all 6 Phi_ij = c_i c_j
phis = [cnt(**{f'c{i}': 1, f'c{j}': 1}) for i in range(1, 5) for j in range(i+1, 5)]
res['scalar_bilinears'] = [str(e) for e in spectrum(sorted({(x, y) for x in phis for y in phis}))]
for Sp in range(0, 7):
    r = {}
    r['twist2_ZZ'] = [str(e) for e in spectrum(fam([(Zl, Zl)], Sp))]
    if Sp <= 5:
        r['fermion_psipsi'] = [str(e) for e in spectrum(fam([(psi, psi)], Sp))]
        r['gluon_FF'] = [str(e) for e in spectrum(fam([(F, F)], Sp))]
        r['F_psi'] = [str(e) for e in spectrum(fam([(F, psi)], Sp))]
        r['psi3psi4_ZF'] = [str(e) for e in spectrum(fam([(psi3, psi4), (Zl, F)], Sp))]
    res[f"S'={Sp}"] = r
    print(Sp, json.dumps(r), flush=True)
Path(__file__).resolve().parent.parent.joinpath('harmonic.json').write_text(json.dumps(res, indent=1))
print(json.dumps(res['scalar_bilinears']))
