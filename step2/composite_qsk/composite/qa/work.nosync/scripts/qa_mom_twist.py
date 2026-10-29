"""QA third method: twist-two family tr(D_n^k Z D_n^{S-k} Z) on the two-field shell at random complex null momenta
(not collinear), P7 (b) and (c) from the scalar-sector chains of qa_mom_chains (F0 = -1/2 fixed on the k=1 chain).
Planar rule (derived in compare notes): adjacent ordered pair in a trace of length >= 3: (-1)^r c/2 u; two-site trace:
sum_w c_w u (both leading colour terms).  Evaluation on legs: sum over the assignments of the output legs."""
import json, os, sys, time, itertools
from math import comb
from fractions import Fraction as Fr
import qa_mom as M, qa_mom_chains as CH
from qa_mom import C, I, cc
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
F0 = M.F0
NV = [1, C(0, 1), 0, 0]
CH.setn(NV)


def tosym(v):
    tot = 0
    for k, x in v.items():
        f = 1 if k == '1' else sp.log(int(k[3:]))
        tot += (sp.Rational(x.r.numerator, x.r.denominator) + sp.I * sp.Rational(x.i.numerator, x.i.denominator)) * f
    return sp.expand(tot)


def dotc(a, b):
    t = C()
    for x, y in zip(a, b):
        t = t + cc(x) * cc(y)
    return t


def csym(x):
    x = cc(x)
    return sp.Rational(x.r.numerator, x.r.denominator) + sp.I * sp.Rational(x.i.numerator, x.i.denominator)


def run(S, q1, q2, kind, cache):
    legs = [q1, q2]
    def zzv(a, b, qa, qb):
        key = ('zz', a, b, tuple(map(str, qa)), tuple(map(str, qb)), str(kind))
        if key not in cache:
            cache[key] = tosym(CH.zz(a, b, qa, qb, kind))
        return cache[key]
    def azv(r, a, q):
        key = ('az', r, a, tuple(map(str, q)), str(kind))
        if key not in cache:
            cache[key] = tosym(CH.az(r, a, q, kind, NV))
        return cache[key]
    def zav(a, r, q):
        key = ('za', a, r, tuple(map(str, q)), str(kind))
        if key not in cache:
            cache[key] = tosym(CH.za(a, r, q, kind, NV))
        return cache[key]
    nq = [csym(dotc(NV, q)) for q in legs]
    def spect(L, d):
        return (sp.I * nq[L]) ** d
    f0 = sp.Rational(F0.numerator, F0.denominator)
    out = {}
    for k in range(S // 2 + 1):
        a, b = k, S - k
        # two-field: sum over chains and assignments of K(X_1, X_2), times u (both leading colour terms)
        E1 = 0
        for (i, j) in ((0, 1), (1, 0)):
            E1 += f0 * (zzv(a, b, legs[i], legs[j]) + zzv(b, a, legs[i], legs[j]))
        # three-field part: one-output kernels, -(c/2) per adjacent ordered pair, sum over assignments
        E3 = 0
        def pairAZ(r, aZ, dsp):      # pair (A, Z) -> c = F0 az
            return sum(-sp.Rational(1, 2) * f0 * azv(r, aZ, legs[L]) * spect(1 - L, dsp) for L in (0, 1))
        def pairZA(aZ, r, dsp):      # pair (Z, A) -> c = -F0 za
            return sum(-sp.Rational(1, 2) * (-f0) * zav(aZ, r, legs[L]) * spect(1 - L, dsp) for L in (0, 1))
        for r in range(a):           # D^a Z = sum_r C(a, r+1) [d^r A, d^{a-1-r} Z]
            cf = comb(a, r + 1); x = a - 1 - r; y = b
            # tr(A X Y): (A,X) spect Y ; (Y,A) spect X
            E3 += cf * (pairAZ(r, x, y) + pairZA(y, r, x))
            # - tr(X A Y): (X,A) spect Y ; (A,Y) spect X
            E3 -= cf * (pairZA(x, r, y) + pairAZ(r, y, x))
        for r in range(b):           # second letter
            cf = comb(b, r + 1); xp = a; yp = b - 1 - r
            # tr(X' A Y'): (X',A) spect Y' ; (A,Y') spect X'
            E3 += cf * (pairZA(xp, r, yp) + pairAZ(r, yp, xp))
            # - tr(X' Y' A): (Y',A) spect X' ; (A,X') spect Y'
            E3 -= cf * (pairZA(yp, r, xp) + pairAZ(r, xp, yp))
        EO = (sp.I * nq[0]) ** a * (sp.I * nq[1]) ** b + (sp.I * nq[1]) ** a * (sp.I * nq[0]) ** b
        tot = sp.expand(E1 + E3)
        if kind == 'D0':
            target = -sp.Rational(2, 3) * (4 + S) * EO
        else:
            bv = kind[1]
            ibq = sp.I * (csym(dotc(bv, q1)) + csym(dotc(bv, q2)))
            target = -sp.Rational(1, 3) * ibq * EO
        out[k] = {'two_field_pair_part': str(sp.expand(E1)), 'three_field_completion_part': str(sp.expand(E3)),
                  'total': str(tot), 'target': str(sp.expand(target)),
                  'ratio_total/target': str(sp.simplify(tot / sp.expand(target))) if sp.expand(target) != 0 else None,
                  'agree': sp.expand(tot - target) == 0}
    return out


if __name__ == '__main__':
    Smax = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    # complex null momenta (generic, not collinear): Pythagorean quadruples with an imaginary last component
    q1 = [1, 2, 2, C(0, 3)]
    q2 = [2, C(0, 3), 6, C(0, 5)]   # 4 - 9 + 36 - 25 = 6? adjust below
    q2 = [2, -3, 6, C(0, 7)]
    assert dotc(q1, q1).iszero() and dotc(q2, q2).iszero()
    b_gen = [1, 2, -1, 3]
    res = {'kinematics': {'n': str(NV), 'q1': str(q1), 'q2': str(q2), 'b_generic': b_gen}, 'by_S': {}}
    cache = {}
    for S in range(Smax + 1):
        t0 = time.time()
        r = {'D0': run(S, q1, q2, 'D0', cache), 'P_generic_b': run(S, q1, q2, ('P', b_gen), cache),
             'P_b=n': run(S, q1, q2, ('P', NV), cache)}
        r['seconds'] = round(time.time() - t0, 1)
        res['by_S'][S] = r
        print('S=%d' % S, {kk: {k: v['agree'] for k, v in vv.items()} for kk, vv in r.items() if kk != 'seconds'}, r['seconds'], 's', flush=True)
        json.dump(res, open(os.path.join(HERE, '..', 'qa_mom_twist.json'), 'w'), indent=1)
