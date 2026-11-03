"""Generators, position injection and the canonical densities of the defect-cycle sums.

S(x) at the vertex point y: S(x) = sum_mu y^mu S_mu + S_0 (position degree <= 1).  A graph with the vertex S_mu carries
the injected momentum r = t e_mu (absorbed by the last leg); its canonical density (Def c11-canonical, last leg eliminated)
    alpha = G_0 - i sum_mu d/dt G_mu ,   beta^mu = -i G_mu(0)     (value = alpha delta + beta^mu d_mu delta).
Generators (Def at-gens; qa_s1 conventions):
  Delta : S_0 = even(dv = 8)  (w/4 div v: 2 on A, 3 on lam, 0 on c),  S_mu = even(a = 2 e_mu)
  P_b   : S_0 = even(a = b)
  M_om  : S_0 = even(om),  S_mu = even(a = om e_mu column)   (v = om x)
  R_r   : S_0 = even(Wperp = r)
  xi    : S_0 = T1 + ... + T5 (xi in Sigma)
  eta   : S_0 = s0(eta) (eta in Sigma^*),  S_mu = S(rho(e_mu) eta) = S(s[mu] eta)"""
from fractions import Fraction
import numpy as np
import sa_core as SA
import sa_fit as SF

I4 = np.eye(4, dtype=int)
MINUS_I = (Fraction(0), Fraction(-1))
_GEN = {}


def gens_of(spec):
    """{'0': Gen or None, mu: Gen} for a generator spec (cached)"""
    key = repr(sorted((k, np.asarray(v).tolist() if not isinstance(v, str) else v) for k, v in spec.items()))
    if key in _GEN:
        return _GEN[key]
    k = spec['kind']
    out = {}
    if k == 'Delta':
        out['0'] = SA.gen_even('D0', dv=8.0)
        for m in range(4):
            out[m] = SA.gen_even('D%d' % m, a=2 * I4[m])
    elif k == 'P':
        out['0'] = SA.gen_even('P', a=np.asarray(spec['b']))
    elif k == 'M':                      # v = Om x;  v.A = L_v A: (v.A)_m = v.dA_m + Om_{nm} A_n  (qa_s1 om = Om^T)
        Om = np.asarray(spec['Om'], dtype=float)
        out['0'] = SA.gen_even('M0', om=Om.T)
        for m in range(4):
            out[m] = SA.gen_even('M%d' % m, a=Om[:, m])
    elif k == 'R':
        out['0'] = SA.gen_even('R0', Wperp=np.asarray(spec['R'], dtype=float))
    elif k == 'xi':
        out['0'] = SA.gen_odd(np.asarray(spec['xi']), name='xi' + str(abs(hash(tuple(np.asarray(spec['xi']).tolist())))),
                              which=spec.get('which', ('T1', 'T2', 'T3', 'T4', 'T5')))
    elif k == 'eta':
        eta = np.asarray(spec['eta'])
        out['0'] = SA.gen_s0(eta, name='s0eta')
        for m in range(4):
            xm = np.rint((SF.SS[m] @ eta).real).astype(int) + 1j * np.rint((SF.SS[m] @ eta).imag).astype(int)
            out[m] = SA.gen_odd(xm, name='eta%d_%s' % (m, abs(hash(tuple(eta.tolist())))))
    else:
        raise ValueError(k)
    _GEN[key] = out
    return out


def reduce_colour(res, nl, w, tt):
    """{colour key: 5-vector}; n=2: '(Y0,Y1)'; n=3: '(Y0,[Y1,Y2])'; n>=4: cyclic words"""
    red = SA.colour_reduce(res, nl)
    out = {}
    for key, lst in red.items():
        acc = SA.vzero()
        for item in lst:
            if nl == 3:
                c, val = item
                acc = SA.vadd(acc, SA.vscale(val[w][tt], c))
            else:
                acc = SA.vadd(acc, item[w][tt])
        out[key] = acc
    return out


def _shift_last(legs, mu):
    L2 = [dict(l) for l in legs]
    L2[-1] = dict(L2[-1]); L2[-1]['p1'] = L2[-1]['p1'] - I4[mu]
    return L2


def density(spec, legs, other=None, what=('fin',), slot_x='x', slot_y='y'):
    """canonical density of the defect-cycle sum along spec (arity one) or along (spec, other) (arity two, other
    without position part).  Returns {w: {'alpha': {colour: 5vec}, 'beta': [ {colour: 5vec} ] * 4}}"""
    g = gens_of(spec)
    nl = len(legs)
    go = None if other is None else gens_of(other)['0']
    def run(gen, L, r1=None):
        gens = [gen] if go is None else [gen, go]
        names = [slot_x] if go is None else [slot_x, slot_y]
        return SA.cycle_sum(L, gens, slotnames=names, r1=r1, what=what)
    out = {}
    res0 = run(g['0'], legs)
    for w in what:
        out[w] = {'alpha': reduce_colour(res0, nl, w, 0), 'beta': [{} for _ in range(4)]}
    for m in range(4):
        if m not in g:
            continue
        resm = run(g[m], _shift_last(legs, m), r1=I4[m])
        for w in what:
            d1 = reduce_colour(resm, nl, w, 1)
            d0 = reduce_colour(resm, nl, w, 0)
            A = out[w]['alpha']
            for ck, v in d1.items():
                A[ck] = SA.vadd(A.get(ck, SA.vzero()), SA.vscale(v, MINUS_I))
            for ck, v in d0.items():
                out[w]['beta'][m][ck] = SA.vscale(v, MINUS_I)
    return out


def sign2(spec_x):
    """(-1)^{|x|} for the arity-two normalisation"""
    return -1 if spec_x['kind'] in ('xi', 'eta') else 1


def xi_generic(chir='P', coef=(1, 2, -1, 1, 0, 1, -2, 1)):
    B = SA.E.BASIS['lamP' if chir == 'P' else 'lamM']
    return sum(c * np.asarray(v) for c, v in zip(coef, B)).astype(complex)


def eta_generic(chir='P', coef=(1, -1, 2, 1, 1, 0, -1, 2)):
    B = SA.E.BASIS['lamsP' if chir == 'P' else 'lamsM']
    return sum(c * np.asarray(v) for c, v in zip(coef, B)).astype(complex)
