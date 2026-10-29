"""Engine A, part C: composite kernels on legs with the general core (Gaussian momenta, exact r-derivative).

kernel(Xa, Xb, legs, kind) -> (letter part, vertex-vertex part), two values in (Q+iQ) (x) span{1,log2,log3,log5,log7};
C^(2) = letter part + kappa * vertex-vertex part (kappa = +1, -1, 0: the three conventions VV_a, VV_b, VV_0; the
mirror-symmetry gate holds for each; see notes.txt).
kinds: ('P', b): C^(2)_{P_b};  'D0': C^(2),0_Delta = s_0 graph + sum_mu (-i) d/dr_mu (s_mu graph at r = t e_mu);
('Ds', mu): s_mu graph at r = 0 (= 2 C^(2)_{P_mu}, the position part);  ('Pspec', b, m): eps^-m of Z_{P_b};
('Plog', b), ('Peps', b, m): plain composite cycles along P_b (log eps, eps^-m) -- letter part only.
Letters (ltype, zeta, deriv); legs (ltype, pol, momentum (Gaussian-integer 4-vector))."""
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_gcore as G, cq_eval as V
import ea_kern as EK

vzero, vadd, vscale, vstr, viszero = V.vzero, V.vadd, V.vscale, V.vstr, V.viszero


def _val(rr, E_, za, zb):
    """[rational, log(E/(E-1))] at one t-order -> value"""
    r0, r1 = rr
    return V.from_fin((r0, r1, E_), za, zb)


def logkernel(Xa, Xb, legs, kd='log'):
    """two-letter kernel K(X_a, X_b) (no composite line, no special edge) at Gaussian momenta: log eps coefficient
    (kd='log') or ('eps', m); intrinsic sign as ea_kern.kernel.  Units hbar/(16 pi^2), stripped colour word."""
    ta, za, fa = Xa; tb, zb, fb = Xb
    tot = vzero()
    for st in V.structures(legs):
        outs = V.blocks_of(st, legs)
        res = G.gconfig(outs, None, None, None, ta, tb, fa=fa, fb=fb, kinds=(kd,), plain=True)
        if not res:
            continue
        tot = vadd(tot, _val(res[kd][0], res['E'], za, zb))
    return tot


def _sum(Xa, Xb, legs, v2kind, v2dir=None, r1=None, kd='fin', tord=0, plain=False):
    ta, za, fa = Xa; tb, zb, fb = Xb
    tl, tv = vzero(), vzero()
    for st in V.structures(legs):
        outs = V.blocks_of(st, legs)
        nb = len(outs)
        for v2pos in range(nb + 1):
            for spos in ([None] if plain else range(nb + 2)):
                res = G.gconfig(outs, v2pos, spos, v2kind, ta, tb, fa=fa, fb=fb, v2dir=v2dir, r1=r1, kinds=(kd,),
                                plain=plain)
                if not res:
                    continue
                v = _val(res[kd][tord], res['E'], za, zb)
                if res['vv']:
                    tv = vadd(tv, v)
                else:
                    tl = vadd(tl, v)
    return tl, tv


INTR = {'on': True}       # intrinsic pair sign (-1)^{|X_a|+|X_b|} as for the two-letter kernels (ea_kern.kernel)


def kernel(Xa, Xb, legs, kind):
    tl, tv = _kernel(Xa, Xb, legs, kind)
    if INTR['on'] and (Q.PARITY_T[Xa[0]] + Q.PARITY_T[Xb[0]]) % 2 == 1:
        m1 = (Fraction(-1), Fraction(0))
        return vscale(tl, m1), vscale(tv, m1)
    return tl, tv


def _kernel(Xa, Xb, legs, kind):
    if isinstance(kind, tuple) and kind[0] == 'P':
        return _sum(Xa, Xb, legs, 'P', v2dir=kind[1])
    if isinstance(kind, tuple) and kind[0] == 'Pv':         # C^(2)_{P_beta}, beta a 4-vector (tuple)
        return _sum(Xa, Xb, legs, 'Pv', v2dir=np.asarray(kind[1]))
    if isinstance(kind, tuple) and kind[0] == 'Ds':
        return _sum(Xa, Xb, legs, 's', v2dir=kind[1])
    if kind == 'D0':
        tl, tv = _sum(Xa, Xb, legs, 's0')
        mi = (Fraction(0), Fraction(-1))
        for mu in range(4):
            r1 = np.zeros(4, dtype=int); r1[mu] = 1
            a, b = _sum(Xa, Xb, legs, 's', v2dir=mu, r1=r1, tord=1)
            tl = vadd(tl, vscale(a, mi)); tv = vadd(tv, vscale(b, mi))
        return tl, tv
    if isinstance(kind, tuple) and kind[0] == 'Pspec':
        return _sum(Xa, Xb, legs, 'P', v2dir=kind[1], kd=('spec', kind[2]))
    if isinstance(kind, tuple) and kind[0] == 'Plog':
        return _sum(Xa, Xb, legs, 'P', v2dir=kind[1], kd='log', plain=True)
    if isinstance(kind, tuple) and kind[0] == 'Peps':
        return _sum(Xa, Xb, legs, 'P', v2dir=kind[1], kd=('eps', kind[2]), plain=True)
    raise ValueError(kind)


def combine(pair, kappa):
    tl, tv = pair
    return vadd(tl, vscale(tv, (Fraction(kappa), Fraction(0))))
