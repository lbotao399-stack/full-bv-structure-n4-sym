"""Step 2, Part D4, engine A, task (1): two-point sector of the arity-one functional d/dxi W(P(eps,infty), I - I^CT + xi S_M)
at one loop, and the one-loop two-point functions of W itself (self-energies), exact eps-expansions.

Pieces (closed forms of the 1PR pieces from d4_signs, see d4_notes.txt):
  internal insertion  Pi_int : 1PI cycles with the S_M vertex on the cycle (S_M's (A, A^*) term joined to an A--A and a
                               c--A^*|_V component of P: a composite A|_V--c line), insertions S_M, X, Y;
  external leg        C_P(S_M, W_1^(2)) : S_M with one leg external, the other joined by P to W_1^(2) = Sigma - I^CT_2.
Output per sector: canonical form alpha delta(Sigma) + beta^mu (d_mu delta)(Sigma) (first leg = field leg, momentum p fixed),
fitted on the invariant structures of the sector, coefficients over the basis of d4_integrals.BASIS, in units
hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y))."""
import os, sys, json, time, itertools
from fractions import Fraction
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import s2_core as C
from r2_core import Q2
import d4_cycles as DC
import d4_integrals as DI
from s2_graphs import Poly

TVALS = [0, 1, -1, 2, -2, 3, -3]


def lagrange_d0(ts):
    """weights w_t with f'(0) = sum w_t f(t) for polynomials of degree < len(ts)"""
    out = []
    for i, ti in enumerate(ts):
        # derivative at 0 of the Lagrange basis polynomial l_i
        others = [t for j, t in enumerate(ts) if j != i]
        denom = Fraction(1)
        for t in others:
            denom *= (ti - t)
        # l_i(x) = prod (x - t)/denom ; l_i'(0) = sum_k prod_{t != t_k} (0 - t) / denom
        s = Fraction(0)
        for k in range(len(others)):
            pr = Fraction(1)
            for j, t in enumerate(others):
                if j != k:
                    pr *= (0 - t)
            s += pr
        out.append(s / denom)
    return out


LW = lagrange_d0(TVALS)


def acc_add(A, B, c=Q2(1)):
    for k, v in B.items():
        A[k] = A[k] + v * c if k in A else v * c
    return A


def groups_of(qs):
    return [i - 1 for i in range(1, len(qs)) if any(int(x) != 0 for x in qs[i])]


def cycle_value(insertions, canon_odd, scale, degree, p2, extra_sign=1):
    """1/2 sum over sequences of the expansion (dict over BASIS) of 16 pi^2 x (cycle value); odd parameters reordered
    to canon_odd; S_M insertions carry the factor -1 (Hessian of xi S_M)"""
    cyc = DC.Cycle(insertions, degree)
    tot = {}
    nseq = 0
    for order, Ls in cyc.sequences():
        N, qs = cyc.numerator_fn(order, Ls)
        poly = DC.poly_from(N, degree)
        if not poly.t:
            continue
        nseq += 1
        sgn = cyc.parity_sign(order, canon_odd) * extra_sign
        terms = DI.terms(poly, qs, degree)
        G2 = groups_of(qs)
        val = DI.expand_two_point(terms, len(qs), G2, p2)
        acc_add(tot, val, Q2(Fraction(sgn, 2) / scale))
    return tot, nseq


def cycle_derivative(build, canon_odd, scale, degree, p2):
    """d/dt at t=0 of the cycle value of insertions build(t) (list), exact: N_t interpolated in k at each t in TVALS,
    Lagrange derivative in t; offsets linear in t"""
    cycs = {t: DC.Cycle(build(t), degree) for t in TVALS}
    c0 = cycs[0]
    tot = {}
    for order, Ls in c0.sequences():
        polys = {}
        qss = {}
        for t in TVALS:
            N, qs = cycs[t].numerator_fn(order, Ls)
            polys[t] = DC.poly_from(N, degree)
            qss[t] = qs
        if all(not polys[t].t for t in TVALS):
            continue
        dN = Poly(4)
        for w, t in zip(LW, TVALS):
            dN = dN + polys[t].scale(Q2(w))
        qs0 = qss[0]
        dqs = [np.array([Fraction(int(qss[1][i][m]) - int(qs0[i][m])) for m in range(4)]) for i in range(len(qs0))]
        # check linearity of offsets
        for i in range(len(qs0)):
            assert all(int(qss[2][i][m]) - int(qs0[i][m]) == 2 * dqs[i][m] for m in range(4))
        sgn = c0.parity_sign(order, canon_odd)
        G2 = groups_of(qs0)
        terms = DI.terms(polys[0], qs0, degree, dN=dN, dqs=dqs)
        val = DI.expand_two_point(terms, len(qs0), G2, p2)
        acc_add(tot, val, Q2(Fraction(sgn, 2) / scale))
    return tot


# ------------------------------------------------------------------ sectors
def pol(rng, t, which=None):
    return C.rpol(rng, t, which)


SECT = {
    # name: (field leg type, which comps, antifield leg type, which)
    'AsP,AP': ('A', 'P', 'As', 'P'),
    'AsV,AV': ('A', 'V', 'As', 'V'),
    'ls,l': ('l', None, 'ls', None),
    'cs,c': ('c', None, 'cs', None),
}


def pi_int(sector, p, a, b, degree=6, dependent='antifield'):
    """internal-insertion bubbles: canonical alpha, beta^mu (dicts over BASIS), in units.
    dependent = 'antifield': legs field (pol a at momentum p, held fixed) and antifield (pol b at -p - r);
    dependent = 'field'    : antifield (pol b at momentum p, held fixed) and field (pol a at -p - r)."""
    tf, wf, ta, wa = SECT[sector]
    p = np.asarray(p, dtype=int); p2 = int(p @ p)
    if dependent == 'antifield':
        fixed = (tf, a); dep = (ta, b)
    else:
        fixed = (ta, b); dep = (tf, a)
    X = DC.Bg([(fixed[0], fixed[1], p)])
    scale = 4
    def ins(qf, rmom):
        Y = DC.Bg([(dep[0], dep[1], -p - rmom)])
        return [qf, X, Y]
    # canonical odd order: (xi, r of the odd external leg); S_M = insertion 0, fixed leg = 1, dependent leg = 2
    oddleg = 1 if C.PARITY[fixed[0]] else 2
    canon = [0, oddleg]
    G0, n0 = cycle_value(ins(DC.QF('s0', np.zeros(4, dtype=int)), np.zeros(4, dtype=int)), canon, scale, degree, p2,
                         extra_sign=-1)
    Gmu = []
    dG = {}
    for mu in range(4):
        e = DC.E4[mu]
        g, _ = cycle_value(ins(DC.QF('smu', np.zeros(4, dtype=int), mu=mu), np.zeros(4, dtype=int)), canon, scale,
                           degree, p2, extra_sign=-1)
        Gmu.append(g)
        def build(t, mu=mu, e=e):
            r = t * e
            return ins(DC.QF('smu', r, mu=mu), r)
        d = cycle_derivative(build, canon, scale, degree, p2)
        acc_add(dG, {k: v * Q2(-1) for k, v in d.items()})     # extra_sign -1 (xi S_M Hessian)
    alpha = dict(G0)
    acc_add(alpha, dG, Q2(0, -1))                               # alpha = G_0 - i sum_mu dG_mu
    beta = [{k: v * Q2(0, -1) for k, v in g.items()} for g in Gmu]   # beta^mu = -i G_mu(0)
    cu = Q2(Fraction(1, C.C_ADJ))
    alpha = {k: v * cu for k, v in alpha.items() if not v.iszero()}
    beta = [{k: v * cu for k, v in bm.items() if not v.iszero()} for bm in beta]
    return alpha, beta, n0


def self_energy(legs, p, degree=4):
    """1PI one-loop two-point function W_1^(2) (plane-wave value, coefficient of odd parameters in leg order) on
    legs [(t1, pol1) at p, (t2, pol2) at -p]: cubic bubbles + quartic seagull; units"""
    p = np.asarray(p, dtype=int); p2 = int(p @ p)
    (t1, a1), (t2, a2) = legs
    X = DC.Bg([(t1, a1, p)]); Y = DC.Bg([(t2, a2, -p)])
    odd = [i for i, (t, _) in enumerate(legs) if C.PARITY[t]]
    tot, n2 = cycle_value([X, Y], odd, 4, degree, p2)
    if t1 == 'A' and t2 == 'A':
        XY = DC.Bg([(t1, a1, p), (t2, a2, -p)])
        s1, n1 = cycle_value([XY], [], 8, degree, p2)
        acc_add(tot, s1)
    cu = Q2(Fraction(1, C.C_ADJ))
    return {k: v * cu for k, v in tot.items()}
