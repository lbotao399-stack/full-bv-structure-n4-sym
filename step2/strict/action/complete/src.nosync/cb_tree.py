"""Route B tree-level operations on two-leg plane-wave densities (eb_qb conventions):
  Q (eb_qb.apply_Q, engine B gate G-sign),
  the bracket {G, F} for a quadratic functional G (density g) as the derivation of F that replaces a leg of the result
  by M_G(p) leg, M_G(p) = -K_0 H_G(p,-p)^T (H_G the Hessian, row leg momentum p), with the Koszul rule of apply_Q for even
  G; for odd G, M_G = Pi(-K_0 H_G^T) and no Koszul sign.  The formula is calibrated (gate GB-tree) by: M_{S_free} = q (Prop c1-Q), M_{T1} = engine B's
  validated slot map of {T1,-}, {N, S_free} = {S_free, N} = Q N for the odd field-counting functional N, and the derivation rule
  Q{N,F} = -{QN,F} - (-1)^{|N|}{N,QF} (bracket of Def c1-bracket) for an odd F (t_tree*.py, gate GB-tree)."""
import itertools
import numpy as np
import cb_common as B
E, C, V, Q = B.E, B.C, B.V, B.Q
K0 = np.array(C.k0_matrix().tolist(), dtype=complex)
SPEC_OF = C.SPEC_OF
PARV = np.array(C.PARITY)
SPECIES = ['A', 'lam', 'c', 'As', 'lams', 'cs']


def unit_leg(J, p):
    s = SPEC_OF[J]
    u = np.zeros(V.DIM[s], complex)
    u[J - E.OFF[s]] = 1
    return dict(sp=s, u=u, p=list(p))


def hessian(g, kr, kc):
    """H_JK = (-1)^{pJ pK} g(leg_J(kr), leg_K(kc))  (vertices.hessian_blocks rule, no background)."""
    H = np.zeros((54, 54), complex)
    for J in range(54):
        lJ = unit_leg(J, kr)
        for K in range(54):
            v = g([lJ, unit_leg(K, kc)])
            if v != 0:
                H[J, K] = (-1) ** (PARV[J] * PARV[K]) * v
    return H


PI = np.diag([(-1) ** int(x) for x in PARV]).astype(complex)


def Mmat(g, p, gpar=0):
    mp = [-x for x in p]
    M = -K0 @ hessian(g, p, mp).T
    return M if gpar == 0 else PI @ M


def split_leg(vec, p):
    out = []
    for s in SPECIES:
        u = vec[E.OFF[s]:E.OFF[s] + V.DIM[s]]
        if np.any(u != 0):
            out.append(dict(sp=s, u=u.copy(), p=p))
    return out


def leg_vec(leg):
    v = np.zeros(54, complex)
    v[E.OFF[leg['sp']]:E.OFF[leg['sp']] + V.DIM[leg['sp']]] = leg['u']
    return v


def bracket(g, gpar, f, cache=None):
    """density of {G, F}: G quadratic with density g and parity gpar, F with density f."""
    cache = {} if cache is None else cache

    def h(legs):
        tot = 0j
        for j, leg in enumerate(legs):
            key = tuple(float(x) for x in leg['p'])
            if key not in cache:
                cache[key] = Mmat(g, leg['p'], gpar)
            img = cache[key] @ leg_vec(leg)
            nb = sum(V.PAR[l['sp']] for l in legs[:j])
            sg = (-1) ** (nb * (1 - gpar))
            for nl in split_leg(img, leg['p']):
                new = list(legs); new[j] = nl
                v = f(new)
                if v != 0:
                    tot += sg * v
        return tot
    return h
