"""
firstorder.py -- the gauge sector of the first-order formulation (abelian, momentum space):
fields c (deg -1, w 0), A|_V (deg 0, w 1), B self-dual two-form (deg 0, w 2, 3 components);
antifields B^* (deg 1, w 2), A^*|_V (deg 1, w 3), c^* (deg 2, w 4).
Free action  int ( (B, d_+A) + (B,B) + (dc, A^*) ):  eliminating B gives -1/4 int (dA,dA) + topological term.
Self-dual components: D^a(A) = eta^a_{mn} d_m A_n with the 't Hooft symbols eta^a_{mn}.
"""
import sympy as sp
from freecx import p, P2, I_

# layout: c 0 | A 1-4 | B 5-7 | Bs 8-10 | As 11-14 | cs 15
DIM = 16
SL1 = {'c': (0, 1, 1, -1, 0), 'A': (1, 4, 0, 0, 1), 'B': (5, 3, 0, 0, 2),
       'Bs': (8, 3, 1, 1, 2), 'As': (11, 4, 1, 1, 3), 'cs': (15, 1, 0, 2, 4)}


def r(name):
    s = SL1[name]
    return range(s[0], s[0] + s[1])


def eta():
    """'t Hooft symbols eta^a_{mn}, a = 0..2, m,n = 0..3 (labels 1..4): self-dual."""
    E = []
    for a in range(3):
        M = sp.zeros(4, 4)
        for b in range(3):
            for c in range(3):
                M[b, c] = sp.LeviCivita(a, b, c)
            M[b, 3] = 1 if a == b else 0
            M[3, b] = -1 if a == b else 0
        E.append(M)
    return E


def Dmat(q):
    """momentum symbol of D : A -> (eta^a_{mn} i q_m A_n)_a, a 3x4 matrix."""
    E = eta()
    return sp.Matrix(3, 4, lambda a, n: sum(E[a][m, n] * I_ * q[m] for m in range(4)))


def parity():
    out = []
    for k in ['c', 'A', 'B', 'Bs', 'As', 'cs']:
        out += [SL1[k][2]] * SL1[k][1]
    return out


def Omega1():
    O = sp.zeros(DIM, DIM)
    for a, b in [('A', 'As'), ('B', 'Bs')]:
        for k, l in zip(r(a), r(b)):
            O[k, l] = 1
            O[l, k] = -1
    O[15, 0] = 1
    O[0, 15] = -1
    return O


def Q1(q=p):
    """Linear differential of the first-order free action (pairing-preserving, square zero)."""
    Q = sp.zeros(DIM, DIM)
    D = Dmat(q)
    for m, k in enumerate(r('A')):
        Q[k, 0] = -I_ * q[m]                      # c -> -dc
    for a, k in enumerate(r('Bs')):
        for n, l in enumerate(r('A')):
            Q[k, l] = D[a, n]                     # A -> d_+ A (slot B^*)
    for a, (k, l) in enumerate(zip(r('Bs'), r('B'))):
        Q[k, l] = 2                               # B -> 2 B (slot B^*)
    for n, k in enumerate(r('As')):               # B -> transpose of d_+ (slot A^*), fixed by the pairing
        for a, l in enumerate(r('B')):
            Q[k, l] = sp.Symbol('u_BA') * D[a, n]
    for m, l in enumerate(r('As')):
        Q[15, l] = -I_ * q[m]                     # A^* -> d^* A^*
    return Q


def neg(q):
    return tuple(-t for t in q)


def ansatz1(q=p):
    """General scale-free Q^GF (degree -1, order <= 1, Lorentz covariant, Weyl weight 2)."""
    G = sp.zeros(DIM, DIM)
    g = sp.symbols('h1:6')
    D = Dmat(q)
    for m, l in enumerate(r('A')):
        G[0, l] = g[0] * I_ * q[m]                # A -> c
    for n, k in enumerate(r('A')):                # B^* -> A  (one derivative, self-dual -> vector)
        for a, l in enumerate(r('Bs')):
            G[k, l] = g[1] * D[a, n]
    for k, l in zip(r('A'), r('As')):
        G[k, l] = g[2]                            # A^* -> A
    for a, k in enumerate(r('B')):                # A^* -> B (one derivative, vector -> self-dual)
        for n, l in enumerate(r('As')):
            G[k, l] = g[3] * D[a, n]
    for m, k in enumerate(r('As')):
        G[k, 15] = g[4] * I_ * q[m]               # c^* -> A^*
    return G, list(g)


def check_first_order():
    out = {}
    Pm = sp.diag(*[(-1) ** a for a in parity()])
    O = Omega1()
    u = sp.Symbol('u_BA')
    Q = Q1(p)
    Qm = Q1(neg(p))
    X = (Q.T * O + Pm * O * Qm).applyfunc(sp.expand)
    su = sp.solve([x for x in X if x != 0], u, dict=True)
    out['pairing_fixes_u_BA'] = str(su)
    Q = Q.subs(su[0])
    Qm = Qm.subs(su[0])
    out['Q1^2=0'] = all(sp.expand(x) == 0 for x in Q * Q)
    G, g = ansatz1(p)
    Gm, _ = ansatz1(neg(p))
    H = (Q * G + G * Q).applyfunc(sp.expand) - P2 * sp.eye(DIM)
    A = (G.T * O - Pm * O * Gm).applyfunc(sp.expand)
    E = []
    for M in (H, A):
        for x in M:
            if x != 0:
                E += sp.Poly(x, *p).coeffs()
    out['scale_free_solutions_H=p2_selfadjoint'] = str(sp.solve(E, g, dict=True))
    # H = p^2 alone (drop adjointness) still has no solution with G^2 = 0:
    E2 = []
    for x in H:
        if x != 0:
            E2 += sp.Poly(x, *p).coeffs()
    sol = sp.solve(E2, g, dict=True)
    res = []
    for s in sol:
        Gs = G.subs(s)
        G2 = (Gs * Gs).applyfunc(sp.expand)
        res.append(all(x == 0 for x in G2))
    out['H=p2_solutions'] = str(sol)
    out['H=p2_solutions_with_QGF^2=0'] = res
    # Hermitian adjoint (Costello [C07, 5.2]) with an L^2 metric of weights kappa_X on the slots, real p
    kB, kBs, kA, kAs, kc, kcs = sp.symbols('k_B k_Bs k_A k_As k_c k_cs', positive=True)
    w = [kc] + [kA] * 4 + [kB] * 3 + [kBs] * 3 + [kAs] * 4 + [kcs]
    h = sp.diag(*w)
    Qc = Q.subs({t: t for t in p})
    Qdag = h.inv() * Qc.H * h
    Hc = (Q * Qdag + Qdag * Q).applyfunc(sp.expand)
    out['costello_adjoint_(Qdag)^2=0'] = all(sp.expand(x) == 0 for x in Qdag * Qdag)
    out['costello_adjoint_H_on_B_slot_entry'] = str(sp.factor(Hc[5, 5]))
    out['costello_adjoint_H_offdiagonal_B_to_A'] = str(sp.factor(Hc[1, 5]))
    return out


if __name__ == '__main__':
    for k, v in check_first_order().items():
        print(k, ':', v)
