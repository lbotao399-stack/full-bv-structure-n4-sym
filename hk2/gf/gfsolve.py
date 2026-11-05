"""
gfsolve.py -- the non-minimal complexes, the canonical transformation of the gauge fermion,
and the general (scale-free, local, order <= 1, Lorentz covariant) ansatz for a gauge-fixing operator.
"""
import sympy as sp
import freecx as F
from freecx import p, P2, I_, rng

n1, n2, al, kap, kap2 = sp.symbols('n1 n2 alpha kappa kappa2')


def idx(name):
    return rng(name)[0]


# ---------------------------------------------------------------- reduced non-minimal pair (cbar, cbar^*)
def Qprime_red(q=p, alpha=al):
    """Costello frame: Q' = Q (+) Q_nm, Q_nm : cbar^* -> alpha cbar^* in the slot of cbar
    (the linear differential of S + (alpha/2) int (cbar^*, cbar^*))."""
    Q = F.Q_min('red', q)
    Q[idx('cb'), idx('cbs')] = alpha
    return Q


def Nmap_red(q=p, a1=n1, a2=n2):
    """N : cbar -> a1 d cbar (slot A^*|_V),  A|_V -> a2 d^* A (slot cbar^*)."""
    n = F.DIMS['red']
    N = sp.zeros(n, n)
    for m, k in enumerate(rng('AsV')):
        N[k, idx('cb')] = a1 * I_ * q[m]
    for m, l in enumerate(rng('AV')):
        N[idx('cbs'), l] = a2 * (-I_ * q[m])
    return N


def symplectic_condition(Phi, Phim, config):
    """Phi(p)^T Omega Phi(-p) = Omega."""
    Om = F.Omega(config)
    return (Phi.T * Om * Phim - Om).applyfunc(sp.expand)


def Q_Psi_red(q=p, alpha=al, a1=n1, a2=n2):
    N = Nmap_red(q, a1, a2)
    n = N.shape[0]
    Phi = sp.eye(n) + N
    Phinv = sp.eye(n) - N
    return (Phinv * Qprime_red(q, alpha) * Phi).applyfunc(sp.expand)


def QGF_Psi_red(q=p, g1=1, g2=-1):
    """Frame Psi gauge-fixing operator: antifields -> fields.
    A^* -> -2 A^* (A slot), lambda^* -> -dslash lambda^*, cbar^* -> g1 cbar^* (c slot), c^* -> g2 c^* (cbar slot)."""
    n = F.DIMS['red']
    G = sp.zeros(n, n)
    Gold = F.QGF_old(q)
    for k, l in zip(rng('AV'), rng('AsV')):
        G[k, l] = -2
    for k, l in zip(rng('AW'), rng('AsW')):
        G[k, l] = -2
    for k in rng('lam'):
        for l in rng('lams'):
            G[k, l] = Gold[k, l]
    G[idx('c'), idx('cbs')] = g1
    G[idx('cb'), idx('cs')] = g2
    return G


# ---------------------------------------------------------------- general ansatz
def ansatz(config, q=p):
    """Most general odd operator of degree -1, order <= 1, Lorentz and SO(6) covariant, homogeneous of
    Weyl weight 2 (no dimensionful constant): every arrow X -> Y with deg Y = deg X - 1 and
    k = w_Y - w_X + 2 derivatives, k in {0,1}, with the covariant structures of that k."""
    n = F.DIMS[config]
    G = sp.zeros(n, n)
    syms = []
    sl = F.slots_of(config)
    # spinor sector: lambda^* -> a * i rho(p) lambda^*  (the unique covariant structure with k = 1)
    Rm, Rs = F.rho_p(q)
    a = sp.Symbol('g_lam')
    syms.append(a)
    for i_, k in enumerate(rng('lam')):
        for j_, l in enumerate(rng('lams')):
            if Rs[i_, j_] != 0:
                G[k, l] = a * I_ * Rs[i_, j_]
    a = sp.Symbol('g_AW')
    syms.append(a)
    for k, l in zip(rng('AW'), rng('AsW')):
        G[k, l] = a
    gauge = [s for s in sl if s[0] in ('c', 'AV', 'AsV', 'cs', 'cb', 'cbs', 'b', 'bs')]
    for X in gauge:
        for Y in gauge:
            if Y[4] != X[4] - 1:
                continue
            k = Y[5] - X[5] + 2
            if k not in (0, 1):
                continue
            dx, dy = X[2], Y[2]
            name = 'g_%s_%s' % (X[0], Y[0])
            if k == 0 and dx == dy == 1:
                a = sp.Symbol(name)
                syms.append(a)
                G[rng(Y[0])[0], rng(X[0])[0]] = a
            elif k == 0 and dx == dy == 4:
                a = sp.Symbol(name)
                syms.append(a)
                for kk, ll in zip(rng(Y[0]), rng(X[0])):
                    G[kk, ll] = a
            elif k == 1 and dx == 1 and dy == 4:
                a = sp.Symbol(name)
                syms.append(a)
                for m, kk in enumerate(rng(Y[0])):
                    G[kk, rng(X[0])[0]] = a * I_ * q[m]
            elif k == 1 and dx == 4 and dy == 1:
                a = sp.Symbol(name)
                syms.append(a)
                for m, ll in enumerate(rng(X[0])):
                    G[rng(Y[0])[0], ll] = a * I_ * q[m]
    return G, syms


def solve_ansatz(Q, Qm, config, extra_eqs=()):
    """Impose H = [Q,G] = |p|^2 and the adjointness of G; return the linear solution set,
    and the residual of G^2 on that set."""
    G, syms = ansatz(config, p)
    Gm, _ = ansatz(config, F.neg(p))
    n = F.DIMS[config]
    H = (Q * G + G * Q).applyfunc(sp.expand) - P2 * sp.eye(n)
    A = F.adj_check_GF(G, Gm, config)
    eqs = []
    for M in (H, A):
        for x in M:
            x = sp.expand(x)
            if x != 0:
                poly = sp.Poly(x, *p)
                eqs += [c for c in poly.coeffs()]
    eqs += list(extra_eqs)
    sol = sp.solve(eqs, syms, dict=True)
    return G, syms, sol, eqs


# ---------------------------------------------------------------- non-minimal pair with b kept
def Qprime_b(q=p, t=None):
    """Costello frame with the trivial pair (cbar, b): S + int (b, cbar^*), Q : b -> b (cbar slot),
    cbar^* -> t cbar^* (b^* slot), t fixed by the pairing."""
    Q = F.Q_min('b', q)
    Q[idx('cb'), idx('b')] = 1
    Q[idx('bs'), idx('cbs')] = sp.Symbol('t_b') if t is None else t
    return Q


def Nmap_b(q=p, a1=n1, a2=n2, k1=kap, k2=kap2):
    n = F.DIMS['b']
    N = sp.zeros(n, n)
    for m, k in enumerate(rng('AsV')):
        N[k, idx('cb')] = a1 * I_ * q[m]
    for m, l in enumerate(rng('AV')):
        N[idx('cbs'), l] = a2 * (-I_ * q[m])
    N[idx('bs'), idx('cb')] = k1
    N[idx('cbs'), idx('b')] = k2
    return N


def Q_Psi_b(q=p, t=1, a1=-1, a2=1, k1=kap, k2=kap2):
    N = Nmap_b(q, a1, a2, k1, k2)
    n = N.shape[0]
    return ((sp.eye(n) - N) * Qprime_b(q, t) * (sp.eye(n) + N)).applyfunc(sp.expand)
