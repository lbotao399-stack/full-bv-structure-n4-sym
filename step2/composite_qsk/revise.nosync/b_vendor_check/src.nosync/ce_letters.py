"""Letters of the families in the engine realisation (code-only bases), as charge eigenstates.

Cartan generators: rotations in the planes (4,5), (6,7), (8,9) of V_10 (R-symmetry, Z = A_4 + i A_5 has charges
(1,0,0)) and (0,1), (2,3) (Lorentz, engine components).  On the lambda slot the generator of the rotation with
vector matrix L (delta A = theta L A) is the spinor matrix S fixed by rho-equivariance
    S rho(e_M) - rho(e_M) S' = sum_N L_{NM} rho(e_N)   (S on Sigma, S' on Sigma^*),
solved numerically-exactly with sympy.  A letter is a linear functional z on a slot; it has charge vector q when
z^T S = i q z^T for every Cartan generator.  Output evaluation vectors: field configurations u in a slot; for the
collinear on-shell kinematics p = t nbar, u is taken in the on-shell space (ker rho(nbar) on lambda, transverse
polarisations on A|_V) with z.u = 1."""
import os, sys, itertools
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
import ce_core as C

I = sp.I
R10 = C.Rm      # rho(e_M): Sigma -> Sigma^*
Rb10 = C.Rbm    # rho(e_M): Sigma^* -> Sigma


def vec_L(a, b):
    L = sp.zeros(10, 10)
    L[b, a] = 1
    L[a, b] = -1          # delta A_a = -theta A_b, delta A_b = theta A_a
    return L


def spinor_gen(a, b):
    """S on Sigma (16x16) and Sp on Sigma^* for the rotation in the plane (a, b): S = c Gamma^a Gamma^b restricted."""
    L = vec_L(a, b)
    # ansatz: S = x * (Rb_a R_b) (Sigma -> Sigma), Sp = x * (R_a Rb_b) (Sigma^* -> Sigma^*)
    x = sp.Symbol('x')
    S = x * Rb10[a] * R10[b]
    Sp = x * R10[a] * Rb10[b]
    eqs = []
    for M in range(10):
        lhs = Sp * R10[M] - R10[M] * S       # acting on Sigma -> Sigma^*: rho(SA) ... covariance in the form below
        rhs = sum((L[N, M] * R10[N] for N in range(10)), sp.zeros(16, 16))
        eqs += list((lhs - rhs).applyfunc(sp.expand))
    sol = sp.solve([e for e in eqs if e != 0], x, dict=True)
    assert sol, 'no spinor generator'
    xv = sol[0][x]
    return (S.subs(x, xv)), (Sp.subs(x, xv))


CARTAN_R = [(4, 5), (6, 7), (8, 9)]
CARTAN_L = [(0, 1), (2, 3)]
_GEN = {}


def gens():
    if not _GEN:
        for pl in CARTAN_R + CARTAN_L:
            _GEN[pl] = spinor_gen(*pl)
    return _GEN


def lam_letters():
    """common left eigenvectors z (16-dim, on the lambda slot) of the five Cartan generators; returns list of
    (z, charges (R1,R2,R3,L1,L2), chirality '+'|'-')."""
    G = gens()
    # simultaneous diagonalisation of the commuting family acting by right multiplication z^T S
    mats = [G[pl][0].T for pl in CARTAN_R + CARTAN_L]
    # generic combination
    comb = sum((sp.Rational(c) * m for c, m in zip((1, 10, 100, 1000, 10000), mats)), sp.zeros(16, 16))
    out = []
    for ev, mult, vecs in comb.eigenvects():
        for v in vecs:
            v = sp.Matrix(v)
            ch = []
            for m in mats:
                w = (m * v).applyfunc(sp.expand)
                # w = i q v
                k = [i for i in range(16) if v[i] != 0][0]
                q = sp.nsimplify(sp.expand(w[k] / v[k] / I))
                assert (w - I * q * v).applyfunc(sp.expand).is_zero_matrix
                ch.append(q)
            chir = '+' if all(v[i] == 0 for i in range(8, 16)) else ('-' if all(v[i] == 0 for i in range(8)) else '?')
            out.append((v, tuple(ch), chir))
    return out


def A_charge(z):
    """charges of a letter z on the A slot (10-vector), z^T L = i q z^T."""
    ch = []
    for pl in CARTAN_R + CARTAN_L:
        L = vec_L(*pl)
        w = (L.T * z).applyfunc(sp.expand)
        nz = [i for i in range(10) if z[i] != 0]
        q = sp.nsimplify(sp.expand(w[nz[0]] / z[nz[0]] / I))
        if not (w - I * q * z).applyfunc(sp.expand).is_zero_matrix:
            return None
        ch.append(q)
    return tuple(ch)


def weyl_kernel(p):
    """on-shell space of lambda for momentum p: ker of rho(p): Sigma -> Sigma^* (16 -> 16)."""
    rho = sum((p[m] * R10[m] for m in range(4)), sp.zeros(16, 16))
    return rho.nullspace()


if __name__ == '__main__':
    Z = sp.Matrix([0, 0, 0, 0, 1, I, 0, 0, 0, 0])
    print('Z charges', A_charge(Z))
    for v, ch, chir in lam_letters():
        print(ch, chir, [i for i in range(16) if v[i] != 0])
