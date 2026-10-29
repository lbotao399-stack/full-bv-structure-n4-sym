"""Finite-dimensional BV model for the audit of D1 (a), (b), (c).

Fibre slots (degree of the field slot, total parity, Weyl weight), as in the manuscript
(DESIGN.md: c in -1, A, lambda in 0, A*, lambda* in +1, c* in +2; lambda in Pi Sigma, so lambda odd,
lambda* even; Weyl weights 0, 1, 3/2, 3, 5/2, 4 of def:c6-sigma):
  c  (-1, odd,  0)    AL, AT (0, even, 1)    l1, l2 (0, odd, 3/2)
  ALs, ATs (1, odd, 3)   ls1, ls2 (1, even, 5/2)   cs (2, even, 4)
AL/AT, ALs/ATs mimic the longitudinal/transverse parts of A|_V, A*|_V.
Each slot is tensored with R[s]/(s^{N+1}) (s plays the role of |p|): basis e_{slot,a} ~ slot (x) s^a.
Pairing (through K_0): slot (x) s^a with partner (x) s^{N-a}; multiplication by s is self-adjoint.

Arrows (momentum space shape of def:c3-qgf and prop:c1-Q):
  Q   : c -> AL (s),  AT -> ATs (-1/2 s^2),  ALs -> cs (-s),  l1 -> ls2 (s),  l2 -> ls1 (-s)
  QGF : ALs -> AL (gF, the Feynman arrow; gF = -2 Feynman, 0 'Landau-type' variant),
        ATs -> AT (-2),  AL -> c (s),  cs -> ALs (-s),  ls2 -> l1 (s),  ls1 -> l2 (-s)
so that H = [Q,QGF] = s^2 on every slot, (Q x 1 + 1 x Q)K_0 = 0, (QGF x 1 - 1 x QGF)K_0 = 0
(prop:c3-adjoint) are checked below, and (QGF)^2 : ALs -> c, cs -> AL is nonzero iff gF != 0.
M (dilatation, v = 2E): on slot (x) s^a the number 2w - 2a + N - 4 (the constant makes M preserve the
pairing; it plays the role of the measure d^4p); then [M,Q] = 0, [M,QGF] = -4 QGF, [M,H] = -4H.
"""
from fractions import Fraction as Fr
from sp import (Alg, mat_mul, mat_add, mat_scale, mat_id, mat_comm, mat_maxabs, mat_inv,
                tensor_X1, tensor_1X, tensor_act, graded_swap)

SLOTS = [  # name, degree, parity, weight, partner
    ('c', -1, 1, Fr(0), 'cs'),
    ('AL', 0, 0, Fr(1), 'ALs'),
    ('AT', 0, 0, Fr(1), 'ATs'),
    ('l1', 0, 1, Fr(3, 2), 'ls1'),
    ('l2', 0, 1, Fr(3, 2), 'ls2'),
    ('ALs', 1, 1, Fr(3), 'AL'),
    ('ATs', 1, 1, Fr(3), 'AT'),
    ('ls1', 1, 0, Fr(5, 2), 'l1'),
    ('ls2', 1, 0, Fr(5, 2), 'l2'),
    ('cs', 2, 0, Fr(4), 'c'),
]

def build(N=4, gF=Fr(-2), eps=Fr(1, 3), T=Fr(7, 4), slots=SLOTS):
    names = [s[0] for s in slots]
    idx = {}
    par, deg_field, wt = [], [], []
    for (nm, dg, p, w, _) in slots:
        for a in range(N + 1):
            idx[(nm, a)] = len(par)
            par.append(p); deg_field.append(dg); wt.append(w)
    n = len(par)
    alg = Alg(par, deg=[-d for d in deg_field])  # functional degree of x^i is -deg(e_i)

    def arrow(table):
        X = {}
        for (src, dst, coef, r) in table:
            for a in range(N + 1 - r):
                X[(idx[(dst, a + r)], idx[(src, a)])] = Fr(coef)
        return X

    Q = arrow([('c', 'AL', 1, 1), ('AT', 'ATs', Fr(-1, 2), 2), ('ALs', 'cs', -1, 1),
               ('l1', 'ls2', 1, 1), ('l2', 'ls1', -1, 1)])
    tabG = [('ATs', 'AT', -2, 0), ('AL', 'c', 1, 1), ('cs', 'ALs', -1, 1),
            ('ls2', 'l1', 1, 1), ('ls1', 'l2', -1, 1)]
    if gF != 0:
        tabG.append(('ALs', 'AL', gF, 0))
    G = arrow(tabG)
    S = arrow([(nm, nm, 1, 1) for nm in names])          # multiplication by s
    Hs = mat_mul(S, S)
    K0 = {}
    for (nm, dg, p, w, prt) in slots:
        for a in range(N + 1):
            K0[(idx[(nm, a)], idx[(prt, N - a)])] = Fr(1)
    Mdiag = {}
    for (nm, dg, p, w, _) in slots:
        for a in range(N + 1):
            Mdiag[(idx[(nm, a)], idx[(nm, a)])] = 2 * w - 2 * a + N - 4
    M = {k: v for k, v in Mdiag.items() if v != 0}

    def poly_of_H(coeffs):
        """sum_k coeffs[k] H^k (exact; H^k = 0 for 2k > N)"""
        out = {}; Hk = mat_id(n)
        for k, c in enumerate(coeffs):
            if 2 * k > N:
                break
            out = mat_add(out, (c, Hk))
            Hk = mat_mul(Hk, Hs)
        return out

    from math import factorial
    def expH(t):   # e^{-tH}
        return poly_of_H([Fr((-t) ** k) / factorial(k) for k in range(N // 2 + 1)])
    def fH(e, Tt):  # int_e^Tt e^{-tH} dt = sum_k (-H)^k (Tt^{k+1} - e^{k+1})/(k+1)!
        return poly_of_H([Fr((-1) ** k) * (Fr(Tt) ** (k + 1) - Fr(e) ** (k + 1)) / factorial(k + 1)
                          for k in range(N // 2 + 1)])
    def ker(A):    # kernel (A (x) 1) K_0
        return mat_mul(A, K0)

    m = dict(N=N, n=n, par=par, deg_field=deg_field, wt=wt, idx=idx, alg=alg, Q=Q, G=G, H=Hs,
             K0=K0, M=M, S=S, eps=Fr(eps), T=Fr(T), gF=gF, expH=expH, fH=fH, ker=ker,
             poly_of_H=poly_of_H, names=names, slots=slots)
    m['K0inv'] = mat_inv(K0, n)
    m['Kt'] = lambda t: ker(expH(t))
    m['Pstar'] = lambda e, Tt: mat_mul(G, fH(e, Tt))            # operator int_e^T QGF e^{-tH} dt
    m['P'] = lambda e, Tt: ker(m['Pstar'](e, Tt))                # (Pstar (x) 1) K_0 = int (QGF x 1) K_t dt
    m['L'] = lambda t: ker(mat_mul(G, expH(t)))                 # L_t = (QGF x 1) K_t
    return m

def kernel_checks(m):
    """kernel-level identities of ch03 / ch07 in the model (all must be 0)"""
    par, n = m['par'], m['n']
    Q, G, H, K0, M = m['Q'], m['G'], m['H'], m['K0'], m['M']
    e, T = m['eps'], m['T']
    out = {}
    out['Q^2'] = mat_maxabs(mat_mul(Q, Q))
    out['[Q,QGF]-H'] = mat_maxabs(mat_add(mat_comm(Q, G, 1, 1), (-1, H)))
    out['[QGF,H]'] = mat_maxabs(mat_comm(G, H, 1, 0))
    out['K0 graded symmetric'] = mat_maxabs(mat_add(graded_swap(K0, par), (-1, K0)))
    out['(Qx1+1xQ)K0'] = mat_maxabs(tensor_act(Q, 1, K0, par))
    out['(QGFx1-1xQGF)K0'] = mat_maxabs(mat_add(tensor_X1(G, K0), (-1, tensor_1X(G, 1, K0, par))))
    out['(Mx1+1xM)K0'] = mat_maxabs(tensor_act(M, 0, K0, par))
    out['(Hx1-1xH)K0'] = mat_maxabs(mat_add(tensor_X1(H, K0), (-1, tensor_1X(H, 0, K0, par))))
    out['[M,Q]'] = mat_maxabs(mat_comm(M, Q, 0, 1))
    out['[M,QGF]+4QGF'] = mat_maxabs(mat_add(mat_comm(M, G, 0, 1), (4, G)))
    out['[M,H]+4H'] = mat_maxabs(mat_add(mat_comm(M, H, 0, 0), (4, H)))
    P = m['P'](e, T)
    out['P graded symmetric'] = mat_maxabs(mat_add(graded_swap(P, par), (-1, P)))
    out['Kt graded symmetric'] = mat_maxabs(mat_add(graded_swap(m['Kt'](e), par), (-1, m['Kt'](e))))
    out['(Qx1+1xQ)P - (K_eps - K_T)'] = mat_maxabs(mat_add(tensor_act(Q, 1, P, par), (-1, m['Kt'](e)), m['Kt'](T)))
    # (Mx1+1xM)P = -4(eps d_eps + T d_T)P,  d_eps P = -L_eps, d_T P = L_T
    rhs = mat_add((4 * e, m['L'](e)), (-4 * T, m['L'](T)))
    out['(Mx1+1xM)P + 4(eps d_eps + T d_T)P'] = mat_maxabs(mat_add(tensor_act(M, 0, P, par), (-1, rhs)))
    G2 = mat_mul(G, G)
    out['|(QGF)^2|'] = mat_maxabs(G2)
    # components of (QGF)^2 : which slots
    comps = sorted({(m['slots'][j // (m['N'] + 1)][0], m['slots'][i // (m['N'] + 1)][0]) for (i, j) in G2})
    out['(QGF)^2 arrows (src->dst)'] = comps
    return out
