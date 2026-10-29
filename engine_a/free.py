"""Free BV data in a two-mode model (momenta +p, -p, one colour, Vol stripped) -- engine A.
Derives: Q = {S_free, -} on coordinates, and the propagator P(eps,infty) = (e^{-eps p^2}/p^2) Phat
from  partial_P(Q(y_a y_b)) = partial_{K_eps}(y_a y_b) = -e^{-eps p^2} Delta_0(y_a y_b)
([C07 Lemma 10.2.1] at T = infty) with the DESIGN ansatz for Q^GF (Feynman type)."""
import warnings
import numpy as np
import sympy as sp
from grass import GA, PAR, newvar
from cliff import s, sb, VDIR
warnings.filterwarnings('ignore')
p = sp.symbols('p0:4', real=True)
p2 = sum(x ** 2 for x in p)
I = sp.I

# component list: (name, index, parity)
COMP = [('c', 0, 1)] + [('A', M, 0) for M in range(10)] + [('l', a, 1) for a in range(16)] + \
       [('As', M, 1) for M in range(10)] + [('ls', a, 0) for a in range(16)] + [('cs', 0, 0)]
VAR = {}
n = 1000
for sg in (+1, -1):
    for (nm, i, par) in COMP:
        VAR[(nm, i, sg)] = newvar(n, par); n += 1

def y(nm, i, sg):
    return GA.var(VAR[(nm, i, sg)])

def mom(sg):
    return [sg * x for x in p]

def sympify_mat(m):
    return sp.Matrix(16, 16, lambda a, b: sp.nsimplify(complex(m[a, b]).real) + I * sp.nsimplify(complex(m[a, b]).imag))

S_ = [sympify_mat(m) for m in s]
SB = [sympify_mat(m) for m in sb]

def S_free():
    S = GA()
    # -1/4 (dA,dA)_Hodge, cross terms of the +p and -p modes:  -1/4 * 2 * sum_{M<N} F+_{MN} F-_{MN}
    def F(sg, M, N):
        k = mom(sg)
        t = GA()
        if M < 4:
            t = t + y('A', N, sg) * (I * k[M])
        if N < 4:
            t = t - y('A', M, sg) * (I * k[N])
        return t
    for M in range(10):
        for N in range(M + 1, 10):
            S = S + F(+1, M, N) * F(-1, M, N) * sp.Rational(-1, 2)
    # 1/2 (lambda, dslash lambda) = 1/2 lambda^T sb_mu d_mu lambda
    for (s1, s2) in ((+1, -1), (-1, +1)):
        k2 = mom(s2)
        for mu in VDIR:
            for a in range(16):
                for b in range(16):
                    c = SB[mu][a, b]
                    if c != 0:
                        S = S + y('l', a, s1) * y('l', b, s2) * (sp.Rational(1, 2) * c * I * k2[mu])
    # (dc, A^*) = d_mu c A*_mu
    for (s1, s2) in ((+1, -1), (-1, +1)):
        k1 = mom(s1)
        for mu in VDIR:
            S = S + y('c', 0, s1) * y('As', mu, s2) * (I * k1[mu])
    return S.map(sp.expand)

PAIRS = [(('A', M), ('As', M)) for M in range(10)] + [(('ls', a), ('l', a)) for a in range(16)] + [(('cs', 0), ('c', 0))]

def Delta0(F):
    out = GA()
    for (xn, xi), (qn, qi) in PAIRS:
        for sg in (+1, -1):
            out = out + F.dleft(VAR[(qn, qi, -sg)]).dleft(VAR[(xn, xi, sg)])
    return out

def bracket(f, g, pf):
    return Delta0(f * g) - Delta0(f) * g - (f * Delta0(g)) * ((-1) ** pf)

if __name__ == '__main__':
    import json, itertools
    Sf = S_free()
    # Q on every coordinate
    Qy = {}
    for key, v in VAR.items():
        Qy[key] = bracket(Sf, GA.var(v), 0).map(sp.expand)
    # ansatz for Phat (first factor at +p); unknowns a_AA, b_cAs, b_Asc, g_ll
    aAA, bc, bA, gl, gl2 = sp.symbols('aAA bcAs bAsc gll gll2')
    def Phat(k1, k2, sgn):
        """Phat^{k1 k2} with first factor at momentum sgn*p"""
        (n1, i1), (n2, i2) = k1, k2
        q = mom(sgn)
        if n1 == 'A' and n2 == 'A':
            return aAA if i1 == i2 else 0
        if n1 == 'c' and n2 == 'As':
            return bc * I * q[i2] if i2 < 4 else 0
        if n1 == 'As' and n2 == 'c':
            return bA * I * q[i1] if i1 < 4 else 0
        if n1 == 'l' and n2 == 'l':
            return gl * I * sum(S_[mu][i1, i2] * q[mu] for mu in VDIR) + gl2 * I * sum(SB[mu][i1, i2] * q[mu] for mu in VDIR)
        return 0
    def dPhat(F):
        """partial_Phat on a quadratic GA: sum over monomials y_a y_b (a at +p or -p)"""
        tot = 0
        inv = {v: k for k, v in VAR.items()}
        for key, c in F.t.items():
            if len(key) != 2:
                continue
            ka, kb = inv[key[0]], inv[key[1]]
            if ka[2] == kb[2]:
                continue
            tot += c * Phat((ka[0], ka[1]), (kb[0], kb[1]), ka[2])
        return sp.expand(tot)
    eqs = []
    comps = [(nm, i) for (nm, i, _) in COMP]
    rng = np.random.default_rng(0)
    sample = [(a, b) for a in comps for b in comps]
    for (a, b) in sample:
        ya = GA.var(VAR[(a[0], a[1], +1)]); yb = GA.var(VAR[(b[0], b[1], -1)])
        pa = PAR[VAR[(a[0], a[1], +1)]]
        Qprod = Qy[(a[0], a[1], +1)] * yb + ya * Qy[(b[0], b[1], -1)] * ((-1) ** pa)
        lhs = dPhat(Qprod.map(sp.expand))
        rhs = -p2 * Delta0(ya * yb).t.get((), 0)
        e = sp.expand(lhs - rhs)
        if e != 0:
            eqs.append(e)
    # collect coefficients in p
    polys = []
    for e in eqs:
        polys += sp.Poly(e, *p).coeffs()
    sol = sp.solve(polys, [aAA, bc, bA, gl, gl2], dict=True)
    print('nontrivial equations:', len(eqs))
    print('solution:', sol)
    # graded symmetry check of Phat: Phat^{ba}(-p) = (-1)^{|a||b|} Phat^{ab}(p)
    so = sol[0]
    bad = 0
    for a in comps:
        for b in comps:
            pa = PAR[VAR[(a[0], a[1], 1)]]; pb = PAR[VAR[(b[0], b[1], 1)]]
            x = sp.expand(sp.sympify(Phat(a, b, +1)).subs(so))
            yv = sp.expand(sp.sympify(Phat(b, a, -1)).subs(so))
            if sp.expand(yv - (-1) ** (pa * pb) * x) != 0:
                bad += 1
    print('graded symmetry violations:', bad)
    # Q on fields printed for the record (A at +p -> coefficient structure)
    print('Q(A_0^+) =', {k: v for k, v in list(Qy[('A', 0, 1)].t.items())[:6]})
    print('Q(c^+)   =', Qy[('c', 0, 1)].t)
    print('Q(As_0^+)=', Qy[('As', 0, 1)].t)
    json.dump({str(k): str(v) for k, v in so.items()}, open('logs/propagator_solution.json', 'w'))
