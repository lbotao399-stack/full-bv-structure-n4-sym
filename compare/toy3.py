"""0-dim brute-force test of the one-loop formula
   W_1 = sum_n (1/(2n)) str[(P o (-1)^{|.|} o d2I)^n],   (d2I)_{jl} = d_j (d_l I)  (right slot first),
against hbar log(e^{hbar d_P} e^{I/hbar}) with d_P(phi^i phi^j) = P^{ij} (graded symmetric P), for an interaction
quadratic in the fluctuation with background-linear coefficients, including ODD background parameters.
Independent of engine_a / engine_b."""
import sympy as sp
from grassmann3 import G

x1, x2 = sp.symbols('x1 x2')           # even fluctuations
T1, T2, R1, R2 = 0, 1, 2, 3             # odd generator ids: theta1, theta2 (fluct), r1, r2 (background)
b, tau = sp.symbols('b tau')            # even background, order-counting parameter
FL = [('e', x1), ('e', x2), ('o', T1), ('o', T2)]
PAR = [0, 0, 1, 1]

def d(F, j):
    k, v = FL[j]
    return F.deven(v) if k == 'e' else F.dodd(v)

# graded-symmetric propagator: P^{x x} symmetric, P^{theta theta} antisymmetric
P = sp.zeros(4, 4)
P[0, 0], P[0, 1], P[1, 0], P[1, 1] = 3, 2, 2, -1
P[2, 3], P[3, 2] = 5, -5

def dP(F):
    out = G()
    for k in range(4):
        for l in range(4):
            if P[k, l] != 0:
                out = out + d(d(F, k), l) * sp.Rational(P[k, l], 2)
    return out

r1, r2 = G.gen(R1), G.gen(R2)
th1, th2 = G.gen(T1), G.gen(T2)
X1, X2 = G.c(x1), G.c(x2)
I = (G.c(b * (2 * x1**2 - 3 * x1 * x2 + x2**2)) + th1 * th2 * (7 * b)
     + r1 * (X1 * th1 * 2 + X2 * th2 * (-1) + X1 * th2 * 3)
     + r2 * (X2 * th1 * 5 + X1 * th1 * (-2) + X2 * th2 * 4)) * tau

def at_zero(F):
    out = {}
    for k, v in F.d.items():
        if T1 in k or T2 in k:
            continue
        v = v.subs({x1: 0, x2: 0})
        out[k] = out.get(k, 0) + v
    return G(out)

# brute force: Z = sum_n dP^n(I^n)/(n!)^2  (I quadratic in the fluctuation), hbar = 1
Z = [G.c(1)]
In = G.c(1)
for n in range(1, 4):
    In = In * I
    F = In
    for _ in range(n):
        F = dP(F)
    Z.append(at_zero(F) * sp.Rational(1, sp.factorial(n) ** 2))
z1, z2, z3 = Z[1], Z[2], Z[3]
W = {1: z1, 2: z2 - z1 * z1 * sp.Rational(1, 2), 3: z3 - z1 * z2 + z1 * z1 * z1 * sp.Rational(1, 3)}

def hessian(conv):
    H = [[None] * 4 for _ in range(4)]
    for j in range(4):
        for l in range(4):
            Hjl = d(d(I, l), j) if conv == 'right_first' else d(d(I, j), l)
            H[j][l] = at_zero(Hjl)
    return H

def matmul(A, B):
    n = len(A)
    return [[sum((A[i][k] * B[k][j] for k in range(n)), G()) for j in range(n)] for i in range(n)]

def formula(conv, parity_side):
    H = hessian(conv)
    M = [[G() for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for l in range(4):
            acc = G()
            for j in range(4):
                if P[i, j] != 0:
                    sgn = (-1) ** PAR[j] if parity_side == 'middle' else 1
                    acc = acc + H[j][l] * (P[i, j] * sgn)
            M[i][l] = acc
    out = {}
    Mn = M
    for n in range(1, 4):
        if n > 1:
            Mn = matmul(Mn, M)
        out[n] = sum((Mn[i][i] * (-1) ** PAR[i] for i in range(4)), G()) * sp.Rational(1, 2 * n)
    return out

def order(F, n):
    return G({k: sp.expand(v).coeff(tau, n) for k, v in F.d.items()})

res = {}
for conv in ('right_first', 'left_first'):
    for side in ('middle', 'none'):
        F = formula(conv, side)
        ok = all((order(W[n], n) - order(F[n], n)).iszero() for n in (1, 2, 3))
        res[(conv, side)] = ok
        print(conv, side, 'matches brute force at n=1,2,3:', ok)
print('odd-background trilinear coefficient (r1 r2 b):', order(W[3], 3).coeff((R1, R2)))
assert res[('right_first', 'middle')]
