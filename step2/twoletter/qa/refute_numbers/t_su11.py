"""Fermion test (new, not run by either engine): su(1|1) sector {Z, psi}.  Beisert: H_12 = 1 - Pi (graded permutation),
so with lc^natural = -2u H and n_F(A)=2, n_F(lambda)=4 the planar pair kernel must be
   K(Z,psi) -> alpha u (Z psi) + beta u (psi Z),  alpha + beta = 3 (protection of tr Z^{L-1} psi),  beta = 2 (swap).
Computed here from scratch: own Clifford realisation, own Wick/Koszul signs, own Schwinger integration."""
import json, itertools
import numpy as np, sympy as sp
from ref_core import *
from ref_cliff import SIG, SIGB
cl = Cliff(SIG, SIGB)
# exact Clifford check
ok = all((cl.sigb[M] * cl.sig[N] + cl.sigb[N] * cl.sig[M] - (2 if M == N else 0) * sp.eye(16)).is_zero_matrix
         for M in range(10) for N in range(10))
print('exact Clifford relations:', ok)
# rotation generators: vectors L_MN, spinors S_MN = (1/2) sigb_M sig_N on Sigma; letter charges via transpose
def Lvec(M, N):
    m = sp.zeros(10, 10); m[M, N] = 1; m[N, M] = -1; return m
def Ssp(M, N):
    return sp.Rational(1, 2) * cl.sigb[M] * cl.sig[N]
# covariance check rho(L v) = [S, rho(v)] on Sigma -> Sigma^*: sig(Lv) = S'^... (check on Sigma: sigb_K -> ...)
okc = True
for (M, N) in [(0, 1), (4, 5), (8, 9), (2, 7)]:
    S = Ssp(M, N)
    Sstar = -S.T          # action on Sigma^* (contragredient)
    for K in range(10):
        Lv = Lvec(M, N)[:, K]
        lhs = sum((Lv[J] * cl.sig[J] for J in range(10)), sp.zeros(16, 16))
        rhs = Sstar * cl.sig[K] - cl.sig[K] * S
        okc &= (lhs - rhs).is_zero_matrix
print('spinor covariance:', okc)
z = [0] * 8 + [1, sp.I]
Lz = [sp.simplify(sum(Lvec(8, 9).T[i, j] * z[j] for j in range(10)) / z[i]) for i in (8, 9)]
print('Z letter charge in plane 89 (L^T z = q z):', Lz)
qZ = Lz[0]
planes = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]
Ts = [Ssp(*pl).T for pl in planes]
found = []
for l1, l2 in itertools.product([sp.I / 2, -sp.I / 2], repeat=2):
    target = [l1, l2, qZ / 2, qZ / 2, qZ / 2]
    A = sp.Matrix.vstack(*[Ts[k] - target[k] * sp.eye(16) for k in range(5)])
    ns = A.nullspace()
    if ns:
        found.append((target, ns))
for t, ns in found:
    print('psi letter charges', t, 'dim', len(ns))
zeta = list(found[0][1][0])
print('psi letter zeta =', zeta)
ysym = list(sp.symbols('y0:10')); xsym = list(sp.symbols('x0:16'))
r1 = [1, 2, -1, 3]; r2 = [2, -1, 3, 1]
def kern(ta, za, tb, zb, chains):
    tot = {}
    for verts, rs in chains:
        res = kernel_chain(ta, za, tb, zb, verts, rs, cl=cl)
        for k, v in res.items():
            tot[k] = sp.expand(tot.get(k, 0) + v)
    return tot
def project(expr):
    """expr bilinear in y (A output coords) and x (lambda output coords) -> coefficient c with expr = c (z.y)(zeta.x)
    as functionals: output letter Z = sum z_M A_M, psi = sum zeta_b lambda^b; check rank one."""
    P = sp.Poly(expr, *(ysym + xsym))
    km = sp.zeros(10, 16)
    for mon, c in P.terms():
        iy = [i for i in range(10) if mon[i]]; ix = [i for i in range(16) if mon[10 + i]]
        assert len(iy) == 1 and len(ix) == 1
        km[iy[0], ix[0]] += c
    zz = sp.Matrix(z) * sp.Matrix(zeta).T
    # c = <km, zz> / <zz, zz>
    num = sum(km[i, j] * sp.conjugate(zz[i, j]) for i in range(10) for j in range(16))
    den = sum(zz[i, j] * sp.conjugate(zz[i, j]) for i in range(10) for j in range(16))
    c = sp.simplify(num / den)
    resid = (km - c * zz).applyfunc(sp.simplify)
    return c, resid.is_zero_matrix, km
out = {}
# K(Z, psi): letter a = Z (A slot, selector z), letter b = psi (lambda slot, selector zeta)
ch_ZA_order = [([('AAA', ysym), ('lAl', 'A', 'l', 'l', xsym)], [r1, r2])]          # outputs (A, lambda)
ch_swap = [([('lAl', 'A', 'l', 'l', xsym), ('lAl', 'l', 'l', 'A', ysym)], [r1, r2])]  # outputs (lambda, A)
kA = kern('A', z, 'l', zeta, ch_ZA_order)
kS = kern('A', z, 'l', zeta, ch_swap)
print('K(Z,psi) outputs (A,lambda): singular keys', list(kA.keys()))
print('K(Z,psi) outputs (lambda,A): singular keys', list(kS.keys()))
cA, okA, _ = project(kA['log']); cS, okS, _ = project(kS['log'])
print('log coefficient, units hbar/16pi^2:  (Z psi):', cA, 'rank-one', okA, '  (psi Z):', cS, 'rank-one', okS)
alpha, beta = cA / 2, cS / 2
print('planar (units u): alpha =', alpha, ' beta =', beta, ' alpha+beta =', sp.simplify(alpha + beta))
# reversed pair K(psi, Z)
kA2 = kern('l', zeta, 'A', z, [([('lAl', 'l', 'A', 'l', xsym), ('AAA', ysym)], [r1, r2])])          # outputs (lambda, A)
kS2 = kern('l', zeta, 'A', z, [([('lAl', 'l', 'l', 'A', ysym), ('lAl', 'l', 'A', 'l', xsym)], [r1, r2])])  # (A, lambda)
cA2, okA2, _ = project(kA2['log']); cS2, okS2, _ = project(kS2['log'])
print('K(psi,Z): (psi Z) identity:', cA2 / 2, ' (Z psi) swap:', cS2 / 2, ' [units u]', okA2, okS2)
out = dict(zeta=[str(x) for x in zeta], z=[str(x) for x in z],
           K_Z_psi=dict(identity_units_hbar16pi2=str(cA), swap_units_hbar16pi2=str(cS), alpha_u=str(alpha), beta_u=str(beta)),
           K_psi_Z=dict(identity_u=str(cA2 / 2), swap_u=str(cS2 / 2)),
           eps_parts={'(A,lambda)': {str(k): str(v) for k, v in kA.items() if k != 'log'},
                      '(lambda,A)': {str(k): str(v) for k, v in kS.items() if k != 'log'}})
json.dump(out, open('out_su11.json', 'w'), indent=1)
