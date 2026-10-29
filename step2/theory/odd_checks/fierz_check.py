# Euclidean Fierz identities used in odd_components_design.tex (dictionary conventions, def:wick):
#   sigma_E^m = (s1,s2,s3, i 1), sigmabar_E^{m adot a} = (-s1,-s2,-s3, i 1), eps^{12} = +1, eps_{12} = -1.
# (1) sum_m sigma^m_{a adot} sigma^m_{b bdot} = -2 eps_{ab} eps_{adot bdot}
# (2) sum_m sigmabar^{m adot a} sigmabar^{m bdot b} = -2 eps^{ab} eps^{adot bdot}
# (3) with Lam_adot Grassmann and commuting u, u' (u'^a u_a = 1):
#     sum_m (i Lam_adot sbar^{m adot a} u_a)(i Lam_bdot sbar^{m bdot b} u'_b) = 2 (Lam_adot Lam^adot)  [F-bar identity]
import numpy as np, itertools
s = [np.array([[0,1],[1,0]],complex), np.array([[0,-1j],[1j,0]]), np.array([[1,0],[0,-1]],complex)]
sig = s + [1j*np.eye(2)]
sbar = [-x for x in s] + [1j*np.eye(2)]
epsU = np.array([[0,1],[-1,0]],complex)   # eps^{ab}
epsL = -epsU                               # eps_{ab}
f1 = sum(np.einsum('ab,cd->acbd', m, m) for m in sig)       # [a,b,adot,bdot] ordering: sigma_{a adot} sigma_{b bdot} -> idx (a,b,adot,bdot)
r1 = np.abs(sum(np.einsum('ax,by->abxy', m, m) for m in sig) + 2*np.einsum('ab,xy->abxy', epsL, epsL)).max()
r2 = np.abs(sum(np.einsum('xa,yb->abxy', m, m) for m in sbar) + 2*np.einsum('ab,xy->abxy', epsU, epsU)).max()
# (3): coefficient tensor T^{adot bdot} = sum_m (i sbar^{m adot a} u_a)(i sbar^{m bdot b} u'_b); Lam_adot Lam_bdot antisymmetric
rng = np.random.default_rng(0)
u = rng.normal(size=2)+1j*rng.normal(size=2); up = rng.normal(size=2)+1j*rng.normal(size=2)
uU = epsU@u; upU = epsU@up                  # upper-index spinors u^a = eps^{ab} u_b
norm = upU@u                               # u'^a u_a
up = up/norm; upU = epsU@up
T = sum(np.outer(1j*(m@u), 1j*(m@up)) for m in sbar)   # T^{adot bdot}
# Lam_adot Lam_bdot = (1/2) eps_{adot bdot}?? use: Lam_adot Lam_bdot = -(1/2) eps_{adot bdot} (Lam Lam), with (Lam Lam) = Lam_adot eps^{adot bdot} Lam_bdot
# contraction T^{xy} Lam_x Lam_y = T^{xy} * (-(1/2) eps_{xy}) (LamLam)
coef = (T*(-0.5*epsL)).sum()
print('Fierz sigma sigma residual', r1, ' sigmabar sigmabar residual', r2)
print('coefficient of (Lam_adot Lam^adot) in sum_m (i Lam sbar u)(i Lam sbar u_prime) with u_prime^a u_a = 1 :', np.round(coef,12))
# check of Lam_x Lam_y = -(1/2) eps_{xy} (Lam Lam) for Grassmann Lam: (LamLam) = eps^{xy} Lam_x Lam_y = eps^{12}Lam_1Lam_2+eps^{21}Lam_2Lam_1 = 2 Lam_1 Lam_2
# so Lam_1 Lam_2 = (LamLam)/2 and eps_{12} = -1 -> -(1/2) eps_{12} = 1/2 : consistent.
