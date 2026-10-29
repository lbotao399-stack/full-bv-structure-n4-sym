"""su(1|1) sector evaluated with ENGINE A's kernel function (its own Clifford basis): does the agreed fermion kernel
reproduce Beisert's H_12 = 1 - Pi?  (psi located in engine A's basis by the same charge conditions.)"""
import sys, json
import numpy as np, sympy as sp
sys.path.insert(0, '../../engine_a/src.nosync')
import ea_kern as EK, ea_core as EC
from cliff import s as Sg, sb as Sgb
def S(M, N):   # (1/2) rho_M rho_N on Sigma = (1/2) s[M] sb[N]
    return 0.5 * Sg[M] @ Sgb[N]
z = np.zeros(10, complex); z[8] = 1; z[9] = 1j
planes = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]
qZ = -1j
res = []
for l1 in (0.5j, -0.5j):
    for l2 in (0.5j, -0.5j):
        target = [l1, l2, qZ / 2, qZ / 2, qZ / 2]
        A = np.vstack([S(*pl).T - t * np.eye(16) for pl, t in zip(planes, target)])
        u, sv, vh = np.linalg.svd(A)
        if sv[-1] < 1e-9:
            res.append((target, vh[-1].conj()))
print('psi solutions in engine A basis:', len(res))
zeta = res[0][1]; zeta = zeta / zeta[np.argmax(abs(zeta))]
r1 = np.array([1, 2, -1, 3]); r2 = np.array([2, -1, 3, 1])
def kval(blocks, ta, tb, za, zb, opos):
    # returns 10x16 matrix over output components (A comp M, lambda comp beta) of the log coefficient
    km = np.zeros((10, 16), complex)
    for M in range(10):
        for b in range(16):
            pols = []
            for kind, mom in blocks:
                if kind == 'A':
                    p = np.zeros(10, int); p[M] = 1; pols.append(('cub', 'AW' if M >= 4 else 'AV', p, mom))
                else:
                    p = np.zeros(16, int); p[b] = 1; pols.append(('cub', 'lamP', p, mom))
            out = EK.kernel(pols, ta, tb)
            if out is None or out.get('log') is None:
                continue
            re, im = out['log']
            K = np.array([[float(re[i, j]) + 1j * float(im[i, j]) for j in range(re.shape[1])] for i in range(re.shape[0])])
            km[M, b] = za @ K @ zb
    return km
def proj(km):
    zz = np.outer(z, zeta)
    c = np.vdot(zz, km) / np.vdot(zz, zz)
    return c, np.abs(km - c * zz).max()
zW = z[4:10]
kA = kval([('A', r1), ('l', r2)], 'AW', 'lamP', zW, zeta, None)
kS = kval([('l', r1), ('A', r2)], 'AW', 'lamP', zW, zeta, None)
cA, eA = proj(kA); cS, eS = proj(kS)
print('engine A, K(Z,psi): identity (Z psi) =', cA, ' resid', eA, ';  swap (psi Z) =', cS, ' resid', eS)
print('planar units u: alpha =', cA / 2, ' beta =', cS / 2)
json.dump(dict(alpha_u=[cA.real / 2, cA.imag / 2], beta_u=[cS.real / 2, cS.imag / 2], resid=[eA, eS]),
          open('out_su11_engA.json', 'w'), indent=1)
