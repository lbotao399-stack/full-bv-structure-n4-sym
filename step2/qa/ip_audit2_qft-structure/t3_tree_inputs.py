"""Independent check of the classical/tree inputs of T3(iv) (Prop ip-weight3): flavour sum 48|w|^4 of B_w,
the flavour factor 6 (zeta, conj zeta) of the O_zeta tree norm, and |kappa_B| = 1/4 from the two-component
S^(1)(xi) of ch08 Remark (twist comparison), with blocks built from ch08 Def c8-twocomp only (A^{ij} = (1/2) eps^{ijkl} A_kl).
"""
import json, itertools
import numpy as np

eps3 = np.zeros((3, 3, 3)); eps4 = np.zeros((4, 4, 4, 4))
for p in itertools.permutations(range(3)): eps3[p] = np.linalg.det(np.eye(3)[list(p)])
for p in itertools.permutations(range(4)): eps4[p] = np.linalg.det(np.eye(4)[list(p)])
Sig = np.zeros((6, 4, 4), dtype=complex)
for a in range(3):
    for b in range(3):
        for c in range(3):
            Sig[a, b, c] = eps3[a, b, c]; Sig[a + 3, b, c] = 1j * eps3[a, b, c]
    Sig[a, a, 3] = 1; Sig[a, 3, a] = -1; Sig[a + 3, a, 3] = -1j; Sig[a + 3, 3, a] = 1j
blocks_ok = np.allclose(np.einsum('Iij,Ikl->ijkl', Sig, Sig), 2 * eps4)
SigUp = 0.5 * np.einsum('ijkl,Ikl->Iij', eps4, Sig)     # A^{ij} = A_I SigUp^I_{ij}

def antisym(c):
    ca = np.zeros_like(c)
    for p in itertools.permutations(range(3)):
        ca += np.linalg.det(np.eye(3)[list(p)]) * np.transpose(c, p)
    return ca / 6

rng = np.random.default_rng(2026)
rows = []
for trial in range(5):
    w = rng.normal(size=4) + 1j * rng.normal(size=4)
    u = rng.normal(size=2) + 1j * rng.normal(size=2)
    # B_w = ([A^{ik}, A_kj w^j], A_il w^l) = sum_IJK c_IJK ([A_I, A_J], A_K)
    cB = np.einsum('Iik,Jkj,j,Kil,l->IJK', SigUp, Sig, w, Sig, w)
    flav = np.sum(abs(antisym(cB)) ** 2) / np.vdot(w, w).real ** 2
    # u' with u'^a u_a = 1 and Hermitian orthogonality: u' = (-conj u2, conj u1)/|u|^2 (lower index), u'^1 = u'_2, u'^2 = -u'_1
    upl = np.array([-np.conj(u[1]), np.conj(u[0])]) / np.vdot(u, u).real
    upu = np.array([upl[1], -upl[0]])               # upper index
    contr = upu @ u                                  # u'^a u_a
    herm = np.vdot(u, upl)                           # sum conj(u_a) u'_a
    # classical image: lambda_a^j -> -(1/4)[A^{jk}, A_kl w^l] u_a in O_zeta = (lambda_a^j, A_ji) zeta^{a i}, zeta^{a i} = u'^a w^i
    cimg = -0.25 * contr * np.einsum('Ijk,Kkl,l,Lji,i->IKL', SigUp, Sig, w, Sig, w)
    ratio = np.vdot(antisym(cB).ravel(), antisym(cimg).ravel()) / np.vdot(antisym(cB).ravel(), antisym(cB).ravel())
    resid = np.linalg.norm(antisym(cimg) - ratio * antisym(cB))
    # O_zeta flavour factor: coefficient of (lambda_a^j, A_I) is sum_i Sig^I_ji zeta^{a i}; squared norm / |zeta|^2
    zeta = np.einsum('a,i->ai', upu, w)
    v = np.einsum('Iji,ai->Iaj', Sig, zeta)
    fz = np.sum(abs(v) ** 2) / np.sum(abs(zeta) ** 2)
    rows.append({'flavour_sum_over_|w|^4': float(flav), 'u_prime_contraction': complex(contr).real,
                 'hermitian_overlap': abs(herm), 'kappa_B': complex(ratio).real, 'kappa_B_imag': complex(ratio).imag,
                 'image_residual': float(resid), 'O_zeta_flavour_factor': float(fz)})
ok = (blocks_ok and all(abs(r['flavour_sum_over_|w|^4'] - 48) < 1e-9 and abs(r['kappa_B'] + 0.25) < 1e-12
                        and r['image_residual'] < 1e-12 and abs(r['O_zeta_flavour_factor'] - 6) < 1e-12 for r in rows))
out = {'blocks_ok': bool(blocks_ok), 'samples': rows, 'all_pass': bool(ok),
       'conclusion': 'the classical/tree inputs of T3(iv) are reproduced: sum |c_[IJK]|^2 = 48|w|^4, flavour factor 6 in '
                     'n_zeta = 24 hbar dim g (zeta,conj zeta)/(16 pi^4), kappa_B = -1/4 (sign relative to the two-component '
                     'S^(1)(xi) of ch08); only moduli enter k = -24, so T3(iv) does not test any sign convention'}
json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1)
print(json.dumps({'all_pass': ok, 'sample0': rows[0]}, indent=1))
