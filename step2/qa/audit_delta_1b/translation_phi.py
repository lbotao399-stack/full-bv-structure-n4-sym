"""Invertible-H model (fixed momentum): checks of Prop d-psi (iii) (the Q-primitive Phi_T of the composite
line, at the level of functionals) and of Prop d-P (translations: Psi_v, Ward identity without scaling term).

Model: two copies of the fibre of model.py at one fixed momentum (N = 0, s a nonzero rational), so H = s^2
on every slot.  f := int_eps^T e^{-t s^2} dt is replaced by an arbitrary positive rational (only the
algebraic relations k = s^2 f, P* = f QGF, (Q x 1 + 1 x Q)P = k K_0 = K_eps - K_T are used).
M = (rotation of the two copies) x 1: even, commutes with Q, QGF, H, preserves the pairing; it plays the
role of a translation v. (prop:d-P): [v.,Q] = [v.,QGF] = [v.,H] = 0.
Both convention systems of d1_checks.py are used.
"""
import json, os, sys, random
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp import Alg, mat_mul, mat_add, mat_scale, mat_comm, mat_maxabs, mat_id, tensor_act, graded_swap
from model import SLOTS
from conventions import test_functionals
from d1_checks import System, kernel_of_second_order

HERE = os.path.dirname(os.path.abspath(__file__))
S = str

def build_fixed(s=Fr(3, 2), f=Fr(2, 5), gF=Fr(-2)):
    copies = 2
    idx = {}; par = []; deg = []
    for cpy in range(copies):
        for (nm, dg, p, w, prt) in SLOTS:
            idx[(nm, cpy)] = len(par); par.append(p); deg.append(dg)
    n = len(par)
    alg = Alg(par, deg=[-d for d in deg])
    def arrow(table):
        X = {}
        for (src, dst, coef) in table:
            for cpy in range(copies):
                X[(idx[(dst, cpy)], idx[(src, cpy)])] = Fr(coef)
        return X
    Q = arrow([('c', 'AL', s), ('AT', 'ATs', Fr(-1, 2) * s * s), ('ALs', 'cs', -s), ('l1', 'ls2', s), ('l2', 'ls1', -s)])
    tab = [('ATs', 'AT', -2), ('AL', 'c', s), ('cs', 'ALs', -s), ('ls2', 'l1', s), ('ls1', 'l2', -s)]
    if gF != 0:
        tab.append(('ALs', 'AL', gF))
    G = arrow(tab)
    K0 = {}
    for cpy in range(copies):
        for (nm, dg, p, w, prt) in SLOTS:
            K0[(idx[(nm, cpy)], idx[(prt, cpy)])] = Fr(1)
    # rotation of the copies on every slot (antisymmetric, so it preserves the copy-diagonal pairing)
    M = {}
    for (nm, dg, p, w, prt) in SLOTS:
        M[(idx[(nm, 1)], idx[(nm, 0)])] = Fr(1)
        M[(idx[(nm, 0)], idx[(nm, 1)])] = Fr(-1)
    H = mat_comm(Q, G, 1, 1)
    Ps = mat_scale(f, G)
    k = s * s * f
    return dict(alg=alg, par=par, n=n, Q=Q, G=G, K0=K0, M=M, H=H, Ps=Ps, f=f, k=k, s=s, idx=idx)

def run(gF=Fr(-2)):
    m = build_fixed(gF=gF)
    alg, par, n = m['alg'], m['par'], m['n']
    Q, G, K0, M, H, Ps, k = m['Q'], m['G'], m['K0'], m['M'], m['H'], m['Ps'], m['k']
    ker = lambda A: mat_mul(A, K0)
    P = ker(Ps)
    out = {'model': 'fixed momentum s = 3/2, f = 2/5, two copies, M = copy rotation, gF = %s' % gF}
    out['H - s^2'] = S(mat_maxabs(mat_add(H, (-m['s'] ** 2, mat_id(n)))))
    out['(Mx1+1xM)K0'] = S(mat_maxabs(tensor_act(M, 0, K0, par)))
    out['[M,Q], [M,QGF]'] = [S(mat_maxabs(mat_comm(M, Q, 0, 1))), S(mat_maxabs(mat_comm(M, G, 0, 1)))]
    out['(Qx1+1xQ)P - k K0'] = S(mat_maxabs(mat_add(tensor_act(Q, 1, P, par), (-k, K0))))
    out['(Mx1+1xM)P (translation: Ward without scaling term)'] = S(mat_maxabs(tensor_act(M, 0, P, par)))
    PMP = mat_mul(mat_mul(Ps, M), Ps)
    out['P*MP* - M (P*)^2'] = S(mat_maxabs(mat_add(PMP, (-1, mat_mul(M, mat_mul(Ps, Ps))))))
    out['[Q, P*MP*]'] = S(mat_maxabs(mat_comm(Q, PMP, 1, 0)))
    Phi = mat_scale(1 / Fr(k), mat_mul(Ps, PMP))            # Phi_T = k^{-1} P* (P*MP*), odd
    out['Q Phi + Phi Q - P*MP* (operators)'] = S(mat_maxabs(mat_add(mat_comm(Q, Phi, 1, 1), (-1, PMP))))
    kPhi = ker(Phi)                                         # (Phi x 1)K_0, even kernel
    out['(Phi x 1)K_0 graded symmetric'] = S(mat_maxabs(mat_add(graded_swap(kPhi, par), (-1, kPhi))))
    out['(Qx1+1xQ)(Phi x 1)K_0 - (P*MP* x 1)K_0'] = S(mat_maxabs(mat_add(tensor_act(Q, 1, kPhi, par), (-1, ker(PMP)))))
    Pi = {(i, i): (Fr(-1) if par[i] else Fr(1)) for i in range(n)}
    kPhi_tw = ker(mat_mul(Pi, Phi))                          # the note's 'kernel of an odd operator' (Pi X x 1)K_0
    rng = random.Random(5)
    fs = test_functionals(alg, rng, n_cubic=700, n_quartic=4, terms=50)
    out['number of test functionals'] = len(fs)
    for name in ('II', 'I'):
        sysm = System({'alg': alg}, name)
        sysm.alg = alg
        r = {}
        Qf = lambda F: sysm.D(Q, 1, F)
        SM = sysm.generator(K0, lambda F: sysm.D(M, 0, F), -1)
        r['S_v exists'] = SM is not None
        Cop = lambda F: alg.add(sysm.dphi(P, alg.mul(SM, F)), (-1, alg.mul(sysm.dphi(P, SM), F)), (-1, alg.mul(SM, sysm.dphi(P, F))))
        Dop = lambda F: alg.add(sysm.dphi(P, Cop(F)), (-1, Cop(sysm.dphi(P, F))))
        r['Cop - D_{P*M}'] = S(max(alg.dist(Cop(F), sysm.D(mat_mul(Ps, M), 1, F)) for F in fs))
        kD = kernel_of_second_order(alg, Dop, sysm.dconv)
        r['Dop - d_{kernel(Dop)}'] = S(max(alg.dist(Dop(F), sysm.dphi(kD, F)) for F in fs))
        for c in (1, 2, -1, -2):
            r['kernel(Dop) - %d (P*MP* x 1)K_0' % c] = S(mat_maxabs(mat_add(kD, (-c, ker(PMP)))))
        r['|kernel(Dop)|'] = S(mat_maxabs(kD))
        r['Q Dop + Dop Q'] = S(max(alg.maxabs(alg.add(Qf(Dop(F)), Dop(Qf(F)))) for F in fs))
        r['Q Dop - Dop Q'] = S(max(alg.maxabs(alg.add(Qf(Dop(F)), (-1, Dop(Qf(F))))) for F in fs))
        # functional form of Psi = Q Phi + Phi Q
        dPhi = lambda F: sysm.dphi(kPhi, F)
        dPhit = lambda F: sysm.dphi(kPhi_tw, F)
        for c in (1, 2, -1, -2):
            r['Dop - %d (d_Phi Q - Q d_Phi)' % c] = S(max(alg.dist(Dop(F), alg.scale(c, alg.add(dPhi(Qf(F)), (-1, Qf(dPhi(F)))))) for F in fs))
            r['Dop - %d (d_Phi Q + Q d_Phi)  [note]' % c] = S(max(alg.dist(Dop(F), alg.scale(c, alg.add(dPhi(Qf(F)), Qf(dPhi(F))))) for F in fs))
            r['Dop - %d (d_Phitw Q + Q d_Phitw)  [note, twisted kernel]' % c] = S(max(alg.dist(Dop(F), alg.scale(c, alg.add(dPhit(Qf(F)), Qf(dPhit(F))))) for F in fs))
            r['Dop - %d (d_Phitw Q - Q d_Phitw)  [twisted kernel]' % c] = S(max(alg.dist(Dop(F), alg.scale(c, alg.add(dPhit(Qf(F)), (-1, Qf(dPhit(F)))))) for F in fs))
        r['d_{(P*MP* x1)K0} - (d_Phi Q - Q d_Phi)'] = S(max(alg.dist(sysm.dphi(ker(PMP), F), alg.add(dPhi(Qf(F)), (-1, Qf(dPhi(F))))) for F in fs))
        out['system ' + name] = r
    return out

if __name__ == '__main__':
    res = {'Feynman': run(Fr(-2)), 'without arrow': run(Fr(0))}
    with open(os.path.join(HERE, 'translation_phi_results.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
