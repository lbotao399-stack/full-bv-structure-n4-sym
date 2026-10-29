"""Numerical evaluation at finite eps of the class-(i) graphs (S_M on the cycle, all edges P(eps,infty)) in canonical form,
and of the one-loop two-point functions (class (ii) input), by direct Schwinger integration (no Taylor expansion):
   graph(eps) = sum_k 2^-k eps^{E-2-k} int_{[1,oo)^E} du U^{-2-k} C_k(beta) e^{-eps U rho(beta)}   (units (16 pi^2)^-1),
tensor Gauss-Legendre in z = log u on [0, Z].  Used for the finite parts of Part 1 and the route-A cross-check of the
arity-one obstruction."""
import json, math, sys, time
import numpy as np
import sympy as sp
import delta_engine as D
import delta_canon as C
import loops as L

q = L.q
BET = D.BET

def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def integrand_parts(Num, curs, ops):
    """returns per op: list of (k, numpy callable of betas) and rho callable."""
    E = len(curs)
    shift = [sum(BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
    Cm = L.gaussian_moments(Num, shift)
    rho = sp.expand(sum(BET[e] * sum(x * x for x in curs[e]) for e in range(E)) - sum(x * x for x in shift))
    out = {}
    for nm, op in ops.items():
        lst = []
        for k, Ck in Cm.items():
            # the eta-linear part of e^{-eps U rho} matters for the alpha_E op: handled by expanding rho = rho0 + eta rho1
            lst.append((k, Ck))
        out[nm] = lst
    return Cm, rho

def evaluate_graphs(recs, ops, eps, n=64, Z=None, with_rho_eta=True):
    """recs: enumerate_terms records (P mode); ops: dict name -> function on sympy expressions in (eta, AD, betas).
    The exponential e^{-eps U rho}, with rho = rho0 + eta rho1 (+O(eta^2)), is expanded to first order in eta:
    e^{-eps U rho} = e^{-eps U rho0} (1 - eta eps U rho1).  Returns {op: complex value} (units (16 pi^2)^-1)."""
    ETA = C.ETA
    tot = {nm: 0j for nm in ops}
    Zm = Z if Z is not None else math.log(1.0 / eps) + 14.0
    for r in recs:
        Num = C.trunc_eta(r['Num']); curs = [[C.trunc_eta(x) for x in c] for c in r['curs']]
        E = len(curs)
        shift = [sum(BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
        Cm = L.gaussian_moments(Num, shift)
        rho = sp.expand(sum(BET[e] * sum(x * x for x in curs[e]) for e in range(E)) - sum(x * x for x in shift))
        rho0 = sp.expand(rho.coeff(ETA, 0)); rho1 = sp.expand(rho.coeff(ETA, 1))
        Us = sp.Symbol('Usym')
        epsS = sp.Symbol('epsS')
        z, wz = gl(n, 0.0, Zm)
        grids = np.meshgrid(*([z] * E), indexing='ij')
        W = np.ones_like(grids[0])
        for g in grids:
            W = W * 1.0
        Wt = np.einsum(','.join('ijk'[:E][i] for i in range(E)) + '->' + 'ijk'[:E], *([wz] * E))
        U_e = [np.exp(g) for g in grids]
        Wt = Wt * np.prod(U_e, axis=0)
        Ut = sum(U_e)
        betas = [u / Ut for u in U_e]
        r0 = sp.lambdify(BET[:E], rho0, 'numpy')
        rho0v = np.broadcast_to(np.asarray(r0(*betas), dtype=complex), Ut.shape)
        expo = np.exp(-eps * Ut * rho0v.real)
        for nm, op in ops.items():
            val = 0j
            for k, Ck in Cm.items():
                # full expression with the eta-expanded exponential factor: Ck * (1 - eta eps U rho1)
                ex_plain = sp.expand(op(C.trunc_eta(Ck)))
                ex_rho = sp.expand(op(C.trunc_eta(Ck * rho1 * ETA))) if (with_rho_eta and rho1 != 0) else 0
                f0 = sp.lambdify(BET[:E], ex_plain, 'numpy') if ex_plain != 0 else None
                f1 = sp.lambdify(BET[:E], ex_rho, 'numpy') if ex_rho != 0 else None
                pref = 2.0 ** (-k) * eps ** (E - 2 - k) * Ut ** (-2.0 - k)
                acc = 0
                if f0 is not None:
                    acc = acc + np.broadcast_to(np.asarray(f0(*betas), dtype=complex), Ut.shape)
                if f1 is not None:
                    acc = acc - eps * Ut * np.broadcast_to(np.asarray(f1(*betas), dtype=complex), Ut.shape)
                if f0 is None and f1 is None:
                    continue
                val += np.sum(Wt * pref * acc * expo)
            tot[nm] += complex(r['fac']) * val
    return tot

def evaluate_special(recs, ops, eps, n=48, n1=24, Z=None):
    """numerical one-special-edge graphs at finite eps (special Schwinger time in [0, eps], numerator times -k_s^2),
    same conventions as evaluate_graphs (units (16 pi^2)^-1)."""
    ETA = C.ETA
    tot = {nm: 0j for nm in ops}
    Zm = Z if Z is not None else math.log(1.0 / eps) + 14.0
    for r in recs:
        Num = C.trunc_eta(r['Num']); curs = [[C.trunc_eta(x) for x in c] for c in r['curs']]
        E = len(curs); s = r['s']
        ks = [q[m] - curs[s][m] for m in range(4)]
        Nt = C.trunc_eta(-Num * sum(x * x for x in ks))
        shift = [sum(BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
        Cm = L.gaussian_moments(Nt, shift)
        rho = sp.expand(sum(BET[e] * sum(x * x for x in curs[e]) for e in range(E)) - sum(x * x for x in shift))
        rho0 = sp.expand(rho.coeff(ETA, 0)); rho1 = sp.expand(rho.coeff(ETA, 1))
        us, ws = gl(n1, 0.0, 1.0)
        z, wz = gl(n, 0.0, Zm)
        axes = []; wts = []
        for e in range(E):
            if e == s:
                axes.append(us); wts.append(ws)
            else:
                axes.append(np.exp(z)); wts.append(wz * np.exp(z))
        grids = np.meshgrid(*axes, indexing='ij')
        Wt = wts[0]
        for w_ in wts[1:]:
            Wt = np.multiply.outer(Wt, w_)
        Ut = sum(grids)
        betas = [g / Ut for g in grids]
        r0 = sp.lambdify(BET[:E], rho0, 'numpy')
        rho0v = np.broadcast_to(np.asarray(r0(*betas), dtype=complex), Ut.shape)
        expo = np.exp(-eps * Ut * rho0v.real)
        for nm, op in ops.items():
            val = 0j
            for k, Ck in Cm.items():
                ex_plain = sp.expand(op(C.trunc_eta(Ck)))
                ex_rho = sp.expand(op(C.trunc_eta(Ck * rho1 * ETA))) if rho1 != 0 else 0
                f0 = sp.lambdify(BET[:E], ex_plain, 'numpy') if ex_plain != 0 else None
                f1 = sp.lambdify(BET[:E], ex_rho, 'numpy') if ex_rho != 0 else None
                if f0 is None and f1 is None:
                    continue
                pref = 2.0 ** (-k) * eps ** (E - 2 - k) * Ut ** (-2.0 - k)
                acc = 0
                if f0 is not None:
                    acc = acc + np.broadcast_to(np.asarray(f0(*betas), dtype=complex), Ut.shape)
                if f1 is not None:
                    acc = acc - eps * Ut * np.broadcast_to(np.asarray(f1(*betas), dtype=complex), Ut.shape)
                val += np.sum(Wt * pref * acc * expo)
            tot[nm] += complex(r['fac']) * val
    return tot
