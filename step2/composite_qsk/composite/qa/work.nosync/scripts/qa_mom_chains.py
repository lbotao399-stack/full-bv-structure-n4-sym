"""QA third method: scalar-sector chains (letters A_perp with n-derivatives, A|_V letter with n-derivatives) built from
qa_mom series.  Values are the chain integrals (units hbar/16pi^2) BEFORE the overall factor F0 and the colour sign.
  zz(aS, aP, qV1, qV2, kind)  : k=2 chain  letter(S) -special- V1([A_perp,c],A*_perp) -K- V2(A_V A_perp^2) -P- letter(P)
                                colour word (out V2, out V1) read from the S-letter; value of kind
  az(r, a, q, kind, avec)     : k=1 chain  A letter (r derivatives, polarisation avec) -K- V1 -special- Z letter (a derivatives)
  za(a, r, q, kind, avec)     : k=1 chain  Z letter -special- V1 -K- A letter
kind: ('P', b) for C_{P_b}; 'D0' for the translation-invariant part along Delta (kernel form, Lemma lem:cc-kernels (i):
K_Delta = (v.(x)1 - 1(x)v.) X with v. = 2((x - x_O).d + w) on each end, positions by injected momenta absorbed by omega;
at a letter end the injection enters only the letter derivative factor)."""
import qa_mom as M
from qa_mom import Mom, S_mul, S_add, S_scale, S_const, lin, comp, expo, invsq, special, integrate, vadd, vscale, C, I, cc
from fractions import Fraction as Fr

NV = [1, 0, 0, 0]


def setn(n):
    global NV
    NV = list(n)


def letter(p, nder, sgn):
    """(sgn * i n.p)^nder"""
    out = S_const(1)
    f = S_scale(lin(NV, p), C(0, sgn))
    for _ in range(nder):
        out = S_mul(out, f)
    return out


def dotS(p1, p2):
    out = {}
    for m in range(4):
        out = S_add(out, S_mul(comp(p1, m), comp(p2, m)))
    return out


def Kcore(p):
    return S_scale(S_mul(expo(p, 2), invsq(p, 2)), 2 * I)


def zz_raw(aS, aP, qV1, qV2, Kfac, h1=None, h2=None):
    z4 = [0, 0, 0, 0]
    h1 = h1 or z4; h2 = h2 or z4
    ps = Mom(1, z4)
    pk = Mom(-1, [-(cc(a)) for a in qV1], [-cc(x) for x in h1])
    pP = Mom(-1, [-(cc(a) + cc(b)) for a, b in zip(qV1, qV2)], [-(cc(a) + cc(b)) for a, b in zip(h1, h2)])
    v2vec = pP.add(QL=[-cc(x) for x in qV2])          # p_P - q_V2
    A = special(ps)
    A = S_mul(A, letter(ps, aS, -1))
    A = S_mul(A, Kcore(pk))
    A = S_mul(A, S_scale(dotS(v2vec, pk), C(0, Fr(-1, 2))))
    A = S_mul(A, S_mul(S_const(2), S_mul(expo(pP, 1), invsq(pP, 1))))
    A = S_mul(A, letter(pP, aP, -1))
    A = S_mul(A, Kfac(pk))
    return A


def zz(aS, aP, qV1, qV2, kind):
    d = aS + aP
    if kind == 'D0':
        M.LMAX[0] = d; M.EMAX[0] = 1
        tot = {}
        for mu in range(4):
            e = [0, 0, 0, 0]; e[mu] = 1
            A1 = zz_raw(aS, aP, qV1, qV2, lambda pk: S_scale(comp(pk, mu), I), h1=e)
            tot = vadd(tot, vscale(integrate(A1, d, 1), 2 * C(0, -1)))
            A2 = zz_raw(aS, aP, qV1, qV2, lambda pk: S_scale(comp(pk, mu), C(0, -1)), h2=e)
            tot = vadd(tot, vscale(integrate(A2, d, 1), -2 * C(0, -1)))
        A3 = zz_raw(aS, aP, qV1, qV2, lambda pk: S_const(1))
        tot = vadd(tot, vscale(integrate(A3, d, 0), -2))
        return tot
    b = kind[1]
    M.LMAX[0] = d + 1; M.EMAX[0] = 0
    A = zz_raw(aS, aP, qV1, qV2, lambda pk: S_scale(lin(b, pk), 2 * I))
    return integrate(A, d + 1, 0)


def az_raw(r, a, q, avec, Kfac, h0=None, h1=None):
    z4 = [0, 0, 0, 0]
    h0 = h0 or z4; h1 = h1 or z4
    p = Mom(1, z4)
    pA = Mom(1, z4, h0)                       # letter derivative sees p + h0
    ps = Mom(1, [-cc(x) for x in q], [-cc(x) for x in h1])
    A = letter(pA, r, +1)
    A = S_mul(A, S_mul(lin(avec, p), Kcore(p)))
    A = S_mul(A, special(ps))
    A = S_mul(A, letter(ps, a, -1))
    A = S_mul(A, Kfac(p))
    return A


def az(r, a, q, kind, avec):
    d = r + a + 1
    if kind == 'D0':
        M.LMAX[0] = d; M.EMAX[0] = 1
        tot = {}
        for mu in range(4):
            e = [0, 0, 0, 0]; e[mu] = 1
            A1 = az_raw(r, a, q, avec, lambda p: S_scale(comp(p, mu), I), h0=e)
            tot = vadd(tot, vscale(integrate(A1, d, 1), 2 * C(0, -1)))
            A3 = az_raw(r, a, q, avec, lambda p: S_scale(comp(p, mu), C(0, -1)), h1=e)
            tot = vadd(tot, vscale(integrate(A3, d, 1), -2 * C(0, -1)))
        A2 = az_raw(r, a, q, avec, lambda p: S_const(1))
        tot = vadd(tot, vscale(integrate(A2, d, 0), 2))
        return tot
    b = kind[1]
    M.LMAX[0] = d + 1; M.EMAX[0] = 0
    A = az_raw(r, a, q, avec, lambda p: S_scale(lin(b, p), 2 * I))
    return integrate(A, d + 1, 0)


def za_raw(a, r, q, avec, Kfac, h0=None, h1=None):
    z4 = [0, 0, 0, 0]
    h0 = h0 or z4; h1 = h1 or z4
    p = Mom(1, z4)
    pA = Mom(1, z4, [-cc(x) for x in h0])      # letter derivative (-i n.p + i n.h0) = (-i n.(p - h0))
    ps = Mom(-1, [-cc(x) for x in q], [-cc(x) for x in h1])
    A = letter(pA, r, -1)
    A = S_mul(A, S_mul(lin(avec, p), Kcore(p)))
    A = S_mul(A, special(ps))
    A = S_mul(A, letter(ps, a, -1))
    A = S_mul(A, Kfac(p))
    return A


def za(a, r, q, kind, avec):
    d = r + a + 1
    if kind == 'D0':
        M.LMAX[0] = d; M.EMAX[0] = 1
        tot = {}
        for mu in range(4):
            e = [0, 0, 0, 0]; e[mu] = 1
            A1 = za_raw(a, r, q, avec, lambda p: S_scale(comp(p, mu), I), h1=e)
            tot = vadd(tot, vscale(integrate(A1, d, 1), 2 * C(0, -1)))
            A2 = za_raw(a, r, q, avec, lambda p: S_scale(comp(p, mu), C(0, -1)), h0=e)
            tot = vadd(tot, vscale(integrate(A2, d, 1), -2 * C(0, -1)))
        A3 = za_raw(a, r, q, avec, lambda p: S_const(1))
        tot = vadd(tot, vscale(integrate(A3, d, 0), -2))
        return tot
    b = kind[1]
    M.LMAX[0] = d + 1; M.EMAX[0] = 0
    A = za_raw(a, r, q, avec, lambda p: S_scale(lin(b, p), 2 * I))
    return integrate(A, d + 1, 0)
