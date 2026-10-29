"""Step 2, engine A: independent numerical check of the cubic sector at finite eps (no eps-expansion of any graph).

For one sector and one generic configuration (integer momenta, integer polarisations), every one-loop graph of
Q W(P(eps,infty), I - I^CT[eps]) with these legs is evaluated at finite eps:
  I-a      : value(eps) = (16pi^2)^{-1} sum_j (4^j j!)^{-1} int dt T^{-2-j} F_j(u) e^{-T q(u)}, special time in [0,eps],
             others in [eps,infty).  In the variables (T, u) the T-integral is done in closed form (exponential integral
             E_1 and recursion), the u-integral numerically (scipy, split along the kinks); for the j-terms whose
             e^{-Tq} -> 1 integral converges, that part is taken exactly from s2_integrals and only the remainder is
             integrated numerically (no large cancellations);
  I-b + III(K') = sum_l (e^{-eps p_l^2} - 1) x (W_1^(2)(eps, p_l) - I^CT_2) contracted with I_3 by Khat,
             W_1^(2)(eps) from the ordinary bubbles at finite eps (s2_graphs.PCycles), 1-dim integrals;
  I-d + II-b = sum_l (special bubbles(eps, p_l) - Q I^CT_2) e^{-eps p_l^2}/p_l^2 contracted with I_3 by Pnum;
  II-a + III(K_0): log eps parts cancel identically (s2_trees), eps^{-1} part 0.
The sum is compared with the predicted limit  alpha . G + (tree of Ob_2 with I_3)  at two values of eps."""
import json, sys, time, itertools, math
from fractions import Fraction
import numpy as np
from scipy import integrate, special
import s2_core as C
from r2_core import Q2
import s2_graphs as SG
import s2_integrals as SI
import s2_ga as GAm
import s2_Ia as IA
from grass import GA

LOG2, LOG3 = math.log(2), math.log(3)
OPT = dict(epsabs=1e-14, epsrel=1e-12, limit=400)


def cnum(z):
    return complex(float(z.re), float(z.im))


def poly_fn(P):
    terms = [(np.array(k), cnum(v)) for k, v in P.t.items()]
    def f(*u):
        s_ = 0j
        for k, c in terms:
            m = c
            for ui, e in zip(u, k):
                if e:
                    m *= ui ** e
            s_ += m
        return s_
    return f


def Tint(m, a, b, c):
    """int_a^b T^m e^{-cT} dT, c >= 0, b may be inf (or huge); boundary terms with e^{-cb} = 0 are dropped"""
    if c == 0:
        if m == -1:
            return math.log(b / a)
        return (b ** (m + 1) - a ** (m + 1)) / (m + 1)
    eb = 0.0 if (b == math.inf or c * b > 700) else math.exp(-c * b)
    ea = math.exp(-c * a)
    def bt(k):          # b^k e^{-cb}
        return 0.0 if eb == 0.0 else (b ** k) * eb
    if m >= 0:
        val = (ea - eb) / c
        for k in range(1, m + 1):
            val = ((a ** k) * ea - bt(k)) / c + k / c * val
        return val
    val = special.exp1(c * a) - (special.exp1(c * b) if eb != 0.0 else 0.0)
    for k in range(-2, m - 1, -1):
        val = (bt(k + 1) - (a ** (k + 1)) * ea) / (k + 1) + c / (k + 1) * val
    return val


def special_value(n, Fsum, q, eps, Jexact):
    """finite-eps value of a group of special-edge graphs with common momentum routing; Fsum = {j: Poly}, includes pref.
    Jexact = {j: exact c->0 value (float)} for the j-terms whose c=0 integral converges"""
    qf = poly_fn(q)
    tot = 0j
    for j, Fj in Fsum.items():
        if not Fj.t:
            continue
        ff = poly_fn(Fj)
        m = n - 3 - j
        pre = eps ** (n - 2 - j) / (4 ** j * math.factorial(j))
        sub = j in Jexact
        if n == 3:
            def g(u3, u2, part):
                u1 = 1 - u2 - u3
                a = 1 / min(u2, u3); b = 1 / u1 if u1 > 1e-300 else 1e300
                if a >= b:
                    return 0.0
                c = eps * qf(u2, u3).real
                val = Tint(m, a, b, c) - (Tint(m, a, b, 0.0) if sub else 0.0)
                w = ff(u2, u3)
                return (w.real if part == 0 else w.imag) * val
            res = 0j
            for part in (0, 1):
                acc = 0.0
                # u2 <= u3 piece and u3 < u2 piece; region u1 <= min(u2,u3)
                acc += integrate.dblquad(lambda u3, u2: g(u3, u2, part), 0, 0.5,
                                         lambda u2: max(u2, 1 - 2 * u2, (1 - u2) / 2), lambda u2: 1 - u2,
                                         epsabs=1e-13, epsrel=1e-11)[0]
                acc += integrate.dblquad(lambda u3, u2: g(u3, u2, part), 1 / 3, 1,
                                         lambda u2: max(1 - 2 * u2, (1 - u2) / 2, 0), lambda u2: min(u2, 1 - u2),
                                         epsabs=1e-13, epsrel=1e-11)[0]
                res += acc if part == 0 else 1j * acc
        elif n == 2:
            def g(u2, part):
                u1 = 1 - u2
                a = 1 / u2; b = 1 / u1 if u1 > 1e-300 else 1e300
                if a >= b:
                    return 0.0
                c = eps * qf(u2).real
                val = Tint(m, a, b, c) - (Tint(m, a, b, 0.0) if sub else 0.0)
                w = ff(u2)
                return (w.real if part == 0 else w.imag) * val
            res = integrate.quad(lambda u: g(u, 0), 0.5, 1, **OPT)[0] + 1j * integrate.quad(lambda u: g(u, 1), 0.5, 1, **OPT)[0]
        else:
            raise ValueError(n)
        tot += pre * (res + (Jexact[j] if sub else 0.0))
    return tot


def P_value(n, Fsum, q, eps):
    """finite-eps value of ordinary bubbles (both edges in [eps, infty)) or quartic tadpoles (n = 1), common routing"""
    if n == 1:          # int_eps^infty dt t^{-2} F_0 = F_0 / eps  (F_0 constant)
        F0 = Fsum.get(0)
        return (cnum(F0.t.get((), Q2(0))) / eps) if F0 is not None else 0j
    assert n == 2
    qf = poly_fn(q)
    tot = 0j
    for j, Fj in Fsum.items():
        if not Fj.t:
            continue
        ff = poly_fn(Fj)
        m = -1 - j
        pre = eps ** (-j) / (4 ** j * math.factorial(j))
        def g(u2, part):
            u1 = 1 - u2
            a = 1 / min(u1, u2)
            c = eps * qf(u2).real
            w = ff(u2)
            return (w.real if part == 0 else w.imag) * Tint(m, a, math.inf, c)
        res = 0j
        for lo, hi in ((0, 0.5), (0.5, 1)):
            res += integrate.quad(lambda u: g(u, 0), lo, hi, **OPT)[0] + 1j * integrate.quad(lambda u: g(u, 1), lo, hi, **OPT)[0]
        tot += pre * res
    return tot


def exact_c0(n, Fj, j):
    """exact int dtau T^{-2-j} F_j(u) over [0,1] x [1,infty)^{n-1} (s2_integrals), as float"""
    tot = 0j
    for beta, v in Fj.t.items():
        Np = 2 + j + sum(beta)
        if Np - sum(beta) <= n - 1:
            return None
        c0, c2, c3 = SI.J(n, beta, Np)
        val = float(c0) + float(c2) * LOG2 + float(c3) * LOG3
        tot += cnum(v) * val
    return tot


def group_special(G, special=True):
    groups = {}
    for (seq, Ls, sgn) in G.enumerate():
        n, pref, Fj, q = SG.graph_polys(G, seq, Ls, sgn, special=special)
        key = tuple(seq)
        g = groups.setdefault(key, {'n': n, 'q': q, 'F': {}})
        for j, P in Fj.items():
            P = P.scale(Q2(pref))
            g['F'][j] = g['F'][j] + P if j in g['F'] else P
    return groups


def Ia_value(L, deg, eps):
    G = SG.OrientedCycles(L, deg)
    groups = group_special(G)
    tot = 0j
    for key, g in groups.items():
        Jex = {}
        for j, Fj in g['F'].items():
            if not Fj.t:
                continue
            e = exact_c0(g['n'], Fj, j)
            if e is not None:
                Jex[j] = e
        tot += special_value(g['n'], g['F'], g['q'], eps, Jex)
    return tot / C.C_ADJ


def twopoint_value(legs, deg, eps, special):
    G = (SG.OrientedCycles if special else SG.PCycles)(legs, deg)
    groups = group_special(G, special=special)
    tot = 0j
    for key, g in groups.items():
        if special:
            Jex = {}
            for j, Fj in g['F'].items():
                if Fj.t:
                    e = exact_c0(g['n'], Fj, j)
                    if e is not None:
                        Jex[j] = e
            tot += special_value(g['n'], g['F'], g['q'], eps, Jex)
        else:
            tot += P_value(g['n'], g['F'], g['q'], eps)
    return tot / C.C_ADJ


def run(sector, seed, eps_list):
    rng = np.random.default_rng(seed)
    cfg = IA.SECTORS[sector]
    while True:
        ps = [rng.integers(-2, 3, size=4) for _ in range(2)]
        ps = [-(ps[0] + ps[1])] + ps
        allm = ps + [-p for p in ps]
        if not (all(p.any() for p in ps) and len({tuple(p) for p in allm}) == 6):
            continue
        L = [(t, w, C.rpol(rng, t, w), p) for (t, w), p in zip(cfg['legs'], ps)]
        if all(abs(complex(x)) > 0 for x in cfg['st'](L, 0)[0][:2]):
            break
    # predicted limit
    alpha = json.load(open('logs/alpha_engineA.json'))['alpha_float']
    Gv = cfg['G'](L)
    pred_local = sum(alpha[g] * complex(v) for g, v in Gv.items())
    # mode model: contractions with I_3
    M = GAm.Model(allm)
    rs = GAm.rvars(3)
    legs = [(t, pol, i, rs[i] if C.PARITY[t] else None) for i, (t, w, pol, p) in enumerate(L)]
    odd_r = [rs[i] for i, (t, w, pol, p) in enumerate(L) if C.PARITY[t]]
    perms = list(itertools.permutations([0, 1, 2]))
    I3 = sum((M.cubic(nm, perms) for nm in ['AAA', 'llA', 'AcAs', 'cllS', 'cccS']), GA())
    neg = {0: 3, 1: 4, 2: 5}
    contr = {}
    for l in range(3):
        i, j = [x for x in range(3) if x != l]
        F = M.subst(I3, [legs[i], legs[j]], None)
        pairs = [(l, neg[l]), (neg[l], l)]
        two = lambda Fq: M.subst(Fq, [legs[l]], None)
        K0 = M.kernel('K', l, neg[l]); P1 = M.kernel('P', l, neg[l])
        for nm in ('AA', 'ldl', 'dcAs', 'dcAV', 'dsAdsA'):
            X = P1 if nm == 'dcAV' else K0
            contr[(l, nm)] = cnum(GAm.coeff(M.C(X, F, two(M.quad(nm, pairs))), odd_r))
    tree = sum(Fraction(-27, 8) * contr[(l, 'dcAV')] for l in range(3))      # Ob_2 = -(27/8) int(d^*dc,d^*A|_V)
    tree = complex(tree)
    out = {'sector': sector, 'algebra': C.ALG, 'momenta': [p.tolist() for p in ps],
           'predicted_limit': {'alpha.G': str(pred_local), 'tree_Ob2_I3': str(tree)}, 'eps': {}}
    # two-point data per leg
    for eps in eps_list:
        t0 = time.time()
        Ia = Ia_value(L, cfg['deg'], eps)
        IbIII = 0j; IdIIb = 0j
        for l in range(3):
            t, w, pol, p = L[l]; p2 = int(p @ p)
            kfac = math.expm1(-eps * p2)
            if t == 'A':            # A self-energy (K' joins A^* of I_3 with A of the two-point function)
                rng2 = np.random.default_rng(11)
                M2 = GAm.Model([p, -p])
                if w == 'V':        # two structures: int(A,A) and int(d^*A|_V, d^*A|_V)
                    rows_ = []; rhs_ = []
                    for rep in range(4):
                        a1 = C.rpol(rng2, 'A', 'V'); a2 = C.rpol(rng2, 'A', 'V')
                        se = twopoint_value([('A', 'V', a1, p), ('A', 'V', a2, -p)], 2, eps, False)
                        lg = [('A', a1, 0, None), ('A', a2, 1, None)]
                        aa = cnum(GAm.coeff(M2.subst(M2.quad('AA', [(0, 1), (1, 0)]), lg, None), []))
                        dd = cnum(GAm.coeff(M2.subst(M2.quad('dsAdsA', [(0, 1), (1, 0)]), lg, None), []))
                        rows_.append([aa, dd]); rhs_.append(se)
                    Am = np.array(rows_); bm = np.array(rhs_)
                    sol = np.linalg.lstsq(Am, bm, rcond=None)[0]
                    assert np.linalg.matrix_rank(Am) == 2
                    out.setdefault('AV_selfenergy_two_structure_fit_residual', {})['eps=%g,l=%d' % (eps, l)] = \
                        float(np.max(np.abs(Am @ sol - bm)) / np.max(np.abs(bm)))
                    # I^CT_2 on A|_V: (9/4)/eps int(A,A) - log eps int(dA,dA), int(dA,dA) = p^2 int(A,A) - int(d^*A,d^*A)
                    ren1 = sol[0] - (9 / 4 / eps - p2 * math.log(eps)); ren2 = sol[1] - math.log(eps)
                    IbIII += kfac * (-(ren1 * contr[(l, 'AA')] + ren2 * contr[(l, 'dsAdsA')]))
                    sigs = None
                else:
                    sigs = []
                for rep in range(2 if w == 'P' else 0):
                    a1 = C.rpol(rng2, 'A', 'P'); a2 = C.rpol(rng2, 'A', 'P')
                    se = twopoint_value([('A', 'P', a1, p), ('A', 'P', a2, -p)], 2, eps, False)
                    aa = cnum(GAm.coeff(M2.subst(M2.quad('AA', [(0, 1), (1, 0)]), [('A', a1, 0, None), ('A', a2, 1, None)], None), []))
                    sigs.append(se / aa)
                if w == 'P':
                    out.setdefault('Wperp_selfenergy_ratio_polarisation_independent', {})['eps=%g,l=%d' % (eps, l)] = abs(sigs[0] - sigs[1])
                    # I^CT_2 on A|_Wperp per unit of the int(A,A) plane wave: (9/4)/eps - p^2 log eps
                    ren = sigs[0] - (9 / 4 / eps - p2 * math.log(eps))
                    IbIII += kfac * ren * (-contr[(l, 'AA')])
            if t == 'l':
                rng2 = np.random.default_rng(12)
                u1 = C.rpol(rng2, 'l'); u2 = C.rpol(rng2, 'l')
                se = twopoint_value([('l', None, u1, p), ('l', None, u2, -p)], 1, eps, False)
                M2 = GAm.Model([p, -p]); rr = GAm.rvars(2)
                st = GAm.coeff(M2.subst(M2.quad('ldl', [(0, 1), (1, 0)]), [('l', u1, 0, rr[0]), ('l', u2, 1, rr[1])], None), rr)
                sig = se / cnum(st)
                u3 = C.rpol(rng2, 'l'); u4 = C.rpol(rng2, 'l')
                se2 = twopoint_value([('l', None, u3, p), ('l', None, u4, -p)], 1, eps, False)
                st2 = GAm.coeff(M2.subst(M2.quad('ldl', [(0, 1), (1, 0)]), [('l', u3, 0, rr[0]), ('l', u4, 1, rr[1])], None), rr)
                out.setdefault('lambda_selfenergy_ratio_polarisation_independent', {})['eps=%g,l=%d' % (eps, l)] = \
                    abs(se2 / cnum(st2) - sig)
                ren = sig - 4 * math.log(eps)                  # I^CT_2 lambda part: 4 log eps int(lambda, dslash lambda)
                IbIII += kfac * ren * (-contr[(l, 'ldl')])
            if t == 'c':
                rng2 = np.random.default_rng(13)
                g1 = C.rpol(rng2, 'c'); al = C.rpol(rng2, 'As', 'V')
                se = twopoint_value([('c', None, g1, p), ('As', 'V', al, -p)], 1, eps, False)
                M2 = GAm.Model([p, -p]); rr = GAm.rvars(2)
                st = GAm.coeff(M2.subst(M2.quad('dcAs', [(0, 1), (1, 0)]), [('c', g1, 0, rr[0]), ('As', al, 1, rr[1])], None), rr)
                sig = se / cnum(st)
                ren = sig + math.log(eps)                      # I^CT_2 ghost part: -log eps int(dc, A^*)
                IbIII += kfac * ren * (-contr[(l, 'dcAs')])
            if abs(contr[(l, 'dcAV')]) > 0:
                # special bubbles (c at p, A|_V at -p) at finite eps, per unit of int(dc, A|_V)
                rng3 = np.random.default_rng(14)
                g2 = C.rpol(rng3, 'c'); a3 = C.rpol(rng3, 'A', 'V')
                sbv = twopoint_value([('c', None, g2, p), ('A', 'V', a3, -p)], 3, eps, True)
                M3 = GAm.Model([p, -p]); rr3 = GAm.rvars(1)
                st3 = GAm.coeff(M3.subst(M3.quad('dcAV', [(0, 1), (1, 0)]), [('c', g2, 0, rr3[0]), ('A', a3, 1, None)], None), rr3)
                sig_sb = sbv / cnum(st3)
                IdIIb += (sig_sb + 4.5 / eps) * math.exp(-eps * p2) / p2 * contr[(l, 'dcAV')]
                out.setdefault('special_bubble_check', {})['eps=%g,l=%d' % (eps, l)] = {
                    'sigma_sb*eps': str(sig_sb * eps), 'sigma_sb + 9/(2eps)': str(sig_sb + 4.5 / eps),
                    'expected_limit -27/8 p^2': -27 / 8 * p2}
        total = Ia + IbIII + IdIIb
        dev = total - (pred_local + tree)
        out['eps']['%g' % eps] = {'I-a - alpha.G': str(Ia - pred_local), 'I-d+II-b - tree': str(IdIIb - tree),
                                 'eps*log(eps)': eps * math.log(eps),
                                 'I-a': str(Ia), 'I-b+III(K\')': str(IbIII), 'I-d+II-b': str(IdIIb), 'total': str(total),
                                 'total - predicted limit': str(dev), '|dev|': abs(dev),
                                 '|predicted|': abs(pred_local + tree), 'seconds': round(time.time() - t0, 1)}
        print(sector, 'eps', eps, out['eps']['%g' % eps], flush=True)
    return out


if __name__ == '__main__':
    sector = sys.argv[1] if len(sys.argv) > 1 else 'c,l,l'
    epsl = [float(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [1e-3, 1e-4]
    res = run(sector, int(sys.argv[2]) if len(sys.argv) > 2 else 5, epsl)
    json.dump(res, open('logs/s2_numeric_%s_%s.json' % (sector.replace(',', ''), C.ALG), 'w'), indent=1)
