"""TQ4 (two-body part) and G-plane (route B): the scalar bilinears O = int omega (iota_w1 A, iota_w2 A) along xi.

Two-body part of Q_1(O) = Upsilon_xi(A_perp, A_perp) at za = w1, zb = w2, contracted with the colour pairing of the
singlet: the single-output word gives tr(ad_Y) = 0, the two-output words (0,1), (1,0) both give tr(ad_Y0 ad_Y1)
(identical colour factor c_adj (Y0,Y1)), so the singlet value is P_(0,1) + P_(1,0) in units hbar/(16 pi^2) times
that factor.  It is decomposed on the leg sets:
  A_perp lam_+ :  (rho(w1) xi, Y1)(w2.Y0) + (w1.Y0)(rho(w2) xi, Y1)  [the classical structure {S^(1)(xi), O}]
                  and (w1.w2)(rho(Y0) xi, Y1) [trace structure]
  A_V lam_-    :  (w1.w2) (Gamma(xi, Y1)|_V . Y0) [trace]  (and a possible traceless structure, tested)
  c lam*_+     :  (w1.w2) (xi, Y1) Y0  [trace]
G-plane (Prop qs-plane): for xi = u (x) e_4 and w_1, w_2 with rho(w_a) xi = 0 the singlet value vanishes on every
leg set (the one-body term vanishes on these letters for every quadratic F_1^(1)(xi), since its arrows out of A_perp
are multiples of iota_w Gamma(xi,lam) = (rho(w) xi, lam) = 0)."""
import os, sys, json, time, random
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.domains import QQ_I
import qb_core as B
import qb_job as QJ

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')
LEGSETS = [['lamM'], ['AV', 'lamM'], ['AW', 'lamP'], ['c', 'lamsP']]


def expr(job, P):
    return sp.expand(P.as_expr(*[sp.Symbol(n) for n in job.names])) if P else sp.Integer(0)


def singlet(job, tot, m=0):
    out = {}
    for w, dm in tot.items():
        if len(w) != 2:
            continue
        for b, P in dm.get(m, {}).items():
            out[b] = sp.expand(out.get(b, 0) + expr(job, P))
    return out


def run():
    t0 = time.time()
    res = {}; gates = []
    # rank-one xi = u (x) e_4 (realisation component 3): determines the kernel for every xi by covariance; the
    # structures below are evaluated at the same xi
    xi = [0] * 16
    xi[3] = 1
    xv = sp.Matrix(xi)
    w1 = [sp.Symbol(f'za_{c}') for c in range(4, 10)]
    w2 = [sp.Symbol(f'zb_{c}') for c in range(4, 10)]
    dot = lambda a, b: sum(x * y for x, y in zip(a, b))

    def rhow(w):
        return sum((w[K] * B.Rm[4 + K] for K in range(6)), sp.zeros(16, 16))
    tabs = {}
    for outs in LEGSETS:
        job = QJ.DJob('AW', 'AW', outs, 'xi', xi, tag='bil')
        r = job.kernel(mlist=(0, 1))
        tabs[' '.join(outs)] = (job, r['total'])
        print('leg set', outs, round(r['seconds'], 1), 's', flush=True)
    # ---- decomposition of the singlet value on the structures
    dec = {}
    # A_perp lam_+
    job, tot = tabs['AW lamP']
    s = singlet(job, tot)
    Y0 = [sp.Symbol(f'o0_{c}') for c in range(4, 10)]
    Y1 = sp.zeros(16, 1)
    for c in range(8):
        Y1[c] = sp.Symbol(f'o1_{c}')
    pair = lambda Q, lam: (lam.T * B.Mp * Q)[0]      # (lam, Q*) with Q* in Sigma^*, lam in Sigma
    cl = sp.expand(pair(rhow(w1) * xv, Y1) * dot(w2, Y0) + dot(w1, Y0) * pair(rhow(w2) * xv, Y1))
    trc = sp.expand(dot(w1, w2) * pair(rhow(Y0) * xv, Y1))
    c1, c2 = sp.symbols('c1 c2')
    out_b = {}
    for b, e in s.items():
        eqs = sp.Poly(sp.expand(e - c1 * cl - c2 * trc), *(w1 + w2 + Y0 + list(Y1[:8]))).coeffs()
        sol = sp.solve(eqs, [c1, c2], dict=True)
        out_b[b] = {str(k): str(v) for k, v in sol[0].items()} if sol else 'NOT IN SPAN'
    dec['A_perp lam_+'] = out_b
    # A_V lam_-
    job, tot = tabs['AV lamM']
    s = singlet(job, tot)
    Y0v = [sp.Symbol(f'o0_{c}') for c in range(4)]
    Y1m = sp.zeros(16, 1)
    for c in range(8, 16):
        Y1m[c] = sp.Symbol(f'o1_{c}')
    G = B.Gam(xi, list(Y1m))
    trv = sp.expand(dot(w1, w2) * sum(G[m] * Y0v[m] for m in range(4)))
    c3 = sp.Symbol('c3')
    out_b = {}
    for b, e in s.items():
        eqs = sp.Poly(sp.expand(e - c3 * trv), *(w1 + w2 + Y0v + list(Y1m[8:]))).coeffs()
        sol = sp.solve(eqs, [c3], dict=True)
        out_b[b] = {str(k): str(v) for k, v in sol[0].items()} if sol else 'NOT IN SPAN'
    dec['A_V lam_-'] = out_b
    # c lam*_+
    job, tot = tabs['c lamsP']
    s = singlet(job, tot)
    Yc = sp.Symbol('o0_0')
    Ys = sp.zeros(16, 1)
    for c in range(8, 16):
        Ys[c] = sp.Symbol(f'o1_{c}')
    trc2 = sp.expand(dot(w1, w2) * Yc * (xv.T * B.Mp * Ys)[0])
    c4 = sp.Symbol('c4')
    out_b = {}
    for b, e in s.items():
        eqs = sp.Poly(sp.expand(e - c4 * trc2), *(w1 + w2 + [Yc] + list(Ys[8:]))).coeffs()
        sol = sp.solve(eqs, [c4], dict=True)
        out_b[b] = {str(k): str(v) for k, v in sol[0].items()} if sol else 'NOT IN SPAN'
    dec['c lam*_+'] = out_b
    # single output lam_-: colour singlet value 0 (tr ad_Y); record that it has only the word (0,)
    job, tot = tabs['lamM']
    dec['lam_- (single output)'] = dict(words=[str(list(w)) for w in tot], singlet='0 (tr ad_Y = 0)')
    res['decomposition_eps0'] = dec
    # poles of the singlet values
    pol = {}
    for k, (job, tot) in tabs.items():
        s1 = singlet(job, tot, 1)
        pol[k] = {b: str(e) for b, e in s1.items()}
    res['singlet_poles'] = pol
    # ---- G-plane at xi = u (x) e_4: realisation components xi_3 = 1 (a = 0, W index 3)
    xin = {f'xi{j}': (1 if j == 3 else 0) for j in range(8)}
    xnum = sp.Matrix([xin.get(f'xi{j}', 0) for j in range(8)] + [0] * 8)
    Wsym = sp.symbols('w0:6')
    cond = sum((Wsym[K] * B.Rm[4 + K] for K in range(6)), sp.zeros(16, 16)) * xnum
    solw = sp.solve(list(cond), Wsym, dict=True)[0]
    free = sorted(set(Wsym) - set(solw.keys()), key=lambda s_: s_.name)
    basis = []
    for f in free:
        vec = [sp.sympify(solw.get(Wk, Wk)).subs({g: (1 if g == f else 0) for g in free}) for Wk in Wsym]
        basis.append(vec)
    res['chiral_w_basis'] = [[str(x) for x in v] for v in basis]
    random.seed(11)
    plane = {}
    okp = True
    for k, (job, tot) in tabs.items():
        for trial in range(2):
            a = [Fr(random.randint(-3, 3), random.randint(1, 2)) for _ in basis]
            b2 = [Fr(random.randint(-3, 3), random.randint(1, 2)) for _ in basis]
            wa = [sum(ai * v[K] for ai, v in zip(a, basis)) for K in range(6)]
            wb = [sum(bi * v[K] for bi, v in zip(b2, basis)) for K in range(6)]
            sub = {}
            sub.update({sp.Symbol(f'za_{4 + K}'): wa[K] for K in range(6)})
            sub.update({sp.Symbol(f'zb_{4 + K}'): wb[K] for K in range(6)})
            for m in (0, 1):
                s_ = singlet(job, tot, m)
                for bkey, e in s_.items():
                    v = sp.expand(e.subs(sub))
                    plane[f'{k} | eps^-{m} | {bkey} | trial {trial}'] = str(v)
                    if v != 0:
                        okp = False
    res['G_plane_values'] = plane
    gates.append(dict(id='G-plane', pass_=okp, value='all singlet values zero' if okp else {k: v for k, v in plane.items() if v != '0'},
                      target='Prop qs-plane: Upsilon_xi(A_perp,A_perp) on chiral w1, w2 (rho(w_a) xi = 0, xi = u (x) e_4), colour singlet, every leg set, eps^0 and poles: 0'))
    res['seconds'] = round(time.time() - t0, 1)
    return res, gates


if __name__ == '__main__':
    res, gates = run()
    json.dump(dict(results=res, gates=gates), open(os.path.join(OUT, 'bilinears_xi.json'), 'w'), indent=1)
    print(json.dumps(res['decomposition_eps0'], indent=1))
    for g in gates:
        print(g['id'], 'PASS' if g['pass_'] else 'FAIL')
