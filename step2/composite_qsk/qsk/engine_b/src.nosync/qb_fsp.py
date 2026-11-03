"""Finite special vertices F_xi and the pole identity G-pole (route B) on the pairs of gauginos at n = 0.

At n = 0 the only nonzero two-letter counterterm poles of field pairs are ct_1(lam lam) with one output A (weight
count, notes.txt), and the counterterm cycles ct1_{xi,m} of field pairs vanish by weight; hence
    F_xi(X_a,X_b) = - sum_m ((-1)^m/m!) {S^(1)(xi), ct_m(X_a X_b)}_{H^m} = {S^(1)(xi), ct_1(X_a X_b)}_H,
with {S^(1)(xi), g}_H = g with its A output replaced by H (xi. e)_A = |p|^2 Gamma(xi, lam_out) (xi. of Prop c7-susy,
E0), only on (lam,lam) pairs and on the leg set of one gaugino output.
Pole identity (eq qs-pole, m = 1) on O = int omega X_a X_b with GENERAL colour polarisations of the letters (the
kernels carry colour words): s ct1_{xi,1}(O) = 0 and ct_2 = 0 by weight; QO = 0 (Q annihilates lambda and A_perp
letters); the terms vacuum-pole x spectator do NOT vanish for general colour: the cubic part of {S^(1)(xi),O} with the
vacuum ct_1(A A) = 2 (za.zb), and ct1_{xi,1}(sO) with the vacuum counterterm cycle ct1_xi(c, lam_+)|_0 = (1/2)(xi.z):
    [eps^-1] D_xi(X_a,X_b)|_lam = {S^(1)(xi), ct_1(X_a X_b)} - ct_1({S^(1)(xi),X_a} X_b) + ct_1(X_a {S^(1)(xi),X_b})
                                  - ct1_{xi,1}(sO)  (gaugino pairs; no cubic vacuum term: no lam-A propagator),
{S^(1)(xi), lam letter z.lam} = z.((1/2) rho(dA) xi) = sum_{M<4,N} c_MN d_M A_N, c_MN = (1/4) z^T (rb_M r_N - rb_N r_M) xi
(odd derivation: the minus sign on the second letter).  ct_1 of the pairs (d A, lam), (lam, d A) from the twoletter
engine B (vendored tl_engine.Job, deriv ('gen',1): coefficient of s_M za_N, resp. u_M zb_N)."""
import os, sys, json, time
import sympy as sp
from sympy.polys.domains import QQ_I
import qb_core as B
import qb_job as QJ
T = QJ.CJ.T          # vendored twoletter engine B

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')
I = sp.I
XI8 = B.XI[:8]
XIV = XI8 + [0] * 8
CASES = [('lamP', 'lamP', 'lamP', 'AW'), ('lamP', 'lamM', 'lamM', 'AV'), ('lamM', 'lamM', 'lamP', 'AW'),
         ('lamM', 'lamP', 'lamM', 'AV')]


def to_expr(job, P, rename=None):
    rename = rename or {}
    syms = [sp.Symbol(rename.get(n, n)) for n in job.names]
    return sp.expand(P.as_expr(*syms)) if P else sp.Integer(0)


def slot_comps(slot):
    return B.SLOT[slot][1]


def cMN(zsyms, slot):
    """c_MN for the classical transformation of the gaugino letter sum_i z_i lam_i (i in the slot components)."""
    z = sp.zeros(16, 1)
    for c, s in zip(slot_comps(slot), zsyms):
        z[c] = s
    xv = sp.Matrix(XIV)
    out = {}
    for M in range(4):
        for N in range(10):
            if N == M:
                continue
            v = (z.T * (B.Rbm[M] * B.Rm[N] - B.Rbm[N] * B.Rm[M]) * xv)[0] / 4
            v = sp.expand(v)
            if v != 0:
                out[(M, N)] = v
    return out


def gamma_sub(outslot):
    """o0_M (A output) -> Gamma(xi, lam_new)_M with lam_new the output polarisation in outslot."""
    lam = sp.zeros(16, 1)
    for c in slot_comps(outslot):
        lam[c] = sp.Symbol(f'o0_{c}')
    xv = sp.Matrix(XIV)
    return {sp.Symbol(f'oA_{M}'): sp.expand((xv.T * B.MR[M] * lam)[0]) for M in range(10)}


def ct1_pair_A(sa, sb, Aslot):
    """ct_1(X_a X_b) on the leg set A (tl_engine Job), as expression with o0_M renamed oA_M."""
    job = T.Job(sa, sb, [Aslot])
    res = job.kernel(keep_sequences=False)
    out = {}
    for w, d in res['total'].items():
        P = d.get(1)
        if P:
            ren = {f'o0_{M}': f'oA_{M}' for M in range(10)}
            out[w] = to_expr(job, P, ren)
    return out, res['seconds']


def ct1_deriv_pair(sa, sb, outslot, which):
    """ct_1 of the pair (d A, lam) [which = 'a': A letter first] or (lam, d A) [which = 'b'] on leg set outslot,
    returned as {word: {(M, N): expr}} = coefficient of s_M zA_N (resp. u_M zA_N)."""
    Aslot = sa if which == 'a' else sb
    job = T.Job(sa, sb, [outslot], deriv=('gen', 1))
    res = job.kernel(keep_sequences=False)
    ren = {}
    for c in slot_comps(Aslot):
        ren[('za_' if which == 'a' else 'zb_') + str(c)] = f'zA_{c}'
    out = {}
    for w, d in res['total'].items():
        P = d.get(1)
        if not P:
            continue
        e = to_expr(job, P, ren)
        dv = 's' if which == 'a' else 'u'
        co = {}
        for M in range(4):
            for N in slot_comps(Aslot):
                c = sp.expand(sp.diff(e, sp.Symbol(f'{dv}{M}'), sp.Symbol(f'zA_{N}')))
                if c != 0:
                    co[(M, N)] = c
        out[w] = co
    return out, res['seconds']


def run():
    t0 = time.time()
    res = {}; gates = []
    for (sa, sb, lout, Aslot) in CASES:
        key = f'{sa}|{sb} -> {lout}'
        # D job (route B defect cycles)
        dj = QJ.DJob(sa, sb, [lout], 'xi', XIV, tag='fsp')
        dk = dj.kernel(mlist=(0, 1))
        Dm1 = {w: to_expr(dj, dm.get(1, {}).get('1')) for w, dm in dk['total'].items()}
        D0 = {w: {b: to_expr(dj, P) for b, P in dm.get(0, {}).items()} for w, dm in dk['total'].items()}
        # ct_1(X_a X_b)|_A and the bracket with S^(1)(xi)
        ctA, s1 = ct1_pair_A(sa, sb, Aslot)
        gs = gamma_sub(lout)
        br = {w: sp.expand(e.subs(gs, simultaneous=True)) for w, e in ctA.items()}
        p2 = sum(sp.Symbol(f'r0_{m}') ** 2 for m in range(4))
        Fx = {w: sp.expand(p2 * e) for w, e in br.items()}
        # ct_1({S,X_a} X_b) and ct_1(X_a {S,X_b})
        Aa = 'AV' if sa == 'lamP' else 'AW'
        Ab = 'AV' if sb == 'lamP' else 'AW'
        ka, s2 = ct1_deriv_pair(Aa, sb, lout, 'a')
        kb, s3 = ct1_deriv_pair(sa, Ab, lout, 'b')
        za = [sp.Symbol(f'za_{c}') for c in slot_comps(sa)]
        zb = [sp.Symbol(f'zb_{c}') for c in slot_comps(sb)]
        ca = cMN(za, sa); cb = cMN(zb, sb)
        # -ct1_{xi,1}(sO): sO = {I,X_a}X_b - X_a{I,X_b}, {I, lam letter} = [c, lam]; the vacuum counterterm cycle
        # ct1_xi(c, lam_+) = Z (xi.z) (-,-)^{-1}, Z = 1/2 (all-P cycle, t_gpole2.py), times the other letter as spectator:
        # -ct1(sO) = Z [(Y.za)(xi.zb) [X_b in lam_+] + (Y.zb)(xi.za) [X_a in lam_+]] on the word (0,)
        Z = sp.Rational(1, 2)
        Yv = {c: sp.Symbol(f'o0_{c}') for c in slot_comps(lout)}
        xia = sum(XIV[c] * z for c, z in zip(slot_comps(sa), za))
        xib = sum(XIV[c] * z for c, z in zip(slot_comps(sb), zb))
        Ya = sum(Yv[c] * z for c, z in zip(slot_comps(sa), za) if c in Yv)
        Yb = sum(Yv[c] * z for c, z in zip(slot_comps(sb), zb) if c in Yv)
        sO = sp.expand(Z * ((Ya * xib if sb == 'lamP' else 0) + (Yb * xia if sa == 'lamP' else 0)))
        rhs = {}
        for w in set(br) | set(ka) | set(kb) | set(Dm1) | {(0,)}:
            v = br.get(w, 0)
            v -= sum((ka.get(w, {}).get(mn, 0) * c for mn, c in ca.items()), 0)
            v += sum((kb.get(w, {}).get(mn, 0) * c for mn, c in cb.items()), 0)
            if w == (0,):
                v += sO
            rhs[w] = sp.expand(v)
        diff = {str(w): str(sp.expand(Dm1.get(w, 0) - rhs.get(w, 0))) for w in rhs}
        diff_alt = {str(w): str(sp.expand(Dm1.get(w, 0) + rhs.get(w, 0))) for w in rhs}
        nontriv = any(e != 0 for e in br.values()) or any(c != 0 for d_ in ka.values() for c in d_.values()) or \
            any(c != 0 for d_ in kb.values() for c in d_.values())
        ok = all(v == '0' for v in diff.values()) and nontriv
        res[key] = dict(D_eps_m1={str(w): str(e) for w, e in Dm1.items()},
                        D_eps0={str(w): {b: str(e) for b, e in d.items()} for w, d in D0.items()},
                        ct1_pair_A={str(w): str(e) for w, e in ctA.items()},
                        bracket_S_ct1={str(w): str(e) for w, e in br.items()},
                        F_xi={str(w): str(e) for w, e in Fx.items()},
                        rhs_pole_identity={str(w): str(e) for w, e in rhs.items()},
                        lhs_minus_rhs=diff, lhs_plus_rhs=diff_alt,
                        Upsilon_eps0_plus_F={str(w): {b: str(sp.expand(e + (Fx.get(w, 0) if b == '1' else 0))) for b, e in d.items()}
                                             for w, d in D0.items()},
                        seconds=dict(D=round(dk['seconds'], 1), ctA=round(s1, 1), ka=round(s2, 1), kb=round(s3, 1)))
        gates.append(dict(id=f'G-pole {key}', value=diff, target='[eps^-1]D = {S,ct_1(O)} - ct_1({S,X_a}X_b) + ct_1(X_a{S,X_b}) (eq qs-pole, m = 1, leg set lam)',
                          pass_=ok, note=('with the overall sign of D reversed: ' + ('holds' if all(v == '0' for v in diff_alt.values()) else 'fails'))))
        print(key, 'PASS' if ok else 'FAIL', 'alt-sign', all(v == '0' for v in diff_alt.values()), flush=True)
    # ---- pole identity on (A_perp, lam_+), leg set A_perp (the pair of Example qs-cubic)
    t1 = time.time()
    dj = QJ.DJob('AW', 'lamP', ['AW'], 'xi', XIV, tag='fsp2')
    dk = dj.kernel(mlist=(0, 1))
    Dm1 = {w: to_expr(dj, dm.get(1, {}).get('1')) for w, dm in dk['total'].items()}
    # {S, X_a} for X_a = sum_K za_{4+K} A_{4+K}: the lam_+ letter with polarisation z'_c = sum_K za_{4+K} (xi^T MR_{4+K})_c
    xv = sp.Matrix(XIV)
    zprime = {c: sp.expand(sum(sp.Symbol(f'za_{4 + K}') * (xv.T * B.MR[4 + K])[c] for K in range(6))) for c in range(8)}
    job = T.Job('lamP', 'lamP', ['AW'])
    r = job.kernel(keep_sequences=False)
    ren = {f'za_{c}': f'zl_{c}' for c in range(8)}
    k1 = {}
    for w, d in r['total'].items():
        P = d.get(1)
        if P:
            e = to_expr(job, P, ren)
            k1[w] = sp.expand(e.subs({sp.Symbol(f'zl_{c}'): zprime[c] for c in range(8)}, simultaneous=True))
    # {S, X_b} for the lam_+ letter zb: sum c_MN d_M A_N (A_V); pair (A_perp, d A_V)
    kb, s3 = ct1_deriv_pair('AW', 'AV', 'AW', 'b')
    zb = [sp.Symbol(f'zb_{c}') for c in range(8)]
    cb = cMN(zb, 'lamP')
    # general colour: -ct_1(X_a {S,X_b}_cubic) = +(V/2) Acub (vacuum ct_1(A_perp A_perp) = V (za.zb), V = 2 from
    # tl_engine), -ct1_{xi,1}(sO) = -Z (xi.zb)(za.Y) (Z = 1/2), both on the word (0,)
    zbv = sp.zeros(16, 1)
    for c in range(8):
        zbv[c] = sp.Symbol(f'zb_{c}')
    Acub = 0
    for M in range(4, 10):
        for N in range(4, 10):
            if M != N:
                c3 = (zbv.T * B.Rbm[M] * B.Rm[N] * xv)[0]; c3r = (zbv.T * B.Rbm[N] * B.Rm[M] * xv)[0]
                Acub += sp.Rational(1, 2) * (c3 - c3r) * sp.Symbol(f'za_{M}') * sp.Symbol(f'o0_{N}')
    s1 = sp.expand(sum(sp.Symbol(f'za_{M}') * sp.Symbol(f'o0_{M}') for M in range(4, 10)) * sum(XIV[c] * sp.Symbol(f'zb_{c}') for c in range(8)))
    vac = T.Job('AW', 'AW', []).kernel(keep_sequences=False)['total'][()][1]
    Vv = sp.cancel(to_expr(T.Job('AW', 'AW', []), vac) / sum(sp.Symbol(f'za_{c}') * sp.Symbol(f'zb_{c}') for c in range(4, 10)))
    rhs = {}
    for w in set(k1) | set(kb) | set(Dm1):
        v = -k1.get(w, 0) - sum((kb.get(w, {}).get(mn, 0) * c for mn, c in cb.items()), 0)
        if w == (0,):
            v += sp.Rational(1, 2) * Vv * Acub - sp.Rational(1, 2) * s1
        rhs[w] = sp.expand(v)
    diff = {str(w): str(sp.expand(Dm1.get(w, 0) - rhs.get(w, 0))) for w in rhs}
    diff_alt = {str(w): str(sp.expand(Dm1.get(w, 0) + rhs.get(w, 0))) for w in rhs}
    ok = all(v == '0' for v in diff.values()) and any(rhs.get(w, 0) != 0 for w in rhs)
    res['AW|lamP -> AW'] = dict(D_eps_m1={str(w): str(e) for w, e in Dm1.items()}, rhs_pole_identity={str(w): str(e) for w, e in rhs.items()},
                                lhs_minus_rhs=diff, lhs_plus_rhs=diff_alt,
                                terms='rhs = -ct_1({S,X_a} X_b) - ct_1(X_a {S,X_b}_quad) - ct_1(X_a {S,X_b}_cubic) - ct1_xi(sO) ({S,O} = {S,X_a}X_b + X_a{S,X_b}, X_a even; ct_1(O) = 0; QO = 0 on A_perp, lam letters); vacuum ct_1(A_perp A_perp) = ' + str(Vv) + ' (za.zb), ct1_xi(c, lam_+)|_0 = (1/2)(xi.zb)',
                                seconds=round(time.time() - t1, 1))
    gates.append(dict(id='G-pole AW|lamP -> AW (Example qs-cubic pair)', value=diff,
                      target='[eps^-1]D(A_perp,lam_+)|_{A_perp} = -ct_1({S,A_perp} lam_+) - ct_1(A_perp {S,lam_+}) (eq qs-pole, m = 1)', pass_=ok,
                      note=('with the overall sign of D reversed: ' + ('holds' if all(v == '0' for v in diff_alt.values()) else 'fails'))))
    print('AW|lamP -> AW', 'PASS' if ok else 'FAIL', 'alt-sign', all(v == '0' for v in diff_alt.values()), flush=True)
    res['seconds'] = round(time.time() - t0, 1)
    return res, gates


if __name__ == '__main__':
    res, gates = run()
    json.dump(dict(results=res, gates=gates), open(os.path.join(OUT, 'fsp_gpole.json'), 'w'), indent=1)
    for g in gates:
        print(g['id'], 'PASS' if g['pass_'] else 'FAIL', g.get('note', ''))
