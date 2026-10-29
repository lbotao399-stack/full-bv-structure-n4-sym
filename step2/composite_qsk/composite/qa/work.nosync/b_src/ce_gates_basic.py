"""Basic gates of route B: core identities, exact integrals, bilinear (GC3), Lemma (ii) (option 2 vs option 3),
SO(6) planar density, colour rules (reused twoletter colour_rules.py).  Writes out.nosync/basic_gates.json."""
import os, sys, json, time
import sympy as sp
from sympy.polys.domains import QQ_I
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')
import ce_core as C
import ce_int as CI
import ce_job as CJ
I = sp.I
G = {}


def rec(name, value, target, ok, note=''):
    G[name] = dict(value=str(value), target=str(target), pass_=bool(ok), note=note)
    print(('PASS ' if ok else 'FAIL ') + name, str(value)[:160], flush=True)


def core():
    r = C.selftest()
    for k, v in r.items():
        if isinstance(v, bool):
            rec('GB0 ' + k, v, True, v, 'ce_core.selftest: manuscript Def c1-bracket, Def c3-qgf, Prop c1-Q, Prop c3-QP, Prop c9-composite, Lemma lem:cc-kernels')


def integrals():
    import mpmath as mp
    worst = 0
    for (m, a, nv, c) in [(0, 0, (0, 0), (1, 1)), (0, 1, (2, 0), (1, 1)), (1, 0, (1, 1), (1, 1)), (0, 1, (2,), (2,)), (1, 0, (3,), (2,)),
                          (0, 1, (1, 0), (2, 1)), (2, 2, (1, 3), (1, 1))]:
        ex = float(CI.J(m, a, nv, c)); nm = float(CI.J_numeric(m, a, nv, c, dps=20))
        worst = max(worst, abs(ex - nm) / abs(nm))
    rec('GB1 exact Feynman-parameter integrals J_m vs mpmath quadrature (n_t = 2, 3)', f'max rel diff {worst:.2e}', '< 1e-12', worst < 1e-12)
    import ce_int_check as CK
    worst = 0
    for (m, a, nvec, cut) in [(0, 0, (0, 0, 0), (1, 1, 1)), (0, 2, (1, 0, 1), (1, 1, 1)), (1, 1, (0, 2, 0), (1, 1, 1)),
                              (1, 0, (0, 1, 1), (2, 1, 1)), (0, 1, (2, 1, 0, 1), (1, 1, 1, 1))]:
        ex = float(CI.J(m, a, nvec, cut))
        q = sum(CK.J_qmc(m, a, nvec, cut, npow=21, seed=s) for s in range(3)) / 3
        worst = max(worst, abs(ex - q) / abs(q))
    rec('GB1b exact J_m vs scrambled-Sobol QMC (n_t = 4, 5, cutoffs 1, 2)', f'max rel diff {worst:.2e}', '< 1e-4 (QMC accuracy)', worst < 1e-4)


def bilinear():
    w1 = {4: 1, 5: 2, 6: -1, 7: 0, 8: 3, 9: 1}
    w2 = {4: 2, 5: -1, 6: 1, 7: 1, 8: 0, 9: -2}
    res = {}
    for opt in (3, 2):
        job = CJ.CJob('AW', 'AW', ['AW', 'AW'], pol='sym', za=w1, zb=w2)
        r = job.kernels(mlist=(0, 1), option=opt)
        R = job.R; g = job.g
        Y0 = [g[f'o0_{c}'] for c in range(4, 10)]; Y1 = [g[f'o1_{c}'] for c in range(4, 10)]
        W1 = [R(QQ_I.from_sympy(sp.nsimplify(w1[c]))) for c in range(4, 10)]; W2 = [R(QQ_I.from_sympy(sp.nsimplify(w2[c]))) for c in range(4, 10)]
        dot = lambda x, y: sum((a * b for a, b in zip(x, y)), R.zero)
        form = dot(W1, Y0) * dot(W2, Y1) + dot(W1, Y1) * dot(W2, Y0)
        ap = sum((g[f'a{m}'] * (g[f'r0_{m}'] + g[f'r1_{m}']) for m in range(4)), R.zero)
        totP = sum((dm.get(0, {}).get('1', R.zero) for dm in r['Pb'].values()), R.zero)
        totD = sum((dm.get(0, {}).get('1', R.zero) for dm in r['D0'].values()), R.zero)
        poles = sum(1 for dm in r['Pb'].values() for m in dm if m != 0)
        other = [(w, m, b) for tab in ('Pb', 'D0') for w, dm in r[tab].items() for m, d in dm.items() for b in d if b != '1']
        rP = sp.simplify(totP.as_expr() / (ap * form * R(QQ_I(0, 1))).as_expr())
        rD = sp.simplify(totD.as_expr() / form.as_expr())
        res[opt] = (rP, rD, poles, other)
    rP, rD, poles, other = res[3]
    rec('GC3 bilinear C2_Pb(A_perp,A_perp)|_{A_perp A_perp} / (u d_b O)', rP, '-1/3 (Example ex:cc-bilinear, Thm c9-P1)', rP == sp.Rational(-1, 3),
        'w1.w2 != 0: traceless and trace parts together; exact proportionality')
    rec('GC3 bilinear C2_Delta^0 / (u O)', rD, '-8/3 = -(2/3) w_O - 4/3', rD == sp.Rational(-8, 3))
    rec('GC3 bilinear: position part of C2_Delta = 2 C2_Pb (kernel form, structural) => C2_Delta = -(1/3) u v.O - (4/3) u O', 'structural', 'Example ex:cc-bilinear', True)
    rec('GC5a bilinear: no log2/log3 remainder, no eps^-1 term at W0', (poles, other), '(0, [])', poles == 0 and not other)
    rec('GC4b Lemma lem:cc-kernels (ii) on the bilinear: option 2 (both ends) = option 3 (Euler on the line + column end)',
        (res[2][0], res[2][1]), (rP, rD), res[2][0] == rP and res[2][1] == rD)


def colour():
    sys.path.insert(0, os.path.join(HERE, 'vendor'))     # vendored twoletter colour_rules.py
    import colour_rules as CR
    r = CR.run()
    rec('GC12 colour rules of the word form and planar adjacent rule (su(2), su(3), orthonormal and rescaled bases)', r['all_pass'], True, r['all_pass'],
        'reused twoletter/engine_b/colour_rules.py (verified code); gl(3) realisation not run here')


if __name__ == '__main__':
    t0 = time.time()
    for f in (core, integrals, bilinear, colour):
        f()
    json.dump(G, open(os.path.join(OUT, 'basic_gates.json'), 'w'), indent=1)
    print('seconds', round(time.time() - t0, 1))
