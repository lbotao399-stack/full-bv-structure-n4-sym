"""QA second evaluation of Ob^(1)(eta_+)|_2 (engine A only so far, task T6) with engine B's cycle code (read-only).
S(eta) = sum_mu x^mu S(rho(e_mu) eta) + int s0, s0 = -(rho(A_perp) eta, lam*) (c11-Seta); rho(e_M): Sigma* -> Sigma is
RB[M] in engine B's realisation.  The sign sigma of s0 relative to the x^mu part is fixed here by Q-closedness of the
quadratic part (Q S_2(eta) = 0), not assumed.  Odd generator: theory = raw (engine B arity-one rule)."""
import sys, os, json, random, time, resource
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
BSRC = os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b', 'src.nosync'))
sys.path.insert(0, BSRC)
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_pos as P, eb_qb as Q
import ce_core as C
QA = os.path.abspath(os.path.join(HERE, '..'))
LOG = os.path.join(QA, 'PROGRESS.txt')
def log(s):
    line = time.strftime('%Y-%m-%d %H:%M:%S') + ' q_eta: ' + s
    print(line, flush=True)
    open(LOG, 'a').write(line + '\n')
R, RB, MP = Q.R, Q.RB, Q.MP
LAMP = list(range(8)); LAMM = list(range(8, 16)); LAMS_M = list(range(8)); LAMS_P = list(range(8, 16))

def s0_density(eta, sigma):
    return Q.dens(-sigma, ('A', 'lams'), lambda s: sum(s[0]['u'][I] * ((RB[I] @ eta) @ MP @ s[1]['u']) for I in range(4, 10)))

def quad_xi(xi):
    T = Q.xi_terms(xi)
    return Q.dsum(*[(1, T[k]) for k in ('T1', 'T2', 'T3', 'T4')])

def closed_test(sigma, ntrial=3):
    worst = 0.0
    for spcs in (['A', 'lam'], ['lam', 'A'], ['c', 'lams'], ['A', 'lams'], ['lam', 'As']):
        for _ in range(ntrial):
            eta = np.zeros(16, complex)
            for c in LAMS_P:
                eta[c] = random.randint(-2, 2)
            legs = [Q.rand_leg(s) for s in spcs]
            p0 = legs[0]['p']
            def legs_at(r, legs=legs, p0=p0):
                l0 = dict(legs[0]); l1 = dict(legs[1]); l1['p'] = [-p0[m] - r[m] for m in range(4)]
                return [l0, l1]
            fn = (s0_density(eta, sigma), [quad_xi(RB[mu] @ eta) for mu in range(4)])
            al, be = Q.canon(Q.Qfn(fn), legs_at)
            worst = max(worst, abs(float(al[0])) + abs(float(al[1])), *[abs(float(b[0])) + abs(float(b[1])) for b in be])
    return worst

def main():
    random.seed(4242)
    res = {}
    res['closedness'] = {str(sg): closed_test(sg) for sg in (1, -1)}
    log(f'Q S_2(eta) residual by sign of s0: {res["closedness"]}')
    good = [sg for sg in (1, -1) if res['closedness'][str(sg)] == 0]
    if len(good) != 1:
        log('no unique closed sign; abort'); json.dump(res, open(os.path.join(QA, 'out.nosync', 'q_eta_raw.json'), 'w')); return
    sigma = good[0]
    ETA = sp.symbols('qeta0:16')
    ev = sp.Matrix(ETA)
    def L_s0(s):
        a, ls = s
        w = sp.zeros(16, 1)
        for I in range(4, 10):
            w += a['u'][I] * (C.Rbm[I] * ev)
        return (w.T * C.Mp * sp.Matrix(ls['u']))[0]
    E.V.TERMS['QA_S0_eta'] = (sp.Integer(-sigma), ('A', 'lams'), L_s0, ('k', (0, 1)))
    t_s0 = E.build_template(['QA_S0_eta'], [], extra_odd=0, gen_syms=ETA, name=None)
    tq, t3 = G.tmpl_xi()
    SECT = {'AV_lam+': ('A', range(4), 'lam', LAMP, 2), 'Ap_lam-': ('A', range(4, 10), 'lam', LAMM, 1),
            'c_lams-': ('c', None, 'lams', LAMS_M, 1)}
    res['sigma'] = sigma
    res['sectors'] = {}
    for sec, (s0, r0, s1, r1, nb) in SECT.items():
        recs = []
        for ic in range(nb + 2):
            eta = np.zeros(16, complex)
            for c in LAMS_P:
                eta[c] = random.randint(-2, 2)
            p = [Fr(random.randint(-3, 3)) for _ in range(4)]
            u0 = np.zeros(E.V.DIM[s0], complex); u1 = np.zeros(E.V.DIM[s1], complex)
            if s0 == 'c':
                u0[0] = 1
            else:
                for c in r0:
                    u0[c] = random.randint(-2, 2)
            for c in r1:
                u1[c] = random.randint(-2, 2)
            def legs_at(r, p=p, u0=u0, u1=u1, s0=s0, s1=s1):
                return [dict(sp=s0, u=u0, p=p), dict(sp=s1, u=u1, p=[-p[m] - r[m] for m in range(4)])]
            def gdir(mu, r, eta=eta):
                return Y.Gen('xi', 0, tq, gvec=RB[mu] @ eta, t3=t3, inj=r)
            g0 = Y.Gen('s0', 0, t_s0, gvec=eta)
            t0 = time.time()
            al, be, top = P.position_split(legs_at, gdir, g0, stats={})
            recs.append(dict(eta=[int(x.real) for x in eta], p=[str(x) for x in p], u0=[int(x.real) for x in u0],
                             u1=[int(x.real) for x in u1],
                             alpha={k: {str(m): v.to_json() for m, v in d.items()} for k, d in al.items()},
                             beta={str(mu): {k: {str(m): v.to_json() for m, v in d.items()} for k, d in dd.items()} for mu, dd in be.items()},
                             top_zero=all(v.is_zero() for d in top.values() for v in d.values())))
            log(f'{sec} config {ic}: {time.time() - t0:.1f} s, RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB')
        res['sectors'][sec] = recs
        json.dump(res, open(os.path.join(QA, 'out.nosync', 'q_eta_raw.json'), 'w'))
    log('done')

if __name__ == '__main__':
    main()
