"""Task T3: Ob^(1)(R_r) (gate G-R0, two and three legs) and Ob^(1)(M_omega) on two legs (intrinsic part, position
part along b = Ob^(1)(P_{omega b}), the spin direction).  Even generators: theory component = -(raw)."""
import json, os, random, time
from fractions import Fraction as Fr
import numpy as np
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_pos as P, eb_qb as Q
import eb_t1

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')


def rv(n, lo=-3, hi=3):
    return [Fr(random.randint(lo, hi)) for _ in range(n)]


def allzero(d):
    return all(v.is_zero() for dm in d.values() for v in dm.values())


def run_R(seed=21, log=print):
    random.seed(seed)
    tR = G.tmpl_R()
    out = {}
    secs2 = ['AVAV', 'ApAp', 'lamlam', 'cAsV']
    for sec in secs2:
        zs = []
        for ic in range(2):
            r = np.array([random.randint(-2, 2) for _ in range(15)], complex)
            legs = eb_t1.make_legs(sec)([Fr(0)] * 4)
            res = Y.evaluate(legs, [Y.Gen('R', 1, tR, gvec=r)], mlist=(0, 1))
            zs.append(allzero(res))
        out[sec] = all(zs)
        log(f'T3 R_r two legs {sec}: all cycle sums zero = {out[sec]}')
    # three legs
    for name, sps in [('c_Ap_Asp', [('c', None), ('A', range(4, 10)), ('As', range(4, 10))]),
                      ('AV_Ap_Ap', [('A', range(4)), ('A', range(4, 10)), ('A', range(4, 10))]),
                      ('Ap_lam_lam', [('A', range(4, 10)), ('lam', None), ('lam', None)])]:
        p1, p2 = rv(4), rv(4)
        ps = [p1, p2, [-(p1[m] + p2[m]) for m in range(4)]]
        legs = []
        for (spc, rg), p in zip(sps, ps):
            u = np.zeros(E.V.DIM[spc], complex)
            for c in (rg if rg is not None else range(E.V.DIM[spc])):
                u[c] = random.randint(-2, 2)
            legs.append(dict(sp=spc, u=u, p=p))
        r = np.array([random.randint(-2, 2) for _ in range(15)], complex)
        t0 = time.time()
        res = Y.evaluate(legs, [Y.Gen('R', 1, tR, gvec=r)], mlist=(0, 1))
        out[name] = allzero(res)
        log(f'T3 R_r three legs {name}: all cycle sums zero = {out[name]} ({time.time() - t0:.0f} s)')
    return out


def omega_fn(om6):
    om = np.zeros((4, 4), complex)
    for v, (i, j) in zip(om6, G.VPAIRS):
        om[i, j] = v; om[j, i] = -v
    return om


def run_M(seed=23, log=print):
    random.seed(seed)
    tdir = G.tmpl_dir()
    tsp = G.tmpl_Mspin()
    B = __import__('eb_fit').basis_delta2()
    out = {}
    for sec in ['AVAV', 'ApAp', 'lamlam', 'cAsV']:
        recs = []
        for ic in range(3):
            om6 = [random.randint(-2, 2) for _ in range(6)]
            om = omega_fn(om6)

            def gdir(rho, r, om=om):
                a = np.array([om[nu, rho] for nu in range(4)], complex)   # v^nu = omega_{nu rho} x^rho
                return Y.Gen('M', 1, tdir, gvec=a / 2, inj=r)
            g0 = Y.Gen('M', 1, tsp, gvec=np.array(om6, complex))
            legs_at = eb_t1.make_legs(sec)
            al, be, top = P.position_split(legs_at, gdir, g0)
            legs0 = legs_at([Fr(0)] * 4)
            # candidate: Q int omega_{mu nu} (A_mu, A*_nu)
            dens = Q.dens(1, ('A', 'As'), lambda s, om=om: s[0]['u'][:4] @ om @ s[1]['u'][:4])
            qd = Q.apply_Q(dens)
            cand = Q.to_fr(qd(legs0))
            v0 = al.get('k', {}).get(0, E.CLog())
            v1 = al.get('k', {}).get(1, E.CLog())
            recs.append(dict(raw_eps0=v0.to_json(), raw_pole=v1.to_json(), cand=[str(cand[0]), str(cand[1])],
                             beta_zero=all(allzero(d) for d in be.values()), top_zero=allzero(top)))
        out[sec] = recs
        log(f'T3 M_omega {sec}: {recs}')
    return out


if __name__ == '__main__':
    r = dict(R=run_R(), M=run_M())
    with open(os.path.join(OUT, 'T3.json'), 'w') as f:
        json.dump(r, f, indent=1)
