"""Task B1: Ob^(2)(Delta,P_b) on two legs (degree -1, weight 5, p <= 1) by the defect cycles of
I + c_D S(Delta) + c_P S(P_b): raw coefficient of c_D c_P (generator order Delta, P_b), intrinsic part (canonical
densities, last leg eliminated) and position part, eps^0 and poles.  Theory normalisation: Ob^(2) = -(raw) (notes.txt)."""
import json, os, random, time, sys
from fractions import Fraction as Fr
import numpy as np
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_pos as P

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')
LAMP = list(range(8)); LAMM = list(range(8, 16)); LAMS_M = list(range(8)); LAMS_P = list(range(8, 16))
SECT = {'c_cs': ('c', 'cs', None, None), 'AV_AsV': ('A', 'As', range(4), range(4)),
        'Ap_Asp': ('A', 'As', range(4, 10), range(4, 10)), 'lam+_lams+': ('lam', 'lams', LAMP, LAMS_P),
        'lam-_lams-': ('lam', 'lams', LAMM, LAMS_M)}


def rv(n, lo=-3, hi=3):
    return [Fr(random.randint(lo, hi)) for _ in range(n)]


def make_config(sec):
    s0, s1, r0, r1 = SECT[sec]
    p = rv(4)
    u0 = np.zeros(E.V.DIM[s0], complex); u1 = np.zeros(E.V.DIM[s1], complex)
    if s0 == 'c':
        u0[0] = 1; u1[0] = 1
    else:
        for c in r0:
            u0[c] = random.randint(-2, 2)
        for c in r1:
            u1[c] = random.randint(-2, 2)
    b = rv(4, -2, 2)

    def legs_at(r):
        return [dict(sp=s0, u=u0, p=p), dict(sp=s1, u=u1, p=[-p[m] - r[m] for m in range(4)])]
    return legs_at, b, dict(u0=[str(complex(x)) for x in u0], u1=[str(complex(x)) for x in u1],
                            p=[str(x) for x in p], b=[str(x) for x in b])


def clog_json(v):
    return v.to_json()


def dump_vals(d):
    return {k: {str(m): clog_json(v) for m, v in dm.items()} for k, dm in d.items()}


def run_sector(sec, nconf, seed, log=print):
    random.seed(seed)
    tdir = G.tmpl_dir()
    tw = G.tmpl_weights(2, 3, 0, 'D')

    def gdir(mu, r):
        e = np.zeros(4, complex)
        e[mu] = 1
        return Y.Gen('Delta', 1, tdir, gvec=e, inj=r)
    g0 = Y.Gen('Delta', 1, tw)
    out = []
    for ic in range(nconf):
        legs_at, b, desc = make_config(sec)
        gP = Y.Gen('P', 1, tdir, gvec=np.array([complex(x) / 2 for x in b]))
        t0 = time.time()
        st = {}
        al, be, top = P.position_split(legs_at, gdir, g0, extra_gens=[gP], stats=st, gen_first=True)
        rec = dict(config=desc, alpha=dump_vals(al), beta={str(mu): dump_vals(d) for mu, d in be.items()},
                   top_zero=all(v.is_zero() for dm in top.values() for v in dm.values()),
                   chains=st.get('chains', 0), seconds=round(time.time() - t0, 1))
        out.append(rec)
        log(f'B1 {sec} config {ic}: chains {rec["chains"]} {rec["seconds"]} s top_zero {rec["top_zero"]}')
    return out


if __name__ == '__main__':
    secs = sys.argv[1:] or list(SECT)
    for sec in secs:
        fn = os.path.join(OUT, f'B1_raw_{sec}.json')
        if os.path.exists(fn):
            print('skip', sec)
            continue
        n = {'c_cs': 4, 'AV_AsV': 8, 'Ap_Asp': 4, 'lam+_lams+': 7, 'lam-_lams-': 7}[sec]
        r = run_sector(sec, n, seed=101 + list(SECT).index(sec))
        with open(fn, 'w') as f:
            json.dump(r, f)
