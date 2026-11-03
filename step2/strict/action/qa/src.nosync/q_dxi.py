"""QA third evaluation of Ob^(2)(Delta,xi)|_2 (engine A only so far, flagged by its G-WZ2):
defect cycles of I + c_D S(Delta) + c_xi S(xi) evaluated with engine B's verified cycle code (read-only import; engine B
never computed this pair), intrinsic part (canonical densities, last leg eliminated) and position part beta_mu
(= 2 Ob^(2)(P_{e_mu},xi)), eps^0 and poles.  Normalisation derived in qa/compare.json: theory = -(raw), since
c_D S(Delta) = -S(Delta) c_D and c_xi is even.  xi in S_+ (x) W."""
import sys, os, json, random, time, resource
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
BSRC = os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b', 'src.nosync'))
sys.path.insert(0, BSRC)
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_pos as P, eb_qb as Q
QA = os.path.abspath(os.path.join(HERE, '..'))
LOG = os.path.join(QA, 'PROGRESS.txt')

def log(s):
    line = time.strftime('%Y-%m-%d %H:%M:%S') + ' q_dxi' + CHIR + ': ' + s
    print(line, flush=True)
    with open(LOG, 'a') as f:
        f.write(line + '\n')

LAMP = list(range(8)); LAMM = list(range(8, 16)); LAMS_M = list(range(8)); LAMS_P = list(range(8, 16))
CHIR = sys.argv[1] if len(sys.argv) > 1 else '+'
SECT = {  # name: (sp0, rng0, sp1, rng1, basis labels)
    'lam-_AsV': ('lam', LAMM, 'As', range(4), ['T1']),
    'lam+_Asp': ('lam', LAMP, 'As', range(4, 10), ['T2']),
    'AV_lams+': ('A', range(4), 'lams', LAMS_P, ['T0', 'T3']),
    'Ap_lams-': ('A', range(4, 10), 'lams', LAMS_M, ['T4']),
} if CHIR == '+' else {
    'lam+_AsV': ('lam', LAMP, 'As', range(4), ['T1']),
    'lam-_Asp': ('lam', LAMM, 'As', range(4, 10), ['T2']),
    'AV_lams-': ('A', range(4), 'lams', LAMS_M, ['T0', 'T3']),
    'Ap_lams+': ('A', range(4, 10), 'lams', LAMS_P, ['T4']),
}

def rv(n, lo=-3, hi=3):
    return [Fr(random.randint(lo, hi)) for _ in range(n)]

def make_config(sec):
    s0, r0, s1, r1, _ = SECT[sec]
    p = rv(4)
    u0 = np.zeros(E.V.DIM[s0], complex); u1 = np.zeros(E.V.DIM[s1], complex)
    for c in r0:
        u0[c] = random.randint(-2, 2)
    for c in r1:
        u1[c] = random.randint(-2, 2)
    def legs_at(r):
        return [dict(sp=s0, u=u0, p=p), dict(sp=s1, u=u1, p=[-p[m] - r[m] for m in range(4)])]
    return legs_at

def cl(v):
    return v.to_json()

def main(nextra=2, seed=777 if CHIR == '+' else 778):
    random.seed(seed)
    tdir = G.tmpl_dir(); tw = G.tmpl_weights(2, 3, 0, 'D'); tq, t3 = G.tmpl_xi()
    def gdir(mu, r):
        e = np.zeros(4, complex); e[mu] = 1
        return Y.Gen('Delta', 1, tdir, gvec=e, inj=r)
    g0 = Y.Gen('Delta', 1, tw)
    out = {}
    fn = os.path.join(QA, 'out.nosync', 'q_dxi_raw.json' if CHIR == '+' else 'q_dxi_raw_minus.json')
    if os.path.exists(fn):
        out = json.load(open(fn))
    for sec, (s0, r0, s1, r1, labels) in SECT.items():
        if sec in out:
            log(f'{sec}: cached'); continue
        recs = []
        for ic in range(len(labels) + nextra):
            xi = np.zeros(16, complex)
            for c in (range(8) if CHIR == '+' else range(8, 16)):
                xi[c] = random.randint(-2, 2)
            legs_at = make_config(sec)
            gx = Y.Gen('xi', 0, tq, gvec=xi, t3=t3)
            t0 = time.time(); st = {}
            al, be, top = P.position_split(legs_at, gdir, g0, extra_gens=[gx], stats=st, gen_first=True)
            legs0 = legs_at([Fr(0)] * 4)
            Tb = Q.xi_terms(xi)
            bvals = {l: [str(x) for x in Q.to_fr(Tb[l](legs0))] for l in labels}
            rec = dict(xi=[int(x.real) for x in xi], legs=[dict(sp=l['sp'], u=[int(z.real) for z in l['u']],
                                                               p=[str(x) for x in l['p']]) for l in legs0],
                       alpha={k: {str(m): cl(v) for m, v in d.items()} for k, d in al.items()},
                       beta={str(mu): {k: {str(m): cl(v) for m, v in d.items()} for k, d in dd.items()} for mu, dd in be.items()},
                       top_zero=all(v.is_zero() for d in top.values() for v in d.values()), basis=bvals,
                       chains=st.get('chains', 0), seconds=round(time.time() - t0, 1))
            recs.append(rec)
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
            log(f'{sec} config {ic}: {rec["seconds"]} s, chains {rec["chains"]}, top_zero {rec["top_zero"]}, RSS {rss:.0f} MB')
        out[sec] = recs
        with open(fn, 'w') as f:
            json.dump(out, f)
    log('done')

if __name__ == '__main__':
    main()
