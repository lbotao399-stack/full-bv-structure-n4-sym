"""Route B: Ob^(2)(xi1, xi2) on two legs for odd pairs (xi1 in S_+ (x) W, xi2 in S_- (x) W^* or S_+ (x) W), defect
cycles of J = I + c1 S(xi1) + c2 S(xi2) - c1 c2 S^(2)(xi1,xi2) (commuting CE coordinates of odd generators; generating
coefficient a_2 = -1/2 of Def c2-cochains / ch08), theory component = -(c1 c2 coefficient) = -raw(two letters) +
raw(S^(2) letter).  eps^0 (m = 0) and poles (m = 1).  Raw values stored per configuration."""
import sys, os, json, random, time, resource
from fractions import Fraction as Fr
import numpy as np
import cb_common as B, cb_s2 as S
E, Y, G = B.E, B.Y, B.G
OUT = os.path.join(B.COMPLETE, 'out.nosync')
SPIN = {'+': B.LAMP, '-': B.LAMM}
SECTORS = {
    'pm': {'c_cs': ('c', None, 'cs', None), 'AV_AsV': ('A', range(4), 'As', range(4)),
           'Ap_Asp': ('A', range(4, 10), 'As', range(4, 10)), 'l+_ls+': ('lam', B.LAMP, 'lams', B.LAMS_P),
           'l-_ls-': ('lam', B.LAMM, 'lams', B.LAMS_M)},
    'pp': {'AV_Asp': ('A', range(4), 'As', range(4, 10)), 'Ap_AsV': ('A', range(4, 10), 'As', range(4)),
           'l+_ls-': ('lam', B.LAMP, 'lams', B.LAMS_M), 'l-_ls+': ('lam', B.LAMM, 'lams', B.LAMS_P)},
}
SECTORS['mm'] = dict(SECTORS['pp'])
CHIR = {'pm': ('+', '-'), 'pp': ('+', '+'), 'mm': ('-', '-')}


def rspin(rng):
    x = np.zeros(16, complex)
    while not np.any(x):
        for c in rng:
            x[c] = random.randint(-2, 2)
    return x


def rpol(spc, rng):
    u = np.zeros(E.V.DIM[spc], complex)
    while not np.any(u):
        for c in (rng if rng is not None else range(E.V.DIM[spc])):
            u[c] = random.randint(-2, 2)
    return u


def cl(d):
    return {k: {str(m): v.to_json() for m, v in dm.items()} for k, dm in d.items()}


def run(pair, nconf=8, seed=101):
    random.seed(seed + {'pm': 0, 'pp': 1000, 'mm': 2000}[pair])
    tq, t3 = G.tmpl_xi()
    fn = os.path.join(OUT, f'B_pair_{pair}_raw.json')
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    c1, c2 = CHIR[pair]
    for sec, (s0, r0, s1, r1) in SECTORS[pair].items():
        recs = out.get(sec, [])
        while len(recs) < nconf:
            xi1, xi2 = rspin(SPIN[c1]), rspin(SPIN[c2])
            p = [Fr(random.randint(-3, 3)) for _ in range(4)]
            while not any(p):
                p = [Fr(random.randint(-3, 3)) for _ in range(4)]
            u0, u1 = rpol(s0, r0), rpol(s1, r1)
            legs = [dict(sp=s0, u=u0, p=p), dict(sp=s1, u=u1, p=[-x for x in p])]
            t0 = time.time(); st = {}
            raw2 = Y.evaluate(legs, [Y.Gen('x1', 0, tq, gvec=xi1, t3=t3), Y.Gen('x2', 0, tq, gvec=xi2, t3=t3)],
                              mlist=(0, 1), stats=st)
            tS2 = E.build_template(S.s2_terms(xi1, xi2, f'{pair}{sec}{len(recs)}'), [], extra_odd=0, name=None)
            raw1 = Y.evaluate(legs, [Y.Gen('S2', 0, tS2)], mlist=(0, 1), stats=st)
            recs.append(dict(xi1=[[z.real, z.imag] for z in xi1], xi2=[[z.real, z.imag] for z in xi2],
                             legs=[dict(sp=l['sp'], u=[[z.real, z.imag] for z in l['u']], p=[str(x) for x in l['p']]) for l in legs],
                             raw2=cl(raw2), raw1=cl(raw1), chains=st.get('chains', 0), sec=round(time.time() - t0, 2)))
            out[sec] = recs
            json.dump(out, open(fn, 'w'))
        B.log(f'pair {pair} sector {sec}: {len(recs)} configs (last {recs[-1]["sec"]} s, chains {recs[-1]["chains"]}), '
              f'RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB')
    return out


if __name__ == '__main__':
    E.start_watchdog(1400)
    for pair in sys.argv[1:] or ['pm', 'pp']:
        run(pair)
