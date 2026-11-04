"""Route A: Ob^(2)(xi1, xi2) on two legs (odd pairs), engine A's closed defect cycles (sa_core.cycle_sum, GSimplex):
two-letter cycles with S(xi1), S(xi2) and the S^(2) block; theory component = -(two letters) + (S^(2) cycles)
(commuting CE coordinates, generating coefficient -1/2 S^(2); engine A rule (-1)^{|x|} for the two-letter part).
Configurations: the SAME legs, momenta and spinors as route B (cb_pairs raw files), mapped to engine A's realisation
is impossible (different Clifford matrices), so route A draws its own random data in its realisation."""
import sys, os, json, time, resource
import numpy as np
import ca_common as A, ca_s2 as S2
SA, SO = A.SA, A.SO
OUT = os.path.join(A.COMPLETE, 'out.nosync')
SECTORS = {
    'pm': {'c_cs': ('c', 'cs'), 'AV_AsV': ('AV', 'AsV'), 'Ap_Asp': ('AW', 'AsW'), 'l+_ls+': ('lamP', 'lamsP'),
           'l-_ls-': ('lamM', 'lamsM')},
    'pp': {'AV_Asp': ('AV', 'AsW'), 'Ap_AsV': ('AW', 'AsV'), 'l+_ls-': ('lamP', 'lamsM'), 'l-_ls+': ('lamM', 'lamsP')},
}
SECTORS['mm'] = dict(SECTORS['pp'])
CHIR = {'pm': ('lamP', 'lamM'), 'pp': ('lamP', 'lamP'), 'mm': ('lamM', 'lamM')}


def rngs(seed):
    rng = np.random.default_rng(seed)

    def rp():
        while True:
            p = rng.integers(-3, 4, size=4)
            if p @ p:
                return p

    def rpol(t):
        Bs = SA.E.BASIS[t]
        while True:
            c = rng.integers(-2, 3, size=len(Bs))
            if c.any():
                return sum(ci * np.asarray(v) for ci, v in zip(c, Bs)).astype(complex)
    return rp, rpol


def v5(d):
    return [[str(x[0]), str(x[1])] for x in d]


def run(pair, nconf=8, seed=2101):
    rp, rpol = rngs(seed + {'pm': 0, 'pp': 1, 'mm': 2}[pair])
    fn = os.path.join(OUT, f'A_pair_{pair}_raw.json')
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    c1, c2 = CHIR[pair]
    for sec, (t0_, t1_) in SECTORS[pair].items():
        recs = out.get(sec, [])
        while len(recs) < nconf:
            xi1, xi2 = rpol(c1), rpol(c2)
            p = rp()
            L = [SA.leg(t0_, rpol(t0_), p), SA.leg(t1_, rpol(t1_), -p)]
            t0 = time.time()
            what = ('fin', ('pole', 1))
            r2 = SA.cycle_sum(L, [SA.gen_odd(xi1, 'x1'), SA.gen_odd(xi2, 'x2')], slotnames=['x1', 'x2'], what=what)
            r1 = SA.cycle_sum(L, [], s2=S2.gen_s2(xi1, xi2), what=what)
            rec = dict(xi1=[[z.real, z.imag] for z in xi1], xi2=[[z.real, z.imag] for z in xi2],
                       legs=[dict(t=l['t'], pol=[[z.real, z.imag] for z in l['pol']], p0=[int(x.real) for x in l['p0']]) for l in L])
            for nm, r in (('raw2', r2), ('raw1', r1)):
                rec[nm] = {}
                for w, wn in (('fin', '0'), (('pole', 1), '1')):
                    rec[nm][wn] = v5(SO.reduce_colour(r, 2, w, 0).get('(Y0,Y1)', SA.vzero()))
            rec['sec'] = round(time.time() - t0, 2)
            recs.append(rec)
            out[sec] = recs
            json.dump(out, open(fn, 'w'))
        A.log(f'pair {pair} sector {sec}: {len(recs)} configs (last {recs[-1]["sec"]} s), RSS '
              f'{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB')


if __name__ == '__main__':
    for pair in sys.argv[1:] or ['pm', 'pp', 'mm']:
        run(pair)
