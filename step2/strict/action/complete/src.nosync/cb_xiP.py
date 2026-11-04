"""Route B: Ob^(2)(xi, P_b) on two legs, defect cycles of I + c_xi S(xi) + c_P S(P_b); theory = (-1)^{n_odd legs}
x raw (c_P in front, c^P on the right: c_P S(P) = -S(P) c^P; engine A rule Ob2(x,y) = (-1)^{|x|} coefficient of
theta.. c^y c^x); n_odd is even in every sector, so theory = raw."""
import sys, os, json, random, time, resource
from fractions import Fraction as Fr
import numpy as np
import cb_common as B, cb_pairs as CP
E, Y, G = B.E, B.Y, B.G
OUT = os.path.join(B.COMPLETE, 'out.nosync')
SECT = {'+': {'AV_ls+': ('A', range(4), 'lams', B.LAMS_P), 'Ap_ls-': ('A', range(4, 10), 'lams', B.LAMS_M),
              'l+_Asp': ('lam', B.LAMP, 'As', range(4, 10)), 'l-_AsV': ('lam', B.LAMM, 'As', range(4))},
        '-': {'AV_ls-': ('A', range(4), 'lams', B.LAMS_M), 'Ap_ls+': ('A', range(4, 10), 'lams', B.LAMS_P),
              'l-_Asp': ('lam', B.LAMM, 'As', range(4, 10)), 'l+_AsV': ('lam', B.LAMP, 'As', range(4))}}


def run(chir, nconf=8, seed=303):
    random.seed(seed + (0 if chir == '+' else 7))
    tq, t3 = G.tmpl_xi(); tdir = G.tmpl_dir()
    fn = os.path.join(OUT, f'B_xiP_{"p" if chir == "+" else "m"}_raw.json')
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    for sec, (s0, r0, s1, r1) in SECT[chir].items():
        recs = out.get(sec, [])
        while len(recs) < nconf:
            xi = CP.rspin(CP.SPIN[chir])
            b = np.array([random.randint(-2, 2) for _ in range(4)], complex)
            while not np.any(b):
                b = np.array([random.randint(-2, 2) for _ in range(4)], complex)
            p = [Fr(random.randint(-3, 3)) for _ in range(4)]
            while not any(p):
                p = [Fr(random.randint(-3, 3)) for _ in range(4)]
            legs = [dict(sp=s0, u=CP.rpol(s0, r0), p=p), dict(sp=s1, u=CP.rpol(s1, r1), p=[-x for x in p])]
            t0 = time.time(); st = {}
            raw = Y.evaluate(legs, [Y.Gen('xi', 0, tq, gvec=xi, t3=t3), Y.Gen('P', 1, tdir, gvec=b / 2)], mlist=(0, 1), stats=st)
            recs.append(dict(xi=[[z.real, z.imag] for z in xi], b=[z.real for z in b],
                             legs=[dict(sp=l['sp'], u=[[z.real, z.imag] for z in l['u']], p=[str(x) for x in l['p']]) for l in legs],
                             raw=CP.cl(raw), chains=st.get('chains', 0), sec=round(time.time() - t0, 2)))
            out[sec] = recs
            json.dump(out, open(fn, 'w'))
        B.log(f'(xi{chir},P_b) sector {sec}: {len(recs)} configs (last {recs[-1]["sec"]} s, chains {recs[-1]["chains"]})')


if __name__ == '__main__':
    for c in sys.argv[1:] or ['+', '-']:
        run(c)
