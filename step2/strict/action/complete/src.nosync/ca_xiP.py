"""Route A: Ob^(2)(xi, P_b) on two legs, engine A cycles with gens [xi, P_b] (ghost c^P on the right, after the odd
legs); theory = (-1)^{|xi|} x coefficient = -(raw)."""
import sys, os, json, time
import numpy as np
import ca_common as A, ca_pairs as CPA
SA, SO = A.SA, A.SO
OUT = os.path.join(A.COMPLETE, 'out.nosync')
SECT = {'+': {'AV_ls+': ('AV', 'lamsP'), 'Ap_ls-': ('AW', 'lamsM'), 'l+_Asp': ('lamP', 'AsW'), 'l-_AsV': ('lamM', 'AsV')},
        '-': {'AV_ls-': ('AV', 'lamsM'), 'Ap_ls+': ('AW', 'lamsP'), 'l-_Asp': ('lamM', 'AsW'), 'l+_AsV': ('lamP', 'AsV')}}


def run(chir, nconf=8, seed=3303):
    rp, rpol = CPA.rngs(seed + (0 if chir == '+' else 1))
    fn = os.path.join(OUT, f'A_xiP_{"p" if chir == "+" else "m"}_raw.json')
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    for sec, (t0_, t1_) in SECT[chir].items():
        recs = out.get(sec, [])
        while len(recs) < nconf:
            xi = rpol('lamP' if chir == '+' else 'lamM')
            b = rp()
            p = rp()
            L = [SA.leg(t0_, rpol(t0_), p), SA.leg(t1_, rpol(t1_), -p)]
            t0 = time.time()
            r = SA.cycle_sum(L, [SA.gen_odd(xi, 'x'), SA.gen_even('P', a=b.astype(float))], slotnames=['x', 'P'],
                             what=('fin', ('pole', 1)))
            rec = dict(xi=[[z.real, z.imag] for z in xi], b=[int(x) for x in b],
                       legs=[dict(t=l['t'], pol=[[z.real, z.imag] for z in l['pol']], p0=[int(x.real) for x in l['p0']]) for l in L],
                       raw={'0': CPA.v5(SO.reduce_colour(r, 2, 'fin', 0).get('(Y0,Y1)', SA.vzero())),
                            '1': CPA.v5(SO.reduce_colour(r, 2, ('pole', 1), 0).get('(Y0,Y1)', SA.vzero()))},
                       sec=round(time.time() - t0, 2))
            recs.append(rec)
            out[sec] = recs
            json.dump(out, open(fn, 'w'))
        A.log(f'(xi{chir},P_b) sector {sec}: {len(recs)} configs (last {recs[-1]["sec"]} s)')


if __name__ == '__main__':
    for c in sys.argv[1:] or ['+', '-']:
        run(c)
