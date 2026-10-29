"""Evaluate engine B's CJob (copied source, read-only use) at engine A's off-shell W0 table points (bosonic pairs).
Usage: python3 -B b_at_a_points.py PAIR LEGSET [maxrows]"""
import os, sys, json, time, resource
HERE = os.path.dirname(os.path.abspath(__file__))
BSRC = os.path.join(HERE, '..', 'b_src')
sys.path.insert(0, BSRC)
os.chdir(BSRC)
import sympy as sp
import ce_job as CJ
from sympy.polys.domains import QQ_I

A_TAB = os.path.join(HERE, '..', '..', '..', 'engine_a', 'kernels_C2.json')
pair, legset = sys.argv[1], sys.argv[2]
maxrows = int(sys.argv[3]) if len(sys.argv) > 3 else 99
tab = json.load(open(A_TAB))[pair]['leg_sets'][legset]
ta, tb = pair.split('|')
COMPS = {'AV': list(range(4)), 'AW': list(range(4, 10)), 'AsV': list(range(4)), 'AsW': list(range(4, 10)), 'c': [0], 'cs': [0]}
def cnum(s):
    return complex(s.replace('(', '').replace(')', ''))
def todict(vec, t):
    return {c: sp.nsimplify(cnum(x).real) + sp.I * sp.nsimplify(cnum(x).imag) for c, x in zip(COMPS[t], vec)}
LOG = {'1': 1, 'log2': sp.log(2), 'log3': sp.log(3), 'log5': sp.log(5), 'log7': sp.log(7)}
out = []
beta = [1, 2, -1, 3]
for ir, row in enumerate(tab['rows'][:maxrows]):
    t0 = time.time()
    outs = [l['type'] for l in row['legs']]
    pol = [todict(l['pol'], l['type']) for l in row['legs']]
    za = todict(row['zeta_a'], ta); zb = todict(row['zeta_b'], tb)
    job = CJ.CJob(ta, tb, outs, pol=pol, za=za, zb=zb)
    res = job.kernels(mlist=(0,))
    g = job.g
    subsd = {}
    for j, l in enumerate(row['legs']):
        for m in range(4):
            subsd[sp.Symbol(f'r{j}_{m}')] = l['mom'][m]
    for m in range(4):
        subsd[sp.Symbol(f'a{m}')] = beta[m]
    vals = {}
    for tabn in ('Pb', 'D0'):
        for w, dm in res[tabn].items():
            tot = 0
            for b, P in dm.get(0, {}).items():
                tot += sp.expand(P.as_expr().subs(subsd)) * LOG.get(b, sp.Symbol('L_' + str(b)))
            vals[f'{tabn} word {w}'] = str(sp.nsimplify(sp.expand(tot)))
    rec = {'row': ir, 'word_A': row['word'], 'outs_B': outs, 'A_Pbeta': row['C2_Pbeta'], 'A_D0': row.get('C2_0_Delta'),
           'B': vals, 'seconds': round(time.time() - t0, 1),
           'rss_MB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20}
    print(json.dumps(rec), flush=True)
    out.append(rec)
json.dump(out, open(os.path.join(HERE, '..', f'b_at_a_{pair.replace("|","_")}_{legset.replace(" ","_")}.json'), 'w'), indent=1)
