"""Derivative tables (route B): K(d^alpha phi_a, d^beta phi_b) from the point-split generating function
sum s^alpha u^beta/(alpha! beta!) K(...), s,u symbolic, truncated at total order n = W - W0.
Selection (resources): every pair type at W = W0 + 1 for the leg sets with at most two outputs, and the scalar and
gauge pairs (AV|AV, AV|AW, AW|AW) at W = W0 + 2 for the leg sets with at most two outputs.
The strings keep s,u: the coefficient of s^alpha u^beta times alpha! beta! is the kernel of the pair with those derivatives.
Writes kernels_deriv.json incrementally."""
import os, sys, json, time, resource
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tl_engine as T

CENSUS = os.path.join(HERE, '..', 'theory', 'census.json')
OUT = os.path.join(HERE, 'kernels_deriv.json')


def run(maxout=2, only=None):
    cen = json.load(open(CENSUS))['output_leg_sets_by_pair_weight']
    try:
        res = json.load(open(OUT))
    except Exception:
        res = {}
    plan = [(p, 1) for p in cen] + [(p, 2) for p in ('AV|AV', 'AV|AW', 'AW|AW')]
    for pair, n in plan:
        key = f'{pair} @ W0+{n}'
        if only is not None and key != only:
            continue
        if key in res and res[key].get('complete'):
            continue
        if only is None:          # one subprocess per (pair, n): bounded memory
            import subprocess
            r = subprocess.run([sys.executable, __file__, '--one', key], cwd=HERE)
            res = json.load(open(OUT)) if os.path.exists(OUT) else {}
            if r.returncode != 0:
                print(f'{key} FAILED (exit {r.returncode})', flush=True)
            continue
        sa, sb = pair.split('|')
        Ws = list(cen[pair].keys())
        W0 = Fraction(Ws[0])
        Wk = str(W0 + n) if (W0 + n).denominator > 1 else str(int(W0 + n))
        if Wk not in cen[pair]:
            continue
        sets = list(cen[pair][Wk]['log'])
        for k in cen[pair][Wk]['eps_max_m']:
            ks = '' if k == '(none)' else k
            if ks not in sets:
                sets.append(ks)
        sets = [s for s in sets if len(s.split()) <= maxout]
        entry = dict(letters=[sa, sb], W=Wk, derivatives_on_letters=n,
                     convention='coefficient of s^alpha u^beta (|alpha|+|beta| = n) times alpha! beta! = K(d^alpha X_a, d^beta X_b); '
                                'units hbar/(16 pi^2); words as in the W0 tables', leg_sets={})
        t0 = time.time()
        for ls in sets:
            outs = ls.split() if ls else []
            t1 = time.time()
            job = T.Job(sa, sb, outs, deriv=('gen', n))
            K = job.kernel(keep_sequences=False)
            sidx = [job.names.index(f's{m}') for m in range(4)] + [job.names.index(f'u{m}') for m in range(4)]
            tot = {}
            for w, d in K['total'].items():
                ws = ','.join(map(str, w))
                for m, p in d.items():
                    p = job.R.from_dict({mon: c for mon, c in p.terms() if sum(mon[i] for i in sidx) == n})
                    if p:
                        tot.setdefault(ws, {})[f'eps^-{m}' if m != 'log' else 'log'] = T.poly_str(job, p)
            entry['leg_sets'][ls if ls else '(none)'] = dict(outputs=outs, total=tot, seconds=round(time.time() - t1, 1))
            print(f'{key:22s} [{ls}] words={len(tot)} {time.time() - t1:.1f}s', flush=True)
        entry['complete'] = True
        res[key] = entry
        json.dump(res, open(OUT, 'w'), indent=0)
        with open(os.path.join(HERE, 'PROGRESS.md'), 'a') as f:
            f.write(f'  kernels_deriv: {key} done ({time.time() - t0:.0f} s), peak RSS '
                    f'{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB\n')
    return res


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--one':
        run(only=sys.argv[2])
    else:
        run()
