"""T1 (route B), pair weight W = W0 (no derivatives on the letters): all 21 pair types, every leg set of the census
at W0 (log and eps^-m), every ordered output sequence (topology), symbolic letter components, output polarisations
and output momenta.  Writes kernels_W0.json incrementally (one pair at a time) and checks the output bound (G7)."""
import os, sys, json, time, resource
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tl_engine as T

CENSUS = os.path.join(HERE, '..', 'theory', 'census.json')
OUT = os.path.join(HERE, 'kernels_W0.json')


def mom_degree(job, poly):
    idx = [i for i, nm in enumerate(job.names) if nm.startswith('r')]
    degs = set(sum(mon[i] for i in idx) for mon in poly.itermonoms())
    return sorted(degs)


PARTS = os.path.join(HERE, 'parts_W0')


def one(pair, ls):
    """compute one leg set in this process and write parts_W0/<pair>__<ls>.json"""
    sa, sb = pair.split('|')
    cen = json.load(open(CENSUS))['output_leg_sets_by_pair_weight'][pair]
    W0 = Fraction(list(cen.keys())[0])
    d_ab = T.DEG[sa] + T.DEG[sb]
    outs = ls.split() if ls else []
    t1 = time.time()
    job = T.Job(sa, sb, outs)
    K = job.kernel()
    w_out = sum(T.WEYL[o] for o in outs)
    D = W0 - w_out
    g7 = dict(n_out=len(outs), bound=str(W0 + d_ab), n_out_ok=len(outs) <= W0 + d_ab)
    deg_ok = True
    tot = {}
    for w, d in K['total'].items():
        ws = ','.join(map(str, w))
        tot[ws] = {}
        for m, p in d.items():
            tot[ws][f'eps^-{m}' if m != 'log' else 'log'] = T.poly_str(job, p)
            want = int(D) if m == 'log' else int(D - 2 * m)
            if mom_degree(job, p) != [want]:
                deg_ok = False
    g7['momentum_degrees_ok'] = deg_ok
    byseq = []
    for seq, sres in K['by_sequence']:
        byseq.append(dict(topology=[list(b) for b in seq],
                          vertices=['quartic' if len(b) == 2 else 'cubic' for b in seq],
                          words={','.join(map(str, w)): {(f'eps^-{m}' if m != 'log' else 'log'): T.poly_str(job, p)
                                                         for m, p in d.items()} for w, d in sres.items()}))
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20
    ent = dict(outputs=outs, w_out=str(w_out), superficial_degree=str(D), total=tot, by_topology=byseq, G7=g7,
               seconds=round(time.time() - t1, 2), peak_rss_MB=round(rss))
    os.makedirs(PARTS, exist_ok=True)
    json.dump(ent, open(os.path.join(PARTS, f'{pair.replace("|", "_")}__{ls.replace(" ", "_") or "none"}.json'), 'w'))
    print(f'{pair:10s} [{ls}] words={len(tot)} G7={g7} {time.time() - t1:.1f}s peak {rss:.0f} MB', flush=True)


def run(pairs=None, verbose=True):
    """one subprocess per leg set (bounded memory); results merged into kernels_W0.json."""
    import subprocess
    cen = json.load(open(CENSUS))
    legs = cen['output_leg_sets_by_pair_weight']
    try:
        res = json.load(open(OUT))
    except Exception:
        res = {}
    for pair, v in legs.items():
        if pairs is not None and pair not in pairs:
            continue
        if pair in res and res[pair].get('complete'):
            continue
        sa, sb = pair.split('|')
        W0s = list(v.keys())[0]
        W0 = Fraction(W0s)
        sets = list(v[W0s]['log'])
        for k in v[W0s]['eps_max_m']:
            ks = '' if k == '(none)' else k
            if ks not in sets:
                sets.append(ks)
        entry = dict(letters=[sa, sb], letter_names=['X_a', 'X_b'], W=str(W0), d_ab=T.DEG[sa] + T.DEG[sb], leg_sets={})
        t0 = time.time()
        failed = []
        for ls in sets:
            fn = os.path.join(PARTS, f'{pair.replace("|", "_")}__{ls.replace(" ", "_") or "none"}.json')
            if not os.path.exists(fn):
                r = subprocess.run([sys.executable, __file__, '--one', pair, ls], cwd=HERE)
                if r.returncode != 0 or not os.path.exists(fn):
                    failed.append(ls)
                    print(f'{pair} [{ls}] FAILED (exit {r.returncode})', flush=True)
                    continue
            entry['leg_sets'][ls if ls else '(none)'] = json.load(open(fn))
        entry['complete'] = not failed
        entry['failed_leg_sets'] = failed
        entry['seconds'] = round(time.time() - t0, 1)
        res[pair] = entry
        json.dump(res, open(OUT, 'w'), indent=0)
        with open(os.path.join(HERE, 'PROGRESS.md'), 'a') as f:
            f.write(f'  kernels_W0: {pair} {"done" if not failed else "INCOMPLETE " + str(failed)} ({entry["seconds"]} s)\n')
    return res


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--one':
        one(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else '')
    else:
        run(sys.argv[1:] or None)
