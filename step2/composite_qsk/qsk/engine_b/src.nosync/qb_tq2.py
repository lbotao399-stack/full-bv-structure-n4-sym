"""Task TQ2 (route B): Upsilon_xi on the pairs of field and ghost letters (Table qs-fieldpairs), n = 0 letter
derivatives, every leg set of census_qsk.json with output weight <= W_0 + 1/2, every chain (all placements of T1..T5),
general colour (words) and the planar adjacent form; xi in S_+ (x) W symbolic (xi0..xi7) unless stated, letters and
outputs symbolic.  Poles (eps^-1) and eps^0 exactly.  Resumable: one json file per leg set in out.nosync/tq2/.
Usage: python3 qb_tq2.py [maxsec] [pair filter]"""
import os, sys, json, time, signal
from fractions import Fraction as F
import sympy as sp
from sympy.polys.domains import QQ_I
import qb_core as B
import qb_job as QJ

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync', 'tq2')
os.makedirs(OUT, exist_ok=True)
CENSUS = os.path.join(HERE, '..', '..', 'theory', 'src.nosync', 'census_qsk.json')
PAIRS = ['AV|AV', 'AV|AW', 'AV|lamP', 'AV|lamM', 'AV|c', 'AW|AW', 'AW|lamP', 'AW|lamM', 'AW|c', 'lamP|lamP',
         'lamP|lamM', 'lamP|c', 'lamM|lamM', 'lamM|c', 'c|c']
PROG = os.path.join(HERE, '..', 'PROGRESS.txt')


def log(msg):
    with open(PROG, 'a') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S') + ' ' + msg + '\n')


def legsets(gen='xi', nder='0'):
    d = json.load(open(CENSUS))
    p = d['generators'][gen]['pairs']
    out = []
    for k in PAIRS:
        e = p.get(k)
        if e is None:
            continue
        b = e['by_letter_derivatives'][nder]
        W = F(b['W'])
        for n, v in b['defect_leg_sets'].items():
            if F(v['w_out']) <= W + F(1, 2):
                outs = [] if n == '(none)' else n.split()
                out.append(dict(pair=k, legset=n, outs=outs, W=str(W), chains_census=v['chains'], k_I_min=v['k_I_min'],
                                w_out=v['w_out'], placements_census=v['placements']))
    return out


def cost(ls):
    return (len(ls['outs']), ls['chains_census'])


class Timeout(Exception):
    pass


def _alarm(sig, frm):
    raise Timeout()


def planar(job, tot, m=0):
    """planar adjacent form: a term a C_w -> (-1)^r a/2 units u, outputs in trace order Y_wr ... Y_w1; odd outputs
    reordered from label order to trace order with the permutation sign; identical outputs: the labelled multilinear
    form is kept (the division by the symmetry factor is applied when a functional is formed, twoletter revision)."""
    R = job.R
    res = {}
    for w, dm in tot.items():
        d = dm.get(m, {})
        order = tuple(reversed(w))
        odd = [j for j in order if job.odd_out[j]]
        sg = V_perm(odd)
        for b, P in d.items():
            key = (order, b)
            res[key] = res.get(key, R.zero) + P * QQ_I(F((-1) ** len(w) * sg, 2), 0)
    return {k: v for k, v in res.items() if v}


def V_perm(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def run_one(ls, xi, tag, maxsec):
    sa, sb = ls['pair'].split('|')
    fn = os.path.join(OUT, f"{ls['pair'].replace('|', '_')}__{ls['legset'].replace(' ', '_').replace('(', '').replace(')', '')}__{tag}.json")
    if os.path.exists(fn):
        return 'cached', fn
    t0 = time.time()
    signal.signal(signal.SIGALRM, _alarm)
    signal.alarm(int(maxsec))
    try:
        # superficial degree D = W + 1/2 - w_out: for D = 0 the coefficients are momentum independent and are computed
        # at zero output momenta (exact: they are polynomials of degree 0 in the momenta)
        D0 = (F(ls['W']) + F(1, 2) - F(ls['w_out'])) == 0
        job = QJ.DJob(sa, sb, ls['outs'], 'xi', xi, tag='tq2' + tag, zero_mom=(D0 and len(ls['outs']) >= 3))
        res = job.kernel(mlist=(0, 1, 2))
    except Timeout:
        signal.alarm(0)
        log(f"TQ2 {ls['pair']} [{ls['legset']}] {tag}: TIMEOUT after {maxsec} s")
        return 'timeout', None
    except MemoryError as e:
        signal.alarm(0)
        log(f"TQ2 {ls['pair']} [{ls['legset']}] {tag}: MEMORY {e}")
        return 'memory', None
    signal.alarm(0)
    tot = res['total']
    nterms = sum(len(P.terms()) for dm in tot.values() for d in dm.values() for P in d.values())
    pl = planar(job, tot, 0)
    out = dict(ls, xi=tag, sign_version=QJ.SIGN_VERSION, zero_output_momenta=job.omom[0][0] == job.R.zero if ls['outs'] else False, seconds=round(time.time() - t0, 1), configs=res['configs'], nterms=nterms,
               peak_MB=round(QJ.peak_mb(), 1), ring=job.names,
               kernel={f'word {list(w)}': {f'eps^-{m}': {b: QJ.poly_str(job, P) for b, P in d.items()} for m, d in dm.items()}
                       for w, dm in tot.items()},
               planar_eps0={f'trace order {list(o)} | {b}': QJ.poly_str(job, P) for (o, b), P in pl.items()},
               poles_present=sorted(set(m for dm in tot.values() for m in dm if m > 0)),
               bases=sorted(set(b for dm in tot.values() for d in dm.values() for b in d)))
    json.dump(out, open(fn, 'w'), indent=0)
    log(f"TQ2 {ls['pair']} [{ls['legset']}] {tag}: {out['seconds']} s, {nterms} terms, poles {out['poles_present']}, bases {out['bases']}, peak {out['peak_MB']} MB")
    return 'done', fn


def run_one_comp(ls, j, maxsec):
    xi = [0] * 16
    xi[j] = 1
    return run_one(ls, xi, f'xiP_e{j}', maxsec)


if __name__ == '__main__':
    import subprocess
    if len(sys.argv) > 1 and sys.argv[1] == '--one':
        pair, lsname, mode, maxsec = sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5])
        ls = [x for x in legsets() if x['pair'] == pair and x['legset'] == lsname][0]
        if mode == 'sym':
            st, fn = run_one(ls, B.XI[:8] + [0] * 8, 'xiP', maxsec)
        else:
            st, fn = run_one_comp(ls, int(mode[1:]), maxsec)
        print(st)
        sys.exit(0 if st in ('done', 'cached') else 3)
    maxsec = float(sys.argv[1]) if len(sys.argv) > 1 else 1800
    filt = sys.argv[2] if len(sys.argv) > 2 else None
    lss = sorted(legsets(), key=cost)
    if filt:
        lss = [x for x in lss if x['pair'] == filt]
    for ls in lss:
        fn_sym = os.path.join(OUT, f"{ls['pair'].replace('|', '_')}__{ls['legset'].replace(' ', '_').replace('(', '').replace(')', '')}__xiP.json")
        fn_e3 = fn_sym.replace('__xiP.json', '__xiP_e3.json')
        if os.path.exists(fn_sym) or os.path.exists(fn_e3):
            print(ls['pair'], ls['legset'], 'cached', flush=True)
            continue
        if len(ls['outs']) >= 2 or ls.get('k_I_min', 0) >= 2:
            r = subprocess.CompletedProcess([], 3, stdout='skip sym', stderr='')
        else:
            r = subprocess.run([sys.executable, __file__, '--one', ls['pair'], ls['legset'], 'sym', str(maxsec)],
                               capture_output=True, text=True)
        st = r.stdout.strip().split('\n')[-1] if r.stdout.strip() else 'error'
        if r.returncode != 0:
            # the kernel at the single rank-one xi = e_3 (u (x) e_4) with symbolic letters and outputs determines it for
            # every xi (Spin(4) x SU(4) covariance, S_+ (x) W irreducible, linearity); run in a fresh process
            sts = []
            for j in [3]:
                rj = subprocess.run([sys.executable, __file__, '--one', ls['pair'], ls['legset'], f'e{j}', str(maxsec)],
                                    capture_output=True, text=True)
                sts.append(rj.stdout.strip().split('\n')[-1] if rj.stdout.strip() else 'error:' + rj.stderr.strip()[-200:])
            st = 'components: ' + ','.join(sts)
        print(ls['pair'], ls['legset'], st, flush=True)
