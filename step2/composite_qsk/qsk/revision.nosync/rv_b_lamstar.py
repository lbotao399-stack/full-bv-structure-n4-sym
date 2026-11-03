"""Revision (Part Q): route B with the T5 placement of a lam* output (finding Q4) on the n = 0 leg sets of TQ2 that
contain a lam* output.  Each stored file of the first production pass (engine_b/out.nosync/tq2/*lams*) is recomputed
with the fixed code (same xi, symbolic letters and outputs) and compared term by term after the sign convention of
qb_assemble.fix_sign (files without sign_version 2 carry -D); the contribution of the new placement alone is recorded.
Writes revision.nosync/out/b_lamstar_n0.json and the recomputed tables to revision.nosync/out/tq2_b_fixed/."""
import os, sys, json, glob, time
import sympy as sp
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_b/src.nosync'))
import qb_core as B
import qb_job as QJ

OUT = os.path.join(QSK, 'revision.nosync/out')
TAB = os.path.join(OUT, 'tq2_b_fixed')
os.makedirs(TAB, exist_ok=True)


def kernel_tab(job, res):
    return {f'word {list(w)}': {f'eps^-{m}': {b: QJ.poly_str(job, P) for b, P in d.items()} for m, d in dm.items()}
            for w, dm in res['total'].items()}


import re
_BR = re.compile(r'(\w+)\[([^\]]*)\]')


def sym(s):
    if s == '0':
        return 0
    s2 = _BR.sub(lambda m: m.group(1) + '__' + re.sub(r'\W', '_', m.group(2)), s)
    return sp.expand(sp.sympify(s2.replace('^', '**')))


def same_P(r1, r2):
    t1, t2 = r1['total'], r2['total']
    for w in set(t1) | set(t2):
        a = t1.get(w, {}); b = t2.get(w, {})
        for m in set(a) | set(b):
            aa = a.get(m, {}); bb = b.get(m, {})
            for base in set(aa) | set(bb):
                x = aa.get(base); y = bb.get(base)
                if x is None or y is None:
                    if (x if x is not None else y):
                        return False
                    continue
                if (x - y.ring.from_dict(dict(y.items())) if x.ring != y.ring else x - y):
                    return False
    return True


def same(t1, t2, sign):
    keys = set(t1) | set(t2)
    for w in keys:
        a = t1.get(w, {}); b = t2.get(w, {})
        for m in set(a) | set(b):
            aa = a.get(m, {}); bb = b.get(m, {})
            for base in set(aa) | set(bb):
                x = sym(aa.get(base, '0')); y = sym(bb.get(base, '0'))
                if sp.expand(x - sign * y) != 0:
                    return False
    return True


def main():
    t0 = time.time()
    rows = []
    for fn in sorted(glob.glob(os.path.join(QSK, 'engine_b/out.nosync/tq2/*lams*.json'))):
        d = json.load(open(fn))
        sa, sb = d['pair'].split('|')
        tag = d['xi']
        if tag == 'xiP':
            xi = B.XI[:8] + [0] * 8
        else:
            xi = [0] * 16; xi[int(tag.split('_e')[1])] = 1
        zero_mom = bool(d.get('zero_output_momenta'))
        res = {}
        for flag in (False, True):
            QJ.T5_LAMS_OUTPUT = flag
            job = QJ.DJob(sa, sb, d['outs'], 'xi', xi, tag='rvl' + tag + str(flag), zero_mom=zero_mom)
            res[flag] = (job, job.kernel(mlist=(0, 1, 2)))
        QJ.T5_LAMS_OUTPUT = True
        tab_old_code = kernel_tab(*res[False]); tab_new = kernel_tab(*res[True])
        stored_sign = 1 if d.get('sign_version') == 2 else -1
        rec = dict(file=os.path.basename(fn), pair=d['pair'], legset=d['legset'], xi=tag,
                   stored_equals_old_code=same(tab_old_code, d['kernel'], stored_sign),
                   new_equals_old_code=same(tab_new, tab_old_code, 1),
                   new_equals_old_code_exact=same_P(res[True][1], res[False][1]) if res[True][0].names == res[False][0].names else None,
                   configs_old=res[False][1]['configs'], configs_new=res[True][1]['configs'],
                   seconds=round(time.time() - t0, 1), peak_MB=round(QJ.peak_mb(), 1))
        json.dump(dict(d, kernel=tab_new, sign_version=2, revision='T5 with lam* output included (finding Q4)',
                       configs=res[True][1]['configs']), open(os.path.join(TAB, os.path.basename(fn)), 'w'), indent=0)
        rows.append(rec)
        print(rec, flush=True)
        json.dump(rows, open(os.path.join(OUT, 'b_lamstar_n0.json'), 'w'), indent=1)
    print('ALL', len(rows), 'unchanged by the new placement:', sum(r['new_equals_old_code'] for r in rows),
          'stored = old code:', sum(r['stored_equals_old_code'] for r in rows))


if __name__ == '__main__':
    main()
