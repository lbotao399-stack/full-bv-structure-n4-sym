"""Revision (Part Q): engine A TQ2 tables after the exact contraction (rerun of run_all.py) against the tables of the first
production run (engine_a/bak.nosync/out_before_revision.nosync).  Writes revision.nosync/out/tq2_regression_A.json."""
import os, json
QSK = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/qsk'
OLD = os.path.join(QSK, 'engine_a/bak.nosync/out_before_revision.nosync')
NEW = os.path.join(QSK, 'engine_a/out.nosync')
out = {'tables': {}}
allsame = True
for n in range(4):
    a = json.load(open(os.path.join(OLD, 'tq2_n%d.json' % n)))
    b = json.load(open(os.path.join(NEW, 'tq2_n%d.json' % n)))
    key = lambda r: (tuple(r['pair']), r['leg_set'], tuple(r['ordering']))
    da = {key(r): (r['eps0'], json.dumps(r['poles'], sort_keys=True)) for r in a}
    db = {key(r): (r['eps0'], json.dumps(r['poles'], sort_keys=True)) for r in b}
    diff = [list(map(str, k)) for k in set(da) | set(db) if da.get(k) != db.get(k)]
    out['tables']['n=%d' % n] = {'rows_old': len(a), 'rows_new': len(b), 'differing_rows': diff[:20], 'n_differing': len(diff)}
    allsame = allsame and not diff
out['summary'] = 'TQ2 n = 0..3: %s' % ('all rows identical (eps^0 and poles)' if allsame else 'DIFFERENCES, see tables')
json.dump(out, open(os.path.join(QSK, 'revision.nosync/out/tq2_regression_A.json'), 'w'), indent=1)
print(out['summary'], {k: v['n_differing'] for k, v in out['tables'].items()})
