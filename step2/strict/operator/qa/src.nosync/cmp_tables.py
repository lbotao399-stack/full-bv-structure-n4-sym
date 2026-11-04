"""QA comparator: engine A vs engine B table rows (operator side), basis-independent row data.
Reads only results/compact tables of both engines; writes qa/out.nosync/cmp_tables.json."""
import json, collections, os, sys
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/operator'
A = json.load(open(f'{ROOT}/engine_a/out.nosync/tables_compact.json'))
rb = json.load(open(f'{ROOT}/engine_b/results.json'))
B = rb['OP1_tables_n0']['rows'] + rb['OP1_tables_n1']['rows']

def norm_legs(s):
    if s in ('(none)', '', None):
        return ()
    return tuple(sorted(s.split()))

def keyA(x):
    return (tuple(sorted(x['pair'])), norm_legs(x['leg_set']), x['n'])

def keyB(x):
    return (tuple(sorted(x['pair'].split('|'))), norm_legs(x['legset']), x['n'])

GEN = {'D0': 'D0', 'P': 'P', 'M0': 'M0', 'R': 'R'}
Aidx = collections.defaultdict(list)
for x in A:
    if x['generator'] in GEN and x['n'] <= 1:
        Aidx[(x['generator'],) + keyA(x)].append(x)

out = {'overlap_rows': 0, 'by_gen': {}, 'D_mismatch': [], 'conj_disagree': [], 'pole_presence_disagree': [],
       'F_nonzero_A_vs_B': [], 'B_rows_without_A': [], 'notes': []}
cnt = collections.Counter()
for xb in B:
    kb = keyB(xb)
    for g in GEN:
        rowsA = Aidx.get((g,) + kb)
        if not rowsA:
            if g == 'D0':
                out['B_rows_without_A'].append(list(map(str, kb)))
            continue
        cnt[g] += 1
        # superficial degree: B gives D for EP (orbital) and ESW blocks; A per row D (max over orderings)
        DA = max(int(r['D']) for r in rowsA)
        DB = xb.get('D', {}).get('EP')
        if g == 'P' and DB is not None:
            # P rows: A's D for P is the EP-type degree
            if int(DB) != DA:
                out['D_mismatch'].append({'gen': g, 'key': list(map(str, kb)), 'A': DA, 'B': DB})
        # conjecture verdicts
        vA = all(r.get('conj_holds', True) for r in rowsA)
        vB = xb.get('G_D1pair', {}).get(g)
        if vB is not None and vB != vA:
            out['conj_disagree'].append({'gen': g, 'key': list(map(str, kb)), 'A': vA, 'B': vB,
                                         'B_has': xb.get('has')})
        # finite special vertex nonzero (A only on single-output rows with F_x_single_output)
        fA = [r.get('F_x_single_output') for r in rowsA if 'F_x_single_output' in r]
        if fA and g == 'D0':
            out['F_nonzero_A_vs_B'].append({'key': list(map(str, kb)),
                                            'A_F_nonzero': any(v not in ('0', 0) for v in fA),
                                            'B_F_nonzero': xb.get('finite_special_vertex_nonzero')})
        # pole presence along D0: A poles nonzero?  B G_pole defined?
        if g == 'D0':
            pA = any(any(v not in ('0', 0) for v in (r.get('poles') or {}).values()) for r in rowsA)
            out.setdefault('D0_poles', []).append({'key': list(map(str, kb)), 'A_nonzero_poles': pA,
                                                    'B_pole_gate': xb.get('G_pole', {}).get('D0')})
out['overlap_rows'] = dict(cnt)
out['n_B_rows'] = len(B)
os.makedirs(f'{ROOT}/qa/out.nosync', exist_ok=True)
json.dump(out, open(f'{ROOT}/qa/out.nosync/cmp_tables.json', 'w'), indent=1)
print('overlap', dict(cnt), 'B rows', len(B))
print('D mismatch', len(out['D_mismatch']), out['D_mismatch'][:5])
print('conj disagree', len(out['conj_disagree']), out['conj_disagree'][:8])
print('F', out['F_nonzero_A_vs_B'])
print('B rows without A (D0)', len(out['B_rows_without_A']), out['B_rows_without_A'][:12])
