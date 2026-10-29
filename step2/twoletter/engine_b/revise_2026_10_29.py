"""Engine B, revision 2026-10-29 (qa/disposition.json N4, N5, N6, F-2, F-3):
 - two_site.json: the D_1 entries become numbers times the unit u (conditional on P7 at theta = -1/3), the classical part
   kept separately; the scheme-independent content gamma = 2 eig(lc^natural) recorded next to them;
 - gates.json: derived W2 gate (entrywise equality with engine A on the compared tables, qa/compare C1/C2, carries
   engine A's W2 gate over to these entries; not an independent run)."""
import json, os
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
ts = json.load(open(os.path.join(HERE, 'two_site.json')))
for S, ent in ts['twist_two'].items():
    old = ent.pop('D1_eigenvalues_units', None)
    if old is None:
        continue
    nums = [x.split(' - ')[0].strip() for x in old]
    ent['D1_eigenvalues_times_u_conditional_on_P7_theta=-1/3'] = nums
    ent['D1_classical_part'] = '-(1/3 + theta) u {S^(1)(Delta), -} (separate operator, not an eigenvalue)'
    ent['gamma_times_u_scheme_independent'] = [str(Fraction(x) / 2) for x in nums]
ts['_revision_2026-10-29'] = ('D_1 strings rewritten (numbers times u; classical part separate; conditional on P7); '
                              'gamma = 2 eig(lc^natural) = D_1/2 recorded')
json.dump(ts, open(os.path.join(HERE, 'two_site.json'), 'w'), indent=1)
g = json.load(open(os.path.join(HERE, 'gates.json')))
g['G3 W2 (derived)'] = dict(
    value='engine B tables equal engine A entrywise on all W_0 tables (389/389 rows, rerun 2026-10-29) and on the W_0+1, '
          'W_0+2 tables (164/164 words); engine A W2 gate 741/741 (232 nontrivial)',
    target='W2 on the compared entries', pass_=True,
    note='derived through qa/compare and results/compare_rerun, not an independent engine B run of W2')
json.dump(g, open(os.path.join(HERE, 'gates.json'), 'w'), indent=1)
print('ok')
