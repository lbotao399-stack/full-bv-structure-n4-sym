"""Remark d-gf: H = [Q, QGF] after deleting (a) the whole arrow A*|_V -> A|_V (here ALs -> AL and ATs -> AT,
the local deletion) and (b) only its longitudinal part ALs -> AL (the variant used by the theory toy and in
d1_checks.py).  Reports H on each slot as a multiple of s^2 (diagonal entries)."""
import json, os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp import mat_mul, mat_add, mat_comm, mat_maxabs
import model

def H_by_slot(drop_T, gF):
    m = model.build(N=4, gF=gF)
    G = dict(m['G'])
    if drop_T:
        idx, N = m['idx'], m['N']
        for a in range(N + 1):
            G.pop((idx[('AT', a)], idx[('ATs', a)]), None)
    H = mat_comm(m['Q'], G, 1, 1)
    G2 = mat_mul(G, G)
    out = {}
    for (nm, *_r) in m['slots']:
        i = m['idx'][(nm, 0)]; j = m['idx'][(nm, 2)]
        out[nm] = str(H.get((j, i), Fr(0)))     # coefficient of s^2 on the slot
    out['|(QGF)^2|'] = str(mat_maxabs(G2))
    return out

res = {'Feynman (all arrows)': H_by_slot(False, Fr(-2)),
       '(a) whole arrow A*|_V -> A|_V deleted (local)': H_by_slot(True, Fr(0)),
       '(b) only longitudinal part ALs -> AL deleted (toy variant; nonlocal in the field theory)': H_by_slot(False, Fr(0))}
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'remark_gf_H_results.json'), 'w') as fh:
    json.dump(res, fh, indent=1)
print(json.dumps(res, indent=1))
