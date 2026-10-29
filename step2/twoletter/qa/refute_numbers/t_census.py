"""Independent census of the log-eps leg sets (Prop prop:tl-census / prop:tl-bound): own chain enumeration from the vertex
list of eq c1-I (chiralities of Prop c1-chirality) and the four edges; leg set allowed iff k >= 1 and w_out <= W.
Compared with theory/census.json output_leg_sets_by_pair_weight for every pair and W <= 7 listed there."""
import json
from fractions import Fraction as F
from collections import Counter
VERTS = [('AV', 'AV', 'AV'), ('AV', 'AW', 'AW'), ('AV', 'AV', 'AV', 'AV'), ('AV', 'AV', 'AW', 'AW'), ('AW', 'AW', 'AW', 'AW'),
         ('lamP', 'lamM', 'AV'), ('lamP', 'lamP', 'AW'), ('lamM', 'lamM', 'AW'), ('AV', 'c', 'AsV'), ('AW', 'c', 'AsW'),
         ('c', 'lamP', 'lamsP'), ('c', 'lamM', 'lamsM'), ('c', 'c', 'cs')]
EDGE = {'AV': 'AV', 'AW': 'AW', 'lamP': 'lamM', 'lamM': 'lamP', 'c': 'AsV', 'AsV': 'c'}
WT = {'AV': F(1), 'AW': F(1), 'lamP': F(3, 2), 'lamM': F(3, 2), 'c': F(0), 'AsV': F(3), 'AsW': F(3), 'lamsP': F(5, 2),
      'lamsM': F(5, 2), 'cs': F(4)}
DEG = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 0, 'lamM': 0, 'AsV': -1, 'AsW': -1, 'lamsP': -1, 'lamsM': -1, 'cs': -2}
C = json.load(open('../../theory/census.json'))['output_leg_sets_by_pair_weight']
def legsets(ta, tb, W, kmax):
    out = set()
    # state: slot type at the end of the current edge (vertex in-slot type), outputs so far, k
    def rec(intype, outs, k):
        if k >= kmax:
            return
        for v in VERTS:
            if intype not in v:
                continue
            rest = list(v); rest.remove(intype)
            for i, ot in enumerate(rest):
                if ot not in EDGE:
                    continue
                o2 = rest[:i] + rest[i + 1:]
                no = outs + o2
                if sum(WT[t] for t in no) > W:
                    continue
                nxt = EDGE[ot]
                if nxt == tb:            # close the chain into letter b
                    out.add(' '.join(sorted(no)))
                rec(nxt, no, k + 1)
    rec(EDGE[ta], [], 0)
    return out
bad = 0; n = 0; rep = {}
for pair, byW in C.items():
    ta, tb = pair.split('|')
    for Ws, d in byW.items():
        W = F(Ws)
        dab = DEG[ta] + DEG[tb]
        mine = legsets(ta, tb, W, kmax=int(W + dab) + 2)
        theirs = set(d['log'])
        n += 1
        if mine != theirs:
            bad += 1
            rep['%s W=%s' % (pair, Ws)] = dict(only_referee=sorted(mine - theirs), only_census=sorted(theirs - mine))
print('pair/weight cells compared:', n, ' differing:', bad)
for k, v in list(rep.items())[:10]:
    print(k, v)
json.dump(dict(cells=n, differing=bad, diffs=rep), open('out_census.json', 'w'), indent=1)
