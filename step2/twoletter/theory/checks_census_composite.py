"""Census of the composite cycles through two letters (the cycles of the dilatation component D_1 that carry the
composite kernel 2(P* v. P* (x) 1)K_0 of Def c9-composite on exactly one edge; Prop c9-composite (ii): it joins a
slot c to a slot A|_V).  Same vertices and P edges as checks_census.py.

Weight rule.  The composite kernel obeys the slot-weight rule of the P edges (an edge between slots of weights w1, w2
has momentum degree w1 + w2 - 4): c-A|_V gives -3, the degree of (Q^GF)^2 |p|^{-4} on that component (Prop c9-composite
(iv)); v. contributes in addition an explicit position (weight -1).  The outputs of a composite cycle carry the degree
d_ab - 1.  Under (H3) the Q-image is local, of Weyl weight W (fields and derivatives, the position at weight -1), at most
linear in the position; Q changes one output of weight w_t into a field of weight w_t - delta_t with delta_t derivatives
(A* -> A: delta 2; lambda* -> lambda, c* -> A*, A|_V -> c: delta 1).  Hence a composite cycle enters C^(2)_Delta at pair
weight W only if w_out - delta_t <= W + 1 for an output t with delta_t > 0, i.e. w_out <= W + 3.
Output: census_composite.json (leg sets by pair type at W_0 with w_out <= W_0 + 3) and a printed summary."""
import json, os
from fractions import Fraction as Fr
from collections import defaultdict
import checks_census as CC

HERE = os.path.dirname(os.path.abspath(__file__))
COMP = {'c': 'AV', 'AV': 'c'}            # the composite edge joins c and A|_V
DELTA = {'AsV': 2, 'AsW': 2, 'lamsP': 1, 'lamsM': 1, 'cs': 1, 'AV': 1}


def comp_chains(t1, t2, wmax):
    """chains letter1 - ... - letter2 with exactly one composite edge (position comp_edge = index of the edge, 0..k)"""
    res = []
    if COMP.get(t1) == t2:
        res.append(dict(k=0, path=[], outputs=(), comp_edge=0))

    def rec(cur, path, outs, wout, ce):
        options = [(CC.PARTNER.get(cur), None)]
        if ce is None and cur in COMP:
            options.append((COMP[cur], len(path)))
        for enter, newce in options:
            if enter is None:
                continue
            ce2 = ce if newce is None else newce
            for (name, ex, o, nd) in CC.vertex_moves(enter):
                w2 = wout + sum(CC.W[x] for x in o)
                if w2 > wmax:
                    continue
                p2 = path + [(name, enter, ex)]
                outs2 = outs + list(o)
                if ce2 is not None and CC.PARTNER.get(ex) == t2:
                    res.append(dict(k=len(p2), path=p2, outputs=tuple(sorted(outs2)), comp_edge=ce2))
                if ce2 is None and COMP.get(ex) == t2:
                    res.append(dict(k=len(p2), path=p2, outputs=tuple(sorted(outs2)), comp_edge=len(p2)))
                if len(p2) < 8:
                    rec(ex, p2, outs2, w2, ce2)
    rec(t1, [], [], Fr(0), None)
    return res


def main():
    pairs = []
    for i, a in enumerate(CC.ACTIVE):
        for b in CC.ACTIVE[i:]:
            pairs.append((a, b))
    out = {}
    lines = []
    for a, b in pairs:
        W0 = CC.W[a] + CC.W[b]
        chs = comp_chains(a, b, W0 + 3)
        sets = defaultdict(lambda: dict(k_min=99, chains=0))
        for ch in chs:
            dout = sum(CC.DEG[x] for x in ch['outputs'])
            assert dout == CC.DEG[a] + CC.DEG[b] - 1, (a, b, ch)
            assert sum(CC.ODD[x] for x in ch["outputs"]) % 2 == (CC.ODD[a] + CC.ODD[b] + 1) % 2   # the composite kernel is odd (degree -1)
            wout = sum(CC.W[x] for x in ch['outputs'])
            act = [x for x in ch['outputs'] if x in DELTA]
            if not act:
                continue
            if min(wout - DELTA[x] for x in act) > W0 + 1:
                continue
            key = ' '.join(ch['outputs']) if ch['outputs'] else '(none)'
            sets[key]['k_min'] = min(sets[key]['k_min'], ch['k'])
            sets[key]['chains'] += 1
            sets[key]['w_out'] = str(wout)
        out[f'{a}|{b}'] = dict(W0=str(W0), leg_sets=dict(sorted(sets.items(), key=lambda kv: (Fr(kv[1]['w_out']), kv[0]))))
        lines.append('%-5s|%-5s W0=%-4s %s' % (a, b, W0, ', '.join('%s[k>=%d]' % (s, v['k_min'])
                                                                  for s, v in out[f'{a}|{b}']['leg_sets'].items())))
    res = dict(description=__doc__, pairs=out)
    json.dump(res, open(os.path.join(HERE, 'census_composite.json'), 'w'), indent=1)
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
