"""Independent re-derivation (referee) of the census of one-loop cycles through two letters.

Written from the manuscript data only (Ch1 eq. c1-I and its component table, Prop c1-chirality, Ch3 Prop c3-P);
it does not import theory/checks_census.py.  Algorithm: breadth-first search over (exit slot, output multiset),
which is a different enumeration from the theory's depth-first path search, followed by an independent count of
the leg sets of Table tab:tl-A / tab:tl-B and of the k = 0 self-contraction powers.

Output: indep_census.json (comparison with theory/census.json) and a printed summary.
"""
import json, os, itertools
from fractions import Fraction as Fr
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
THEORY = os.path.join(HERE, '..', '..', 'theory')

# slot data (Ch5 Def c5-family weights; degrees of coordinate functionals = minus field degree)
WT = dict(c=Fr(0), AV=Fr(1), AW=Fr(1), lamP=Fr(3, 2), lamM=Fr(3, 2), AsV=Fr(3), AsW=Fr(3), lamsP=Fr(5, 2),
          lamsM=Fr(5, 2), cs=Fr(4))
DG = dict(c=1, AV=0, AW=0, lamP=0, lamM=0, AsV=-1, AsW=-1, lamsP=-1, lamsM=-1, cs=-2)
# vertices of I with their slot multisets (Ch1 table after eq c1-I; chiralities from Prop c1-chirality:
# rho(v) exchanges chirality -> lambda lambda A|_V pairs lamP with lamM; rho(w) keeps it -> lamP lamP A_perp,
# lamM lamM A_perp; ([c,lambda],lambda^*) pairs lambda in S_+ (x) W with lambda^* in S_+ (x) W^* (= lamsP)).
VERTS = [('AV', 'AV', 'AV'), ('AV', 'AW', 'AW'), ('lamP', 'lamM', 'AV'), ('lamP', 'lamP', 'AW'), ('lamM', 'lamM', 'AW'),
         ('AV', 'c', 'AsV'), ('AW', 'c', 'AsW'), ('c', 'lamP', 'lamsP'), ('c', 'lamM', 'lamsM'), ('c', 'c', 'cs'),
         ('AV',) * 4, ('AV', 'AV', 'AW', 'AW'), ('AW',) * 4]
# nonzero components of P (Prop c3-P): which slot a coordinate functional of slot t is contracted with
EDGE = dict(AV='AV', AW='AW', lamP='lamM', lamM='lamP', c='AsV', AsV='c')
ACTIVE = ['AV', 'AW', 'lamP', 'lamM', 'c', 'AsV']
WMAX = Fr(7)


def moves(enter):
    """(exit slot, outputs) for each way to traverse a vertex entered through `enter`"""
    res = set()
    for v in VERTS:
        if enter not in v:
            continue
        rest = list(v); rest.remove(enter)
        for i, ex in enumerate(rest):
            if ex not in EDGE:
                continue
            outs = tuple(sorted(rest[:i] + rest[i + 1:]))
            res.add((ex, outs))
    return sorted(res)


def legsets(a, b, wmax):
    """set of output multisets (with k>=1) of chains a -> b, output weight <= wmax; BFS on states"""
    out = set()
    start = (a, ())
    frontier = {start}
    seen = set()
    while frontier:
        new = set()
        for (slot, outs) in frontier:
            if (slot, outs) in seen:
                continue
            seen.add((slot, outs))
            enter = EDGE[slot]
            for ex, o in moves(enter):
                outs2 = tuple(sorted(outs + o))
                if sum(WT[x] for x in outs2) > wmax:
                    continue
                if EDGE.get(ex) == b:
                    out.add(outs2)
                new.add((ex, outs2))
        frontier = new
    return out


def main():
    cen = json.load(open(os.path.join(THEORY, 'census.json')))['output_leg_sets_by_pair_weight']
    report = dict(pairs={}, mismatches=[], k0_eps_rows=[])
    for i, a in enumerate(ACTIVE):
        for b in ACTIVE[i:]:
            key = f'{a}|{b}'
            L = legsets(a, b, WMAX)
            W0 = WT[a] + WT[b]
            rows = {}
            n = 0
            while W0 + n <= WMAX:
                Wp = W0 + n
                logsets = sorted(' '.join(s) for s in L if sum(WT[x] for x in s) <= Wp)
                th = cen[key][str(Wp)]['log']
                if sorted(th) != logsets:
                    report['mismatches'].append(dict(pair=key, W=str(Wp), indep_only=sorted(set(logsets) - set(th)),
                                                     theory_only=sorted(set(th) - set(logsets))))
                # bound n_out <= W + d_ab
                d = DG[a] + DG[b]
                nmax = max([len(s.split()) for s in logsets] + [0])
                # k = 0 self-contraction: exists iff EDGE[a] == b; it gives eps^{-W/2} only, zero for odd W
                k0 = EDGE.get(a) == b
                k0_m = (Wp / 2) if (k0 and (Wp.denominator == 1 and Wp.numerator % 2 == 0)) else None
                th_eps = cen[key][str(Wp)]['eps_max_m']
                if k0 and '(none)' in th_eps:
                    report['k0_eps_rows'].append(dict(pair=key, W=str(Wp), theory_lists_m_up_to=th_eps['(none)'],
                                                      actual_powers=('eps^-%s' % k0_m) if k0_m else 'none (W odd)',
                                                      theory_lists_eps_minus_1=(th_eps['(none)'] >= 1),
                                                      eps_minus_1_actually_present=(k0_m == 1)))
                rows[str(Wp)] = dict(n_log=len(logsets), max_out=nmax, bound=str(Wp + d), bound_ok=nmax <= Wp + d)
                n += 1
            report['pairs'][key] = rows
    # k = 1 existence (Prop census (iii))
    k1 = {}
    for i, a in enumerate(ACTIVE):
        for b in ACTIVE[i:]:
            k1[f'{a}|{b}'] = any(EDGE.get(ex) == b for ex, o in moves(EDGE[a]))
    report['k1_exists'] = k1
    report['k1_absent'] = [k for k, v in k1.items() if not v]
    json.dump(report, open(os.path.join(HERE, 'indep_census.json'), 'w'), indent=1)
    print('leg-set mismatches with theory census.json:', len(report['mismatches']))
    for m in report['mismatches'][:10]:
        print('  ', m)
    print('k=1 absent for:', report['k1_absent'])
    bad = [r for r in report['k0_eps_rows'] if r['theory_lists_eps_minus_1'] and not r['eps_minus_1_actually_present']]
    print('rows where census/Table lists the empty leg set with eps^-1 but the self-contraction has no eps^-1:')
    for r in bad:
        print('  ', r)
    allb = all(v['bound_ok'] for rows in report['pairs'].values() for v in rows.values())
    print('output bound n_out <= W + d_ab holds on all rows:', allb)
    # Table B as printed
    for key in ['AV|AV', 'AV|c', 'c|c', 'lamP|lamM', 'AsV|AsV']:
        print(key, {w: (r['n_log'], r['max_out']) for w, r in report['pairs'][key].items()})


if __name__ == '__main__':
    main()
