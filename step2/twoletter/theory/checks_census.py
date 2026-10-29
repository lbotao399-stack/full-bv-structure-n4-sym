"""Census of the one-loop cycles through two letters of a local operator (two-letter kernels).

Data (manuscript Ch1, Ch3):
  slot types: c, AV (=A|_V), AW (=A|_{W_R^perp}), lamP (lambda in S_+ (x) W), lamM (lambda in S_- (x) W^*),
              AsV (=A^*|_V), AsW (=A^*|_{W_R^perp}), lamsP (lambda^* in S_+ (x) W^*), lamsM (lambda^* in S_- (x) W), cs (=c^*)
  Weyl weights (Def c6-sigma): c 0, A 1, lambda 3/2, A^* 3, lambda^* 5/2, c^* 4
  functional degrees (coordinate functionals, Lemma c4-enumeration): c +1, A 0, lambda 0, A^* -1, lambda^* -1, c^* -2
  vertices of I (eq c1-I, table of component types, Prop c1-chirality for the chiral structure):
     (AV,AV,AV) 1 derivative; (AV,AW,AW) 1 derivative; (lamP,lamM,AV); (lamP,lamP,AW); (lamM,lamM,AW);
     (AV,c,AsV); (AW,c,AsW); (c,lamP,lamsP); (c,lamM,lamsM); (c,c,cs); (AV)^4; (AV)^2(AW)^2; (AW)^4
  nonzero components of P(eps,T) (Prop c3-P): AV-AV, AW-AW, lamP-lamM, c-AsV (numerator degrees 0,0,1,1).
A cycle through two letters (t1,t2) of an operator is a chain  letter1 -P- V_1 -P- V_2 ... V_k -P- letter2:
the half-edge of letter 1 is joined by a component of P to a slot of V_1, V_1 is left through a second slot, ...,
the last slot of V_k is joined to letter 2.  The remaining slots of the V_i are the output legs.
Superficial degree = w_1 + w_2 - (sum of the Weyl weights of the output legs), w_i = Weyl weight of letter i
including its derivatives (proof in twoletter_theory.tex, Lemma lem:tl-degree); it is checked here against the
direct count  deg(numerator) + 4 - 2E.
Output: census.json (all chains with output weight <= WMAX) and a printed summary."""
import json, itertools, os
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
W = {'c': Fr(0), 'AV': Fr(1), 'AW': Fr(1), 'lamP': Fr(3, 2), 'lamM': Fr(3, 2), 'AsV': Fr(3), 'AsW': Fr(3),
     'lamsP': Fr(5, 2), 'lamsM': Fr(5, 2), 'cs': Fr(4)}
DEG = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 0, 'lamM': 0, 'AsV': -1, 'AsW': -1, 'lamsP': -1, 'lamsM': -1, 'cs': -2}
ODD = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 1, 'lamM': 1, 'AsV': 1, 'AsW': 1, 'lamsP': 0, 'lamsM': 0, 'cs': 0}
Z4 = {'c': 0, 'AV': 0, 'AW': 2, 'lamP': 1, 'lamM': 3, 'AsV': 0, 'AsW': 2, 'lamsP': 3, 'lamsM': 1, 'cs': 0}
VERTS = {  # name: (slots, number of derivatives)
    'AV3': (('AV', 'AV', 'AV'), 1), 'AVAW2': (('AV', 'AW', 'AW'), 1),
    'lamlamAV': (('lamP', 'lamM', 'AV'), 0), 'lamPlamPAW': (('lamP', 'lamP', 'AW'), 0),
    'lamMlamMAW': (('lamM', 'lamM', 'AW'), 0),
    'AVcAsV': (('AV', 'c', 'AsV'), 0), 'AWcAsW': (('AW', 'c', 'AsW'), 0),
    'clamPlamsP': (('c', 'lamP', 'lamsP'), 0), 'clamMlamsM': (('c', 'lamM', 'lamsM'), 0),
    'ccCs': (('c', 'c', 'cs'), 0),
    'AV4': (('AV', 'AV', 'AV', 'AV'), 0), 'AV2AW2': (('AV', 'AV', 'AW', 'AW'), 0), 'AW4': (('AW', 'AW', 'AW', 'AW'), 0)}
PARTNER = {'AV': 'AV', 'AW': 'AW', 'lamP': 'lamM', 'lamM': 'lamP', 'c': 'AsV', 'AsV': 'c'}
EDGE_NUMDEG = {frozenset(['AV']): 0, frozenset(['AW']): 0, frozenset(['lamP', 'lamM']): 1, frozenset(['c', 'AsV']): 1}
ACTIVE = ['AV', 'AW', 'lamP', 'lamM', 'c', 'AsV']          # letter types with a component of P
INERT = ['AsW', 'lamsP', 'lamsM', 'cs']                     # letter types without a component of P
WMAX = Fr(7)


def edge_deg(a, b):
    return EDGE_NUMDEG[frozenset([a, b])]


def vertex_moves(enter):
    """all (vertex, exit slot type, outputs) with a slot of type `enter` used to enter"""
    out = []
    for name, (slots, nd) in VERTS.items():
        sl = list(slots)
        if enter not in sl:
            continue
        rest = list(sl); rest.remove(enter)
        seen = set()
        for ex in rest:
            if ex in seen:
                continue
            seen.add(ex)
            if ex not in PARTNER:
                continue
            outs = list(rest); outs.remove(ex)
            out.append((name, ex, tuple(sorted(outs)), nd))
    return out


def chains(t1, t2, wmax):
    """DFS over chains from letter t1 to letter t2 with output weight <= wmax"""
    res = []
    # k = 0
    if PARTNER.get(t1) == t2:
        res.append(dict(k=0, path=[], outputs=(), numdeg=edge_deg(t1, t2)))

    def rec(cur_slot_type, path, outs, wout, numdeg):
        # cur_slot_type: type of the slot through which we left the last vertex (or the letter 1 type at start)
        enter = PARTNER[cur_slot_type]
        ndeg_edge = edge_deg(cur_slot_type, enter)
        for (name, ex, o, nd) in vertex_moves(enter):
            w2 = wout + sum(W[x] for x in o)
            if w2 > wmax:
                continue
            p2 = path + [(name, enter, ex)]
            outs2 = outs + list(o)
            nd2 = numdeg + ndeg_edge + nd
            # close the chain at letter 2
            if PARTNER.get(ex) == t2:
                res.append(dict(k=len(p2), path=p2, outputs=tuple(sorted(outs2)), numdeg=nd2 + edge_deg(ex, t2)))
            if len(p2) < 10:
                rec(ex, p2, outs2, w2, nd2)
    if t1 in PARTNER:
        rec(t1, [], [], Fr(0), 0)
    return res


def canonical(ch, symmetric=False):
    p = tuple(ch['path'])
    if symmetric:
        rev = tuple((name, ex, en) for (name, en, ex) in reversed(ch['path']))
        p = min(p, rev)
    return (p, ch['outputs'])


def main():
    pairs = []
    for i, a in enumerate(ACTIVE):
        for b in ACTIVE[i:]:
            pairs.append((a, b))
    census = {}
    summary = []
    deg_check_fail = 0
    for (a, b) in pairs:
        chs = chains(a, b, WMAX)
        # the reversed chain (b -> a) is the same cycle: keep one orientation (letter a first)
        uniq = {}
        for ch in chs:
            uniq[canonical(ch, symmetric=(a == b))] = ch
        rows = []
        for ch in uniq.values():
            wout = sum(W[x] for x in ch['outputs'])
            E = ch['k'] + 1
            base_sup = W[a] + W[b] - wout              # superficial degree without derivatives on the letters
            # direct count: deg(numerator) + 4 - 2E with no derivatives on the letters
            direct = ch['numdeg'] + 4 - 2 * E
            if Fr(direct) != base_sup:
                deg_check_fail += 1
            dout = sum(DEG[x] for x in ch['outputs'])
            assert dout == DEG[a] + DEG[b], (a, b, ch)
            assert sum(Z4[x] for x in ch['outputs']) % 4 == (Z4[a] + Z4[b]) % 4
            assert sum(ODD[x] for x in ch['outputs']) % 2 == (ODD[a] + ODD[b]) % 2
            nmin = max(Fr(0), -base_sup)               # derivatives needed on the letters for degree >= 0
            rows.append(dict(k=ch['k'], vertices=[p[0] for p in ch['path']],
                             cycle_slots=[[p[1], p[2]] for p in ch['path']],
                             outputs=list(ch['outputs']), n_out=len(ch['outputs']),
                             output_weight=str(wout), base_superficial_degree=str(base_sup),
                             min_letter_derivatives_for_log=str(nmin)))
        rows.sort(key=lambda r: (r['k'], Fr(r['output_weight']), r['vertices']))
        census[f'{a}|{b}'] = rows
        byk = Counter(r['k'] for r in rows)
        summary.append((a, b, len(rows), dict(sorted(byk.items()))))
    # output leg sets per pair type and pair weight Wp (log: output weight = Wp; eps^-m: Wp - 2m)
    legsets = {}
    for key, rows in census.items():
        a, b = key.split('|')
        tab = defaultdict(set)
        for r in rows:
            wout = Fr(r['output_weight'])
            tab[str(wout)].add(' '.join(r['outputs']) if r['outputs'] else '(none)')
        legsets[key] = {wk: sorted(v) for wk, v in sorted(tab.items(), key=lambda kv: Fr(kv[0]))}
    # output leg sets per pair type and number n of derivatives on the two letters (pair weight Wp = w(t1)+w(t2)+n):
    # log eps part: output field weight <= Wp; eps^{-m} part: output field weight <= Wp - 2m
    legsets_n = {}
    for key, rows in census.items():
        a, b = key.split('|')
        d = {}
        for n in range(0, 5):
            Wp = W[a] + W[b] + n
            sets = sorted(set((' '.join(r['outputs']) if r['outputs'] else '(none)', r['output_weight'])
                              for r in rows if Fr(r['output_weight']) <= Wp))
            d[str(n)] = dict(pair_weight=str(Wp),
                             log=[s_ for (s_, w_) in sets],
                             max_m={s_: str((Wp - Fr(w_)) // 2) for (s_, w_) in sets})
        legsets_n[key] = d
    # bound check: n_out <= W_pair + D for the log part (W_pair = output weight), all rows
    bound_viol = 0
    for key, rows in census.items():
        a, b = key.split('|')
        D = DEG[a] + DEG[b]
        for r in rows:
            if r['n_out'] > Fr(r['output_weight']) + D:
                bound_viol += 1
            if r['k'] > r['n_out']:
                bound_viol += 1
    WDET = Fr(5)
    compact = {}
    for key, rows in census.items():
        compact[key] = ['k=%d | %s | out: %s | w_out=%s | base_sup=%s' % (
            r['k'], ', '.join('%s[%s>%s]' % (v, cs[0], cs[1]) for v, cs in zip(r['vertices'], r['cycle_slots'])),
            ' '.join(r['outputs']) if r['outputs'] else '-', r['output_weight'], r['base_superficial_degree'])
            for r in rows if Fr(r['output_weight']) <= WDET]
    counts = {}
    for key, rows in census.items():
        c = Counter((r['k'], r['output_weight']) for r in rows)
        counts[key] = {'k=%d,w_out=%s' % kw: n for kw, n in sorted(c.items(), key=lambda t: (t[0][0], Fr(t[0][1])))}
    out = dict(description=__doc__, WMAX=str(WMAX), WMAX_listed=str(WDET), active_letter_types=ACTIVE,
               inert_letter_types=INERT,
               degree_formula_mismatches=deg_check_fail, bound_violations=bound_viol,
               chains_listed=compact, chain_counts=counts, output_leg_sets_by_letter_derivatives=legsets_n,
               summary=[dict(pair=f'{a}|{b}', n_chains=n, by_k=bk) for (a, b, n, bk) in summary])
    with open(os.path.join(HERE, 'census.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print('degree formula mismatches:', deg_check_fail, ' bound violations:', bound_viol)
    for (a, b, n, bk) in summary:
        print(f'{a:5s}|{b:5s}  chains(out weight<={WMAX}): {n:4d}  by k: {bk}')
    # print the minimal chains (output weight <= 3) for the main pairs
    for key in ['AW|AW', 'AV|AV', 'AV|AW', 'lamP|lamM', 'lamP|lamP', 'AV|c', 'c|AsV', 'AW|lamP', 'AV|lamP', 'c|c']:
        if key not in census:
            continue
        print('==', key)
        for r in census[key]:
            if Fr(r['output_weight']) <= 2 or r['k'] <= 1:
                print('   k=%d  %-40s outputs=%-28s w_out=%s base_sup=%s' % (
                    r['k'], ','.join(r['vertices']), ' '.join(r['outputs']), r['output_weight'],
                    r['base_superficial_degree']))


if __name__ == '__main__':
    main()
