"""Census of the special-edge composite cycles through two letters (composite_theory.tex, Sec. 1).

A special-edge composite cycle through the letters (X_a, X_b) is a chain
    X_a - e_0 - V_1 - e_1 - ... - V_k - e_k - X_b
with k vertices of I, exactly one edge equal to the composite line (components c (x) A|_V, Prop c9-composite (ii);
along Delta it carries S^(1)(Delta) between an A-A and a c-A^*|_V edge, along P_b the same with S^(1)(P_b)),
exactly one edge equal to the special edge K_eps - K_0 (components of K_0: field (x) its antifield,
Prop c4-defect1), and all other edges P(eps, infty) (components A-A, lambda_+ - lambda_-, c - A^*|_V, Prop c3-P).

Slot rule.  Every edge between slots of Weyl weights w1, w2 has momentum degree w1 + w2 - 4 (P: numerator w1+w2-2
over |p|^2; special edge: K_0 has degree 0; composite: c - A|_V gives -3).  Hence the superficial degree of the
chain is D = w_a + w_b - w_out (Lemma lem:c10-degree, same proof).  The eps^0 coefficient of the cycle is a
polynomial of degree D in the momenta for the translation-invariant part along Delta and of degree D + 1 for the
part along P_b (the factor i p_b of S^(1)(P_b)) and for the part of C_Delta linear in the position (twice the P_b
part, Prop prop:cc-position); it vanishes in the limit for negative degree (Lemma lem:cc-special).
Charges: composite line odd of degree -1, special edge odd of degree +1, so the outputs carry the pair degree and the
parity of the two letters; the centre Z/4 of SL(4) is conserved.

Output: census_sp.json, census_sp_summary.txt (written next to this file and copied by the caller).
"""
import json, os, sys
from fractions import Fraction as Fr
from collections import defaultdict, Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tl_census_base as CC   # verbatim copy of step2/twoletter/theory/checks_census.py (vertices, weights, P edges)

COMP = {'c': 'AV', 'AV': 'c'}                                     # composite line c (x) A|_V
SPEC = {'AV': 'AsV', 'AsV': 'AV', 'AW': 'AsW', 'AsW': 'AW', 'lamP': 'lamsP', 'lamsP': 'lamP',
        'lamM': 'lamsM', 'lamsM': 'lamM', 'c': 'cs', 'cs': 'c'}     # K_0: field (x) antifield
SLOTS = ['AV', 'AW', 'lamP', 'lamM', 'c', 'AsV', 'AsW', 'lamsP', 'lamsM', 'cs']
FIELDS = ['AV', 'AW', 'lamP', 'lamM']
ACTIVE = CC.ACTIVE
KMAX = 6


def all_moves(enter):
    """(vertex, exit slot, outputs) for a vertex entered through a slot of type `enter` (any exit slot)."""
    out = []
    for name, (slots, nd) in CC.VERTS.items():
        sl = list(slots)
        if enter not in sl:
            continue
        rest = list(sl); rest.remove(enter)
        seen = set()
        for ex in rest:
            if ex in seen:
                continue
            seen.add(ex)
            outs = list(rest); outs.remove(ex)
            out.append((name, ex, tuple(sorted(outs)), nd))
    return out


def edge_options(slot, used_comp, used_spec):
    """edges leaving a slot: (kind, partner slot)."""
    opts = []
    if slot in CC.PARTNER:
        opts.append(('P', CC.PARTNER[slot]))
    if not used_comp and slot in COMP:
        opts.append(('C', COMP[slot]))
    if not used_spec and slot in SPEC:
        opts.append(('S', SPEC[slot]))
    return opts


def chains(t1, t2, wmax):
    res = []

    def close(slot, path, outs, kinds):
        for kind, partner in edge_options(slot, 'C' in kinds, 'S' in kinds):
            if partner != t2:
                continue
            ks = kinds + [kind]
            if ks.count('C') == 1 and ks.count('S') == 1:
                res.append(dict(k=len(path), path=list(path), outputs=tuple(sorted(outs)), edges=ks))

    def rec(slot, path, outs, wout, kinds):
        close(slot, path, outs, kinds)
        if len(path) >= KMAX:
            return
        for kind, enter in edge_options(slot, 'C' in kinds, 'S' in kinds):
            for (name, ex, o, nd) in all_moves(enter):
                w2 = wout + sum(CC.W[x] for x in o)
                if w2 > wmax:
                    continue
                rec(ex, path + [(name, enter, ex)], outs + list(o), w2, kinds + [kind])
    rec(t1, [], [], Fr(0), [])
    return res


def canon(ch, sym):
    p = (tuple(ch['path']), tuple(ch['edges']))
    if sym:
        rev = (tuple((n, ex, en) for (n, en, ex) in reversed(ch['path'])), tuple(reversed(ch['edges'])))
        p = min(p, rev)
    return (p, ch['outputs'])


def main():
    pairs = [(a, b) for i, a in enumerate(SLOTS) for b in SLOTS[i:]]
    out = {}
    lines = []
    NDER = 3
    for a, b in pairs:
        W0 = CC.W[a] + CC.W[b]
        d_ab = CC.DEG[a] + CC.DEG[b]
        chs = chains(a, b, W0 + NDER + 1)
        uniq = {}
        for ch in chs:
            uniq[canon(ch, a == b)] = ch
        rows = []
        for ch in uniq.values():
            wout = sum(CC.W[x] for x in ch['outputs'])
            assert sum(CC.DEG[x] for x in ch['outputs']) == d_ab, (a, b, ch)
            assert sum(CC.ODD[x] for x in ch['outputs']) % 2 == (CC.ODD[a] + CC.ODD[b]) % 2, (a, b, ch)
            assert sum(CC.Z4[x] for x in ch['outputs']) % 4 == (CC.Z4[a] + CC.Z4[b]) % 4, (a, b, ch)
            nout = len(ch['outputs'])
            rows.append(dict(k=ch['k'], wout=wout, nout=nout, outputs=ch['outputs'],
                             path=['%s[%s>%s]' % p for p in ch['path']], edges=''.join(ch['edges'])))
        if not rows:
            continue
        per_n = {}
        for n in range(NDER + 1):
            W = W0 + n
            sets_P = defaultdict(lambda: dict(chains=0, k_min=99))
            sets_0 = defaultdict(lambda: dict(chains=0, k_min=99))
            for r in rows:
                key = ' '.join(r['outputs']) if r['outputs'] else '(none)'
                if r['wout'] <= W + 1:           # P_b part / position part of C_Delta: degree D + 1 >= 0
                    sets_P[key]['chains'] += 1
                    sets_P[key]['k_min'] = min(sets_P[key]['k_min'], r['k'])
                    sets_P[key]['w_out'] = str(r['wout'])
                if r['wout'] <= W:               # translation-invariant part along Delta: degree D >= 0
                    sets_0[key]['chains'] += 1
                    sets_0[key]['k_min'] = min(sets_0[key]['k_min'], r['k'])
                    sets_0[key]['w_out'] = str(r['wout'])
            # output bound n_out <= W + 1 + d_ab (Delta), W + 1 + d_ab (P_b)
            viol = [r for r in rows if r['wout'] <= W + 1 and r['nout'] > W + 1 + d_ab]
            per_n[str(n)] = dict(W=str(W), n_out_bound=str(W + 1 + d_ab),
                                 leg_sets_Pb_and_position=dict(sorted(sets_P.items(), key=lambda kv: (Fr(kv[1]['w_out']), kv[0]))),
                                 leg_sets_translation_invariant=dict(sorted(sets_0.items(), key=lambda kv: (Fr(kv[1]['w_out']), kv[0]))),
                                 bound_violations=len(viol))
        out[f'{a}|{b}'] = dict(W0=str(W0), pair_degree=d_ab, field_pair=(a in FIELDS and b in FIELDS),
                               active_pair=(a in ACTIVE and b in ACTIVE), by_letter_derivatives=per_n)
        s0 = per_n['0']['leg_sets_Pb_and_position']
        lines.append('%-5s|%-5s W0=%-4s d=%2d  n=0: %s' % (a, b, W0, d_ab,
                     ', '.join('%s[k>=%d,%d]' % (s, v['k_min'], v['chains']) for s, v in s0.items())))
    summary = dict(n_pair_types_with_cycles=len(out),
                   field_pairs=sorted(k for k, v in out.items() if v['field_pair']),
                   counts_by_n={k: {n: len(v['by_letter_derivatives'][n]['leg_sets_Pb_and_position'])
                                    for n in v['by_letter_derivatives']} for k, v in out.items()},
                   bound_violations=sum(v['by_letter_derivatives'][n]['bound_violations']
                                        for v in out.values() for n in v['by_letter_derivatives']))
    res = dict(description=__doc__, kmax=KMAX, pairs=out, summary=summary)
    json.dump(res, open(os.path.join(HERE, 'census_sp.json'), 'w'), indent=1)
    open(os.path.join(HERE, 'census_sp_summary.txt'), 'w').write('\n'.join(lines) + '\n' + json.dumps(summary, indent=1))
    print('\n'.join(lines))
    print(json.dumps(summary['counts_by_n'], indent=0)[:3000])
    print('bound violations', summary['bound_violations'])


if __name__ == '__main__':
    main()
