#!/usr/bin/env python3
"""cubic_census.py -- graph census for the cubic sector of the one-loop obstruction (Step 2, Part C).

Classes (DESIGN_STEP2.md): I-a (special edge K_eps - K_0 on the loop), I-b (special edge = bridge to a one-loop
two-point subgraph), II-a (-Q I^CT_3), III (-{I_3, I^CT_2}_eps), I-d (Step-1 special bubble joined by P to I_3),
II-b (-Q I^CT_2 joined by P to I_3).  Legs are labelled; graphs related by permuting identical legs are grouped
into one row whose 'number' is the size of the orbit.  Weight = 1/|Aut| with the special edge marked (the 1/2 of
1/2({I,I}_eps - {I,I}) cancels the two orderings of its two vertices)."""
import itertools
import json

CUBIC = {
    '(A|_V)^3': ('AV', 'AV', 'AV'),
    'A|_V (A|_Wperp)^2': ('AV', 'AW', 'AW'),
    'lambda lambda A|_V': ('lam', 'lam', 'AV'),
    'lambda lambda A|_Wperp': ('lam', 'lam', 'AW'),
    'A|_V c A^*|_V': ('AV', 'c', 'AsV'),
    'A|_Wperp c A^*|_Wperp': ('AW', 'c', 'AsW'),
    'c lambda lambda^*': ('c', 'lam', 'lams'),
    'c c c^*': ('c', 'c', 'cs'),
}
QUARTIC = {
    '(A|_V)^4': ('AV', 'AV', 'AV', 'AV'),
    '(A|_V)^2 (A|_Wperp)^2': ('AV', 'AV', 'AW', 'AW'),
    '(A|_Wperp)^4': ('AW', 'AW', 'AW', 'AW'),
}
NAME = {'AV': 'A|_V', 'AW': 'A|_Wperp', 'lam': 'lambda', 'c': 'c', 'AsV': 'A^*|_V', 'AsW': 'A^*|_Wperp',
        'lams': 'lambda^*', 'cs': 'c^*'}
PLINE = {('AV', 'AV'): 'A|_V', ('AW', 'AW'): 'A|_Wperp', ('lam', 'lam'): 'lambda', ('AsV', 'c'): 'c--A^*'}
SPECIAL = {('AV', 'AsV'): 'A|_V--A^*|_V', ('AW', 'AsW'): 'A|_Wperp--A^*|_Wperp', ('lam', 'lams'): 'lambda--lambda^*',
           ('c', 'cs'): 'c--c^*'}
PARTNER = {'AV': 'AsV', 'AsV': 'AV', 'AW': 'AsW', 'AsW': 'AW', 'lam': 'lams', 'lams': 'lam', 'c': 'cs', 'cs': 'c'}
SECTORS = {
    '(c, A|_V, A|_V)': ('c', 'AV', 'AV'),
    '(c, A|_Wperp, A|_Wperp)': ('c', 'AW', 'AW'),
    '(c, lambda, lambda)': ('c', 'lam', 'lam'),
    '(c, c, A^*|_V)': ('c', 'c', 'AsV'),
}


def kind(a, b):
    k = tuple(sorted((a, b)))
    if k in PLINE:
        return 'P', PLINE[k]
    if k in SPECIAL:
        return 'S', SPECIAL[k]
    return None, None


def remove_one(slots, s):
    sl = list(slots)
    if s not in sl:
        return None
    sl.remove(s)
    return tuple(sl)


def ordered_pairs(two):
    return sorted(set(itertools.permutations(two)))


def leg_perms(legs):
    """permutations of leg labels preserving leg types."""
    n = len(legs)
    out = []
    for p in itertools.permutations(range(n)):
        if all(legs[p[i]] == legs[i] for i in range(n)):
            out.append(p)
    return out


def canon(graph, perms):
    """graph = (vertices: tuple of (vertex id, type, frozenset legs)), edges: tuple of ((v,slot),(v,slot),kind,label))
    vertex ids are taken as the sorted tuple of their legs (all vertices here carry legs)."""
    best = None
    for p in perms:
        vs = tuple(sorted((tuple(sorted(p[l] for l in legs)), t) for (vid, t, legs) in graph[0]))
        vmap = {vid: tuple(sorted(p[l] for l in legs)) for (vid, t, legs) in graph[0]}
        es = []
        for (a, sa), (b, sb), k, lab in graph[1]:
            ea, eb = (vmap[a], sa), (vmap[b], sb)
            if eb < ea:
                ea, eb = eb, ea
            es.append((ea, eb, k, lab))
        rep = (vs, tuple(sorted(es)))
        if best is None or rep < best:
            best = rep
    return best


def describe(rep, legs):
    vs, es = rep
    vtxt = []
    for legset, t in vs:
        vtxt.append(f"{t} [external {', '.join(NAME[legs[l]] + '#' + str(l) for l in legset)}]")
    etxt = []
    special = None
    for (va, sa), (vb, sb), k, lab in es:
        s = f"{lab}: {NAME[sa]} at vertex{list(va)} -- {NAME[sb]} at vertex{list(vb)}"
        if k == 'S':
            special = s
        else:
            etxt.append(s)
    return vtxt, etxt, special


def class_Ia(legs):
    perms = leg_perms(legs)
    orbits = {}
    # triangles: vertex j carries leg j
    for types in itertools.product(CUBIC.items(), repeat=3):
        rems = []
        ok = True
        for j, (t, sl) in enumerate(types):
            r = remove_one(sl, legs[j])
            if r is None:
                ok = False
                break
            rems.append(r)
        if not ok:
            continue
        for ords in itertools.product(*[ordered_pairs(r) for r in rems]):
            # ords[j] = (slot toward j+1, slot toward j-1)
            edges = []
            nS = 0
            good = True
            for j in range(3):
                jn = (j + 1) % 3
                a, b = ords[j][0], ords[jn][1]
                k, lab = kind(a, b)
                if k is None:
                    good = False
                    break
                nS += k == 'S'
                edges.append(((j, a), (jn, b), k, lab))
            if not good or nS != 1:
                continue
            g = (tuple((j, types[j][0], frozenset([j])) for j in range(3)), tuple(edges))
            rep = canon(g, perms)
            orbits.setdefault(('triangle', rep), set()).add(canon(g, [tuple(range(3))]))
    # bubbles: quartic q carries legs {i,j}, cubic v carries leg l
    for l in range(3):
        qi = [i for i in range(3) if i != l]
        for (tq, sq), (tv, sv) in itertools.product(QUARTIC.items(), CUBIC.items()):
            rq = sq
            for i in qi:
                rq = remove_one(rq, legs[i]) if rq is not None else None
            rv = remove_one(sv, legs[l])
            if rq is None or rv is None:
                continue
            seen = set()
            for perm in set(itertools.permutations(rv)):
                pairs = tuple(sorted(zip(rq, perm)))
                if pairs in seen:
                    continue
                seen.add(pairs)
                ks = [kind(a, b) for a, b in pairs]
                if any(k is None for k, _ in ks) or sum(k == 'S' for k, _ in ks) != 1:
                    continue
                edges = tuple(((0, a), (1, b), k, lab) for (a, b), (k, lab) in zip(pairs, ks))
                g = (((0, tq, frozenset(qi)), (1, tv, frozenset([l]))), edges)
                rep = canon(g, perms)
                orbits.setdefault(('bubble', rep), set()).add(canon(g, [tuple(range(3))]))
    rows = []
    for (topo, rep), members in sorted(orbits.items(), key=lambda x: (x[0][0], str(x[0][1]))):
        vtxt, etxt, special = describe(rep, legs)
        rows.append({'topology': topo, 'vertices': vtxt, 'P_edges': etxt, 'special_edge': special,
                     'number': len(members), 'weight': '1'})
    return rows


def two_point_graphs(a, b):
    """one-loop 1PI two-point graphs with legs (a, b), all edges P(eps,infty); returns list of (description, weight)."""
    out = []
    for (ta, sa), (tb, sb) in itertools.product(CUBIC.items(), repeat=2):
        ra, rb = remove_one(sa, a), remove_one(sb, b)
        if ra is None or rb is None:
            continue
        seen = set()
        for perm in set(itertools.permutations(rb)):
            pairs = tuple(sorted(zip(ra, perm)))
            if pairs in seen:
                continue
            seen.add(pairs)
            ks = [kind(x, y) for x, y in pairs]
            if any(k != 'P' for k, _ in ks):
                continue
            aut = 2 if pairs[0] == pairs[1] else 1
            out.append({'type': 'bubble', 'vertices': [f'{ta} [leg {NAME[a]}]', f'{tb} [leg {NAME[b]}]'],
                        'lines': [lab for _, lab in ks], 'weight': f'1/{aut}' if aut > 1 else '1'})
    # remove double counting of the unordered pair of vertices when a == b
    if a == b:
        uniq = {}
        for g in out:
            key = (tuple(sorted(g['vertices'])), tuple(sorted(g['lines'])))
            uniq[key] = g
        out = list(uniq.values())
    for tq, sq in QUARTIC.items():
        r = remove_one(sq, a)
        r = remove_one(r, b) if r is not None else None
        if r is None:
            continue
        k, lab = kind(*r)
        if k != 'P':
            continue
        aut = 2 if r[0] == r[1] else 1
        out.append({'type': 'seagull', 'vertices': [f'{tq} [legs {NAME[a]}, {NAME[b]}]'], 'lines': [lab],
                    'weight': f'1/{aut}' if aut > 1 else '1'})
    return out


ICT2 = {  # Theorem thm:c4-quadratic, unit = hbar tr_g(ad ad)/(16 pi^2 (X,Y))
    ('AV', 'AV'): ['(9/4) eps^{-1} int(A|_V,A|_V)', '-log(eps) int(dA|_V,dA|_V)'],
    ('AW', 'AW'): ['(9/4) eps^{-1} int(A|_Wperp,A|_Wperp)', '-log(eps) int(dA|_Wperp,dA|_Wperp)'],
    ('lam', 'lam'): ['4 log(eps) int(lambda, dslash lambda)'],
    ('AsV', 'c'): ['-log(eps) int(dc, A^*)'],
}


def bridge_rows(legs, with_sigma):
    perms = leg_perms(legs)
    rows = {}
    for l in range(3):
        ul = [i for i in range(3) if i != l]
        for tu, su in CUBIC.items():
            r = su
            for i in ul:
                r = remove_one(r, legs[i]) if r is not None else None
            if r is None:
                continue
            s = r[0]
            sbar = PARTNER[s]
            key2 = tuple(sorted((sbar, legs[l])))
            # canonical row key: the orbit of l under identical-leg permutations
            orb = tuple(sorted(set(p[l] for p in perms)))
            rk = (tu, s, orb)
            entry = {'cubic_vertex': f"{tu} [external {', '.join(NAME[legs[i]] for i in ul)}]",
                     'special_edge': f"{NAME[s]} (at the cubic vertex) -- {NAME[sbar]}",
                     'other_side_legs': [NAME[sbar], NAME[legs[l]]], 'number': len(orb)}
            if with_sigma:
                entry['two_point_subgraphs'] = two_point_graphs(sbar, legs[l])
            else:
                entry['I^CT_2_terms'] = ICT2.get(key2, [])
            if (with_sigma and entry['two_point_subgraphs']) or ((not with_sigma) and entry['I^CT_2_terms']):
                rows[rk] = entry
    return list(rows.values())


II_A = {  # -Q I^CT_3, log part: multiples of the tree vertices (Step 1) and the substitution of Q used
    '(c, A|_V, A|_V)': [('(A|_V)^3', 6, 'A|_V <- -dc'), ('A|_V c A^*|_V', 1, 'A^*|_V <- -1/2 d^*dA|_V')],
    '(c, A|_Wperp, A|_Wperp)': [('A|_V (A|_Wperp)^2', 6, 'A|_V <- -dc'),
                                ('A|_Wperp c A^*|_Wperp', 1, 'A^*|_Wperp <- -1/2 d^*dA|_Wperp')],
    '(c, lambda, lambda)': [('lambda lambda A|_V', 10, 'A|_V <- -dc'), ('c lambda lambda^*', 1, 'lambda^* <- dslash lambda')],
    '(c, c, A^*|_V)': [('A|_V c A^*|_V', 1, 'A|_V <- -dc'), ('c c c^*', 1, 'c^* <- d^*A^*|_V')],
}
CONNECT = {  # I-d / II-b: two-leg defect (legs c, A|_V) joined by P(eps,infty) to a cubic vertex
    '(c, A|_V, A|_V)': [('c external, A|_V joined by the A|_V line', '(A|_V)^3 [external A|_V, A|_V]', 1),
                        ('A|_V external, c joined by the c--A^* line', 'A|_V c A^*|_V [external A|_V, c]', 2)],
    '(c, A|_Wperp, A|_Wperp)': [('c external, A|_V joined by the A|_V line',
                                 'A|_V (A|_Wperp)^2 [external A|_Wperp, A|_Wperp]', 1)],
    '(c, lambda, lambda)': [('c external, A|_V joined by the A|_V line', 'lambda lambda A|_V [external lambda, lambda]', 1)],
    '(c, c, A^*|_V)': [('c external, A|_V joined by the A|_V line', 'A|_V c A^*|_V [external c, A^*|_V]', 2)],
}
STEP1_SPECIAL = ['(a) special A|_V--A^*|_V, line A|_V', '(b) special A|_V--A^*|_V, line c--A^*',
                 '(c) special A|_Wperp--A^*|_Wperp, line A|_Wperp', '(d) special lambda^*--lambda, line lambda',
                 '(e) special c^*--c, line c--A^*']


def census():
    out = {}
    for sec, legs in SECTORS.items():
        ia = class_Ia(legs)
        out[sec] = {
            'I-a': ia,
            'I-a_count': {'triangle_rows': sum(r['topology'] == 'triangle' for r in ia),
                          'triangle_graphs': sum(r['number'] for r in ia if r['topology'] == 'triangle'),
                          'bubble_rows': sum(r['topology'] == 'bubble' for r in ia),
                          'bubble_graphs': sum(r['number'] for r in ia if r['topology'] == 'bubble')},
            'I-b': bridge_rows(legs, True),
            'III': bridge_rows(legs, False),
            'II-a': [{'tree_vertex': v, 'log_multiple': m, 'Q_substitution': s} for v, m, s in II_A[sec]],
            'I-d': [{'defect_bubble': b, 'connection': c, 'cubic_vertex': v, 'number': n, 'weight': '1'}
                    for (c, v, n) in CONNECT[sec] for b in STEP1_SPECIAL],
            'II-b': [{'two_leg_vertex': '-Q I^CT_2 = (9/2) eps^{-1} unit int(dc, A|_V)', 'connection': c,
                      'cubic_vertex': v, 'number': n, 'weight': '1'} for (c, v, n) in CONNECT[sec]],
        }
    return out


if __name__ == '__main__':
    c = census()
    for sec, d in c.items():
        print(sec, d['I-a_count'], 'I-b rows', len(d['I-b']), 'III rows', len(d['III']))
        for r in d['I-a']:
            print('   ', r['topology'], r['number'], '|', '; '.join(r['vertices']), '|| S:', r['special_edge'],
                  '|| P:', '; '.join(r['P_edges']))
        for r in d['I-b']:
            print('   I-b', r['cubic_vertex'], '|', r['special_edge'], '|', r['other_side_legs'], 'n=', r['number'])
            for g in r['two_point_subgraphs']:
                print('        ', g)
        for r in d['III']:
            print('   III', r['cubic_vertex'], '|', r['special_edge'], '|', r['I^CT_2_terms'], 'n=', r['number'])
    json.dump(c, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory/cubic_census.json', 'w'), indent=1)
