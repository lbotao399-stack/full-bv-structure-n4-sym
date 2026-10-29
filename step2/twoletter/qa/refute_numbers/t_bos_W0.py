import json, itertools, sympy as sp
from ref_core import *
from inv_eval import ev
KA = json.load(open('../../engine_a/kernels.json'))
forms = KA['invariant_forms_bosonic_W0']['forms']
r1 = list(sp.symbols('r1_0:4')); r2 = list(sp.symbols('r2_0:4'))
Vc = list(range(4)); Wc = list(range(4, 10))
comps = {'AV': Vc, 'AW': Wc}
def e10(M): return unit(10, M)
results = []
def check(pair, key, ta_t, tb_t, verts_builder, outs_t, kinds):
    fa = forms[pair].get(key)
    if fa is None:
        print('missing', pair, key); return
    nbad = 0; ntot = 0
    for A in comps[ta_t]:
        for B in comps[tb_t]:
            for os_ in itertools.product(*[comps[t] for t in outs_t]):
                ys = [e10(M) for M in os_]
                verts, rs = verts_builder(ys)
                res = kernel_chain('A', e10(A), 'A', e10(B), verts, rs)
                mine = res.get(kinds, 0)
                slots = {'a': e10(A), 'b': e10(B)}
                for i, y in enumerate(ys):
                    slots['o%d' % (i + 1)] = y
                theirs = ev(fa['invariant_form'], slots, [r1, r2])
                ntot += 1
                if sp.expand(mine - theirs) != 0:
                    nbad += 1
                    if nbad < 3: print('  mismatch', pair, key, A, B, os_, mine, '|', theirs)
    print(pair, key, 'checked', ntot, 'mismatches', nbad)
    results.append((pair, key, ntot, nbad))
C1 = lambda ys: ([('AAA', ys[0])], [r1])
CC = lambda ys: ([('AAA', ys[0]), ('AAA', ys[1])], [r1, r2])
Q = lambda ys: ([('AAAA', ys[0], ys[1])], [[r1[m] + r2[m] for m in range(4)]])
for pair, ta, tb in [('AV|AV', 'AV', 'AV'), ('AW|AW', 'AW', 'AW'), ('AV|AW', 'AV', 'AW')]:
    for key in forms[pair]:
        legs, topo, kind = [x.strip() for x in key.split('|')]
        if legs == '(none)':
            continue
        outs = legs.split()
        b = {'C': C1, 'CC': CC, 'Q': Q}[topo]
        check(pair, key, ta, tb, b, outs, 'log' if kind == 'log' else int(kind.split('^-')[1]))
# self-contraction
for pair, t in [('AV|AV', 'AV'), ('AW|AW', 'AW')]:
    f = forms[pair]['(none) | 0 | eps^-1']['invariant_form']
    bad = 0
    for A in comps[t]:
        for B in comps[t]:
            mine = kernel_k0('A', e10(A), 'A', e10(B)).get(1, 0)
            theirs = ev(f, {'a': e10(A), 'b': e10(B)}, [])
            bad += (sp.expand(mine - theirs) != 0)
    print(pair, 'k=0 eps^-1 mismatches', bad)
json.dump([list(map(str, r)) for r in results], open('out_bos_W0.json', 'w'), indent=1)
