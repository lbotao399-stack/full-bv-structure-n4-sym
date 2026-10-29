"""Ghost / antifield pairs (Grassmann signs by explicit Koszul reordering): referee R vs engine A invariant forms."""
import json, itertools, sympy as sp
from ref_core import *
from inv_eval import ev
KA = json.load(open('../../engine_a/kernels.json'))
forms = KA['invariant_forms_bosonic_W0']['forms']
r1 = list(sp.symbols('r1_0:4')); r2 = list(sp.symbols('r2_0:4')); r3 = list(sp.symbols('r3_0:4'))
V = list(range(4)); W = list(range(4, 10))
def e10(M): return unit(10, M)
SL = {'AV': ('A', V), 'AW': ('A', W), 'AsV': ('As', V), 'AsW': ('As', W), 'c': ('c', [None]), 'cs': ('cs', [None])}
def sel(t, M):
    return [1] if M is None else e10(M)
rep = []
def check(pair, key, la, lb, outs, vb, kind='log', rsl=None):
    f = forms[pair][key]['invariant_form']
    ta, Va = SL[la]; tb, Vb = SL[lb]
    nb = nt = nz = 0
    for A in Va:
        for B in Vb:
            for os_ in itertools.product(*[SL[o][1] for o in outs]):
                ys = [sel(SL[o][0], M) for o, M in zip(outs, os_)]
                verts = vb(ys)
                rs = rsl or [r1, r2, r3][:len(verts)]
                res = kernel_chain(ta, sel(ta, A), tb, sel(tb, B), verts, rs)
                mine = res.get(kind, 0)
                slots = {'a': sel('A', A) if A is not None else None, 'b': sel('A', B) if B is not None else None}
                for i, (o, M) in enumerate(zip(outs, os_)):
                    slots['o%d' % (i + 1)] = e10(M) if M is not None else None
                theirs = ev(f, slots, [r1, r2, r3])
                nt += 1; nz += (theirs != 0)
                if sp.expand(mine - theirs) != 0:
                    nb += 1
                    if nb <= 3: print('  mismatch', pair, key, A, B, os_, mine, '|', theirs)
    print('%-8s %-22s checked %4d nonzero %4d mismatches %d' % (pair, key, nt, nz, nb))
    rep.append(dict(pair=pair, key=key, checked=nt, nonzero=nz, mismatches=nb))
# c | A*_V
check('c|AsV', 'AV | C | log', 'c', 'AsV', ['AV'], lambda y: [('AcAs', 'As', 'c', 'A', y[0])])
check('c|AsV', 'AV | C | eps^-1', 'c', 'AsV', ['AV'], lambda y: [('AcAs', 'As', 'c', 'A', y[0])], kind=1)
check('c|AsV', 'c AsV | CC | log', 'c', 'AsV', ['c', 'AsV'],
      lambda y: [('AcAs', 'As', 'A', 'c', y[0]), ('AcAs', 'A', 'c', 'As', y[1])])
check('c|AsV', 'AV AV | CC | log', 'c', 'AsV', ['AV', 'AV'],
      lambda y: [('AcAs', 'As', 'c', 'A', y[0]), ('AcAs', 'As', 'c', 'A', y[1])])
check('c|AsV', 'AV AV AV | CCC | log', 'c', 'AsV', ['AV', 'AV', 'AV'],
      lambda y: [('AcAs', 'As', 'c', 'A', y[0]), ('AcAs', 'As', 'c', 'A', y[1]), ('AcAs', 'As', 'c', 'A', y[2])])
# c | c
check('c|c', 'c c | CC | log', 'c', 'c', ['c', 'c'],
      lambda y: [('AcAs', 'As', 'A', 'c', y[0]), ('AcAs', 'A', 'As', 'c', y[1])])
# A_V | c
check('AV|c', 'c | C | log', 'AV', 'c', ['c'], lambda y: [('AcAs', 'A', 'As', 'c', y[0])])
check('AV|c', 'AV c | CC | log', 'AV', 'c', ['AV', 'c'], lambda y: [('AAA', y[0]), ('AcAs', 'A', 'As', 'c', y[1])])
check('AV|c', 'c AV | CC | log', 'AV', 'c', ['c', 'AV'],
      lambda y: [('AcAs', 'A', 'As', 'c', y[0]), ('AcAs', 'c', 'As', 'A', y[1])])
# A_perp | c
check('AW|c', 'AW c | CC | log', 'AW', 'c', ['AW', 'c'], lambda y: [('AAA', y[0]), ('AcAs', 'A', 'As', 'c', y[1])])
# A_V | A*_V
check('AV|AsV', 'AsV | C | log', 'AV', 'AsV', ['AsV'], lambda y: [('AcAs', 'A', 'c', 'As', y[0])])
check('AV|AsV', 'c cs | CC | log', 'AV', 'AsV', ['c', 'cs'],
      lambda y: [('AcAs', 'A', 'As', 'c', y[0]), ('cccs', 'c', 'c', 'cs', y[1])])
check('AV|AsV', 'AV AsV | CC | log', 'AV', 'AsV', ['AV', 'AsV'], lambda y: [('AAA', y[0]), ('AcAs', 'A', 'c', 'As', y[1])])
check('AV|AsV', 'AsV AV | CC | log', 'AV', 'AsV', ['AsV', 'AV'],
      lambda y: [('AcAs', 'A', 'c', 'As', y[0]), ('AcAs', 'As', 'c', 'A', y[1])])
check('AV|AsV', 'AW AsW | CC | log', 'AV', 'AsV', ['AW', 'AsW'], lambda y: [('AAA', y[0]), ('AcAs', 'A', 'c', 'As', y[1])])
# A_perp | A*_V
check('AW|AsV', 'AsW | C | log', 'AW', 'AsV', ['AsW'], lambda y: [('AcAs', 'A', 'c', 'As', y[0])])
check('AW|AsV', 'AV AsW | CC | log', 'AW', 'AsV', ['AV', 'AsW'], lambda y: [('AAA', y[0]), ('AcAs', 'A', 'c', 'As', y[1])])
check('AW|AsV', 'AW AsV | CC | log', 'AW', 'AsV', ['AW', 'AsV'], lambda y: [('AAA', y[0]), ('AcAs', 'A', 'c', 'As', y[1])])
check('AW|AsV', 'AsW AV | CC | log', 'AW', 'AsV', ['AsW', 'AV'],
      lambda y: [('AcAs', 'A', 'c', 'As', y[0]), ('AcAs', 'As', 'c', 'A', y[1])])
# A*_V | A*_V
check('AsV|AsV', 'cs | C | log', 'AsV', 'AsV', ['cs'], lambda y: [('cccs', 'c', 'c', 'cs', y[0])])
check('AsV|AsV', 'cs | C | eps^-1', 'AsV', 'AsV', ['cs'], lambda y: [('cccs', 'c', 'c', 'cs', y[0])], kind=1)
check('AsV|AsV', 'AsV AsV | CC | log', 'AsV', 'AsV', ['AsV', 'AsV'],
      lambda y: [('AcAs', 'c', 'A', 'As', y[0]), ('AcAs', 'A', 'c', 'As', y[1])])
check('AsV|AsV', 'AV cs | CC | log', 'AsV', 'AsV', ['AV', 'cs'],
      lambda y: [('AcAs', 'c', 'As', 'A', y[0]), ('cccs', 'c', 'c', 'cs', y[1])])
json.dump(rep, open('out_ghost.json', 'w'), indent=1)
