"""Engine A: index-free (invariant) form of the bosonic W_0 kernels.  For pair types and words whose letters and outputs
are all bosonic with vector indices (AV, AsV: SO(4) index m; AW, AsW: SO(6) index I; c, cs: none), every table entry
polynomial is fitted exactly onto the basis of O(4) x O(6) invariant structures (Kronecker deltas among the slots,
slots contracted with output momenta r_j, scalar products r_i.r_j) and the fit is verified on every entry.  Structures
are printed as d(x,y) (delta), x.rj, rj.rk with slots a, b (letters) and o1, o2, ... (outputs in chain order)."""
import itertools, json
import sympy as sp

VT = {'AV': 'V', 'AsV': 'V', 'AW': 'W', 'AsW': 'W', 'c': None, 'cs': None}
OFFS = {'AV': 0, 'AsV': 0, 'AW': 0, 'AsW': 0}


def perfect_matchings(items):
    if not items:
        yield []; return
    a = items[0]
    for i in range(1, len(items)):
        rest = items[1:i] + items[i + 1:]
        for m in perfect_matchings(rest):
            yield [(a, items[i])] + m


def structures(vslots, wslots, n, D):
    """list of (label, function(comp dict, rsyms) -> sympy expr)"""
    out = []
    for k in range(0, len(vslots) + 1):
        if (len(vslots) - k) % 2:
            continue
        for mslots in itertools.combinations(vslots, k):        # slots contracted with momenta
            rest = [s for s in vslots if s not in mslots]
            if (D - k) < 0 or (D - k) % 2:
                continue
            ndot = (D - k) // 2
            for mv in perfect_matchings(rest):
                for js in itertools.product(range(n), repeat=k):
                    for dots in itertools.combinations_with_replacement(itertools.combinations_with_replacement(range(n), 2), ndot):
                        for mw in perfect_matchings(wslots):
                            lab = ' '.join(['d(%s,%s)' % p for p in mv] + ['%s.r%d' % (s, j + 1) for s, j in zip(mslots, js)] +
                                           ['r%d.r%d' % (i + 1, j + 1) for i, j in dots] + ['dW(%s,%s)' % p for p in mw]) or '1'
                            out.append((lab, (mv, list(zip(mslots, js)), dots, mw)))
    # remove duplicates
    seen = {}
    for lab, s in out:
        seen.setdefault(lab, s)
    return list(seen.items())


def evaluate(s, comp, r):
    mv, mj, dots, mw = s
    v = sp.Integer(1)
    for x, y in mv:
        v *= 1 if comp[x] == comp[y] else 0
    for x, j in mj:
        v *= r[j][comp[x]]
    for i, j in dots:
        v *= sum(r[i][m] * r[j][m] for m in range(4))
    for x, y in mw:
        v *= 1 if comp[x] == comp[y] else 0
    return v


def fit_word(ta, tb, word, entries):
    """entries: {'ia,ib,io..': poly string}; returns (formula string, ok)"""
    types = [ta, tb] + list(word)
    if any(t not in VT for t in types):
        return None, None
    names = ['a', 'b'] + ['o%d' % (j + 1) for j in range(len(word))]
    vslots = [nm for nm, t in zip(names, types) if VT[t] == 'V']
    wslots = [nm for nm, t in zip(names, types) if VT[t] == 'W']
    n = len(word)
    r = [[sp.Symbol('r%d_%d' % (j + 1, m)) for m in range(4)] for j in range(n)]
    rs = [x for row in r for x in row]
    ent = {}
    D = None
    for key, ps in entries.items():
        idx = [int(x) for x in key.split(',') if x != '']
        comp = {}
        for nm, t, i in zip(names, types, idx + [0] * (len(types) - len(idx))):
            comp[nm] = i
        e = sp.sympify(ps.replace('^', '**'), locals={str(x): x for x in rs} | {'i': sp.I})
        ent[tuple(idx)] = e
        if e != 0:
            deg = sp.Poly(e, *rs).total_degree() if rs else 0
            D = deg if D is None else D
    if D is None:
        return '0', True
    structs = structures(vslots, wslots, n, D)
    dims = [4 if VT[t] == 'V' else (6 if VT[t] == 'W' else 1) for t in types]
    unknowns = sp.symbols('x0:%d' % len(structs))
    eqs = []
    for idx in itertools.product(*[range(d) for d in dims]):
        comp = dict(zip(names, idx))
        lhs = ent.get(tuple(idx), sp.Integer(0))
        rhs = sum(u * evaluate(s, comp, r) for u, (lab, s) in zip(unknowns, structs))
        diff = sp.expand(rhs - lhs)
        if diff == 0:
            continue
        P = sp.Poly(diff, *rs) if rs else None
        eqs += list(P.coeffs()) if P is not None else [diff]
    if not structs:
        return None, False
    sol = sp.solve(eqs, unknowns, dict=True)
    if not sol:
        return None, False
    sol = sol[0]
    terms = []
    for u, (lab, s) in zip(unknowns, structs):
        c = sol.get(u, sp.Integer(0))
        if c.free_symbols:
            c = c.subs({x: 0 for x in c.free_symbols})       # remaining freedom: linear relations among structures
        if c != 0:
            terms.append('(%s)*%s' % (sp.nsimplify(c), lab))
    return ' + '.join(terms) if terms else '0', True


def run(raw):
    out = {}
    for pname, ent in raw.items():
        ta, tb = ent['letters']
        for legset, words in ent['leg_sets'].items():
            for wkey, d in words.items():
                word = [] if wkey.startswith('(none)') else wkey.split('|')[0].split()
                for kind, entries in d.items():
                    f, ok = fit_word(ta, tb, word, entries)
                    if f is None and ok is None:
                        continue
                    out.setdefault(pname, {})['%s | %s' % (wkey, kind)] = {'invariant_form': f, 'fit_verified': ok}
    return out


if __name__ == '__main__':
    import os
    raw = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'out.nosync', 'kernels_W0_raw.json')))
    res = run({k: v for k, v in raw.items() if k in ('AV|AV', 'AW|AW', 'AV|c')})
    print(json.dumps(res, indent=1))
