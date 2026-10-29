"""Colour: the abstract reduction used by the I-a engine (r2_engine.colour_of_word + reduce_traces:
tr_g(ad ad ad) = (1/2) tr_g(ad_X ad_Y)/(X,Y) ([c_0,c_1],c_2) with the word-order sign) against explicit traces in su(2) and
su(3) (non-orthonormal bases, rescaled pairing), for every word type of the three-leg I-a graphs: three cubic letters in all
orders, and one quartic letter (keys K12, K13, K14) with one cubic letter in both orders.  Also the tree-level contractions
of fa.py ('f' x 'k' and 'k' x 'f' through (-,-)^{-1})."""
import itertools, random, json
import sympy as sp
import fa
import realization as Rz
import r2_engine as R
import fa

out = {}
for gname, L in (('su2', Rz.su2_data()), ('su3', Rz.su3_data())):
    n = L['n']; kap = L['kap']; kinv = L['kinv']; f = L['f']; cadj = L['cadj']
    def br(x, y):
        return [sum(x[a] * y[b] * f[a][b][c] for a in range(n) for b in range(n)) for c in range(n)]
    def pr(x, y):
        return sum(x[a] * kap[a, b] * y[b] for a in range(n) for b in range(n))
    E = [[1 if i == d else 0 for i in range(n)] for d in range(n)]
    def Cmat(kind, cols, key):
        if kind == 'cubic':
            c = cols[0]
            return sp.Matrix(n, n, lambda d, e: pr(br(c, E[d]), E[e]))
        c1, c2 = cols
        if key == 'K12':
            return sp.Matrix(n, n, lambda d, e: pr(br(c1, c2), br(E[d], E[e])))
        if key == 'K13':
            return sp.Matrix(n, n, lambda d, e: pr(br(c1, E[d]), br(c2, E[e])))
        if key == 'K14':
            return sp.Matrix(n, n, lambda d, e: pr(br(c1, E[e]), br(c2, E[d])))
    rng = random.Random(4)
    cs = [[sp.Integer(rng.randint(-2, 2)) for _ in range(n)] for _ in range(3)]
    fval = pr(br(cs[0], cs[1]), cs[2])
    ok = True; ratios = set()
    words = [[('cubic', (i,)) for i in perm] for perm in itertools.permutations(range(3))]
    for (i, j) in itertools.combinations(range(3), 2):
        k = [x for x in range(3) if x not in (i, j)][0]
        words += [[('quartic', (i, j)), ('cubic', (k,))], [('cubic', (k,)), ('quartic', (i, j))]]
    for word in words:
        keysets = [['f'] if kd == 'cubic' else ['K12', 'K13', 'K14'] for kd, _ in word]
        for keys in itertools.product(*keysets):
            mats = [kinv * Cmat(kd, [cs[x] for x in ids], key) for (kd, ids), key in zip(word, keys)]
            Mx = mats[0]
            for X in mats[1:]:
                Mx = Mx * X
            explicit = sp.nsimplify(sp.simplify(Mx.trace()))
            red = R.reduce_traces(R.colour_of_word([R.letter_M(kd, ids, key) for (kd, ids), key in zip(word, keys)]), [0, 1, 2])
            c = red.get('cadj f', 0)
            abstract = sp.Rational(c.numerator, c.denominator) * cadj * fval
            good = sp.simplify(explicit - abstract) == 0
            ok &= good
            if c != 0 and fval != 0:
                ratios.add(str(sp.nsimplify(explicit / (cadj * fval))))
    # tree contractions: ([c0,c1], t_a) kinv^{ab} (t_b, c2) = ([c0,c1],c2); (c0, t_a) kinv ([t_b, c1], c2) = ([c0,c1],c2)
    t1 = sum(pr(br(cs[0], cs[1]), E[a]) * kinv[a, b] * pr(E[b], cs[2]) for a in range(n) for b in range(n))
    t2 = sum(pr(cs[0], E[a]) * kinv[a, b] * pr(br(E[b], cs[1]), cs[2]) for a in range(n) for b in range(n))
    out[gname] = dict(dim=n, cadj=str(cadj), all_words_agree=bool(ok), multiples_of_cadj_f=sorted(ratios),
                      tree_f_k=bool(sp.simplify(t1 - fval) == 0), tree_k_f=bool(sp.simplify(t2 - fval) == 0))
    print(gname, out[gname], flush=True)
json.dump(out, open('colour_check.json', 'w'), indent=1)
