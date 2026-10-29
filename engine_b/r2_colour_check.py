"""Check colour(word) = tr((-,-)^{-1} C_1 ... (-,-)^{-1} C_n) (C from the vertices.py colour keys) against the abstract
reduction r2_engine.reduce_traces(colour_of_word(...)) on su(3) in a non-orthonormal basis with rescaled pairing."""
import itertools, random
from fractions import Fraction as Fr
import sympy as sp
from realization import su3_data
import r2_engine as R

L = su3_data(); n = L['n']; kap = L['kap']; kinv = L['kinv']; f = L['f']
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
def adm(c):   # matrix of ad_c in basis t_a
    return sp.Matrix(n, n, lambda a, b: br(c, E[b])[a])
random.seed(2)
cs = [[sp.Integer(random.randint(-2, 2)) for _ in range(n)] for _ in range(4)]
cadj = L['cadj']
def basisval(k, legs):
    c = [cs[i] for i in legs]
    if k == 'cadj k': return cadj * pr(c[0], c[1])
    if k == 'cadj f': return cadj * pr(br(c[0], c[1]), c[2])
    if k == 'cadj K12': return cadj * pr(br(c[0], c[1]), br(c[2], c[3]))
    if k == 'cadj K14': return cadj * pr(br(c[0], c[3]), br(c[1], c[2]))
    if k == 'str4':
        A = [adm(x) for x in c]
        return sum(((A[s[0]] * A[s[1]] * A[s[2]] * A[s[3]]).trace() for s in itertools.permutations(range(4)))) / 24
ok = True
words = [[('cubic', (0,)), ('cubic', (1,)), ('cubic', (2,))], [('cubic', (2,)), ('cubic', (0,)), ('cubic', (1,))],
         [('cubic', (0,)), ('cubic', (1,)), ('cubic', (2,)), ('cubic', (3,))], [('cubic', (1,)), ('cubic', (3,)), ('cubic', (0,)), ('cubic', (2,))],
         [('quartic', (0, 1)), ('cubic', (2,))], [('cubic', (2,)), ('quartic', (0, 1))]]
for word in words:
    legs = sorted(i for _, ids in word for i in ids)
    keysets = [['f'] if k == 'cubic' else ['K12', 'K13', 'K14'] for k, _ in word]
    for keys in itertools.product(*keysets):
        mats = [kinv * Cmat(k, [cs[i] for i in ids], key) for (k, ids), key in zip(word, keys)]
        M = mats[0]
        for X in mats[1:]:
            M = M * X
        explicit = sp.nsimplify(sp.simplify(M.trace()))
        red = R.reduce_traces(R.colour_of_word([R.letter_M(k, ids, key) for (k, ids), key in zip(word, keys)]), legs)
        abstract = sp.nsimplify(sp.simplify(sum(sp.Rational(v.numerator, v.denominator) * basisval(k, legs) for k, v in red.items())))
        good = sp.simplify(explicit - abstract) == 0
        ok &= good
        print(word, keys, 'explicit', explicit, 'abstract', abstract, good)
# four-letter word with two quartic letters
for keys in itertools.product(['K12', 'K13', 'K14'], repeat=2):
    word = [('quartic', (0, 2)), ('quartic', (1, 3))]
    mats = [kinv * Cmat('quartic', [cs[i] for i in ids], key) for (k, ids), key in zip(word, keys)]
    explicit = sp.simplify((mats[0] * mats[1]).trace())
    red = R.reduce_traces(R.colour_of_word([R.letter_M(k, ids, key) for (k, ids), key in zip(word, keys)]), [0, 1, 2, 3])
    abstract = sum(sp.Rational(v.numerator, v.denominator) * basisval(k, [0, 1, 2, 3]) for k, v in red.items())
    good = sp.simplify(explicit - abstract) == 0; ok &= good
    print(word, keys, good)
# Jacobi: K13 = K12 + K14
K = lambda a, b, c, d: pr(br(cs[a], cs[b]), br(cs[c], cs[d]))
print('Jacobi K13 = K12 + K14:', sp.simplify(K(0, 2, 1, 3) - K(0, 1, 2, 3) - K(0, 3, 1, 2)) == 0)
print('ALL OK' if ok else 'MISMATCH')
