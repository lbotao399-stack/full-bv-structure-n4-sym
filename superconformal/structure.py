"""
structure.py -- identification of the realization in algebra.py with psl(4|4) in the [SW (17)] block form.

Subspaces (named by [SW (8),(17)]):
  P       translations                                 block (adot, alpha)
  K       special conformal |x|^2 b - 2 (b.x) x        block (alpha, adot)
  M, Mbar rotations acting only on S_+ / only on S_-   blocks (alpha,alpha) / (adot,adot)
  Delta   2E                                           diagonal (alpha,alpha)+(adot,adot)
  R       so(W_R^perp)                                 block (i,i)
  Q       S_+ (x) W   in Sigma                         block (i, alpha)
  Qbar    S_- (x) W^* in Sigma                         block (adot, i)
  Sbar    S_- (x) W   in Sigma^*                       block (i, adot)
  S       S_+ (x) W^* in Sigma^*                       block (alpha, i)
with S_+ (x) W := {+1 eigenspace of rho(q(dvol)) = rho(e_0)rho(e_1)rho(e_2)rho(e_3) on Sigma} (dvol = e^0^e^1^e^2^e^3)
and on Sigma^* the image of S_+ (x) W under rho(V) is S_- (x) W, etc.
Block rule of supermatrices [SW (17),(19)]: [block(a,b), block(c,d)] in block(a,d) [if b=c] + block(c,b) [if d=a].
"""
import itertools
import sympy as sp
import algebra as Al
import model as Mo

X = Al.X


def chirality_ops():
    S, Sb = Al.SIGM, Al.SIGBM
    Om = Sb[0] * S[1] * Sb[2] * S[3]          # on Sigma: rho(e0)rho(e1)rho(e2)rho(e3)
    Oms = S[0] * Sb[1] * S[2] * Sb[3]         # on Sigma^*
    return Om, Oms


def eigbasis(Mx, ev):
    ns = (Mx - ev * sp.eye(16)).nullspace()
    return [list(v) for v in ns]


def basis():
    """dict name -> list of algebra elements spanning the subspace."""
    Om, Oms = chirality_ops()
    B = {}
    B['P'] = [Al.P([1 if m == mu else 0 for m in range(4)]) for mu in range(4)]
    B['K'] = [Al.Kb([1 if m == mu else 0 for m in range(4)]) for mu in range(4)]
    rots = []
    for a in range(4):
        for b in range(a + 1, 4):
            om = [[0] * 4 for _ in range(4)]
            om[a][b], om[b][a] = 1, -1
            rots.append(Al.Mrot(om))
    # self-dual / anti-self-dual split by action on S_+ / S_- (subspaces of Sigma)
    Qp = eigbasis(Om, 1)
    Qm = eigbasis(Om, -1)
    def acts_on(rot, sub):
        for s in sub:
            r = Al.eo_bracket(rot, Al.Qgen(s))
            if not r.is_zero():
                return True
        return False
    # find combinations: solve for rotations killing S_- (M) and killing S_+ (Mbar)
    cs = sp.symbols('m0:6')
    gen = rots[0].scale(cs[0])
    for k in range(1, 6):
        gen = gen + rots[k].scale(cs[k])
    def kill_space(sub):
        eqs = []
        for s in sub:
            r = Al.eo_bracket(gen, Al.Qgen(s))
            eqs += [e for e in list(r.Q) + list(r.eta) if e != 0]
        sol = sp.linsolve(eqs, cs)
        (tup,) = list(sol)
        free = sorted(set().union(*[sp.sympify(t).free_symbols for t in tup]), key=str)
        out = []
        for f in free:
            vals = [sp.sympify(t).subs({g: (1 if g == f else 0) for g in free}) for t in tup]
            el = rots[0].scale(vals[0])
            for k in range(1, 6):
                el = el + rots[k].scale(vals[k])
            out.append(el)
        return out
    B['M'] = kill_space(Qm)       # acts only on S_+
    B['Mbar'] = kill_space(Qp)    # acts only on S_-
    B['Delta'] = [Al.Dil()]
    Rl = []
    for i in Mo.WP:
        for j in Mo.WP:
            if i < j:
                Rl.append(Al.Rsym({(i, j): 1, (j, i): -1}))
    B['R'] = Rl
    B['Q'] = [Al.Qgen(s) for s in Qp]
    B['Qbar'] = [Al.Qgen(s) for s in Qm]
    B['Sbar'] = [Al.Sgen(s) for s in eigbasis(Oms, -1)]
    B['S'] = [Al.Sgen(s) for s in eigbasis(Oms, 1)]
    return B


BLOCK = {'P': [(2, 1)], 'K': [(1, 2)], 'M': [(1, 1)], 'Mbar': [(2, 2)], 'Delta': [(1, 1), (2, 2)], 'R': [(3, 3)],
         'Q': [(3, 1)], 'Qbar': [(2, 3)], 'Sbar': [(3, 2)], 'S': [(1, 3)]}
DIAG = {(1, 1): ['M', 'Delta'], (2, 2): ['Mbar', 'Delta'], (3, 3): ['R', 'Delta']}
OFFD = {(2, 1): 'P', (1, 2): 'K', (3, 1): 'Q', (2, 3): 'Qbar', (3, 2): 'Sbar', (1, 3): 'S'}


def allowed(n1, n2):
    out = set()
    for (a, b) in BLOCK[n1]:
        for (c, d) in BLOCK[n2]:
            tgt = []
            if b == c:
                tgt.append((a, d))
            if d == a:
                tgt.append((c, b))
            for t in tgt:
                if t in DIAG:
                    out.update(DIAG[t])
                else:
                    out.add(OFFD[t])
    return out


class Coords:
    """coordinates in the 62-element basis assembled from basis()."""
    def __init__(self, B):
        self.names = []
        self.even, self.odd = [], []
        for n in ['P', 'K', 'M', 'Mbar', 'Delta', 'R']:
            for e in B[n]:
                self.even.append(e)
                self.names.append(n)
        for n in ['Q', 'Qbar', 'Sbar', 'S']:
            for e in B[n]:
                self.odd.append(e)
                self.names.append(n)
        self.ne = len(self.even)
        mons = set()
        for e in self.even:
            for M in range(10):
                p = sp.Poly(e.c[M], *(list(X) + [Al.Y[w] for w in Mo.WP]))
                for m, _ in p.terms():
                    mons.add((M, m))
        self.mons = sorted(mons)
        self.Be = sp.Matrix([[self._feat(e)[k] for e in self.even] for k in range(len(self.mons))])
        self.Le = (self.Be.T * self.Be).inv() * self.Be.T
        self.Bo = sp.Matrix([[*(list(o.Q) + list(o.eta))] for o in self.odd]).T
        self.Lo = self.Bo.inv()

    def _feat(self, e):
        f = [0] * len(self.mons)
        idx = {m: k for k, m in enumerate(self.mons)}
        for M in range(10):
            if e.c[M] == 0:
                continue
            p = sp.Poly(e.c[M], *(list(X) + [Al.Y[w] for w in Mo.WP]))
            for m, cf in p.terms():
                if (M, m) not in idx:
                    raise ValueError('not in even span')
                f[idx[(M, m)]] = cf
        return f

    def coords(self, el):
        if isinstance(el, Al.Even):
            v = self.Le * sp.Matrix(self._feat(el))
            # verify exact reconstruction
            assert sp.simplify(self.Be * v - sp.Matrix(self._feat(el))) == sp.zeros(len(self.mons), 1)
            return list(v) + [0] * len(self.odd)
        v = self.Lo * sp.Matrix(list(el.Q) + list(el.eta))
        return [0] * self.ne + list(v)

    def support(self, el):
        c = self.coords(el)
        return set(self.names[k] for k, x in enumerate(c) if sp.simplify(x) != 0)

    def element(self, k):
        return self.even[k] if k < self.ne else self.odd[k - self.ne]

    def parity(self, k):
        return 0 if k < self.ne else 1


def block_rule_check(B, C):
    """for every pair of named subspaces and every pair of basis elements, the bracket lies in the allowed blocks."""
    bad = []
    names = list(B.keys())
    for n1, n2 in itertools.combinations_with_replacement(names, 2):
        al = allowed(n1, n2) | allowed(n2, n1)
        for x in B[n1]:
            for y in B[n2]:
                br = Al.bracket(x, y)
                if br.is_zero():
                    continue
                sup = C.support(br)
                if not sup <= al:
                    bad.append((n1, n2, sorted(sup), sorted(al)))
                    break
    return bad


def ad_matrix(C, x):
    cols = []
    for k in range(C.ne + len(C.odd)):
        cols.append(C.coords(Al.bracket(x, C.element(k))))
    return sp.Matrix(cols).T


def killing(C, x, y):
    Ax, Ay = ad_matrix(C, x), ad_matrix(C, y)
    P = Ax * Ay
    return sp.simplify(sum(((-1) ** C.parity(k)) * P[k, k] for k in range(P.shape[0])))
