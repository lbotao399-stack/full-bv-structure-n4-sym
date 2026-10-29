"""
tower.py -- component form of  d_CE S_g + 1/2 {S_g, S_g} = 0  [ESW Def 1.26, 1.27(c)] on actual functionals.

Convention (locked in conventions.py and against ESW Lemma 3.6 by machine, see verify.json 'conventions'):
  * {A,B}: the [C07 Lemma 4.3.2] bracket of classical_spec Sec 3 (= ESW's bracket: reproduces ESW Lemma 3.6 item 3);
  * cochain bracket: ESW Def 1.26 formula with |g| := degree of g in C^bullet(g, O_loc) (= 0 for S_g);
  * d_CE: ESW Def 1.26 formula with |f| := degree of f in C^bullet(g, O_loc) (= 0 for S_g) and ONE extra overall
    sign, i.e. (-1)^{... + i + j + |f| + 1}; this is the sign for which ESW Lemma 3.6 item 2,
    (d_CE S^(1))(Q1,Q2) = S^(1)(Gamma(Q1,Q2)), holds.  (With the printed sign and |f| = 0 one gets the CE
    differential of the opposite algebra; the printed formula with |f| := map degree is not consistent at
    arity 3 -- conventions.py.)
"""
import itertools
from fractions import Fraction as Fr
import engine as E
import algebra as Al


def sgn_perm(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def koszul(pars, perm):
    s = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                s += pars[perm[i]] * pars[perm[j]]
    return -1 if s % 2 else 1


class Tower:
    def __init__(self, V, S, cochains):
        """cochains: dict k -> function(list of algebra elements) -> polynomial density (S^(k) value)."""
        self.V = V
        self.S = S
        self.ch = cochains
        self.cache = {}

    def val(self, k, args):
        if k == 0:
            return self.S
        f = self.ch.get(k)
        if f is None:
            return {}
        return f(args)

    def dCE(self, k, args):
        """(d_CE S^(k))(x_1..x_{k+1})."""
        n = len(args)
        pars = [Al.parity(a) for a in args]
        out = {}
        for i in range(n):
            for j in range(i + 1, n):
                e = pars[i] * sum(pars[p] for p in range(i))
                e += pars[j] * sum(pars[p] for p in range(j) if p != i)
                e += (i + 1) + (j + 1) + 0 + 1          # |f| = 0 (degree in C^bullet), plus the overall sign
                br = Al.bracket(args[i], args[j])
                if br.is_zero():
                    continue
                rest = [args[p] for p in range(n) if p not in (i, j)]
                v = self.val(k, [br] + rest)
                if v:
                    E.padd_into(out, v, (Fr(-1 if e % 2 else 1), Fr(0)))
        return out

    def cbr(self, k, l, args):
        """{S^(k), S^(l)}(x_1..x_{k+l})  with |g| = 0."""
        n = len(args)
        pars = [Al.parity(a) for a in args]
        out = {}
        for left in itertools.combinations(range(n), k):
            right = [p for p in range(n) if p not in left]
            perm = list(left) + right
            e1 = 0 * k + sum(pars[p] for p in left) * (l + 0)
            sign = sgn_perm(perm) * koszul(pars, perm) * (-1 if e1 % 2 else 1)
            # evaluate the higher-arity factor first (cheap, often identically zero)
            if l >= k:
                B = self.val(l, [args[p] for p in right])
                if not B:
                    continue
                A = self.val(k, [args[p] for p in left])
            else:
                A = self.val(k, [args[p] for p in left])
                if not A:
                    continue
                B = self.val(l, [args[p] for p in right])
            if not A or not B:
                continue
            E.padd_into(out, E.bracket(self.V, A, B), (Fr(sign), Fr(0)))
        return out

    def MC(self, args, split=False):
        """component of d_CE S_g + 1/2 {S_g,S_g} at (x_1..x_n)."""
        n = len(args)
        parts = {'dCE': self.dCE(n - 1, args) if n >= 2 else {}}
        for k in range(0, n + 1):
            parts['br%d%d' % (k, n - k)] = E.pscale(self.cbr(k, n - k, args), Fr(1, 2))
        tot = {}
        for p in parts.values():
            E.padd_into(tot, p)
        if split:
            return tot, parts
        return tot

    def check(self, args):
        tot = self.MC(args)
        ok, bad = E.total_derivative_test(self.V, tot)
        return ok, len(tot)
