"""
ops.py -- linear operators on sections induced by quadratic functionals, in position space.
For a quadratic functional F the operator X_F is defined by {F, l} = l o X_F on every coordinate
functional l (for F = S_0 this is Q, Proposition c1-Q); on coordinates
   (X_F phi)^a = -kappa_{a b} (-1)^{|a||F|} E_b F ,  b the partner of a.
Sections: dicts {field index: sympy expression in x}.
"""
import sympy as sp
import bvjet

X = bvjet.X


def eval_linear(M, L, e):
    """evaluate a superpolynomial linear in the jets on a section e."""
    out = 0
    for (ev, od), c in L.items():
        if ev and not od:
            assert len(ev) == 1 and ev[0][1] == 1
            (a, al), _ = ev[0]
        elif od and not ev:
            assert len(od) == 1
            a, al = od[0]
        else:
            raise ValueError('not linear')
        f = e.get(a, 0)
        if f == 0:
            continue
        for m in range(4):
            for _ in range(al[m]):
                f = sp.diff(f, X[m])
        out += c * f
    return sp.expand(out)


class Op:
    def __init__(self, M, rows, parity):
        self.M = M
        self.rows = rows          # dict a -> linear superpolynomial (or callable)
        self.parity = parity

    def __call__(self, e):
        out = {}
        for a, L in self.rows.items():
            v = eval_linear(self.M, L, e) if isinstance(L, dict) else L(e)
            if v != 0:
                out[a] = v
        return out


def partner_of(M, a):
    for (x, y) in M.kappa:
        if x == a:
            return y
    return None


def op_from_functional(M, F):
    pF = M.parity_of(F)
    rows = {}
    for a in range(len(M.fields)):
        b = partner_of(M, a)
        if b is None:
            continue
        Eb = M.simplify(M.EL(F, b))
        if Eb:
            rows[a] = M.scale(Eb, -M.kappa[(a, b)] * (-1) ** (M.par(a) * pF))
    return Op(M, rows, (pF + 1) % 2)


def op_from_momentum(M, Gmat, p, parity):
    """operator with matrix G(p) in momentum space, d -> i p, i.e. p_m -> -i d_m."""
    rows = {}
    for a in range(Gmat.shape[0]):
        L = {}
        for b in range(Gmat.shape[1]):
            q = sp.expand(Gmat[a, b])
            if q == 0:
                continue
            for mon, cf in sp.Poly(q, *p).terms():
                L = M.plus(L, M.jet(M.fields[b].name, mon, cf * (-sp.I) ** sum(mon)))
        if L:
            rows[a] = L
    return Op(M, rows, parity)


def compose(A, B):
    return lambda e: A(B(e))


def sec_sub(e1, e2, s=1):
    out = dict(e1)
    for k, v in e2.items():
        out[k] = sp.expand(out.get(k, 0) - s * v)
        if out[k] == 0:
            del out[k]
    return out


def sec_scale(e, c):
    return {k: sp.expand(c * v) for k, v in e.items() if sp.expand(c * v) != 0}


def gcomm(A, B, e):
    """graded commutator [A,B] e = A(B e) - (-1)^{|A||B|} B(A e)."""
    s = (-1) ** (A.parity * B.parity)
    return sec_sub(A(B(e)), B(A(e)), s)


def generic_section(M, slot_indices, tag='f'):
    fs = {}
    for a in slot_indices:
        fs[a] = sp.Function('%s%d' % (tag, a))(*X)
    return fs


def is_zero_section(e):
    return all(sp.simplify(v) == 0 for v in e.values())
