"""
freecx.py -- the free BV complex of the reduced theory in momentum space, one colour component
(the free theory is diagonal in g; the pairing (-,-) of g is set to 1 here and restored as (-,-)^{-1}
in the propagator components).

A section of momentum p is u e^{i(p,x)}; d -> i p, d^* on one-forms -> -i (p, -).
Operators are matrices acting on the column of fibre components.

Slot layout ('min' = minimal complex of Chapter 1; 'red' = minimal + reduced non-minimal pair (cbar, cbar^*);
'b' = minimal + (cbar, b, cbar^*, b^*)):
  c 0 | A|_V 1-4 | A|_Wperp 5-10 | lambda 11-26 | A^*|_V 27-30 | A^*|_Wperp 31-36 | lambda^* 37-52 | c^* 53
  | cbar 54 | cbar^* 55 | b 56 | b^* 57
"""
import sympy as sp
import realization as R

p = sp.symbols('p1:5', real=True)
P2 = sum(t * t for t in p)
I_ = sp.I

SLOTS = [  # name, start, dim, total parity, Z-degree, Weyl weight
    ('c', 0, 1, 1, -1, 0),
    ('AV', 1, 4, 0, 0, 1),
    ('AW', 5, 6, 0, 0, 1),
    ('lam', 11, 16, 1, 0, sp.Rational(3, 2)),
    ('AsV', 27, 4, 1, 1, 3),
    ('AsW', 31, 6, 1, 1, 3),
    ('lams', 37, 16, 0, 1, sp.Rational(5, 2)),
    ('cs', 53, 1, 0, 2, 4),
    ('cb', 54, 1, 1, 1, 2),
    ('cbs', 55, 1, 0, 0, 2),
    ('b', 56, 1, 0, 0, 2),
    ('bs', 57, 1, 1, 1, 2),
]
SL = {s[0]: s for s in SLOTS}
DIMS = {'min': 54, 'red': 56, 'b': 58}
FIELDS = {'c', 'AV', 'AW', 'lam', 'cb', 'b'}
PARTNER = {'c': 'cs', 'cs': 'c', 'AV': 'AsV', 'AsV': 'AV', 'AW': 'AsW', 'AsW': 'AW',
           'lam': 'lams', 'lams': 'lam', 'cb': 'cbs', 'cbs': 'cb', 'b': 'bs', 'bs': 'b'}


def rng(name):
    s = SL[name]
    return range(s[1], s[1] + s[2])


def slots_of(config):
    n = DIMS[config]
    return [s for s in SLOTS if s[1] < n]


def parity_vector(config):
    out = []
    for s in slots_of(config):
        out += [s[3]] * s[2]
    return out


def Pi(config):
    return sp.diag(*[(-1) ** a for a in parity_vector(config)])


def Omega(config):
    """Fibre matrix of the odd pairing: <u e^{ipx}, w e^{-ipx}> = u^T Omega w (volume factored out).
    Even members A, lambda^*, c^*, cbar^*, b; their partners A^*, lambda, c, cbar, b^*."""
    n = DIMS[config]
    O = sp.zeros(n, n)
    C = R.Cpair()
    for a, b in [('AV', 'AsV'), ('AW', 'AsW')]:
        for k, l in zip(rng(a), rng(b)):
            O[k, l] = 1
            O[l, k] = -1
    # (lambda^*, lambda): <e0 in lambda^* slot, e1 in lambda slot> = (e1, e0)_{Sigma x Sigma^*} = e1^T C e0
    for ia, a in enumerate(rng('lam')):
        for ib, b in enumerate(rng('lams')):
            if C[ia, ib] != 0:
                O[b, a] = C[ia, ib]
                O[a, b] = -C[ia, ib]
    pairs = [('cs', 'c')]
    if config in ('red', 'b'):
        pairs.append(('cbs', 'cb'))
    if config == 'b':
        pairs.append(('b', 'bs'))
    for even, odd in pairs:
        k, l = rng(even)[0], rng(odd)[0]
        O[k, l] = 1
        O[l, k] = -1
    return O


def rho_p(q):
    Rm, Rs = R.rho_vec(list(q))
    return Rm, Rs


def Q_min(config='min', q=p):
    """The operator Q of Proposition c1-Q, extended by zero to the non-minimal slots."""
    n = DIMS[config]
    Q = sp.zeros(n, n)
    q2 = sum(t * t for t in q)
    c0 = rng('c')[0]
    for m, k in enumerate(rng('AV')):
        Q[k, c0] = -I_ * q[m]                          # c -> -dc
    for m, k in enumerate(rng('AsV')):                 # A|_V -> -1/2 d^*d A
        for n_, l in enumerate(rng('AV')):
            Q[k, l] = -sp.Rational(1, 2) * ((q2 if m == n_ else 0) - q[m] * q[n_])
    for k, l in zip(rng('AsW'), rng('AW')):
        Q[k, l] = -sp.Rational(1, 2) * q2
    cs0 = rng('cs')[0]
    for m, l in enumerate(rng('AsV')):
        Q[cs0, l] = -I_ * q[m]                         # A^* -> d^*(A^*|_V)
    Rm, Rs = rho_p(q)
    for a, k in enumerate(rng('lams')):                # lambda -> dslash lambda = i rho(p) lambda
        for b, l in enumerate(rng('lam')):
            if Rm[a, b] != 0:
                Q[k, l] = I_ * Rm[a, b]
    return Q


def QGF_old(q=p):
    """The Feynman-gauge operator of the former Chapter 3 (minimal complex)."""
    n = DIMS['min']
    G = sp.zeros(n, n)
    for k, l in zip(rng('AV'), rng('AsV')):
        G[k, l] = -2
    for k, l in zip(rng('AW'), rng('AsW')):
        G[k, l] = -2
    c0 = rng('c')[0]
    for m, l in enumerate(rng('AV')):
        G[c0, l] = I_ * q[m]                           # A -> -d^*(A|_V) = i iota_p A
    cs0 = rng('cs')[0]
    for m, k in enumerate(rng('AsV')):
        G[k, cs0] = I_ * q[m]                          # c^* -> d c^*
    Rm, Rs = rho_p(q)
    for a, k in enumerate(rng('lam')):                 # lambda^* -> - dslash lambda^* = -i rho(p) lambda^*
        for b, l in enumerate(rng('lams')):
            if Rs[a, b] != 0:
                G[k, l] = -I_ * Rs[a, b]
    return G


def neg(q):
    return tuple(-t for t in q)


def adj_check_Q(Q, Qm, config):
    """<Qe0,e1> + (-1)^{|e0|}<e0,Qe1> = 0  <=>  Q(p)^T Om + Pi Om Q(-p) = 0."""
    Om = Omega(config)
    X = Q.T * Om + Pi(config) * Om * Qm
    return X.applyfunc(sp.expand)


def adj_check_GF(G, Gm, config):
    """<Ge0,e1> = (-1)^{|e0|}<e0,Ge1>  <=>  G(p)^T Om - Pi Om G(-p) = 0."""
    Om = Omega(config)
    X = G.T * Om - Pi(config) * Om * Gm
    return X.applyfunc(sp.expand)


def is_zero(M):
    return all(sp.simplify(x) == 0 for x in M)


def embed(M, n_from, n_to):
    out = sp.zeros(n_to, n_to)
    out[0:n_from, 0:n_from] = M
    return out


def K0hat(config):
    """Components of the kernel of the identity K_0 (Definition c1-bracket): -1 times the inverse fibre pairing
    on A (x) A^*, A^* (x) A, lambda (x) lambda^*, lambda^* (x) lambda, c (x) c^*, c^* (x) c (and cbar, b pairs)."""
    n = DIMS[config]
    K = sp.zeros(n, n)
    C = R.Cpair()
    for a, b in [('AV', 'AsV'), ('AW', 'AsW')]:
        for k, l in zip(rng(a), rng(b)):
            K[k, l] = -1
            K[l, k] = -1
    for ia, a in enumerate(rng('lam')):
        b = rng('lams')[ia]
        K[a, b] = -C[ia, ia]
        K[b, a] = -C[ia, ia]
    pairs = [('c', 'cs')]
    if config in ('red', 'b'):
        pairs.append(('cb', 'cbs'))
    if config == 'b':
        pairs.append(('b', 'bs'))
    for x, y in pairs:
        k, l = rng(x)[0], rng(y)[0]
        K[k, l] = -1
        K[l, k] = -1
    return K
