"""Engine B, Step 2: the four three-leg sectors of degree 1, random plane-wave configurations, the target basis G1..G5
and a complete basis of Spin(4) x Spin(6)-invariant multilinear forms (colour ([c_0,c_1],c_2)) for the residual test.
Leg order in every sector: leg 0 = c; legs 1, 2 as named.  Vertex functions: coefficient of the odd r's in leg order."""
import itertools, random, sys, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
EB = os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b'))
if EB not in sys.path:
    sys.path.insert(0, EB)
import vertices as V

I = sp.I
Z = sp.Integer(0)
PLUS = list(range(8)); MINUS = list(range(8, 16))

SECTORS = ['cAVAV', 'cApAp', 'clamlam', 'ccAsV']
SECTOR_LABEL = {'cAVAV': '(c, A|_V, A|_V)', 'cApAp': '(c, A|_Wperp, A|_Wperp)', 'clamlam': '(c, lambda, lambda)', 'ccAsV': '(c, c, A*|_V)'}
OMEGA = {'cAVAV': 2, 'cApAp': 2, 'clamlam': 1, 'ccAsV': 1}

def ri(rng):
    return sp.Integer(rng.choice([-3, -2, -1, 1, 2, 3]))

def config(sector, seed):
    rng = random.Random(1000 + 37 * seed)
    p0 = [ri(rng) for _ in range(4)]; p1 = [ri(rng) for _ in range(4)]
    p2 = [-x - y for x, y in zip(p0, p1)]
    c0 = dict(sp='c', u=[sp.Integer(1)], p=p0)
    if sector == 'cAVAV':
        l1 = dict(sp='A', u=[ri(rng) for _ in range(4)] + [Z] * 6, p=p1)
        l2 = dict(sp='A', u=[ri(rng) for _ in range(4)] + [Z] * 6, p=p2)
    elif sector == 'cApAp':
        l1 = dict(sp='A', u=[Z] * 4 + [ri(rng) for _ in range(6)], p=p1)
        l2 = dict(sp='A', u=[Z] * 4 + [ri(rng) for _ in range(6)], p=p2)
    elif sector == 'clamlam':
        l1 = dict(sp='lam', u=[ri(rng) for _ in range(16)], p=p1)
        l2 = dict(sp='lam', u=[ri(rng) for _ in range(16)], p=p2)
    elif sector == 'ccAsV':
        l1 = dict(sp='c', u=[sp.Integer(1)], p=p1)
        l2 = dict(sp='As', u=[ri(rng) for _ in range(4)] + [Z] * 6, p=p2)
    return [c0, l1, l2]

def dot(a, b, comps=range(4)):
    return sum(a[i] * b[i] for i in comps)

def eps4(a, b, c, d):
    tot = 0
    for perm in itertools.permutations(range(4)):
        s = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if perm[i] > perm[j]:
                    s = -s
        tot += s * a[perm[0]] * b[perm[1]] * c[perm[2]] * d[perm[3]]
    return tot

def lam_bil(u1, u2, mom, P1=None, P2=None):
    """(lambda_1, rho(mom) lambda_2) = lambda_1^T Mp R_mu lambda_2 mom_mu (code basis), optional chiral projections."""
    v1 = sp.Matrix([u1[i] if (P1 is None or i in P1) else 0 for i in range(16)])
    v2 = sp.Matrix([u2[i] if (P2 is None or i in P2) else 0 for i in range(16)])
    return sp.expand(sum((v1.T * V.MR[m] * v2)[0] * mom[m] for m in range(4)))

def structures(sector, legs, n):
    """complete list of invariant multilinear forms of momentum degree omega + 2n (colour ([c0,c1],c2) implicit)."""
    p0, p1, p2 = legs[0]['p'], legs[1]['p'], legs[2]['p']
    u1, u2 = legs[1]['u'], legs[2]['u']
    S = {}
    deg = OMEGA[sector] + 2 * n
    if sector == 'cAVAV':
        if deg == 0:
            S['(a1.a2)'] = dot(u1, u2)
        elif deg == 2:
            S['(a1.a2)(p1.p1)'] = dot(u1, u2) * dot(p1, p1); S['(a1.a2)(p2.p2)'] = dot(u1, u2) * dot(p2, p2)
            S['(a1.a2)(p1.p2)'] = dot(u1, u2) * dot(p1, p2)
            S['(a1.p1)(a2.p1)'] = dot(u1, p1) * dot(u2, p1); S['(a1.p1)(a2.p2)'] = dot(u1, p1) * dot(u2, p2)
            S['(a1.p2)(a2.p1)'] = dot(u1, p2) * dot(u2, p1); S['(a1.p2)(a2.p2)'] = dot(u1, p2) * dot(u2, p2)
            S['eps(a1,a2,p1,p2)'] = eps4(u1[:4], u2[:4], p1, p2)
    elif sector == 'cApAp':
        P = range(4, 10)
        if deg == 0:
            S['(b1.b2)'] = dot(u1, u2, P)
        elif deg == 2:
            S['(b1.b2)(p1.p1)'] = dot(u1, u2, P) * dot(p1, p1); S['(b1.b2)(p2.p2)'] = dot(u1, u2, P) * dot(p2, p2)
            S['(b1.b2)(p1.p2)'] = dot(u1, u2, P) * dot(p1, p2)
    elif sector == 'clamlam':
        if deg == 1:
            for nm, mom in (('p1', p1), ('p2', p2)):
                S[f'(lam1+, rho({nm}) lam2-)'] = lam_bil(u1, u2, mom, PLUS, MINUS)
                S[f'(lam1-, rho({nm}) lam2+)'] = lam_bil(u1, u2, mom, MINUS, PLUS)
    elif sector == 'ccAsV':
        if deg == 1:
            S['(a*.p0)'] = dot(u2, p0); S['(a*.p1)'] = dot(u2, p1)
    return S

# ------------------------------------------------------------------ target basis G1..G5 as local functionals
def _G1(s):   # ([c, A|_V], d^*d A|_V):  (d^*d A)_mu -> |p|^2 a_mu - p_mu (p.a)
    c, a1, a2 = s
    p = a2['p']
    return c['u'][0] * sum(a1['u'][m] * (dot(p, p) * a2['u'][m] - p[m] * dot(p, a2['u'])) for m in range(4))
def _G2(s):   # ([c, A|_V], d d^* A|_V):  (d d^* A)_mu -> p_mu (p.a)
    c, a1, a2 = s
    p = a2['p']
    return c['u'][0] * sum(a1['u'][m] * p[m] * dot(p, a2['u']) for m in range(4))
def _G3(s):   # ([c, A|_Wperp], d^*d A|_Wperp)
    c, a1, a2 = s
    p = a2['p']
    return c['u'][0] * dot(p, p) * sum(a1['u'][n] * a2['u'][n] for n in range(4, 10))
def _G4(s):   # ([c, dc], A^*), pairing of the action term (dc, A^*) = sum_mu (d_mu c, A^*_mu)
    c1, c2, As = s
    return c1['u'][0] * c2['u'][0] * sum(I * c2['p'][m] * As['u'][m] for m in range(4))
def _G5(s):   # ([c, lambda], dslash lambda), pairing of 1/2 (lambda, dslash lambda) = 1/2 sum_mu lam^T MR_mu d_mu lam
    c, l1, l2 = s
    u1 = sp.Matrix(l1['u']); u2 = sp.Matrix(l2['u'])
    return c['u'][0] * sum((u1.T * V.MR[m] * u2)[0] * I * l2['p'][m] for m in range(4))
GTERMS = {
    'G1': (sp.Integer(1), ('c', 'A', 'A'), _G1, ('f', (0, 1, 2))),
    'G2': (sp.Integer(1), ('c', 'A', 'A'), _G2, ('f', (0, 1, 2))),
    'G3': (sp.Integer(1), ('c', 'A', 'A'), _G3, ('f', (0, 1, 2))),
    'G4': (sp.Integer(1), ('c', 'c', 'As'), _G4, ('f', (0, 1, 2))),
    'G5': (sp.Integer(1), ('c', 'lam', 'lam'), _G5, ('f', (0, 1, 2))),
}
SECTOR_G = {'cAVAV': ['G1', 'G2'], 'cApAp': ['G3'], 'clamlam': ['G5'], 'ccAsV': ['G4']}

def G_value(name, legs):
    V.TERMS['_G'] = GTERMS[name]
    try:
        v = V.vertex('_G', legs).get('f', 0)
    finally:
        del V.TERMS['_G']
    return sp.expand(v)

def solve_exact(rows, names):
    """rows: list of (target value, {name: value}).  Picks rows greedily until the structure matrix has full column rank,
    solves exactly, and returns (solution dict, residual on every row).  Raises if rank deficient."""
    sel = []
    for i, r in enumerate(rows):
        trial = sel + [i]
        A = sp.Matrix([[rows[j][1][nm] for nm in names] for j in trial])
        if A.rank() == len(trial):
            sel = trial
        if len(sel) == len(names):
            break
    if len(sel) < len(names):
        raise ValueError('rank deficient structure system')
    A = sp.Matrix([[rows[j][1][nm] for nm in names] for j in sel])
    b = sp.Matrix([rows[j][0] for j in sel])
    sol = A.LUsolve(b).applyfunc(lambda x: sp.expand(sp.radsimp(sp.together(sp.expand(x)))))
    res = [sp.expand(r[0] - sum(sol[k] * r[1][nm] for k, nm in enumerate(names))) for r in rows]
    return {nm: sol[k] for k, nm in enumerate(names)}, res
