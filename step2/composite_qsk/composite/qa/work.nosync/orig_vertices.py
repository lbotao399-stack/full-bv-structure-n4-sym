"""Engine B: action of [ESW (31)] reduced to R^4, evaluated on plane waves (code-only component bases).
Vertex function convention: I evaluated on sum_j r_j u_j e^{i p_j.x}; V = coefficient of the product of the odd r_j in leg order
(multilinear part, (2pi)^4 delta stripped).  Colour is stripped: every term carries one of the colour tensors
  'f'   : ftil(c1,c2,c3) = ([c1,c2],c3)  (legs in leg order),
  'K12','K13','K14' : ([c1,c2],[c3,c4]), ([c1,c3],[c2,c4]), ([c1,c4],[c2,c3]),
  'k'   : (c1,c2)  (quadratic terms),
and V is returned as {colour key: Lorentz value}."""
import itertools
import sympy as sp
from realization import R, Rb, pairing

I = sp.I
Mp = pairing()
Mpinv = Mp.inv()
MR = [Mp * R[M] for M in range(10)]        # symmetric 16x16 : (Q1, rho(e_M) Q2) = Q1^T MR_M Q2 for Q1,Q2 in Sigma
DIM = {'A': 10, 'lam': 16, 'c': 1, 'As': 10, 'lams': 16, 'cs': 1}
PAR = {'A': 0, 'lam': 1, 'c': 1, 'As': 1, 'lams': 0, 'cs': 0}
V4D = range(4)

def dot4(a, b):
    return sum(a[m] * b[m] for m in range(4))

# ---- action terms: (name, coefficient, slot species, Lorentz function(slot legs), colour slot-structure) ----
def L_AAA(s):   # -1/2 sum_{mu,N} (d_mu A_N , [A_mu , A_N]) ; slot0 = differentiated A_N, slot1 = A_mu, slot2 = A_N
    a, b, c = s
    return sum(I * a['p'][m] * a['u'][N] * b['u'][m] * c['u'][N] for m in V4D for N in range(10))
def L_AAAA(s):  # -1/8 sum_{M,N} ([A_M,A_N],[A_M,A_N])
    a, b, c, d = s
    return sum(a['u'][M] * b['u'][N] * c['u'][M] * d['u'][N] for M in range(10) for N in range(10))
def L_lamlamA(s):  # 1/2 (lam, rho([A,lam])) = 1/2 sum_M lam1^T MR_M lam2 A_M ; colour (c_lam1,[c_A,c_lam2])
    l1, A, l2 = s
    u1 = sp.Matrix(l1['u']); u2 = sp.Matrix(l2['u'])
    return sum(A['u'][M] * (u1.T * MR[M] * u2)[0] for M in range(10))
def L_AcAs(s):  # ([A_M, c], A*_M)
    A, c, As = s
    return sum(A['u'][M] * As['u'][M] for M in range(10)) * c['u'][0]
def L_clamlams(s):  # ([c,lam],lam*) = c * lam^T Mp lam*
    c, l, ls = s
    return c['u'][0] * (sp.Matrix(l['u']).T * Mp * sp.Matrix(ls['u']))[0]
def L_cccs(s):
    return s[0]['u'][0] * s[1]['u'][0] * s[2]['u'][0]
# quadratic (free) terms
def L_kinA(s):  # -1/4 (dA,dA) = -1/4 sum_{M<N} (d_M A_N - d_N A_M)^2  (d only along R^4)
    a, b = s
    def dA(x, M, N):
        v = 0
        if M < 4: v += I * x['p'][M] * x['u'][N]
        if N < 4: v -= I * x['p'][N] * x['u'][M]
        return v
    return sum(dA(a, M, N) * dA(b, M, N) for M in range(10) for N in range(M + 1, 10))
def L_kinlam(s):  # 1/2 (lam, dslash lam) = 1/2 sum_mu lam1^T MR_mu d_mu lam2
    a, b = s
    return sum((sp.Matrix(a['u']).T * MR[m] * sp.Matrix(b['u']))[0] * I * b['p'][m] for m in V4D)
def L_kingh(s):  # (dc, A*)
    c, As = s
    return sum(I * c['p'][m] * c['u'][0] * As['u'][m] for m in V4D)

TERMS = {
    'AAA': (sp.Rational(-1, 2), ('A', 'A', 'A'), L_AAA, ('f', (0, 1, 2))),
    'AAAA': (sp.Rational(-1, 8), ('A', 'A', 'A', 'A'), L_AAAA, ('bb', (0, 1, 2, 3))),
    'lamlamA': (sp.Rational(1, 2), ('lam', 'A', 'lam'), L_lamlamA, ('f', (0, 1, 2))),
    'AcAs': (sp.Integer(1), ('A', 'c', 'As'), L_AcAs, ('f', (0, 1, 2))),
    'clamlams': (sp.Integer(1), ('c', 'lam', 'lams'), L_clamlams, ('f', (0, 1, 2))),
    'cccs': (sp.Rational(1, 2), ('c', 'c', 'cs'), L_cccs, ('f', (0, 1, 2))),
    'kinA': (sp.Rational(-1, 4), ('A', 'A'), L_kinA, ('k', (0, 1))),
    'kinlam': (sp.Rational(1, 2), ('lam', 'lam'), L_kinlam, ('k', (0, 1))),
    'kingh': (sp.Integer(1), ('c', 'As'), L_kingh, ('k', (0, 1))),
}
CUBIC = ['AAA', 'lamlamA', 'AcAs', 'clamlams', 'cccs']
QUARTIC = ['AAAA']
FREE = ['kinA', 'kinlam', 'kingh']

def perm_sign(seq):
    seq = list(seq); s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s

def colour_key(struct, sigma):
    kind, sl = struct
    legs = [sigma[x] for x in sl]
    if kind == 'f':
        return 'f', perm_sign(legs)
    if kind == 'k':
        return 'k', 1
    # ([a,b],[c,d]) canonical
    a, b, c, d = legs; s = 1
    if 0 not in (a, b):
        a, b, c, d = c, d, a, b
    if a != 0:
        a, b = b, a; s = -s
    if c > d:
        c, d = d, c; s = -s
    return {1: 'K12', 2: 'K13', 3: 'K14'}[b], s

def vertex(term, legs):
    """legs: list of dict(sp=species, u=list, p=list of 4).  Returns {colour key: value}."""
    coef, slots, Lf, struct = TERMS[term]
    out = {}
    n = len(slots)
    if len(legs) != n:
        return out
    for sigma in itertools.permutations(range(n)):      # slot s gets leg sigma[s]
        if any(legs[sigma[s]]['sp'] != slots[s] for s in range(n)):
            continue
        odd_in_slot_order = [sigma[s] for s in range(n) if PAR[slots[s]]]
        rs = perm_sign(odd_in_slot_order)
        key, cs = colour_key(struct, sigma)
        val = Lf([legs[sigma[s]] for s in range(n)])
        out[key] = out.get(key, 0) + coef * rs * cs * val
    return {k: sp.expand(v) for k, v in out.items()}

def total_vertex(terms, legs):
    out = {}
    for t in terms:
        for k, v in vertex(t, legs).items():
            out[k] = sp.expand(out.get(k, 0) + v)
    return out

def sym_leg(spc, name, p):
    return dict(sp=spc, u=list(sp.symbols(f'{name}0:{DIM[spc]}')), p=list(p))

def hessian_blocks(terms, ext_legs, k, kp, species_pairs):
    """Colour-stripped Hessian h_{JK} = D_{x_J} D_{y_K} I(ext + x e^{ik.} + y e^{ik'.}) with the r's of the external legs
    stripped (they stand in front, in leg order).  Returns {(sJ,sK): {colour key: matrix}}."""
    out = {}
    epar = sum(PAR[l['sp']] for l in ext_legs)
    for sJ, sK in species_pairs:
        X = sym_leg(sJ, 'xx', k); Y = sym_leg(sK, 'yy', kp)
        V = total_vertex(terms, ext_legs + [X, Y])
        pJ, pK = PAR[sJ], PAR[sK]
        sign = (-1) ** (pK * (epar + pJ) + pJ * epar)
        blk = {}
        for key, val in V.items():
            M = sp.zeros(DIM[sJ], DIM[sK])
            poly = sp.Poly(val, *(X['u'] + Y['u']))
            for mon, cf in poly.terms():
                iJ = [i for i in range(DIM[sJ]) if mon[i] == 1]
                iK = [i for i in range(DIM[sK]) if mon[DIM[sJ] + i] == 1]
                assert len(iJ) == 1 and len(iK) == 1 and sum(mon) == 2
                M[iJ[0], iK[0]] += sign * cf
            blk[key] = M
        out[(sJ, sK)] = blk
    return out
