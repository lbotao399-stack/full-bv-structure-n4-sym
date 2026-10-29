"""Engine B round 2, T2: one-loop QME defect at the UV boundary, quadratic sector (external c, A|_V).

Q W = (connected graphs with one D_eps vertex),  D_eps = 1/2({I,I}_eps - {I,I}_0) + hbar Delta_eps I - Q I^CT[eps] - {I,I^CT[eps]}_eps.
In the two-point sector at one loop:
  * hbar Delta_eps I = 0 (tadpole colour factor ([X,t_a],t^a)-type traces vanish; K has no A-A component);
  * {I, I^CT}_eps has >= 3 legs;
  * 1/2({I,I}_eps - {I,I}_0) = -1/2 [d_{K_eps - K_0}(I I)]_cross, closed by one P(eps,infty) edge:
       D-graph = -1/2 hbar [d_P d_Khat (I I)]_cc = - d/d delta W_1[P + delta Khat]   (delta odd, Khat = K_eps - K_0),
    evaluated with the toy-verified rule (toy_r2.py): kernel entries carry the parity operator on their row index.
  * -Q I^CT[eps]: from the agreed quadratic counterterm.
Special edge in momentum space: Khat(k) = (e^{-eps k^2} - 1) * F, F^{XY} = -(fibre identity) for the field-antifield pairs
(A,A*), (lambda,lambda*), (c,c*) [Delta_0 = sum d_x d_xi, d_{K_0} = -Delta_0, contraction rule d_K(f_X f_Y) = K^{XY}].
Output: exact eps^{-2}, eps^{-1}, log eps, eps^0 parts of every graph class (units u = hbar c_adj/(16 pi^2), colour (c_0,c_1))."""
import itertools, json, sys, time
from fractions import Fraction as Fr
import sympy as sp
import vertices as V
import loops as L

I = sp.I
eps = sp.Symbol('epsilon', positive=True)
q = L.q
p = sp.symbols('p0:4')
a = sp.symbols('a0:4')
INT6 = ['A', 'lam', 'c', 'As', 'lams', 'cs']
OFF = {}; _o = 0
for s in INT6:
    OFF[s] = _o; _o += V.DIM[s]
NT = _o
PARD = []
for s in INT6:
    PARD += [V.PAR[s]] * V.DIM[s]

def blocks_of(kind):
    if kind == 'P':
        return {'AV': ('A', range(4), 'A', range(4)), 'Ap': ('A', range(4, 10), 'A', range(4, 10)), 'lam': ('lam', range(16), 'lam', range(16)),
                'cAs': ('c', [0], 'As', range(4)), 'Asc': ('As', range(4), 'c', [0])}
    return {'A-A*': ('A', range(10), 'As', range(10)), 'A*-A': ('As', range(10), 'A', range(10)),
            'lam-lam*': ('lam', range(16), 'lams', range(16)), 'lam*-lam': ('lams', range(16), 'lam', range(16)),
            'c-c*': ('c', [0], 'cs', [0]), 'c*-c': ('cs', [0], 'c', [0])}

def P_block(name, mom):
    """P(eps,infty) numerator block (P = N e^{-eps k^2}/k^2), in the 54-dim space."""
    N = sp.zeros(NT, NT)
    rs, rr, cs, cc = blocks_of('P')[name]
    full = L.prop_numerator(list(mom), {('AV', 'AV')} if name == 'AV' else {('Ap', 'Ap')} if name == 'Ap' else
                            {('lam', 'lam')} if name == 'lam' else {('c', 'As')} if name == 'cAs' else {('As', 'c')})
    o4 = L.OFF
    for i in rr:
        for j in cc:
            N[OFF[rs] + i, OFF[cs] + j] = full[o4[rs] + i, o4[cs] + j]
    return N

def K_block(name):
    """fibre of K_eps - K_0 without the scalar (e^{-eps k^2} - 1): F^{XY} = -(identity of the fibre)."""
    F = sp.zeros(NT, NT)
    Minv = V.Mpinv
    if name in ('A-A*', 'A*-A'):
        x, y = ('A', 'As') if name == 'A-A*' else ('As', 'A')
        for M in range(10):
            F[OFF[x] + M, OFF[y] + M] = -1
    elif name == 'lam-lam*':      # K^{u_i v_j} = -(Mp^{-1})_{ji}
        for i in range(16):
            for j in range(16):
                F[OFF['lam'] + i, OFF['lams'] + j] = -Minv[j, i]
    elif name == 'lam*-lam':      # K^{v_j u_i} = -(Mp^{-1})_{ji}
        for j in range(16):
            for i in range(16):
                F[OFF['lams'] + j, OFF['lam'] + i] = -Minv[j, i]
    elif name == 'c-c*':
        F[OFF['c'], OFF['cs']] = -1
    elif name == 'c*-c':
        F[OFF['cs'], OFF['c']] = -1
    return F

def hessian(leg, k, kp):
    pairs = [(x, y) for x in INT6 for y in INT6]
    bl = V.hessian_blocks(V.CUBIC, [leg], list(k), list(kp), pairs)
    H = sp.zeros(NT, NT)
    for (x, y), d in bl.items():
        for key, M in d.items():
            assert key == 'f'
            for i in range(M.rows):
                for j in range(M.cols):
                    if M[i, j] != 0:
                        H[OFF[x] + i, OFF[y] + j] += M[i, j]
    return H.applyfunc(sp.expand)

def PiRows(H):
    return sp.Matrix(H.rows, H.cols, lambda i, j: (-1) ** PARD[i] * H[i, j])

def STr(mats):
    M = mats[0]
    for X in mats[1:]:
        M = (M * X).applyfunc(sp.expand)
    return sp.expand(sum((-1) ** PARD[i] * M[i, i] for i in range(NT) if M[i, i] != 0))

_uint_cache = {}
def uint(a_, b_):
    """int_1^infty u^a (1+u)^{-b} du (exact)."""
    key = (a_, b_)
    if key not in _uint_cache:
        u = sp.Symbol('u', positive=True)
        _uint_cache[key] = sp.nsimplify(sp.integrate(u**a_ * (1 + u)**(-b_), (u, 1, sp.oo)))
    return _uint_cache[key]

def special_bubble(Num, cK, cP):
    """int d^4q/(2pi)^4 Num(q) (e^{-eps (q-cK)^2} - 1) e^{-eps (q-cP)^2}/(q-cP)^2  times 16 pi^2;
    returns {power of eps: coefficient} for powers <= 0 (exact)."""
    bK, bP = sp.symbols('bK bP')
    shift = [bK * cK[m] + bP * cP[m] for m in range(4)]
    C = L.gaussian_moments(Num, shift)
    D2 = sp.expand(sum((cK[m] - cP[m]) ** 2 for m in range(4)))
    out = {}
    # A-term: sum_k 2^{-k} eps^{-1-k} int_1^oo du (1+u)^{-2-k} C_k(bK = 1/(1+u), bP = u/(1+u)) sum_j (-eps D2 u/(1+u))^j / j!
    for k, Ck in C.items():
        Ck = sp.expand(Ck)
        for j in range(0, k + 2):
            pw = -1 - k + j
            if pw > 0:
                continue
            # integrand (1+u)^{-2-k} bK^m bP^n (u/(1+u))^j  -> u^{n+j} (1+u)^{-2-k-m-n-j}
            poly = sp.Poly(Ck, bK, bP)
            val = 0
            for (m_, n_), cf in poly.terms():
                val += cf * uint(n_ + j, 2 + k + m_ + n_ + j)
            term = sp.Rational(1, 2**k) * val * (-D2) ** j / sp.factorial(j)
            out[pw] = sp.expand(out.get(pw, 0) + term)
    # B-term: -sum_k 2^{-k} C_k(bK=0, bP=1) eps^{-1-k}/(1+k)
    for k, Ck in C.items():
        val = sp.expand(Ck.subs({bK: 0, bP: 1}))
        out[-1 - k] = sp.expand(out.get(-1 - k, 0) - sp.Rational(1, 2**k) * val / (1 + k))
    return {k: v for k, v in out.items() if v != 0}

def leg_c():
    return dict(sp='c', u=[sp.Integer(1)], p=list(p))
def leg_A():
    return dict(sp='A', u=list(a) + [0] * 6, p=[-x for x in p])

def sub(Mfull, rs, rr, cs, cc):
    return sp.Matrix(len(list(rr)), len(list(cc)), lambda i, j: Mfull[OFF[rs] + list(rr)[i], OFF[cs] + list(cc)[j]])

def run():
    t0 = time.time()
    legs = {'c': leg_c(), 'A': leg_A()}
    inj = {'c': list(p), 'A': [-x for x in p]}
    results = {}
    for word in (('c', 'A'), ('A', 'c')):
        c1 = inj[word[0]]
        lines = [[0] * 4, c1]                          # line 0: q, line 1: q - c1
        k1, kp1 = [-x for x in q], [q[m] - c1[m] for m in range(4)]
        k2, kp2 = [-(q[m] - c1[m]) for m in range(4)], list(q)
        h1 = PiRows(hessian(legs[word[0]], k1, kp1)); h2 = PiRows(hessian(legs[word[1]], k2, kp2))
        for pos in (0, 1):
            odd_before = 1 if (pos == 1 and word[0] == 'c') else 0
            for kname, kb in blocks_of('K').items():
                Kf = K_block(kname)
                for pname, pb in blocks_of('P').items():
                    Pm = P_block(pname, [q[m] - lines[1 - pos][m] for m in range(4)])
                    b0, b1 = (kb, pb) if pos == 0 else (pb, kb)
                    X0 = sub(Kf if pos == 0 else Pm, *b0); X1 = sub(Pm if pos == 0 else Kf, *b1)
                    H1 = sub(h1, b0[2], b0[3], b1[0], b1[1]); H2 = sub(h2, b1[2], b1[3], b0[0], b0[1])
                    if H1.is_zero_matrix or H2.is_zero_matrix:
                        continue
                    M = (X0 * H1).applyfunc(sp.expand)
                    M = (M * X1).applyfunc(sp.expand)
                    M = (M * H2).applyfunc(sp.expand)
                    Num = sp.expand((-1) ** V.PAR[b0[0]] * M.trace())
                    if Num == 0:
                        continue
                    res = special_bubble(Num, lines[pos], lines[1 - pos])
                    w = sp.Rational(1, 4) * (-1) ** odd_before
                    res = {k: sp.expand(w * v) for k, v in res.items()}
                    key = (word, pos, kname, pname)
                    results[key] = res
                    print(key, {k: str(v) for k, v in res.items()}, round(time.time() - t0, 1), flush=True)
    return results

def classify(key):
    word, pos, kname, pname = key
    # internal legs at the c-vertex: one end of the K line and one end of the P line
    # letter order: word[0] is between line 0 (row side) and line 1; word[1] between line 1 and line 0.
    kspec = blocks_of('K')[kname]; pspec = blocks_of('P')[pname]
    # line 'pos' is K: rows species = col side of the previous letter... record species at the c-vertex
    if word[0] == 'c':
        # c-letter rows = col species of line 0, cols = row species of line 1
        rows_line, cols_line = (kspec, pspec) if pos == 0 else (pspec, kspec)
        c_int = (rows_line[2], cols_line[0])
        a_int = (cols_line[2], rows_line[0])
    else:
        rows_line, cols_line = (kspec, pspec) if pos == 0 else (pspec, kspec)
        a_int = (rows_line[2], cols_line[0])
        c_int = (cols_line[2], rows_line[0])
    kl = {'A-A*': 'A-A*', 'A*-A': 'A-A*', 'lam-lam*': 'lambda-lambda*', 'lam*-lam': 'lambda-lambda*', 'c-c*': 'c-c*', 'c*-c': 'c-c*'}[kname]
    pl = {'AV': 'A|_V', 'Ap': 'A|_Wperp', 'lam': 'lambda', 'cAs': 'c-A*', 'Asc': 'c-A*'}[pname]
    def vname(ext, ints):
        s = set(ints)
        if ext == 'c':
            if s == {'A', 'As'}: return 'AcA*'
            if s == {'lam', 'lams'}: return 'c lambda lambda*'
            if s == {'c', 'cs'}: return 'ccc*'
        else:
            if s == {'A'}: return 'AAA'
            if s == {'c', 'As'}: return 'AcA*'
            if s == {'lam'}: return 'lambda lambda A'
        return str(ints)
    return f'v1={vname("c", c_int)}[c] v2={vname("A", a_int)}[A|_V] special edge {kl} (K_eps - K_0), P edge {pl}'

if __name__ == '__main__':
    res = run()
    classes = {}
    for key, r in res.items():
        cl = classify(key)
        d = classes.setdefault(cl, {})
        for k, v in r.items():
            d[k] = sp.expand(d.get(k, 0) + v)
    tot = {}
    for cl, d in classes.items():
        for k, v in d.items():
            tot[k] = sp.expand(tot.get(k, 0) + v)
    pa = sum(p[m] * a[m] for m in range(4)); p2 = sum(x * x for x in p)
    def as_struct(v):
        v = sp.expand(v)
        c1 = sp.simplify(v / (I * pa)) if v != 0 else 0
        return str(sp.factor(c1))
    out = dict(classes={cl: {str(k): str(v) for k, v in d.items()} for cl, d in classes.items()},
               classes_over_i_pa={cl: {str(k): as_struct(v) for k, v in d.items()} for cl, d in classes.items()},
               total={str(k): str(v) for k, v in tot.items()}, total_over_i_pa={str(k): as_struct(v) for k, v in tot.items()})
    with open('r2_T2_raw.json', 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
