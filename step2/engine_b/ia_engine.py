"""Engine B, Step 2, Part C: graphs of type I-a (special edge K_eps - K_0 ON the loop, 1PI), exact eps-expansion.

D-graph formula (toy-verified in engine_b/toy_r2.py, used in engine_b/r2_T2.py for the two-point sector):
   graphs with one special edge = - d/d delta W_1[P + delta Khat],  Khat = K_eps - K_0,
   coefficient of (r's in leg order) = sum_{words} sum_{pos} (1/(2n)) * rsign(word) * (-1)^{odd legs before pos}
                                       * colour(word) * STr[X_0 Pi h_0 X_1 Pi h_1 ... X_{n-1} Pi h_{n-1}],
   X_pos = fibre of K_0 (F^{XY} = -identity of the field-antifield fibre) times the scalar (e^{-eps k^2} - 1),
   X_e   = N(k_e) e^{-eps k_e^2}/k_e^2 (P(eps,infty)) for e != pos;  h = letter Hessians (cubic: one leg, quartic: two A legs).
Schwinger form: (e^{-eps k^2} - 1) = - int_0^eps k^2 e^{-t k^2} dt,  e^{-eps k^2}/k^2 = int_eps^infty e^{-t k^2} dt.
Loop momentum: Gaussian (loops.gaussian_moments), int d^4q/(2pi)^4 -> (16 pi^2)^{-1} T^{-2} sum_k C_k (2T)^{-k} e^{-T rho},
rho = sum_e beta_e |P_e|^2 - |sum_e beta_e P_e|^2,  beta_e = t_e/T.  With t = eps u (u_special in [0,1], others in [1,oo)):
   eps^{E-2-k+j} * 2^{-k} (-1)^j/j! * [C_k rho^j](beta) * T^{-2-k+j},
C_k is homogeneous of degree (deg Num - 2k) in the leg momenta, so the power of eps fixes the momentum degree:
   eps^n  <->  momentum degree omega + 2n  (omega = superficial degree).  No half-integer powers, no log eps:
   the eps^n coefficient (n = -1, 0) is a finite sum of convergent master integrals (masters.py; integrand ~ T^{-E+n}
   over E-1 non-compact directions), and the remainder is O(eps log eps) (IR region k ~ p: the special edge supplies
   an extra k^2, so the non-analytic part of g(p) = graph at eps = 1 starts at degree omega + 2).
Explicit component bases live only in code (engine_b/realization.py)."""
import itertools, sys, os, time, json, random
from fractions import Fraction as Fr
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
EB = os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b'))
if EB not in sys.path:
    sys.path.insert(0, EB)
import vertices as V
import loops as L
import r2_engine as R
import masters as M

I = sp.I
q = L.q
INT6 = ['A', 'lam', 'c', 'As', 'lams', 'cs']
PARD = {s: V.PAR[s] for s in INT6}

# ------------------------------------------------------------------ edges
AVc = list(range(4)); APc = list(range(4, 10)); L16 = list(range(16))
PEDGE = {'AV': ('A', AVc, 'A', AVc), 'Ap': ('A', APc, 'A', APc), 'lam': ('lam', L16, 'lam', L16),
         'cAs': ('c', [0], 'As', AVc), 'Asc': ('As', AVc, 'c', [0])}
KEDGE = {'AV-AsV': ('A', AVc, 'As', AVc), 'AsV-AV': ('As', AVc, 'A', AVc),
         'Ap-Asp': ('A', APc, 'As', APc), 'Asp-Ap': ('As', APc, 'A', APc),
         'lam-lams': ('lam', L16, 'lams', L16), 'lams-lam': ('lams', L16, 'lam', L16),
         'c-cs': ('c', [0], 'cs', [0]), 'cs-c': ('cs', [0], 'c', [0])}
ELABEL = {'AV': 'A|_V', 'Ap': 'A|_Wperp', 'lam': 'lambda', 'cAs': 'c-A*', 'Asc': 'c-A*',
          'AV-AsV': 'K[A|_V-A*|_V]', 'AsV-AV': 'K[A|_V-A*|_V]', 'Ap-Asp': 'K[A|_Wperp-A*|_Wperp]',
          'Asp-Ap': 'K[A|_Wperp-A*|_Wperp]', 'lam-lams': 'K[lambda-lambda*]', 'lams-lam': 'K[lambda-lambda*]',
          'c-cs': 'K[c-c*]', 'cs-c': 'K[c-c*]'}

_PN = None
def P_numerator(name, mom):
    """numerator block of P(eps,infty) at line momentum mom (P = N e^{-eps k^2}/k^2), rows x cols of the edge."""
    global _PN
    rs, rr, cs, cc = PEDGE[name]
    mask = {'AV': {('AV', 'AV')}, 'Ap': {('Ap', 'Ap')}, 'lam': {('lam', 'lam')}, 'cAs': {('c', 'As')}, 'Asc': {('As', 'c')}}[name]
    full = L.prop_numerator(list(mom), mask)
    o4 = L.OFF
    return sp.Matrix(len(rr), len(cc), lambda i, j: full[o4[rs] + rr[i], o4[cs] + cc[j]])

_MPINV = V.Mpinv
def K_fibre(name):
    """fibre of K_0: F^{XY} = -(identity of the fibre); lambda: K^{u_i v_j} = -(Mp^{-1})_{ji} (r2_T2.K_block)."""
    rs, rr, cs, cc = KEDGE[name]
    if rs == 'lam':
        return sp.Matrix(16, 16, lambda i, j: -_MPINV[j, i])
    if rs == 'lams':
        return sp.Matrix(16, 16, lambda j, i: -_MPINV[j, i])
    return -sp.eye(len(rr))

# ------------------------------------------------------------------ letters
KS = sp.symbols('hk0:4'); KPS = sp.symbols('hkp0:4')

class Letter:
    """Hessian of the cubic (one external leg) or quartic (two external A legs) part of I, colour stripped.
    blocks[(sJ, sK)][key] = sympy matrix (DIM sJ x DIM sK), affine in (KS, KPS); row slot plane wave e^{i KS.x},
    column slot e^{i KPS.x} (vertices.hessian_blocks)."""
    def __init__(self, kind, legs, leg_ids):
        self.kind = kind; self.leg_ids = tuple(leg_ids)
        terms = V.CUBIC if kind == 'cubic' else V.QUARTIC
        pairs = [(a, b) for a in INT6 for b in INT6]
        bl = V.hessian_blocks(terms, [legs[i] for i in leg_ids], list(KS), list(KPS), pairs)
        self.blocks = {}
        for (a, b), d in bl.items():
            for key, Mx in d.items():
                if not Mx.is_zero_matrix:
                    self.blocks.setdefault((a, b), {})[key] = Mx
        self.ext_species = [legs[i]['sp'] for i in leg_ids]

    def block(self, a, rr, b, cc, key, kv, kpv):
        Mx = self.blocks[(a, b)][key]
        sub = Mx.extract(list(rr), list(cc))
        rep = {KS[m]: kv[m] for m in range(4)}; rep.update({KPS[m]: kpv[m] for m in range(4)})
        sub = sub.xreplace(rep).applyfunc(sp.expand)
        sg = (-1) ** PARD[a]                     # parity operator on the row index (Pi h)
        return sub * sg

def vertex_name(kind, ext_sp, a, b):
    if kind == 'quartic':
        return 'AAAA'
    s = {a, b}
    if ext_sp == 'c':
        return {frozenset({'A', 'As'}): 'AcA*', frozenset({'lam', 'lams'}): 'c lambda lambda*', frozenset({'c', 'cs'}): 'ccc*'}.get(frozenset(s), '?')
    if ext_sp == 'A':
        return {frozenset({'A'}): 'AAA', frozenset({'lam'}): 'lambda lambda A', frozenset({'c', 'As'}): 'AcA*'}.get(frozenset(s), '?')
    if ext_sp == 'lam':
        return {frozenset({'lam', 'A'}): 'lambda lambda A', frozenset({'c', 'lams'}): 'c lambda lambda*'}.get(frozenset(s), '?')
    if ext_sp == 'As':
        return {frozenset({'A', 'c'}): 'AcA*'}.get(frozenset(s), '?')
    return '?'

# ------------------------------------------------------------------ Schwinger integration of one numerator
BET = sp.symbols('bb0:4')

def schwinger_coeffs(Num, curs, s):
    """Num: polynomial in q (loop momentum of line 0), lines e with momentum q - curs[e] (numeric), special line s.
    Returns {n: masters-number} for n = -1, 0 (coefficients of eps^n), units (16 pi^2)^{-1}; includes the special-edge
    factor -k_s^2."""
    E = len(curs)
    ks = [q[m] - curs[s][m] for m in range(4)]
    Nt = sp.expand(-Num * sum(x * x for x in ks))
    shift = [sum(BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
    C = L.gaussian_moments(Nt, shift)
    P2 = lambda v: sum(x * x for x in v)
    rho = sp.expand(sum(BET[e] * P2(curs[e]) for e in range(E)) - P2(shift))
    order = [s] + [e for e in range(E) if e != s]        # u1 = special
    out = {}
    for n in (-1, 0):
        acc = {}
        for k, Ck in C.items():
            j = k + 2 - E + n
            if j < 0:
                continue
            expr = sp.expand(Ck * rho ** j) * sp.Rational((-1) ** j, 2 ** k * sp.factorial(j))
            if expr == 0:
                continue
            poly = sp.Poly(expr, *BET[:E])
            for mon, cf in poly.terms():
                al = [mon[e] for e in order]
                mm = E - n + sum(al)
                J = M.J3(al[0], al[1], al[2], mm) if E == 3 else M.J2(al[0], al[1], mm)
                for kk, v in J.items():
                    acc[kk] = acc.get(kk, 0) + cf * sp.Rational(v.numerator, v.denominator)
        out[n] = {kk: sp.expand(v) for kk, v in acc.items() if sp.expand(v) != 0}
    return out

# ------------------------------------------------------------------ main enumeration
def partitions(legs):
    idx = list(range(len(legs)))
    out = [[('cubic', (i,)) for i in idx]]
    Aidx = [i for i in idx if legs[i]['sp'] == 'A']
    for pr in itertools.combinations(Aidx, 2):
        out.append([('quartic', pr)] + [('cubic', (i,)) for i in idx if i not in pr])
    return out

def canon_class(seq):
    """seq: list over letters k of (vertex name, ext label, label of the line AFTER letter k). Canonical up to rotation
    and reflection."""
    n = len(seq)
    cands = []
    for rot in range(n):
        s = seq[rot:] + seq[:rot]
        cands.append(tuple(s))
        verts = [x[:2] for x in s][::-1]
        lines = [x[2] for x in s]
        lines_ref = [lines[(n - 2 - j) % n] for j in range(n)]
        cands.append(tuple((verts[j][0], verts[j][1], lines_ref[j]) for j in range(n)))
    return min(cands)

def ext_label(leg):
    sp_ = leg['sp']
    if sp_ in ('A', 'As'):
        nz = [i for i, x in enumerate(leg['u']) if x != 0]
        tag = 'V' if all(i < 4 for i in nz) else ('Wperp' if all(i >= 4 for i in nz) else 'all')
        return ('A|_' if sp_ == 'A' else 'A*|_') + tag
    return {'c': 'c', 'lam': 'lambda', 'cs': 'c*', 'lams': 'lambda*'}[sp_]

def evaluate(legs, verbose=False, only_classes=None, collect=None):
    """All I-a graphs for the given external legs (numeric polarisations and momenta).  Returns {class: {n: number}}
    with number = {'1','log2','log3'} -> complex rational, coefficient of unit * ([c_0,c_1],c_2) at eps^n,
    unit = hbar tr_g(ad ad)/(16 pi^2 (X,Y))."""
    nleg = len(legs)
    odd = [V.PAR[l['sp']] for l in legs]
    letters_cache = {}
    def letter(kind, ids):
        key = (kind, ids)
        if key not in letters_cache:
            letters_cache[key] = Letter(kind, legs, ids)
        return letters_cache[key]
    res = {}
    nterms = 0
    allE = list(PEDGE) + list(KEDGE)
    for part in partitions(legs):
        n = len(part)
        for perm in itertools.permutations(range(n)):
            word = [part[i] for i in perm]
            curs = [[sp.Integer(0)] * 4]
            for (kind, ids) in word:
                inj = [sum(legs[i]['p'][m] for i in ids) for m in range(4)]
                curs.append([curs[-1][m] + inj[m] for m in range(4)])
            assert all(x == 0 for x in curs[-1])
            curs = curs[:n]
            odd_seq = [i for (kind, ids) in word for i in ids if odd[i]]
            rsign = R.perm_parity(odd_seq)
            lets = [letter(kind, ids) for (kind, ids) in word]
            for pos in range(n):
                odd_before = sum(odd[i] for (kind, ids) in word[:pos] for i in ids)
                w = sp.Rational(1, 2 * n) * rsign * (-1) ** odd_before
                choices = [list(KEDGE) if e == pos else list(PEDGE) for e in range(n)]
                for edges in itertools.product(*choices):
                    spec = [KEDGE[edges[e]] if e == pos else PEDGE[edges[e]] for e in range(n)]
                    ok = True
                    for kk in range(n):
                        a = spec[kk][2]; b = spec[(kk + 1) % n][0]
                        if (a, b) not in lets[kk].blocks:
                            ok = False; break
                    if not ok:
                        continue
                    keysets = [list(lets[kk].blocks[(spec[kk][2], spec[(kk + 1) % n][0])].keys()) for kk in range(n)]
                    for keys in itertools.product(*keysets):
                        Ms = [R.letter_M(word[kk][0], word[kk][1], keys[kk]) for kk in range(n)]
                        col = R.reduce_traces(R.colour_of_word(Ms), list(range(nleg)))
                        ckey = 'cadj f' if nleg == 3 else 'cadj k'
                        cf = col.get(ckey, 0)
                        assert set(col) <= {ckey}, col
                        if cf == 0:
                            continue
                        mats = []
                        for kk in range(n):
                            line_mom = [q[m] - curs[kk][m] for m in range(4)]
                            if kk == pos:
                                mats.append(K_fibre(edges[kk]))
                            else:
                                mats.append(P_numerator(edges[kk], line_mom))
                            kv = [-(q[m] - curs[kk][m]) for m in range(4)]
                            nxt = curs[(kk + 1) % n] if kk + 1 < n else [0] * 4
                            kpv = [q[m] - nxt[m] for m in range(4)]
                            a, rr = spec[kk][2], spec[kk][3]
                            b, cc = spec[(kk + 1) % n][0], spec[(kk + 1) % n][1]
                            mats.append(lets[kk].block(a, rr, b, cc, keys[kk], kv, kpv))
                        Mx = mats[0]
                        for X in mats[1:]:
                            Mx = (Mx * X).applyfunc(sp.expand)
                            if Mx.is_zero_matrix:
                                break
                        if Mx.is_zero_matrix:
                            continue
                        Num = sp.expand((-1) ** PARD[spec[0][0]] * Mx.trace())
                        if Num == 0:
                            continue
                        fac = w * sp.Rational(cf.numerator, cf.denominator)
                        if collect is not None:
                            collect.append(dict(fac=fac, Num=Num, curs=[list(c) for c in curs], s=pos, word=word, edges=edges))
                            val = {}
                        else:
                            val = schwinger_coeffs(Num, curs, pos)
                        seq = []
                        for kk in range(n):
                            kind, ids = word[kk]
                            vn = vertex_name(kind, legs[ids[0]]['sp'], spec[kk][2], spec[(kk + 1) % n][0])
                            el = ','.join(sorted(ext_label(legs[i]) for i in ids))
                            seq.append((vn, el, ELABEL[edges[(kk + 1) % n]]))
                        cl = canon_class(seq)
                        nterms += 1
                        if collect is not None:
                            collect[-1]['class'] = class_name(cl)
                        d = res.setdefault(cl, {})
                        for nn, num in val.items():
                            dd = d.setdefault(nn, {})
                            for kk2, v in num.items():
                                dd[kk2] = sp.expand(dd.get(kk2, 0) + fac * v)
                        if verbose:
                            print(' ', [w_[1] for w_ in word], 'pos', pos, edges, keys, 'col', cf, {nn: {k: str(sp.expand(fac * v)) for k, v in num.items()} for nn, num in val.items()}, flush=True)
    return res, nterms

def class_name(cl):
    return ' | '.join(f'{v}[{e}] -{ln}-' for v, e, ln in cl)
