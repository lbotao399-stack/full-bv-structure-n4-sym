"""Engine B, round 2: exact one-loop log eps coefficients of cubic and quartic vertex functions.

Formula (round-1 verified one-loop formula, expanded multilinearly in the external legs):
  W_1 = (hbar/2) sum_n (1/n) STr[(P Pi H)^n],  H = sum over letters (single-leg cubic Hessians, two-leg quartic Hessians)
  => coefficient of r_1 r_2 ... (odd legs, canonical order) =
     sum_{words}  (hbar/(2n)) * sign(r-reordering) * colour(word) * int d^4q/(2pi)^4 Num_word(q) prod_f e^{-eps(q-P_f)^2}/(q-P_f)^2 .
UV boundary, method B: eps d/deps of a graph = -sum_e (graph with t_e pinned at eps); the pinned graphs are Taylor-expanded in the
external momenta (lambda-scaling, all external momenta -> lambda P), and the loop momentum integral is done exactly
(radial: int_0^infty q^{3+2m-2b} e^{-E eps q^2} dq; angular: exact average over the 24 short roots of F4, a spherical 5-design).
Summing the pinned edges gives, for a graph of superficial degree w,
   log eps coefficient  = -(16 pi^2)^{-1} < [lambda^w] Num(qhat, lambda P) prod_f |qhat - lambda P_f|^{-2} >_{|qhat|=1},
   eps^{-1/2} coefficient (w = 1): proportional to < Num(qhat, 0) >, an odd polynomial: zero.
All arithmetic exact (Gaussian integers with a tracked denominator).  Explicit component bases live only here."""
import itertools, math, os, pickle, random, time
from fractions import Fraction as Fr
from functools import reduce
import numpy as np
import sympy as sp
import vertices as V
from realization import Rb

HERE = os.path.dirname(os.path.abspath(__file__))
PAR = dict(V.PAR)
DIM = dict(V.DIM)

# ------------------------------------------------------------------ exact Gaussian-integer matrices
def _lcm(a, b):
    return a * b // math.gcd(a, b)

def to_fr(c):
    c = sp.nsimplify(c)
    re, im = sp.re(c), sp.im(c)
    return Fr(int(sp.fraction(re)[0]), int(sp.fraction(re)[1])), Fr(int(sp.fraction(im)[0]), int(sp.fraction(im)[1]))

class GZ:
    """matrix (re + i im)/den with re, im object arrays of python ints."""
    __slots__ = ('re', 'im', 'den')
    def __init__(self, re, im, den=1):
        self.re = re; self.im = im; self.den = den
    @staticmethod
    def zeros(r, c):
        return GZ(np.zeros((r, c), dtype=object) * 0, np.zeros((r, c), dtype=object) * 0, 1)
    @staticmethod
    def from_frac(Mre, Mim):
        """Mre, Mim: nested lists / arrays of Fractions."""
        Mre = np.array(Mre, dtype=object); Mim = np.array(Mim, dtype=object)
        den = 1
        for x in list(Mre.flat) + list(Mim.flat):
            den = _lcm(den, Fr(x).denominator)
        f = np.vectorize(lambda x: int(Fr(x) * den), otypes=[object])
        return GZ(f(Mre), f(Mim), den)
    def mm(self, o):
        return GZ(self.re.dot(o.re) - self.im.dot(o.im), self.re.dot(o.im) + self.im.dot(o.re), self.den * o.den)
    def add(self, o):
        d = _lcm(self.den, o.den); a = d // self.den; b = d // o.den
        return GZ(self.re * a + o.re * b, self.im * a + o.im * b, d)
    def scale_int(self, s):          # s python int
        return GZ(self.re * s, self.im * s, self.den)
    def mul_i(self):
        return GZ(-self.im, self.re, self.den)
    def rowsign(self, signs):
        s = np.array(signs, dtype=object).reshape(-1, 1)
        return GZ(self.re * s, self.im * s, self.den)
    def trace(self):
        return Fr(int(sum(np.diag(self.re))), self.den), Fr(int(sum(np.diag(self.im))), self.den)
    def is_zero(self):
        return not (np.any(self.re != 0) or np.any(self.im != 0))
    def sub(self, rows, cols):
        return GZ(self.re[np.ix_(rows, cols)], self.im[np.ix_(rows, cols)], self.den)

class Dual:
    """M0 + lambda M1 (mod lambda^2)."""
    __slots__ = ('a', 'b')
    def __init__(self, a, b=None):
        self.a = a; self.b = b
    def mm(self, o):
        a = self.a.mm(o.a)
        b = None
        if self.b is not None:
            b = self.b.mm(o.a)
        if o.b is not None:
            t = self.a.mm(o.b)
            b = t if b is None else b.add(t)
        return Dual(a, b)

# ------------------------------------------------------------------ propagators (closed forms, round-1 derived and checked)
MPINV = np.array(V.Mpinv.tolist(), dtype=object)
RB = [np.array(Rb[m].tolist(), dtype=object) for m in range(4)]

def _np_sym_to_GZ(M):
    """M: sympy-number matrix (list of lists) -> GZ."""
    re = [[to_fr(x)[0] for x in row] for row in M]
    im = [[to_fr(x)[1] for x in row] for row in M]
    return GZ.from_frac(re, im)

_NLAM_BASIS = None
def nlam_basis():
    """N_lam(q) = i sum_m q_m Rb[m] Mpinv  (16x16): returns list of 4 GZ matrices B_m with N = i sum q_m B_m."""
    global _NLAM_BASIS
    if _NLAM_BASIS is None:
        out = []
        for m in range(4):
            M = (Rb[m] * V.Mpinv)
            out.append(_np_sym_to_GZ(M.tolist()))
        _NLAM_BASIS = out
    return _NLAM_BASIS

EDGE = {   # edge type: (row species, row comps, col species, col comps)
    'AV': ('A', list(range(4)), 'A', list(range(4))),
    'Ap': ('A', list(range(4, 10)), 'A', list(range(4, 10))),
    'lam': ('lam', list(range(16)), 'lam', list(range(16))),
    'cAs': ('c', [0], 'As', list(range(4))),
    'Asc': ('As', list(range(4)), 'c', [0]),
}
LINE_LABEL = {'AV': 'A|_V', 'Ap': 'A|_Wperp', 'lam': 'lambda', 'cAs': 'c-A*', 'Asc': 'c-A*'}

def prop_block(etype, qv, cur, lam_dual):
    """P numerator block at line momentum q - lambda*cur, as Dual (lambda part only if lam_dual)."""
    if etype in ('AV', 'Ap'):
        n = 4 if etype == 'AV' else 6
        re = np.zeros((n, n), dtype=object) * 0
        for i in range(n):
            re[i, i] = 2
        return Dual(GZ(re, np.zeros((n, n), dtype=object) * 0, 1), None)
    if etype == 'lam':
        B = nlam_basis()
        def lin(v):
            acc = None
            for m in range(4):
                if v[m] == 0:
                    continue
                t = B[m].scale_int(int(v[m]))
                acc = t if acc is None else acc.add(t)
            if acc is None:
                acc = GZ.zeros(16, 16)
            return acc.mul_i()
        a = lin(qv)
        b = lin([-x for x in cur]) if lam_dual else None
        return Dual(a, b)
    # ghost line: N_{c,As_m} = N_{As_m,c} = -i (q - lambda cur)_m
    def ghost(v):
        im = np.array([-int(v[m]) for m in range(4)], dtype=object)
        re = np.zeros(4, dtype=object) * 0
        if etype == 'cAs':
            return GZ(re.reshape(1, 4), im.reshape(1, 4), 1)
        return GZ(re.reshape(4, 1), im.reshape(4, 1), 1)
    a = ghost(qv)
    b = ghost([-x for x in cur]) if lam_dual else None
    return Dual(a, b)

# ------------------------------------------------------------------ letters (Hessians with stripped r's, colour stripped by key)
def _affine_coeffs(expr, syms):
    expr = sp.expand(expr)
    if expr == 0:
        return None
    poly = sp.Poly(expr, *syms)
    out = {}
    for mon, cf in poly.terms():
        assert sum(mon) <= 1, 'letter not affine in momenta'
        key = None if sum(mon) == 0 else syms[mon.index(1)]
        out[key] = out.get(key, 0) + cf
    return out

class Letter:
    """kind 'cubic' (one external leg) or 'quartic' (two external legs).  ext: list of (leg index, species, numeric u).
    Holds, for every internal species pair and colour key, the affine data of h_{JK}(k, kp, p_ext)."""
    def __init__(self, kind, ext, int_species, terms=None):
        self.kind = kind; self.ext = ext
        ks = list(sp.symbols('k0:4')); kps = list(sp.symbols('kp0:4'))
        pes = [list(sp.symbols(f'pe{j}_0:4')) for j in range(len(ext))]
        legs = [dict(sp=s, u=[sp.Integer(x) for x in u], p=pes[j]) for j, (_, s, u) in enumerate(ext)]
        if terms is None:
            terms = V.CUBIC if kind == 'cubic' else V.QUARTIC
        pairs = [(a, b) for a in int_species for b in int_species]
        blocks = V.hessian_blocks(terms, legs, ks, kps, pairs)
        syms = ks + kps + [x for pe in pes for x in pe]
        self.data = {}   # (sJ,sK) -> key -> dict(var -> (re Fr matrix, im Fr matrix))
        for (a, b), bl in blocks.items():
            for key, M in bl.items():
                ent = {}
                nz = False
                for i in range(M.rows):
                    for j in range(M.cols):
                        if M[i, j] == 0:
                            continue
                        co = _affine_coeffs(M[i, j], syms)
                        if co is None:
                            continue
                        nz = True
                        for var, cf in co.items():
                            ent.setdefault(var, {})[(i, j)] = to_fr(cf)
                if nz:
                    self.data.setdefault((a, b), {})[key] = (M.rows, M.cols, ent)
        self.ks = ks; self.kps = kps; self.pes = pes

    def keys(self, a, b):
        return list(self.data.get((a, b), {}).keys())

    def _decomp(self, a, b, key, rows, cols):
        ck = (a, b, key, tuple(rows), tuple(cols))
        if not hasattr(self, '_cache'):
            self._cache = {}
        if ck in self._cache:
            return self._cache[ck]
        if (a, b) not in self.data or key not in self.data[(a, b)]:
            self._cache[ck] = None
            return None
        R, C, ent = self.data[(a, b)][key]
        ridx = {r: i for i, r in enumerate(rows)}; cidx = {c: i for i, c in enumerate(cols)}
        mats = {}
        for var, mat in ent.items():
            re = [[Fr(0)] * len(cols) for _ in rows]; im = [[Fr(0)] * len(cols) for _ in rows]
            nz = False
            for (i, j), (cr, ci) in mat.items():
                if i in ridx and j in cidx:
                    re[ridx[i]][cidx[j]] += cr; im[ridx[i]][cidx[j]] += ci; nz = True
            if nz:
                mats[var] = GZ.from_frac(re, im)
        if not mats:
            self._cache[ck] = None
            return None
        zero = GZ.zeros(len(rows), len(cols))
        const = mats.get(None, zero)
        cq = []
        for m in range(4):
            t = zero
            if self.kps[m] in mats:
                t = t.add(mats[self.kps[m]])
            if self.ks[m] in mats:
                t = t.add(mats[self.ks[m]].scale_int(-1))
            cq.append(t)
        d = dict(const=const, cq=cq, ck=[mats.get(self.ks[m], zero) for m in range(4)],
                 ckp=[mats.get(self.kps[m], zero) for m in range(4)],
                 cpe=[[mats.get(pe[m], zero) for m in range(4)] for pe in self.pes])
        self._cache[ck] = d
        return d

    def block(self, a, b, key, rows, cols, qv, cur_in, cur_out, injs, lam_dual):
        """Dual GZ of the sub-block (rows, cols) at k = -(q - lam cur_in), kp = q - lam cur_out, pe_j = lam inj_j."""
        d = self._decomp(a, b, key, rows, cols)
        if d is None:
            return None
        A0 = d['const']
        for m in range(4):
            if qv[m]:
                A0 = A0.add(d['cq'][m].scale_int(int(qv[m])))
        B1 = None
        if lam_dual:
            B1 = GZ.zeros(len(rows), len(cols))
            for m in range(4):
                if cur_in[m]:
                    B1 = B1.add(d['ck'][m].scale_int(int(cur_in[m])))
                if cur_out[m]:
                    B1 = B1.add(d['ckp'][m].scale_int(-int(cur_out[m])))
                for j, inj in enumerate(injs):
                    if inj[m]:
                        B1 = B1.add(d['cpe'][j][m].scale_int(int(inj[m])))
        return Dual(A0, B1)

# ------------------------------------------------------------------ colour algebra (abstract, exact)
def colour_of_word(Ms):
    """Ms: list of noncommutative polynomials {tuple of leg indices: coeff} (the endomorphism M of each letter, in word order).
    colour = tr_g(M_n ... M_1).  Returns {tuple: coeff} of ordered traces T(tuple) = tr_g(ad_{c_t1} ... ad_{c_tk})."""
    acc = {(): Fr(1)}
    for M in Ms:                 # left-multiply: acc <- M acc  (so final = M_n ... M_1)
        new = {}
        for w1, c1 in M.items():
            for w2, c2 in acc.items():
                w = tuple(w1) + tuple(w2)
                new[w] = new.get(w, 0) + c1 * c2
        acc = new
    return {w: c for w, c in acc.items() if c != 0}

def perm_parity(seq):
    s = 1; seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s

def reduce_traces(T, legs):
    """T: {tuple: coeff}; legs: sorted list of the external leg indices.  Basis:
       2 legs: 'cadj k'  = c_adj (c_0,c_1)
       3 legs: 'cadj f'  = c_adj ([c_0,c_1],c_2)
       4 legs: 'cadj K12' = c_adj ([c_0,c_1],[c_2,c_3]),  'cadj K14' = c_adj ([c_0,c_3],[c_1,c_2]),
               'str4' = (1/24) sum_sigma tr_g(ad ad ad ad)  (symmetrised adjoint trace, not proportional to c_adj).
       c_adj := tr_g(ad_X ad_Y)/(X,Y)."""
    pos = {l: i for i, l in enumerate(legs)}
    out = {}
    def add(k, v):
        out[k] = out.get(k, 0) + v
    for w, c in T.items():
        assert sorted(w) == sorted(legs) and len(set(w)) == len(w), w
        s = [pos[x] for x in w]
        if len(s) == 2:
            add('cadj k', c)
        elif len(s) == 3:
            add('cadj f', c * Fr(1, 2) * perm_parity(s))
        elif len(s) == 4:
            i0 = s.index(0); r = s[i0:] + s[:i0]
            if r[1] > r[3]:
                r = [r[0], r[3], r[2], r[1]]
            t = tuple(r)
            if t == (0, 1, 2, 3):      # t1 = d + (1/6)(K12 - K14)
                add('str4', c); add('cadj K12', c * Fr(1, 6)); add('cadj K14', -c * Fr(1, 6))
            elif t == (0, 1, 3, 2):    # t2 = d - (1/3)K12 - (1/6)K14
                add('str4', c); add('cadj K12', -c * Fr(1, 3)); add('cadj K14', -c * Fr(1, 6))
            elif t == (0, 2, 1, 3):    # t3 = d + (1/6)K12 + (1/3)K14
                add('str4', c); add('cadj K12', c * Fr(1, 6)); add('cadj K14', c * Fr(1, 3))
            else:
                raise ValueError(t)
        else:
            raise ValueError(len(s))
    return {k: v for k, v in out.items() if v != 0}

def letter_M(kind, legs_of_letter, key):
    if kind == 'cubic':
        assert key == 'f'
        return {(legs_of_letter[0],): Fr(1)}
    i, j = legs_of_letter
    if key == 'K12':
        return {(i, j): Fr(1), (j, i): Fr(-1)}
    if key == 'K13':
        return {(j, i): Fr(-1)}
    if key == 'K14':
        return {(i, j): Fr(-1)}
    raise ValueError(key)

# ------------------------------------------------------------------ angular data: the 24 short roots of F4 (scaled by 2)
QPTS = [tuple(2 * s if i == m else 0 for i in range(4)) for m in range(4) for s in (1, -1)] + \
       [tuple(x) for x in itertools.product((1, -1), repeat=4)]
assert len(QPTS) == 24 and all(sum(x * x for x in v) == 4 for v in QPTS)

def check_design():
    """exact check: averages of monomials up to degree 5 over QPTS/2 equal the S^3 averages."""
    def sphere_avg(mon):
        if any(e % 2 for e in mon):
            return Fr(0)
        k = sum(mon) // 2
        num = 1
        for e in mon:
            num *= math.prod(range(e - 1, 0, -2)) if e > 0 else 1
        den = math.prod(range(4, 4 + 2 * k, 2)) if k > 0 else 1
        return Fr(num, den)
    for deg in range(1, 6):
        for mon in itertools.product(range(deg + 1), repeat=4):
            if sum(mon) != deg:
                continue
            avg = sum(Fr(math.prod(v[m] ** mon[m] for m in range(4)), 2 ** deg) for v in QPTS) / 24
            if avg != sphere_avg(mon):
                return False, mon
    return True, None

# ------------------------------------------------------------------ external sets, words, evaluation
class Leg:
    def __init__(self, idx, species, comps, u, p):
        self.idx = idx; self.species = species; self.comps = comps; self.u = u; self.p = p
    @property
    def odd(self):
        return PAR[self.species] == 1

def letters_partitions(legs):
    """all ways of grouping the legs into letters: singles (cubic) and pairs of A legs (quartic)."""
    idx = [l.idx for l in legs]
    Aidx = [l.idx for l in legs if l.species == 'A']
    out = [[('cubic', (i,)) for i in idx]]
    pairs = list(itertools.combinations(Aidx, 2))
    for pr in pairs:
        rest = [i for i in idx if i not in pr]
        out.append([('quartic', pr)] + [('cubic', (i,)) for i in rest])
    if len(idx) == 4 and len(Aidx) == 4:
        seen = set()
        for pr in pairs:
            rest = tuple(i for i in idx if i not in pr)
            key = frozenset([pr, rest])
            if key in seen:
                continue
            seen.add(key)
            out.append([('quartic', pr), ('quartic', rest)])
    return out

def vertex_label(kind, ext_species, row_e, col_e):
    """name of the vertex term for a letter with given external species and internal edge types (row/col edges)."""
    if kind == 'quartic':
        return 'AAAA'
    gh = row_e in ('cAs', 'Asc') or col_e in ('cAs', 'Asc')
    lam = 'lam' in (row_e, col_e)
    if ext_species == 'A':
        return 'AcA*' if gh else ('lambda lambda A' if lam else 'AAA')
    if ext_species == 'lam':
        return 'lambda lambda A' if ('AV' in (row_e, col_e) or 'Ap' in (row_e, col_e)) else 'c lambda lambda*'
    if ext_species == 'As':
        return 'AcA*'
    if ext_species == 'lams':
        return 'c lambda lambda*'
    if ext_species == 'cs':
        return 'ccc*'
    if ext_species == 'c':
        return 'AcA*' if gh else 'c lambda lambda*'
    raise ValueError(ext_species)

class ExternalSet:
    def __init__(self, name, legs, omega, int_species=('A', 'lam', 'c', 'As')):
        self.name = name; self.legs = legs; self.omega = omega; self.int_species = list(int_species)
        self._letters = {}

    def letter(self, kind, leg_ids):
        key = (kind, leg_ids)
        if key not in self._letters:
            ext = [(i, self.legs[i].species, self.legs[i].u) for i in leg_ids]
            self._letters[key] = Letter(kind, ext, self.int_species)
        return self._letters[key]

    def edge_types_for(self, species_row):
        return [e for e, d in EDGE.items() if d[0] == species_row]

    def evaluate(self, use_cyclic=False, verbose=False):
        """returns list of term records: dict(word, edges, keys, colour(basis), sign, n, value(avg), graph_key)"""
        recs = []
        lam_dual = self.omega >= 1
        for part in letters_partitions(self.legs):
            n = len(part)
            perms = list(itertools.permutations(range(n)))
            if use_cyclic:
                perms = [p for p in perms if p[0] == 0]
            for perm in perms:
                word = [part[i] for i in perm]
                # injections, cumulative momenta
                curs = [[0] * 4]
                for (kind, ids) in word:
                    inj = [sum(self.legs[i].p[m] for i in ids) for m in range(4)]
                    curs.append([curs[-1][m] + inj[m] for m in range(4)])
                assert all(x == 0 for x in curs[-1]), 'momentum not conserved'
                # r sign
                odd_seq = [i for (kind, ids) in word for i in ids if self.legs[i].odd]
                rsign = perm_parity(odd_seq)
                letters = [self.letter(kind, ids) for (kind, ids) in word]
                # edges: N_k before letter k (k = 0..n-1), letter k rows = col species of N_k, cols = row species of N_{k+1}
                for edges in itertools.product(list(EDGE.keys()), repeat=n):
                    ok = True
                    for kk in range(n):
                        a = EDGE[edges[kk]][2]; b = EDGE[edges[(kk + 1) % n]][0]
                        if (a, b) not in letters[kk].data:
                            ok = False; break
                    if not ok:
                        continue
                    keysets = [letters[kk].keys(EDGE[edges[kk]][2], EDGE[edges[(kk + 1) % n]][0]) for kk in range(n)]
                    for keys in itertools.product(*keysets):
                        val = self._word_value(word, letters, edges, keys, curs, lam_dual)
                        if val == (0, 0):
                            continue
                        Ms = [letter_M(word[kk][0], word[kk][1], keys[kk]) for kk in range(n)]
                        col = reduce_traces(colour_of_word(Ms), sorted(l.idx for l in self.legs))
                        if not col:
                            continue
                        recs.append(dict(word=word, edges=edges, keys=keys, colour=col, rsign=rsign, n=n,
                                         value=val, E=n, mult=(n if use_cyclic else 1)))
                        if verbose:
                            print(word, edges, keys, col, val, flush=True)
        return recs

    def _word_value(self, word, letters, edges, keys, curs, lam_dual):
        """returns (log-coefficient integrand average) = <[lambda^omega] Num(qhat,lambda P) prod |qhat - lambda P_f|^-2> as complex Fr pair,
        and also the eps^{-1/2}-type average for omega=1 (stored in self._last_lam0)."""
        n = len(word); E = n
        D = 2 * E - 4 + self.omega
        tot_re = Fr(0); tot_im = Fr(0)
        lam0_re = Fr(0); lam0_im = Fr(0)
        for qv in QPTS:
            mats = []
            for kk in range(n):
                et = edges[kk]
                mats.append(prop_block(et, qv, curs[kk], lam_dual))
                a = EDGE[et][2]; acomps = EDGE[et][3]
                nb = EDGE[edges[(kk + 1) % n]]
                b = nb[0]; bcomps = nb[1]
                injs = [self.legs[i].p for i in word[kk][1]]
                blk = letters[kk].block(a, b, keys[kk], acomps, bcomps, qv, curs[kk], curs[kk + 1], injs, lam_dual)
                if blk is None:
                    return (Fr(0), Fr(0))
                sg = [(-1) ** PAR[a]] * len(acomps)
                blk = Dual(blk.a.rowsign(sg), blk.b.rowsign(sg) if blk.b is not None else None)
                mats.append(blk)
            prod = reduce(lambda x, y: x.mm(y), mats)
            ssg = (-1) ** PAR[EDGE[edges[0]][0]]
            t0 = prod.a.trace()
            if self.omega == 0:
                tot_re += ssg * t0[0]; tot_im += ssg * t0[1]
            else:
                t1 = prod.b.trace() if prod.b is not None else (Fr(0), Fr(0))
                vdot = sum(sum(qv[m] * curs[f][m] for m in range(4)) for f in range(n))   # sum_f v.P_f, line f momentum q - P_f
                # 2^{-D}[2 N1(v) + N0(v) sum_f v.P_f]
                tot_re += ssg * (2 * t1[0] + t0[0] * vdot); tot_im += ssg * (2 * t1[1] + t0[1] * vdot)
                lam0_re += ssg * t0[0]; lam0_im += ssg * t0[1]
        scale = Fr(1, 24 * 2 ** D)
        self._last_lam0 = (lam0_re * Fr(1, 24 * 2 ** (D)), lam0_im * Fr(1, 24 * 2 ** (D)))
        if self.omega == 1 and (lam0_re != 0 or lam0_im != 0):
            raise AssertionError('eps^{-1/2} average nonzero')
        return (tot_re * scale, tot_im * scale)

def assemble(recs, legs):
    """log eps coefficient of the r-coefficient of W_1, colour-basis decomposed.
    Units: 'cadj ...' keys in u = hbar c_adj/(16 pi^2); 'str4' in hbar/(16 pi^2).
    log coefficient per word = (1/(2n)) * rsign * colour * (-avg)."""
    out = {}
    for r in recs:
        w = Fr(r['mult'], 2 * r['n']) * r['rsign']
        for k, c in r['colour'].items():
            re = -w * c * r['value'][0]; im = -w * c * r['value'][1]
            o = out.setdefault(k, [Fr(0), Fr(0)])
            o[0] += re; o[1] += im
    return {k: (v[0], v[1]) for k, v in out.items() if v != [0, 0]}

def graph_key(rec, legs):
    word = rec['word']; edges = rec['edges']; n = rec['n']
    seq = []
    for kk in range(n):
        kind, ids = word[kk]
        ext = [(legs[i].species, 'V' if legs[i].comps == list(range(4)) else ('perp' if legs[i].comps == list(range(4, 10)) else 'all')) for i in ids]
        vl = vertex_label(kind, legs[ids[0]].species, edges[kk], edges[(kk + 1) % n])
        ln = LINE_LABEL[edges[(kk + 1) % n]]
        seq.append((vl, tuple(ext), ln))
    # canonical cyclic sequence up to rotation and reflection (lines: element k = (vertex k, line after k))
    cands = []
    for rot in range(n):
        s = seq[rot:] + seq[:rot]
        cands.append(tuple(s))
        # reflection: vertices reversed, line after vertex k becomes line before
        verts = [x[:2] for x in s][::-1]
        lines = [x[2] for x in s]
        lines_ref = [lines[(n - 2 - j) % n] for j in range(n)]
        cands.append(tuple((verts[j][0], verts[j][1], lines_ref[j]) for j in range(n)))
    return min(cands)

def fmt_fr(x):
    re, im = x
    if im == 0:
        return str(re)
    if re == 0:
        return f'{im}*i'
    return f'{re} + {im}*i'
