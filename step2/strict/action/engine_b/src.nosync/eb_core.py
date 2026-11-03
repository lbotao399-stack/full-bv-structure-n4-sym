"""Engine B (route B), action side, strict Costello route: core numerics.

Reused verified code (read-only copies in vendor/, identical to step2/composite_qsk/composite/engine_b/src.nosync and its
vendor/): vertices.py (action terms on plane waves, Hessian blocks, Def c1-vertex convention: coefficient of the odd
leg parameters in leg order), realization.py (Clifford realisation), loops.py (propagator numerator N), ce_core.py
(54-component slot space, N(k), K_0), ce_int.py (exact Feynman-parameter integrals J(m; a; n; cut) with one special
Schwinger parameter in [0,1] and ordinary parameters >= 1, values in Q + sum_p Q log p).

Cycle formula (step2/engine_b ia_engine word formula, extended to letters of the cochain I + c_x S(x) + ...):
   value = sum_{letter partitions} sum_{words} sum_{pos} (1/(2n)) rsign (-1)^{odd before pos} colour(word)
           * STr[E_0 Pi h_0 E_1 Pi h_1 ... E_{n-1} Pi h_{n-1}],
E_pos = K_0 (e^{-eps k^2} - 1) (special edge, Schwinger time in [0, eps]), E_e = N(k_e) e^{-eps k_e^2}/k_e^2 otherwise,
edge e carries k_e = l - curs[e]; letter kk has row slot momentum curs[kk] - l and column slot momentum l - curs[kk+1].
Loop integral: Gaussian in l, Schwinger t = eps u (special u in [0,1], others u >= 1): coefficient of eps^{-m} =
(16 pi^2)^{-1} sum_k 2^{-k} C_k (-rho)^j / j! (j = k + 2 - n_t - m) integrated by ce_int.J, C_k = (Lap^k g)(B)/(2^k k!),
B = sum_e beta_e curs[e], rho = sum_e beta_e |curs_e|^2 - |B|^2, g = -|k_pos|^2 STr[...].
Numerics: the STr is evaluated on the simplex grid of integer loop momenta in complex double precision (block-sparse
batched products; every value is a dyadic rational far below 2^53, checked), converted exactly to Gaussian rationals and
interpolated exactly (forward differences); everything after the STr is exact (sympy ring over QQ_I, ce_int)."""
import os, sys, itertools, math, time, json, pickle, resource, threading, subprocess
from fractions import Fraction as Fr
from functools import lru_cache
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VEND = os.path.join(HERE, 'vendor')
if VEND not in sys.path:
    sys.path.insert(0, VEND)
import sympy as sp
import ce_core as C
import ce_int as CI
from sympy.polys.rings import ring
from sympy.polys.domains import QQ_I

V = C.V
NT = C.NT                      # 54
OFF = C.OFF
INT6 = C.INT6
PAR = C.PAR
PARITY = np.array(C.PARITY)
CACHE = os.path.join(os.path.dirname(HERE), 'out.nosync', 'cache')
os.makedirs(CACHE, exist_ok=True)

# ------------------------------------------------------------------ memory watchdog
MEM_LIMIT_MB = 1400


def start_watchdog(limit_mb=MEM_LIMIT_MB, period=0.5):
    pid = os.getpid()

    def run():
        while True:
            try:
                out = subprocess.run(['ps', '-o', 'rss=', '-p', str(pid)], capture_output=True, text=True).stdout.strip()
                if out and int(out) / 1024 > limit_mb:
                    os.write(2, f'WATCHDOG: RSS {int(out) / 1024:.0f} MB > {limit_mb} MB, exiting\n'.encode())
                    os._exit(7)
            except Exception:
                pass
            time.sleep(period)
    th = threading.Thread(target=run, daemon=True)
    th.start()


def peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20


# ------------------------------------------------------------------ block groups of the 54-dim slot space
GROUPS = [('AV', 0, 4), ('Ap', 4, 10), ('lamP', 10, 18), ('lamM', 18, 26), ('c', 26, 27),
          ('AsV', 27, 31), ('Asp', 31, 37), ('lamsM', 37, 45), ('lamsP', 45, 53), ('cs', 53, 54)]
GPAR = {'AV': 0, 'Ap': 0, 'lamP': 1, 'lamM': 1, 'c': 1, 'AsV': 1, 'Asp': 1, 'lamsM': 0, 'lamsP': 0, 'cs': 0}
assert OFF == {'A': 0, 'lam': 10, 'c': 26, 'As': 27, 'lams': 37, 'cs': 53}


def to_blocks(Fc, Fl=None, tol=0.0):
    """dense const (54,54) and l-coefficients (4,54,54) -> {(gi,gj): (const block, l-coef block)} nonzero only."""
    out = {}
    for gi, a0, a1 in GROUPS:
        for gj, b0, b1 in GROUPS:
            c = Fc[a0:a1, b0:b1]
            lcoef = Fl[:, a0:a1, b0:b1] if Fl is not None else None
            if np.any(c != 0) or (lcoef is not None and np.any(lcoef != 0)):
                out[(gi, gj)] = (c.copy(), None if lcoef is None or not np.any(lcoef != 0) else lcoef.copy())
    return out


def eval_blocks(bl, pts):
    """blocks at the loop-momentum points pts (npts,4) -> {(gi,gj): (npts,ni,nj)}."""
    out = {}
    for key, (c, lc) in bl.items():
        if lc is None:
            out[key] = np.broadcast_to(c, (pts.shape[0],) + c.shape)
        else:
            out[key] = c[None, :, :] + np.einsum('pm,mij->pij', pts, lc)
    return out


def bmul(A, B):
    out = {}
    Bby = {}
    for (gk, gj), X in B.items():
        Bby.setdefault(gk, []).append((gj, X))
    for (gi, gk), X in A.items():
        for gj, Y in Bby.get(gk, []):
            Z = np.matmul(X, Y)
            if (gi, gj) in out:
                out[(gi, gj)] = out[(gi, gj)] + Z
            else:
                out[(gi, gj)] = Z
    return out


def bstr(A, npts):
    tot = np.zeros(npts, dtype=complex)
    for (gi, gj), X in A.items():
        if gi == gj:
            tot += (-1) ** GPAR[gi] * np.trace(X, axis1=1, axis2=2)
    return tot


# ------------------------------------------------------------------ templates
KS = sp.symbols('hk0:4')
KPS = sp.symbols('hkp0:4')
PAIRS = [(a, b) for a in INT6 for b in INT6]


def _sym_leg(spc, name):
    return dict(sp=spc, u=list(sp.symbols(f'{name}u0:{V.DIM[spc]}')), p=list(sp.symbols(f'{name}p0:4')))


def build_template(terms, ext_species, extra_odd=0, gen_syms=None, name=None):
    """Hessian template of the functional sum(terms) with external legs of the given species (symbolic), row slot
    momentum KS, column slot momentum KPS.  extra_odd = 1 adds an odd parameter in front (CE coordinate of an even
    generator): factor (-1)^{p_J + p_K}.  gen_syms: symbols of the generator data (entries linear in them) -> leading
    index of the template.  Returns {colour key: array T[gen comp?, ext comps..., nmom, 54, 54]}, momentum slot index:
    0 constant, 1+4*i+m external leg i momentum m, then KS (4), KPS (4)."""
    if name is not None:
        fn = os.path.join(CACHE, f'tmpl_{name}.pkl')
        if os.path.exists(fn):
            with open(fn, 'rb') as f:
                return pickle.load(f)
    ext = [_sym_leg(s, f'e{i}') for i, s in enumerate(ext_species)]
    bl = V.hessian_blocks(terms, ext, list(KS), list(KPS), PAIRS)
    nmom = 1 + 4 * len(ext) + 8
    momsyms = [s for e in ext for s in e['p']] + list(KS) + list(KPS)
    dims = [V.DIM[s] for s in ext_species]
    gsy = list(gen_syms) if gen_syms else []
    lead = [len(gsy)] if gsy else []
    usyms = [s for leg in ext for s in leg['u']]
    gens = gsy + usyms + momsyms
    out = {}
    for (sJ, sK), d in bl.items():
        pJ, pK = PAR[sJ], PAR[sK]
        sgn = (-1) ** (extra_odd * (pJ + pK))
        for key, Mx in d.items():
            if Mx.is_zero_matrix:
                continue
            T = out.setdefault(key, np.zeros(tuple(lead + dims) + (nmom, NT, NT), dtype=complex))
            for i in range(Mx.rows):
                for j in range(Mx.cols):
                    e = Mx[i, j]
                    if e == 0:
                        continue
                    P = sp.Poly(sp.expand(e), *gens)
                    for mon, cf in P.terms():
                        idx = []
                        off = 0
                        if gsy:
                            hit = [c for c in range(len(gsy)) if mon[c] == 1]
                            assert len(hit) == 1 and sum(mon[:len(gsy)]) == 1, (mon, key)
                            idx.append(hit[0])
                            off = len(gsy)
                        for dd in dims:
                            hit = [c for c in range(dd) if mon[off + c] == 1]
                            assert len(hit) == 1 and sum(mon[off:off + dd]) == 1, (mon, key)
                            idx.append(hit[0])
                            off += dd
                        assert sum(mon[off:]) <= 1, (mon, key)
                        mm = [t for t in range(len(momsyms)) if mon[off + t]]
                        slot = 0 if not mm else 1 + mm[0]
                        T[tuple(idx) + (slot, OFF[sJ] + i, OFF[sK] + j)] += sgn * complex(cf)
    if name is not None:
        with open(os.path.join(CACHE, f'tmpl_{name}.pkl'), 'wb') as f:
            pickle.dump(out, f)
    return out


_EDGE = None


def edge_templates():
    global _EDGE
    if _EDGE is None:
        fn = os.path.join(CACHE, 'edges.pkl')
        if os.path.exists(fn):
            with open(fn, 'rb') as f:
                _EDGE = pickle.load(f)
        else:
            k = sp.symbols('ek0:4')
            N = C.N_matrix(k)
            N0 = np.array(N.subs({x: 0 for x in k}).tolist(), dtype=complex)
            Nk = np.zeros((4, NT, NT), dtype=complex)
            for m in range(4):
                Dm = N.applyfunc(lambda e: sp.diff(e, k[m]))
                assert all(sp.diff(e, x) == 0 for e in Dm for x in k)
                Nk[m] = np.array(Dm.tolist(), dtype=complex)
            K0 = np.array(C.k0_matrix().tolist(), dtype=complex)
            _EDGE = dict(N0=N0, Nk=Nk, K0=K0)
            with open(fn, 'wb') as f:
                pickle.dump(_EDGE, f)
    return _EDGE


def edge_factor(curs_e, special):
    """edge with momentum k = l - curs_e: P edge numerator N(k) (const, l-coef) or K_0 (special)."""
    E = edge_templates()
    if special:
        return E['K0'].copy(), None
    c = E['N0'] - np.einsum('m,mij->ij', np.asarray(curs_e, dtype=complex), E['Nk'])
    return c, E['Nk'].copy()


def letter_factor(T, ext_us, ext_ps, curs_in, curs_out, gen=None):
    """numeric Hessian (const, l-coef) of a letter: T template array for one colour key, ext_us/ext_ps lists of numeric
    polarisations / momenta (complex arrays), row momentum k = curs_in - l, column momentum kp = l - curs_out;
    multiplied by Pi on the row (parity of the row component)."""
    X = T
    if gen is not None:
        X = np.tensordot(np.asarray(gen, dtype=complex), X, axes=([0], [0]))
    for u in ext_us:
        X = np.tensordot(np.asarray(u, dtype=complex), X, axes=([0], [0]))
    # X: (nmom, 54, 54)
    nE = len(ext_us)
    mvec = np.zeros(X.shape[0], dtype=complex)
    mvec[0] = 1
    for i, p in enumerate(ext_ps):
        mvec[1 + 4 * i:1 + 4 * i + 4] = p
    o = 1 + 4 * nE
    mvec[o:o + 4] = curs_in
    mvec[o + 4:o + 8] = -np.asarray(curs_out, dtype=complex)
    c = np.tensordot(mvec, X, axes=([0], [0]))
    lc = -X[o:o + 4] + X[o + 4:o + 8]
    sg = np.array([(-1) ** p for p in PARITY], dtype=float)[:, None]
    c = sg * c
    lc = sg[None, :, :] * lc
    if not np.any(lc != 0):
        lc = None
    return c, lc


# ------------------------------------------------------------------ simplex grid and exact interpolation
@lru_cache(maxsize=None)
def simplex_grid(d):
    pts = [a for a in itertools.product(range(d + 1), repeat=4) if sum(a) <= d]
    return pts


@lru_cache(maxsize=None)
def binom_poly(a):
    """coefficients of binom(x, a) as Fractions, index = power of x."""
    poly = [Fr(1)]
    for t in range(a):
        # multiply by (x - t)
        new = [Fr(0)] * (len(poly) + 1)
        for i, cf in enumerate(poly):
            new[i + 1] += cf
            new[i] -= t * cf
        poly = new
    f = math.factorial(a)
    return tuple(cf / f for cf in poly)


SCALE_BITS = 30


def exact_values(vals):
    """complex doubles that are dyadic rationals -> pairs of exact ints scaled by 2^SCALE_BITS; checks exactness."""
    out = []
    S = 2 ** SCALE_BITS
    for z in vals:
        re, im = z.real * S, z.imag * S
        ire, iim = round(re), round(im)
        if abs(re - ire) > 1e-6 * max(1, abs(re)) * 0 + 1e-3 or abs(im - iim) > 1e-3:
            raise ValueError(f'non-dyadic value {z!r}')
        if max(abs(z.real), abs(z.imag)) > 2.0 ** 45:
            raise ValueError(f'value too large for exactness {z!r}')
        out.append((ire, iim))
    return out


def interpolate(vals, d):
    """exact polynomial in l (dict exponent tuple -> (Fraction re, Fraction im)) from values on simplex_grid(d)."""
    pts = simplex_grid(d)
    ev = exact_values(vals)
    g = {a: [v[0], v[1]] for a, v in zip(pts, ev)}
    for i in range(4):
        new = {}
        for a in pts:
            ai = a[i]
            accr = acci = 0
            for b in range(ai + 1):
                cb = math.comb(ai, b) * (-1) ** (ai - b)
                aa = a[:i] + (b,) + a[i + 1:]
                accr += cb * g[aa][0]
                acci += cb * g[aa][1]
            new[a] = [accr, acci]
        g = new
    S = 2 ** SCALE_BITS
    poly = {}
    for a, (cr, ci) in g.items():
        if cr == 0 and ci == 0:
            continue
        bps = [binom_poly(x) for x in a]
        for e in itertools.product(*[range(len(bp)) for bp in bps]):
            cf = bps[0][e[0]] * bps[1][e[1]] * bps[2][e[2]] * bps[3][e[3]]
            if cf == 0:
                continue
            cur = poly.setdefault(e, [Fr(0), Fr(0)])
            cur[0] += cf * Fr(cr, S)
            cur[1] += cf * Fr(ci, S)
    return {e: v for e, v in poly.items() if v[0] != 0 or v[1] != 0}


# ------------------------------------------------------------------ exact Schwinger / Feynman-parameter integration
_RINGS = {}


def get_ring(nt):
    if nt not in _RINGS:
        names = [f'l{m}' for m in range(4)] + [f'b{e}' for e in range(nt)]
        R, *gens = ring(names, QQ_I)
        _RINGS[nt] = (R, gens)
    return _RINGS[nt]


def to_ring(poly, R):
    d = {}
    nt = R.ngens - 4
    for e, (re, im) in poly.items():
        d[e + (0,) * nt] = QQ_I(re, im)
    return R.from_dict(d) if d else R.zero


class CLog:
    """complex number in Q + sum_p Q log p (real and imaginary LogNum)."""
    __slots__ = ('re', 'im')

    def __init__(self, re=None, im=None):
        self.re = re if re is not None else CI.LogNum(0)
        self.im = im if im is not None else CI.LogNum(0)

    def __add__(self, o):
        return CLog(self.re + o.re, self.im + o.im)

    def scale(self, z):
        """multiply by a Gaussian rational z = (re, im) Fractions."""
        a, b = Fr(z[0]), Fr(z[1])
        return CLog(self.re * a - self.im * b, self.re * b + self.im * a)

    def is_zero(self):
        return self.re.is_zero() and self.im.is_zero()

    def to_json(self):
        return {'re': {k: str(v) for k, v in self.re.items().items()}, 'im': {k: str(v) for k, v in self.im.items().items()}}

    def __repr__(self):
        return f'({self.re!r}) + i({self.im!r})'

    def __complex__(self):
        return complex(float(self.re), float(self.im))


def schwinger(poly_l, curs, pos, mlist=(0, 1)):
    """poly_l: exact STr polynomial in the loop momentum l (dict). curs: list of numeric (Fraction) 4-vectors (edge e
    momentum l - curs[e]); pos: special edge.  Returns {m: CLog} coefficient of eps^{-m} in units (16 pi^2)^{-1}."""
    nt = len(curs)
    R, gens = get_ring(nt)
    lv = gens[:4]
    bv = gens[4:]
    g = to_ring(poly_l, R)
    ks = [lv[m] - R(QQ_I(Fr(curs[pos][m]), 0)) for m in range(4)]
    g = -g * sum(x * x for x in ks)
    # Laplacian powers
    def lap(f):
        return sum((f.diff(lv[m]).diff(lv[m]) for m in range(4)), R.zero)
    B = [sum((bv[e] * R(QQ_I(Fr(curs[e][m]), 0)) for e in range(nt)), R.zero) for m in range(4)]
    c2 = [sum(Fr(x) ** 2 for x in curs[e]) for e in range(nt)]
    rho = sum((bv[e] * R(QQ_I(c2[e], 0)) for e in range(nt)), R.zero) - sum((x * x for x in B), R.zero)
    out = {m: CLog() for m in mlist}
    k = 0
    Lk = g
    order = [pos] + [e for e in range(nt) if e != pos]
    while Lk != 0:
        # C_k = Lap^k g (B) / (2^k k!)
        Ck = Lk.compose(list(zip(lv, B))) if Lk != 0 else R.zero
        Ck = Ck * R(QQ_I(Fr(1, 2 ** k * math.factorial(k)), 0))
        for m in mlist:
            j = k + 2 - nt - m
            if j < 0:
                continue
            term = Ck * (R.one if j == 0 else rho ** j) * R(QQ_I(Fr((-1) ** j, 2 ** k * math.factorial(j)), 0))
            for mon, cf in term.terms():
                assert all(x == 0 for x in mon[:4])
                al = [mon[4 + e] for e in order]
                Jv = CI.J(m, al[0], tuple(al[1:]), (1,) * (nt - 1))
                cre = Fr(int(cf.x.numerator), int(cf.x.denominator)); cim = Fr(int(cf.y.numerator), int(cf.y.denominator))
                out[m] = out[m] + CLog(Jv * cre, Jv * cim)
        Lk = lap(Lk)
        k += 1
    return out
