"""Route B engine for the two-letter kernels (twoletter_theory.tex, Def def:tl-dressed, Thm thm:tl-decomp).

Definition used (engine-B conventions, verified on Ch5 by engine_b/o_graphs.py):
    one loop  W_1 = (hbar/2) sum_n (1/n) STr[(P Pi H)^n],  H = Hessian of I at the outputs,
    delta-linear part with a two-letter vertex M = int omega (zeta_a.phi_a)(zeta_b.phi_b) (slot order a, b):
        ct(M) = (hbar/2) sum_{ordered vertex sequences} STr[P Pi h_1 P Pi h_2 ... P Pi H_M]  (singular part),
    i.e. K(X_a, X_b) = ct(int omega X_a X_b) (sigma_ab = 1 for an adjacent pair).  For odd M a parameter delta of
    degree -|M| is put in front.  The sign of the odd-M Hessian relative to vertices.hessian_blocks is fixed by the
    one-body Ward identity W1 (Prop prop:tl-W1): the naive left-derivative argument gives -1, but then the pairs
    (A_perp,c) and (A_V,c) give lc(gauge rotation) = -u (gauge rotation) and lc(int omega(zeta,[A_V,c])) with a one-field
    term +u int omega (zeta_m, d_m c), i.e. all three W1 statements reversed; with +1 all three hold (gates.json W1).
    ODD_M_SIGN = +1 is therefore used; it affects only the pairs with exactly one odd letter.
Momentum space: propagator P(q) = N(q) e^{-eps q^2}/q^2 (Feynman gauge, ghost sector on L), numerators N from
loops.prop_numerator (inversion of the free Hessian).  Edge j carries q - cur_j, cur_j = sum of the output momenta
injected before it.  Singular part: pinned-edge formula (tl_integrals.T) after the Gaussian loop integral,
    num_j(alpha) = (Lap_q^j N)(q = c(alpha)) / (j! 4^j),  c = sum_j alpha_j cur_j,
    Phi = sum alpha_j |cur_j|^2 - |c|^2,  m_j = j + 2 - E,
    [eps^-m] = (1/m) T(E,m,F_m),  [log eps] = -T(E,0,F_0),  F_m = sum_{m_j>=m} num_j (-Phi)^{m_j-m}/(m_j-m)!
all in units hbar/(16 pi^2) (the factor 1/2 of hbar/2 included).
Colour: every vertex block carries a colour key; cubic ([E,X],Y) -> W = ad_E; quartic K12 -> ad_E1 ad_E2 - ad_E2 ad_E1,
K13 -> -ad_E2 ad_E1, K14 -> -ad_E1 ad_E2 (E1,E2 the two outputs in template order).  The chain colour is
(1 (x) W_k...W_1)(-,-)^{-1} contracted with (H_M column leg, H_M row leg); the kernel is reported as
    K(X_a,X_b) = sum_w L_w  C_w,   C_w = (1 (x) ad_{Y_{w_1}} ... ad_{Y_{w_r}}) (-,-)^{-1},
X_a paired with the first factor, X_b with the second (Lemma lem:tl-colour (i)).
Odd outputs: the multilinear form has the odd parameters of the outputs in front in increasing label order."""
import os, sys, time, itertools, math
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
# vendored copy (step2/twoletter/engine_b/tl_engine.py, 2026-10-29): Step-1 modules vendored next to it
sys.path.insert(0, HERE)
_cwd = os.getcwd()
os.chdir(HERE)
import sympy as sp
import vertices as V
import loops as L
os.chdir(_cwd)
from sympy.polys.rings import ring
from sympy.polys.domains import QQ_I
import tl_integrals as TI
import resource

MEM_CAP_MB = 1400


def check_mem():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20   # bytes on macOS
    if rss > MEM_CAP_MB:
        raise MemoryError(f'peak RSS {rss:.0f} MB exceeds the cap {MEM_CAP_MB} MB')

# ---------------------------------------------------------------- colour key for the two-letter vertex M
_orig_ck = V.colour_key


def _ck(struct, sigma):
    kind, sl = struct
    if kind == 'm':
        legs = [sigma[x] for x in sl]
        return 'm%d%d' % tuple(legs), 1
    return _orig_ck(struct, sigma)


V.colour_key = _ck
ODD_M_SIGN = +1

SLOT = {'AV': ('A', list(range(4))), 'AW': ('A', list(range(4, 10))),
        'lamP': ('lam', list(range(8))), 'lamM': ('lam', list(range(8, 16))),
        'c': ('c', [0]), 'AsV': ('As', list(range(4))), 'AsW': ('As', list(range(4, 10))),
        'lamsM': ('lams', list(range(8))), 'lamsP': ('lams', list(range(8, 16))), 'cs': ('cs', [0])}
WEYL = {'c': Fraction(0), 'AV': Fraction(1), 'AW': Fraction(1), 'lamP': Fraction(3, 2), 'lamM': Fraction(3, 2),
        'AsV': Fraction(3), 'AsW': Fraction(3), 'lamsP': Fraction(5, 2), 'lamsM': Fraction(5, 2), 'cs': Fraction(4)}
DEG = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 0, 'lamM': 0, 'AsV': -1, 'AsW': -1, 'lamsP': -1, 'lamsM': -1, 'cs': -2}
INT = L.INT
OFF = L.OFF
NTOT = L.NTOT
PAR = V.PAR
SPEC_OF = []
for s in INT:
    SPEC_OF += [s] * V.DIM[s]
PARITY = [PAR[s] for s in SPEC_OF]


def comp_name(slot, i):
    """two-component name of component i of a slot (realisation of engine_b/realization.py)."""
    sp_, comps = SLOT[slot]
    if sp_ in ('A', 'As'):
        nm = 'A' if sp_ == 'A' else 'As'
        return f'{nm}_{i}' if i < 4 else f'{nm}_I{i - 3}'
    if sp_ == 'lam':
        if i < 8:
            return f'lam_a{i // 4}^{i % 4 + 1}'
        return f'lambar_ad{(i - 8) // 4}_{i % 4 + 1}'
    if sp_ == 'lams':
        if i < 8:
            return f'lambars^ad{i // 4}{i % 4 + 1}'
        return f'lams^a{(i - 8) // 4}_{i % 4 + 1}'
    return sp_


# ---------------------------------------------------------------- templates (sympy -> term lists)
_K = sp.symbols('kk0:4')
_KP = sp.symbols('kq0:4')
_QS = sp.symbols('Qs0:4')
_TEMPL = {}


def _terms(expr, syms):
    expr = sp.expand(expr)
    if expr == 0:
        return []
    P = sp.Poly(expr, *syms)
    return [(mon, QQ_I.from_sympy(c)) for mon, c in P.terms()]


def _block_template(terms, ext, extsyms, extra=()):
    pairs = [(a, b) for a in INT for b in INT]
    blocks = V.hessian_blocks(terms, ext, list(_K), list(_KP), pairs)
    syms = list(extsyms) + list(_K) + list(_KP) + list(extra)
    out = {}
    for (sJ, sK), bl in blocks.items():
        for key, M in bl.items():
            for i in range(M.rows):
                for j in range(M.cols):
                    if M[i, j] != 0:
                        out.setdefault(key, []).append((OFF[sJ] + i, OFF[sK] + j, _terms(M[i, j], syms)))
    return syms, out


def cubic_template(spc):
    key = ('cubic', spc)
    if key not in _TEMPL:
        e = sp.symbols(f'e0_0:{V.DIM[spc]}'); r = sp.symbols('r0_0:4')
        ext = [dict(sp=spc, u=list(e), p=list(r))]
        _TEMPL[key] = _block_template(V.CUBIC, ext, list(e) + list(r))
    return _TEMPL[key]


def quartic_template():
    key = ('quartic',)
    if key not in _TEMPL:
        e0 = sp.symbols('e0_0:10'); r0 = sp.symbols('r0_0:4'); e1 = sp.symbols('e1_0:10'); r1 = sp.symbols('r1_0:4')
        ext = [dict(sp='A', u=list(e0), p=list(r0)), dict(sp='A', u=list(e1), p=list(r1))]
        _TEMPL[key] = _block_template(V.QUARTIC, ext, list(e0) + list(r0) + list(e1) + list(r1))
    return _TEMPL[key]


def M_template(spa, spb):
    key = ('M', spa, spb)
    if key not in _TEMPL:
        za = sp.symbols(f'za_0:{V.DIM[spa]}'); zb = sp.symbols(f'zb_0:{V.DIM[spb]}')

        def Lf(s):
            return sum(za[i] * s[0]['u'][i] for i in range(V.DIM[spa])) * \
                sum(zb[j] * s[1]['u'][j] for j in range(V.DIM[spb]))
        V.TERMS['_M'] = (sp.Integer(1), (spa, spb), Lf, ('m', (0, 1)))
        try:
            _TEMPL[key] = _block_template(['_M'], [], [], list(za) + list(zb))
        finally:
            del V.TERMS['_M']
    return _TEMPL[key]


def prop_template():
    key = ('P',)
    if key not in _TEMPL:
        N = L.prop_numerator(list(_QS))
        ent = []
        for i in range(NTOT):
            for j in range(NTOT):
                if N[i, j] != 0:
                    ent.append((i, j, _terms(N[i, j], _QS)))
        _TEMPL[key] = (list(_QS), ent)
    return _TEMPL[key]


# ---------------------------------------------------------------- substitution of term lists into a ring
def subst(terms, images, R):
    """terms: [(exponent tuple over template symbols, coeff)]; images: list of ring elements (one per symbol)."""
    out = R.zero
    for mon, c in terms:
        t = R(c)
        for img, e in zip(images, mon):
            if e:
                if img == 0:
                    t = R.zero
                    break
                t = t * img ** e
        out += t
    return out


# ---------------------------------------------------------------- job
class Job:
    """One pair type with one output leg list.
    slot_a, slot_b: letter slot names; outs: list of output slot names (labels 0..n-1);
    pol: 'sym' (symbolic output polarisations) or a list of dicts comp -> number (fixed projections);
    za, zb: None (symbolic letter components) or dict comp -> number;
    deriv: None, ('gen', N) (point splitting s,u to order N), ('null', N, n, nbar) (letters (n.d)^a; outputs momenta
    r_j = t_j nbar + w_j e3)."""

    def __init__(self, slot_a, slot_b, outs, pol='sym', za=None, zb=None, deriv=None, Emax=None):
        self.sa, self.sb, self.outs = slot_a, slot_b, list(outs)
        self.deriv = deriv
        n = len(outs)
        self.Emax = Emax if Emax is not None else n + 1
        names = [f'q{m}' for m in range(4)] + [f'al{e}' for e in range(self.Emax)]
        self.pol_mode = pol
        for j, o in enumerate(outs):
            spc, comps = SLOT[o]
            if pol == 'sym':
                names += [f'o{j}_{c}' for c in comps]
            if deriv is not None and deriv[0] == 'null':
                names += [f't{j}', f'w{j}']
            else:
                names += [f'r{j}_{m}' for m in range(4)]
        spa, ca = SLOT[slot_a]; spb, cb = SLOT[slot_b]
        if za is None:
            names += [f'za_{c}' for c in ca]
        if zb is None:
            names += [f'zb_{c}' for c in cb]
        if deriv is not None:
            if deriv[0] == 'gen':
                names += [f's{m}' for m in range(4)] + [f'u{m}' for m in range(4)]
            else:
                names += ['ss', 'su']
        self.R = ring(names, QQ_I)[0]
        R = self.R
        g = dict(zip(names, R.gens))
        self.g = g
        self.q = [g[f'q{m}'] for m in range(4)]
        self.al = [g[f'al{e}'] for e in range(self.Emax)]
        self.al_idx = [names.index(f'al{e}') for e in range(self.Emax)]
        self.q_idx = [names.index(f'q{m}') for m in range(4)]
        self.names = names
        # output images
        self.opol = []
        self.omom = []
        for j, o in enumerate(outs):
            spc, comps = SLOT[o]
            vec = [R.zero] * V.DIM[spc]
            for c in comps:
                if pol == 'sym':
                    vec[c] = g[f'o{j}_{c}']
                else:
                    vec[c] = R(QQ_I.from_sympy(sp.nsimplify(pol[j].get(c, 0))))
            self.opol.append(vec)
            if deriv is not None and deriv[0] == 'null':
                nb = deriv[3]
                self.omom.append([g[f't{j}'] * R(QQ_I.from_sympy(sp.sympify(nb[m]))) + (g[f'w{j}'] if m == 3 else R.zero)
                                  for m in range(4)])
            else:
                self.omom.append([g[f'r{j}_{m}'] for m in range(4)])

        def zvec(z, slot, nm):
            spc, comps = SLOT[slot]
            vec = [R.zero] * V.DIM[spc]
            for c in comps:
                vec[c] = g[f'{nm}_{c}'] if z is None else R(QQ_I.from_sympy(sp.nsimplify(z.get(c, 0))))
            return vec
        self.zav = zvec(za, slot_a, 'za'); self.zbv = zvec(zb, slot_b, 'zb')
        self.odd_M = (PAR[spa] + PAR[spb]) % 2 == 1

    # derivative factors on the letters (momentum p incoming at the letter)
    def dfac(self, p, which):
        R = self.R
        if self.deriv is None:
            return R.one
        if self.deriv[0] == 'gen':
            N = self.deriv[1]
            sv = [self.g[(f's{m}' if which == 'a' else f'u{m}')] for m in range(4)]
            x = R(QQ_I(0, 1)) * sum((sv[m] * p[m] for m in range(4)), R.zero)
        else:
            N = self.deriv[1]; nv = self.deriv[2]
            x = R(QQ_I(0, 1)) * self.g['ss' if which == 'a' else 'su'] * \
                sum((R(QQ_I.from_sympy(sp.sympify(nv[m]))) * p[m] for m in range(4)), R.zero)
        out = R.one; t = R.one
        for k in range(1, N + 1):
            t = t * x * R(QQ_I(Fraction(1, k), 0))
            out += t
        return out

    def truncate(self, P):
        """drop the terms of total order > N in the point-splitting variables."""
        if self.deriv is None:
            return P
        N = self.deriv[1]
        names = ('s0', 's1', 's2', 's3', 'u0', 'u1', 'u2', 'u3') if self.deriv[0] == 'gen' else ('ss', 'su')
        idx = [self.names.index(n) for n in names]
        return self.R.from_dict({m: c for m, c in P.terms() if sum(m[i] for i in idx) <= N})

    def _pmat(self, mom):
        syms, ent = prop_template()
        R = self.R
        rows = {}
        for (i, j, terms) in ent:
            v = subst(terms, mom, R)
            if v:
                rows.setdefault(i, []).append((j, v))
        return rows

    def _hmat(self, templ, ext_images, k, kp):
        syms, blocks = templ
        images = list(ext_images) + list(k) + list(kp)
        R = self.R
        out = {}
        for key, ent in blocks.items():
            rows = {}
            for (i, j, terms) in ent:
                v = subst(terms, images, R)
                if v:
                    v = v * (-1 if PARITY[i] else 1)        # Pi H
                    rows.setdefault(i, []).append((j, v))
            if rows:
                out[key] = rows
        return out

    def _Mmat(self, k, kp):
        spa, _ = SLOT[self.sa]; spb, _ = SLOT[self.sb]
        syms, blocks = M_template(spa, spb)
        R = self.R
        da = {}; db = {}
        out = {}
        for key, ent in blocks.items():
            # key m01: slot a on the row leg (momentum k), m10: slot a on the column leg (momentum kp)
            if key == 'm01':
                fac = self.truncate(self.dfac(k, 'a') * self.dfac(kp, 'b'))
            else:
                fac = self.truncate(self.dfac(kp, 'a') * self.dfac(k, 'b'))
            images = [R.zero] * 8 + list(k) + list(kp)  # no ext symbols
            images = list(k) + list(kp) + list(self.zav) + list(self.zbv)
            rows = {}
            for (i, j, terms) in ent:
                v = subst(terms, images, R)
                if v:
                    v = v * fac * (-1 if PARITY[i] else 1)
                    if self.odd_M and ODD_M_SIGN == -1:
                        v = -v
                    rows.setdefault(i, []).append((j, v))
            if rows:
                out[key] = rows
        return out

    # ------------------------------------------------------------ one ordered sequence of output blocks
    def sequence(self, seq, split=False):
        """seq: list of tuples of output labels (size 1: cubic vertex; size 2: quartic, both outputs A).
        Returns (E, curs, {word: Num}) with Num the integrand numerator (polynomial in q)."""
        R = self.R
        q = self.q
        zero4 = [R.zero] * 4
        cur = list(zero4)
        mats = []      # list of (kind, {key: rows}) in product order: P, H, P, H, ..., P, M
        curs = []
        vinfo = []
        for blk in seq:
            curs.append(list(cur))
            mats.append(('P', {'': self._pmat([q[m] - cur[m] for m in range(4)])}))
            inj = [sum((self.omom[j][m] for j in blk), R.zero) for m in range(4)]
            nxt = [cur[m] + inj[m] for m in range(4)]
            kk = [cur[m] - q[m] for m in range(4)]
            kq = [q[m] - nxt[m] for m in range(4)]
            if len(blk) == 1:
                j = blk[0]
                spc = SLOT[self.outs[j]][0]
                templ = cubic_template(spc)
                ext = list(self.opol[j]) + list(self.omom[j])
            else:
                j1, j2 = blk
                templ = quartic_template()
                ext = list(self.opol[j1]) + list(self.omom[j1]) + list(self.opol[j2]) + list(self.omom[j2])
            mats.append(('H', self._hmat(templ, ext, kk, kq)))
            vinfo.append(blk)
            cur = nxt
        curs.append(list(cur))
        mats.append(('P', {'': self._pmat([q[m] - cur[m] for m in range(4)])}))
        kk = [cur[m] - q[m] for m in range(4)]
        kq = list(q)
        Mm = self._Mmat(kk, kq)
        E = len(seq) + 1
        # odd-output reordering sign (sequence order -> label order)
        odd_seq = [j for blk in seq for j in blk if PAR[SLOT[self.outs[j]][0]]]
        sgn = V.perm_sign(odd_seq)
        # column support of M
        result = {}
        pieces = {}
        sing_cache = {}
        keylists = [list(m[1].keys()) for m in mats]
        for combo in itertools.product(*keylists):
            allcols = sorted(set(j for Mrows in Mm.values() for rws in Mrows.values() for (j, v) in rws))
            NumK = {mk: R.zero for mk in Mm}
            plistK = {mk: [] for mk in Mm}
            for i0 in allcols:              # one start row at a time (bounded memory)
                vec = {i0: R.one}
                for (kind, km), key in zip(mats, combo):
                    rows = km[key]
                    new = {}
                    for i, vi in vec.items():
                        for (j, mij) in rows.get(i, ()):
                            new[j] = new.get(j, R.zero) + vi * mij
                    vec = {j: v for j, v in new.items() if v}
                    if not vec:
                        break
                if not vec:
                    continue
                check_mem()
                ps = -1 if PARITY[i0] else 1
                for mkey, Mrows in Mm.items():
                    for i, vi in vec.items():
                        for (j, mij) in Mrows.get(i, ()):
                            if j == i0:
                                if split and self.deriv is not None:
                                    # the letter factor depends on q: reduce each piece separately (bounded memory)
                                    plistK[mkey].append((R.one, self.singular(vi * (mij * ps), curs, E)))
                                elif split:
                                    ck = (i0, i)
                                    if ck not in sing_cache:
                                        sing_cache[ck] = self.singular(vi, curs, E)
                                    plistK[mkey].append((mij * ps, sing_cache[ck]))
                                else:
                                    NumK[mkey] += vi * (mij * ps)
                sing_cache.clear()
                del vec
            for mkey in Mm:
                Num = NumK[mkey]; plist = plistK[mkey]
                if not Num and not plist:
                    continue
                # colour word
                words = {(): 1}
                hk = [key for (kind, km), key in zip(mats, combo) if kind == 'H']
                for blk, key in zip(vinfo, hk):
                    if len(blk) == 1:
                        W = {(blk[0],): 1}
                    else:
                        j1, j2 = blk
                        W = {'K12': {(j1, j2): 1, (j2, j1): -1}, 'K13': {(j2, j1): -1}, 'K14': {(j1, j2): -1}}[key]
                    nw = {}
                    for w1, c1 in W.items():
                        for w0, c0 in words.items():
                            nw[w1 + w0] = nw.get(w1 + w0, 0) + c1 * c0
                    words = nw
                for w, c in words.items():
                    if c == 0:
                        continue
                    if mkey == 'm01':
                        w2 = tuple(reversed(w)); c2 = c * (-1) ** len(w)
                    else:
                        w2 = w; c2 = c
                    if split:
                        pieces.setdefault(w2, []).extend((m_ * (c2 * sgn), v_) for m_, v_ in plist)
                    else:
                        result[w2] = result.get(w2, R.zero) + Num * (c2 * sgn)
        if split:
            return E, curs, pieces
        result = {w: v for w, v in result.items() if v}
        return E, curs, result

    # ------------------------------------------------------------ singular part
    def singular(self, Num, curs, E):
        """Returns {m: poly (m>=1), 'log': poly}, polynomials free of q and alpha, units hbar/(16 pi^2).
        Num = sum_mu q^mu C_mu; the Gaussian/pinned-edge reduction is done on each loop monomial q^mu (small
        polynomials in alpha and the output momenta) and multiplied by C_mu at the end."""
        R = self.R
        q = self.q
        al = self.al[:E]
        c = [sum((al[j] * curs[j][m] for j in range(E)), R.zero) for m in range(4)]
        Phi = sum((al[j] * sum((curs[j][m] ** 2 for m in range(4)), R.zero) for j in range(E)), R.zero) - \
            sum((c[m] ** 2 for m in range(4)), R.zero)
        qi = self.q_idx
        groups = {}
        for mon, cf in Num.terms():
            mu = tuple(mon[i] for i in qi)
            rest = list(mon)
            for i in qi:
                rest[i] = 0
            groups.setdefault(mu, {})[tuple(rest)] = cf
        subs = list(zip(q, c))
        negPhi_pows = [R.one]
        out = {}
        for mu, cdict in groups.items():
            Cmu = R.from_dict(cdict)
            qm = R.one
            for m_, e_ in zip(q, mu):
                if e_:
                    qm = qm * m_ ** e_
            laps = [qm]
            for j in range(sum(mu) // 2):
                Dj = laps[-1]
                laps.append(sum((Dj.diff(q[m]).diff(q[m]) for m in range(4)), R.zero))
            nums = []
            for j, Dj in enumerate(laps):
                v = Dj.compose(subs) if Dj else R.zero
                nums.append(v * R(QQ_I(Fraction(1, math.factorial(j) * 4 ** j), 0)))
            mmax = max(len(nums) - 1 + 2 - E, 0)
            while len(negPhi_pows) <= mmax + 1:
                negPhi_pows.append(negPhi_pows[-1] * (-Phi))
            for m in range(0, mmax + 1):
                F = R.zero
                for j, nj in enumerate(nums):
                    mj = j + 2 - E
                    if mj >= m and nj:
                        F += nj * negPhi_pows[mj - m] * R(QQ_I(Fraction(1, math.factorial(mj - m)), 0))
                if not F:
                    continue
                G = self._Tapply(F, E, m)
                if not G:
                    continue
                if m == 0:
                    G = -G
                else:
                    G = G * R(QQ_I(Fraction(1, m), 0))
                key = 'log' if m == 0 else m
                out[key] = out.get(key, R.zero) + G * Cmu
        half = R(QQ_I(Fraction(1, 2), 0))
        out = {k: v * half for k, v in out.items() if v}
        return out

    def _Tapply(self, F, E, m):
        R = self.R
        acc = {}
        ai = self.al_idx
        for mon, cf in F.terms():
            n = tuple(mon[ai[e]] for e in range(E))
            assert all(mon[ai[e]] == 0 for e in range(E, self.Emax))
            rest = list(mon)
            for e in range(E):
                rest[ai[e]] = 0
            rest = tuple(rest)
            t = TI.T(E, m, n)
            if t:
                acc[rest] = acc.get(rest, QQ_I(0, 0)) + cf * QQ_I(t, 0)
        return R.from_dict({k: v for k, v in acc.items() if v})

    # ------------------------------------------------------------ all sequences for the leg list
    def sequences(self):
        n = len(self.outs)
        labels = list(range(n))
        isA = [SLOT[o][0] == 'A' for o in self.outs]
        res = []

        def rec(remaining, cur):
            if not remaining:
                res.append(list(cur))
                return
            for j in remaining:
                rec([x for x in remaining if x != j], cur + [(j,)])
            for a, b in itertools.combinations(remaining, 2):
                if isA[a] and isA[b]:
                    rec([x for x in remaining if x not in (a, b)], cur + [(a, b)])
        rec(labels, [])
        return res

    def kernel(self, keep_sequences=True, verbose=False):
        """Sum over all ordered sequences; returns dict with per-sequence and total singular coefficients:
        {'total': {word: {m/log: poly}}, 'by_sequence': [(seq, {word: {...}})], 'seconds': ...}."""
        t0 = time.time()
        total = {}
        byseq = []
        for seq in self.sequences():
            split = True
            E, curs, words = self.sequence(seq, split=split)
            sres = {}
            cache = {}
            for w, Num in words.items():
                if split:
                    s = {}
                    for (mfac, sv) in Num:          # mfac: q-independent letter factor (no derivatives)
                        for k_, v_ in sv.items():
                            s[k_] = s.get(k_, self.R.zero) + v_ * mfac
                    s = {k_: v_ for k_, v_ in s.items() if v_}
                else:
                    s = self.singular(Num, curs, E)
                if s:
                    sres[w] = s
                    tw = total.setdefault(w, {})
                    for k, v in s.items():
                        tw[k] = tw.get(k, self.R.zero) + v
            if keep_sequences:
                byseq.append((seq, sres))
            if verbose:
                print('   seq', seq, 'words', len(words), 'nonzero', len(sres), round(time.time() - t0, 1), flush=True)
        total = {w: {k: v for k, v in d.items() if v} for w, d in total.items()}
        total = {w: d for w, d in total.items() if d}
        return dict(total=total, by_sequence=byseq, seconds=time.time() - t0)


# ---------------------------------------------------------------- formatting
def poly_str(job, p):
    """sympy-readable string; variables renamed to two-component letter/output names."""
    if not p:
        return '0'
    ren = {}
    for nm in job.names:
        if nm.startswith('o'):
            j, c = nm[1:].split('_')
            ren[nm] = f'Y{j}[{comp_name(job.outs[int(j)], int(c))}]'
        elif nm.startswith('za_'):
            ren[nm] = f'Xa[{comp_name(job.sa, int(nm[3:]))}]'
        elif nm.startswith('zb_'):
            ren[nm] = f'Xb[{comp_name(job.sb, int(nm[3:]))}]'
        elif nm.startswith('r'):
            j, m = nm[1:].split('_')
            ren[nm] = f'p{j}_{m}'
    terms = []
    for mon, cf in sorted(p.terms()):
        fac = []
        for nm, e in zip(job.names, mon):
            if e:
                s = ren.get(nm, nm)
                fac.append(s if e == 1 else f'{s}**{e}')
        a, b = cf.x, cf.y
        if b == 0:
            cs = f'({a})'
        elif a == 0:
            cs = f'({b})*I'
        else:
            cs = f'({a}+({b})*I)'
        terms.append(cs + ('*' + '*'.join(fac) if fac else ''))
    return ' + '.join(terms)
