"""Route B chains for the special-edge composite cycles Z_x(X_a,X_b;eps) (composite_theory.tex Def def:cc-sp).

Formula used (derived in notes.txt, manuscript conventions Delta_0 = -d_{K_0}, {f,g}_K = -(tree with edge K)):
    Z_x(X_a,X_b) = Q-image part = (hbar/2) d_{K_x}(chain trees with one edge K_eps - K_0)
                 = (hbar/4) sum_{vertex sequences} sum_{c != sigma} s(c,sigma) STr[E_0 Pi h_1 E_1 ... E_k Pi H_M],
E_c = composite kernel K_x, E_sigma = K_eps - K_0 = (e^{-eps k^2}-1) k0, other E_j = P(eps,oo) = N e^{-eps k^2}/k^2;
(hbar/2) d_K closes a tree with weight (hbar/4) STr per position (from W_1 = (hbar/2) sum (1/n) STr[(P Pi H)^n]);
s(c,sigma) = (-1)^{n(<sigma) + [c<sigma] + n(<c)}, n(<j) = number of odd outputs at the vertices before edge j
(Koszul sign of the two odd kernels: coefficient of delta_sigma delta_c in front).
Composite kernel in kernel form (route B):
    K_{P_b} = (d_b (x) 1 - 1 (x) d_b) X = 2 i (b.k) MX(k) e^{-2 eps k^2}/k^4      (two segments t', t'' >= eps)
    K_Delta = (v.(x)1 - 1(x)v.) X = [2 i k_mu (x^mu + y^mu) + 2(w_row - w_col)] MX(k) e^{-2eps k^2}/k^4,
x (y) the position of the vertex at the row (column) end; x^mu = -i d/dh_mu of a momentum h injected there and
absorbed by omega (canonical form, Def def:cc-canonical).  Option 1 (check): K_Delta = 2(v.(x)1)X + 4 eps d_eps X with
v. on the ghost end.
Ring variables: q (loop momentum), al (Schwinger alphas), outputs, letter components, derivative variables,
a0..a3 (direction of the line factor), h0..h3 (injected momentum), hm (marker of the constant line factor)."""
import os, sys, time, itertools, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
TL = os.path.join(HERE, 'vendor')     # vendored step2/twoletter/engine_b/tl_engine.py, tl_integrals.py
if TL not in sys.path:
    sys.path.insert(0, TL)
import ce_core as C
import tl_engine as T          # vertices/loops (Step-1 engine_b), colour key of M, verified twoletter conventions
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ_I
import ce_int as CI
import resource

V = C.V
INT6 = C.INT6
OFF = C.OFF
NT = C.NT
PARITY = C.PARITY
SLOT = C.SLOT
MEM_CAP_MB = 1400


def check_mem():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20
    if rss > MEM_CAP_MB:
        raise MemoryError(f'peak RSS {rss:.0f} MB exceeds the cap {MEM_CAP_MB} MB')


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
    pairs = [(a, b) for a in INT6 for b in INT6]
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


def mat_template(name):
    """entries of N (propagator numerator), MX (square kernel numerator) as term lists in _QS; k0 constant."""
    key = ('E', name)
    if key not in _TEMPL:
        if name == 'P':
            Mx = C.N_matrix(list(_QS))
        elif name == 'X':
            Mx = C.MX_matrix(list(_QS))
        elif name == 'S':
            Mx = C.k0_matrix()
        ent = []
        for i in range(NT):
            for j in range(NT):
                if Mx[i, j] != 0:
                    ent.append((i, j, _terms(Mx[i, j], _QS)))
        _TEMPL[key] = ent
    return _TEMPL[key]


subst = T.subst
QUARTER = QQ_I(Fr(1, 4), 0)   # (hbar/2) d_K closing a tree: (hbar/4) STr per position


class CJob:
    """One pair (slot_a, slot_b) with one output leg list; route B special-edge composite cycles.
    pol: 'sym' or list of dicts comp -> number (output projections); za, zb: None (symbolic) or dict;
    deriv: None, ('gen', N) or ('null', N, n, nbar) as in the twoletter engine."""

    def __init__(self, slot_a, slot_b, outs, pol='sym', za=None, zb=None, deriv=None, ntmax=None, xdir=None, nullw=True):
        self.sa, self.sb, self.outs = slot_a, slot_b, list(outs)
        self.nullw = nullw
        self.deriv = deriv
        n = len(outs)
        self.ntmax = ntmax if ntmax is not None else n + 3
        names = [f'q{m}' for m in range(4)] + [f'al{e}' for e in range(self.ntmax)]
        for j, o in enumerate(outs):
            spc, comps = SLOT[o]
            if pol == 'sym':
                names += [f'o{j}_{c}' for c in comps]
            if deriv is not None and deriv[0] == 'null':
                names += [f't{j}'] + ([f'w{j}'] if nullw else [])
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
        names += [f'a{m}' for m in range(4)] + [f'h{m}' for m in range(4)] + ['hm']
        self.R = ring(names, QQ_I)[0]
        R = self.R
        g = dict(zip(names, R.gens))
        self.g = g
        self.names = names
        self.q = [g[f'q{m}'] for m in range(4)]
        self.al = [g[f'al{e}'] for e in range(self.ntmax)]
        self.al_idx = [names.index(f'al{e}') for e in range(self.ntmax)]
        self.q_idx = [names.index(f'q{m}') for m in range(4)]
        self.h_idx = [names.index(f'h{m}') for m in range(4)]
        self.a_idx = [names.index(f'a{m}') for m in range(4)]
        self.hm_idx = names.index('hm')
        self.hv = [g[f'h{m}'] for m in range(4)]
        self.av = [g[f'a{m}'] for m in range(4)]
        self.opol, self.omom = [], []
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
                self.omom.append([g[f't{j}'] * R(QQ_I.from_sympy(sp.sympify(nb[m]))) + (g[f'w{j}'] if (m == 3 and nullw) else R.zero)
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
        self.odd_M = (V.PAR[spa] + V.PAR[spb]) % 2 == 1
        self.odd_out = [V.PAR[SLOT[o][0]] for o in outs]

    # ---------------------------------------------------------------- letter derivative factors (as tl_engine)
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
            t = t * x * R(QQ_I(Fr(1, k), 0))
            out += t
        return out

    def trunc(self, P):
        """drop terms of degree > 1 in h, > 1 in a, > 1 in hm and above the derivative order."""
        if not P:
            return P
        hi, ai, mi = self.h_idx, self.a_idx, self.hm_idx
        if self.deriv is not None:
            N = self.deriv[1]
            dn = ('s0', 's1', 's2', 's3', 'u0', 'u1', 'u2', 'u3') if self.deriv[0] == 'gen' else ('ss', 'su')
            di = [self.names.index(x) for x in dn]
        else:
            N = None; di = []
        keep = {}
        for m, c in P.items():
            if m[hi[0]] + m[hi[1]] + m[hi[2]] + m[hi[3]] > 1:
                continue
            if m[ai[0]] + m[ai[1]] + m[ai[2]] + m[ai[3]] > 1 or m[mi] > 1:
                continue
            if N is not None and sum(m[i] for i in di) > N:
                continue
            keep[m] = c
        return self.R.from_dict(keep) if len(keep) != len(P) else P

    # ---------------------------------------------------------------- matrices
    def _emat(self, name, mom, factor=None):
        ent = mat_template(name)
        R = self.R
        rows = {}
        for (i, j, terms) in ent:
            v = subst(terms, mom, R)
            if v:
                if factor is not None:
                    v = v * factor(i, j)
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
                    v = v * (-1 if PARITY[i] else 1)
                    rows.setdefault(i, []).append((j, v))
            if rows:
                out[key] = rows
        return out

    def _Mmat(self, k, kp, injM=None):
        """injM: 'row' / 'col': the position factor of the composite line sits at the letter on that leg of M; its
        derivative acts on it, i.e. the letter's derivative factor sees the injected momentum h (absorbed by omega)."""
        spa, _ = SLOT[self.sa]; spb, _ = SLOT[self.sb]
        syms, blocks = M_template(spa, spb)
        R = self.R
        kd = [k[m] + (self.hv[m] if injM == 'row' else R.zero) for m in range(4)]
        kpd = [kp[m] + (self.hv[m] if injM == 'col' else R.zero) for m in range(4)]
        out = {}
        for key, ent in blocks.items():
            if key == 'm01':
                fac = self.trunc(self.dfac(kd, 'a') * self.dfac(kpd, 'b'))
            else:
                fac = self.trunc(self.dfac(kpd, 'a') * self.dfac(kd, 'b'))
            images = list(k) + list(kp) + list(self.zav) + list(self.zbv)
            rows = {}
            for (i, j, terms) in ent:
                v = subst(terms, images, R)
                if v:
                    v = v * fac * (-1 if PARITY[i] else 1)
                    rows.setdefault(i, []).append((j, v))
            if rows:
                out[key] = rows
        return out

    def line_factor(self, mode):
        """factor(i, j) multiplying MX entries on the composite line.
        mode 'Pb': 2 i (a.k); 'D2': 2 i (a.k) + hm 2 (w_row - w_col) [Option 2]."""
        R = self.R
        k = self._kline

        def f(i, j):
            ak = sum((self.av[m] * k[m] for m in range(4)), R.zero)
            v = R(QQ_I(0, 2)) * ak
            wr = 1 if C.SPEC_OF[i] == 'A' else 0
            wc = 1 if C.SPEC_OF[j] == 'A' else 0
            if mode == 'D2':
                v = v + self.g['hm'] * R(QQ_I(2 * (wr - wc), 0))
            elif mode == 'D3':
                v = v + self.g['hm'] * R(QQ_I(2 * (wr - wc) - 2, 0))
            elif mode == 'X':
                v = R.one
            return v
        return f

    # ---------------------------------------------------------------- one sequence, one (c, sigma), one injection
    def config(self, seq, c, s, inj, mode, injM=None):
        """inj: vertex index (1..k) at which h is injected, or None.  Returns (types, curs, {word: Num}) with Num the
        STr numerator (polynomial in q), the Koszul sign and the (hbar/4) included, colour words as tl_engine."""
        R = self.R
        q = self.q
        k = len(seq)
        cur = [R.zero] * 4
        curs = []
        mats = []
        vinfo = []
        hsh = lambda j: (inj is not None and j >= inj)
        for b, blk in enumerate(seq):
            j = b                          # edge index before vertex b+1
            cj = [cur[m] + (self.hv[m] if hsh(j) else R.zero) for m in range(4)]
            curs.append(cj)
            mats.append(('E', j, cj))
            inj_m = [sum((self.omom[jj][m] for jj in blk), R.zero) for m in range(4)]
            nxt = [cur[m] + inj_m[m] for m in range(4)]
            cn = [nxt[m] + (self.hv[m] if hsh(j + 1) else R.zero) for m in range(4)]
            kk = [cj[m] - q[m] for m in range(4)]
            kq = [q[m] - cn[m] for m in range(4)]
            if len(blk) == 1:
                jj = blk[0]
                templ = cubic_template(SLOT[self.outs[jj]][0])
                ext = list(self.opol[jj]) + list(self.omom[jj])
            else:
                j1, j2 = blk
                templ = quartic_template()
                ext = list(self.opol[j1]) + list(self.omom[j1]) + list(self.opol[j2]) + list(self.omom[j2])
            mats.append(('H', self._hmat(templ, ext, kk, kq)))
            vinfo.append(blk)
            cur = nxt
        ck = [cur[m] + (self.hv[m] if hsh(k) else R.zero) for m in range(4)]
        curs.append(ck)
        mats.append(('E', k, ck))
        kk = [ck[m] - q[m] for m in range(4)]
        kq = list(q)
        Mm = self._Mmat(kk, kq, injM)
        # edge matrices
        built = []
        types = []
        for item in mats:
            if item[0] == 'H':
                built.append(('H', item[1]))
                continue
            _, j, cj = item
            mom = [q[m] - cj[m] for m in range(4)]
            if j == c:
                self._kline = mom
                built.append(('E', {'': self._emat('X', mom, self.line_factor(mode))}))
                types.append('K')
            elif j == s:
                built.append(('E', {'': self._emat('S', mom)}))
                types.append('S')
            else:
                built.append(('E', {'': self._emat('P', mom)}))
                types.append('P')
        # signs
        nodd = [0]
        for blk in seq:
            nodd.append(nodd[-1] + sum(self.odd_out[jj] for jj in blk))
        ksign = (-1) ** (nodd[s] + (1 if c < s else 0) + nodd[c])
        odd_seq = [jj for blk in seq for jj in blk if self.odd_out[jj]]
        sgn = V.perm_sign(odd_seq) * ksign
        result = {}
        keylists = [list(m[1].keys()) for m in built]
        allcols = sorted(set(j for Mrows in Mm.values() for rws in Mrows.values() for (j, v) in rws))
        for combo in itertools.product(*keylists):
            NumK = {mk: R.zero for mk in Mm}
            for i0 in allcols:
                vec = {i0: R.one}
                for (kind, km), key in zip(built, combo):
                    rows = km[key]
                    new = {}
                    for i, vi in vec.items():
                        for (j, mij) in rows.get(i, ()):
                            new[j] = new.get(j, R.zero) + vi * mij
                    vec = {j: self.trunc(v) for j, v in new.items() if v}
                    vec = {j: v for j, v in vec.items() if v}
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
                                NumK[mkey] += self.trunc(vi * (mij * ps))
            for mkey in Mm:
                Num = NumK[mkey]
                if not Num:
                    continue
                words = {(): 1}
                hk = [key for (kind, km), key in zip(built, combo) if kind == 'H']
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
                for w, cc in words.items():
                    if cc == 0:
                        continue
                    if mkey == 'm01':
                        w2 = tuple(reversed(w)); c2 = cc * (-1) ** len(w)
                    else:
                        w2 = w; c2 = cc
                    result[w2] = result.get(w2, R.zero) + Num * (c2 * sgn)
        result = {w: v for w, v in result.items() if v}
        return types, curs, result

    # ---------------------------------------------------------------- Gaussian + exact integrals
    def integrate(self, Num, curs, types, mlist=(0,), cut_override=None):
        """Num: polynomial in q (includes the (hbar/4) and signs but not the special-edge factor -|k_s|^2).
        Returns {m: {basis: ring element}} (coefficients of eps^{-m}), units hbar/(16 pi^2)."""
        R = self.R
        q = self.q
        # parameters
        params = []    # (edge index, special?, cutoff)
        for j, t in enumerate(types):
            if t == 'S':
                params.append((j, True, 1))
            elif t == 'K':
                if cut_override == 'Ksingle2':
                    params.append((j, False, 2))
                else:
                    params.append((j, False, 1)); params.append((j, False, 1))
            else:
                params.append((j, False, 1))
        nt = len(params)
        assert nt <= self.ntmax
        # order: special first
        sidx = [i for i, p in enumerate(params) if p[1]][0]
        order = [sidx] + [i for i in range(nt) if i != sidx]
        params = [params[i] for i in order]
        al = self.al[:nt]
        js = types.index('S')
        ks = [q[m] - curs[js][m] for m in range(4)]
        Num = Num * (-sum((x * x for x in ks), R.zero))
        cent = [sum((al[p] * curs[params[p][0]][m] for p in range(nt)), R.zero) for m in range(4)]
        Phi = sum((al[p] * sum((curs[params[p][0]][m] ** 2 for m in range(4)), R.zero) for p in range(nt)), R.zero) - \
            sum((cent[m] ** 2 for m in range(4)), R.zero)
        Phi = self.trunc(Phi)
        qi = self.q_idx
        groups = {}
        for mon, cf in Num.items():
            mu = tuple(mon[i] for i in qi)
            rest = list(mon)
            for i in qi:
                rest[i] = 0
            groups.setdefault(mu, {})[tuple(rest)] = cf
        subs = list(zip(q, cent))
        negPhi = [R.one]
        cuts = tuple(p[2] for p in params[1:])
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
                v = self.trunc(Dj.compose(subs)) if Dj else R.zero
                nums.append(v * R(QQ_I(Fr(1, math.factorial(j) * 4 ** j), 0)))
            for m in mlist:
                F = R.zero
                for j, nj in enumerate(nums):
                    i = j + 2 - nt - m
                    if i >= 0 and nj:
                        while len(negPhi) <= i:
                            negPhi.append(self.trunc(negPhi[-1] * (-Phi)))
                        F += self.trunc(nj * negPhi[i] * R(QQ_I(Fr(1, math.factorial(i)), 0)))
                if not F:
                    continue
                G = self._Japply(F, nt, m, cuts)
                for bkey, gv in G.items():
                    d = out.setdefault(m, {})
                    d[bkey] = d.get(bkey, R.zero) + self.trunc(gv * Cmu) * QUARTER
        return {m: {b: v for b, v in d.items() if v} for m, d in out.items()}

    def _Japply(self, F, nt, m, cuts):
        R = self.R
        ai = self.al_idx
        acc = {}
        for mon, cf in F.items():
            n = tuple(mon[ai[e]] for e in range(nt))
            rest = list(mon)
            for e in range(nt):
                rest[ai[e]] = 0
            rest = tuple(rest)
            val = CI.J(m, n[0], n[1:], cuts)
            for bkey, fr in val.items().items():
                d = acc.setdefault(bkey, {})
                d[rest] = d.get(rest, QQ_I(0, 0)) + cf * QQ_I(fr, 0)
        return {b: R.from_dict({k: v for k, v in d.items() if v}) for b, d in acc.items()}

    # ---------------------------------------------------------------- sequences
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

    # ---------------------------------------------------------------- the kernels
    def kernels(self, mlist=(0,), want=('Pb', 'D'), option=3, verbose=False, seqs=None):
        """Returns dict with
           'Pb' : {word: {m: {basis: poly}}}  C_{P_b} (m = 0) and poles along P_b, direction a (linear in a0..a3)
           'D0' : {word: {m: {basis: poly}}}  translation-invariant part along Delta (m = 0 only; poles via 'D0p')
           position part along Delta = 2 'Pb' (structural in the kernel form).
        option 3 (default): x + y - 2x_O = (x - y) + 2(y - x_O): line-local Euler term (-2 X + 8 eps |k|^2 X) and one
        injection at the column end (run B); option 2: injections at both ends (runs A and B), check of Lemma (ii)."""
        t0 = time.time()
        R = self.R
        res = {'Pb': {}, 'D0': {}}
        nconf = 0
        hi, ai, mi = self.h_idx, self.a_idx, self.hm_idx

        def add(tab, w, m, b, v):
            tab.setdefault(w, {}).setdefault(m, {})
            tab[w][m][b] = tab[w][m].get(b, R.zero) + v

        def split(P):
            parts = {}
            for mon, cf in P.items():
                da = sum(mon[i] for i in ai); dh = sum(mon[i] for i in hi); dm = mon[mi]
                parts.setdefault((da, dh, dm), {})[mon] = cf
            return {k: R.from_dict(v) for k, v in parts.items()}

        def contract(Pp):
            red = R.zero
            for mu in range(4):
                red += Pp.diff(self.av[mu]).diff(self.hv[mu])
            return red * R(QQ_I(0, -1))

        D = 'D' in want
        for seq in (seqs if seqs is not None else self.sequences()):
            k = len(seq)
            for c in range(k + 1):
                for s in range(k + 1):
                    if s == c:
                        continue
                    runs = []
                    if option == 3 or not D:
                        if D and c + 1 <= k:
                            runs.append(('B', c + 1, 'D3', True, None))
                        elif D:
                            runs.append(('B', None, 'D3', True, 'row'))
                        else:
                            runs.append(('B', None, 'Pb', True, None))
                        if D:
                            runs.append(('X', None, 'X', False, None))
                    else:
                        if c >= 1:
                            runs.append(('A', c, 'D2', True, None))
                        else:
                            runs.append(('A', None, 'D2', True, 'col'))
                        if c + 1 <= k:
                            runs.append(('B', c + 1, 'Pb', False, None))
                        else:
                            runs.append(('B', None, 'Pb', False, 'row'))
                    for rn, inj, mode, primary, injM in runs:
                        types, curs, words = self.config(seq, c, s, inj, mode, injM)
                        nconf += 1
                        for w, Num in words.items():
                            if rn == 'X':
                                co = self.integrate(Num, curs, types, (1,), cut_override='Ksingle2')
                                for b, P in co.get(1, {}).items():
                                    add(res['D0'], w, 0, b, P * R(QQ_I(8, 0)))
                                continue
                            co = self.integrate(Num, curs, types, mlist if primary else (0,))
                            for m, d in co.items():
                                for b, P in d.items():
                                    for (da, dh, dm), Pp in split(P).items():
                                        if primary and da == 1 and dh == 0 and dm == 0:
                                            add(res['Pb'], w, m, b, Pp)
                                        if D and m == 0 and da == 1 and dh == 1 and dm == 0:
                                            # option 3: 2 (y - x_O) from the run B injection; option 2: x - x_O, y - x_O
                                            fac = 2 if (option == 3 and rn == 'B') else 1
                                            add(res['D0'], w, m, b, contract(Pp) * R(QQ_I(fac, 0)))
                                        if D and m == 0 and primary and da == 0 and dh == 0 and dm == 1:
                                            add(res['D0'], w, m, b, Pp.diff(self.g['hm']))
            if verbose:
                print('   seq', seq, round(time.time() - t0, 1), 's', flush=True)
        for tab in ('Pb', 'D0'):
            res[tab] = {w: {m: {b: v for b, v in d.items() if v} for m, d in dm.items()} for w, dm in res[tab].items()}
            res[tab] = {w: {m: d for m, d in dm.items() if d} for w, dm in res[tab].items()}
            res[tab] = {w: dm for w, dm in res[tab].items() if dm}
        res['seconds'] = time.time() - t0
        res['configs'] = nconf
        return res
