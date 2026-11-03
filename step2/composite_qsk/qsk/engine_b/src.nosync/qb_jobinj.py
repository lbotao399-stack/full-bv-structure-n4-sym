"""Route B, intrinsic parts along eta and K_b (Prop qs-posdec, Thm qs-simplex (ii), (iii)).

A momentum h is injected at the vertex S^(1)(x) and carried by the edges after it to the vertex M, where it is
absorbed by omega (letters without derivatives here, n = 0).  The factor z^mu = (y - x_O)^mu of the Killing data is
-i d/dh_mu at h = 0, z^mu z^nu is -d^2/dh_mu dh_nu.
  eta:  Upsilon^0_eta = sum_mu (-i d/dh_mu) D[S^(1)(xi_mu)]|_{h=0} + D[s_0],  xi_mu = rho(e_mu) eta;
        position part along b: D[S^(1)(rho(b) eta)] (the h = 0 part, graph by graph: gate G-pos).
  K_b:  v(z) = z^2 b - 2 (b.z) z, omega_{rho sigma}(z) = 2(z_rho b_sigma - z_sigma b_rho), div v(z) = -8 (b.z):
        S^(1)(K_{b,x_O}) = sum_{mu nu} z^mu z^nu V_{mu nu} + sum_mu z^mu W_mu with
        V_{mu nu} f = (delta_{mu nu} b.d - b_mu d_nu - b_nu d_mu) f / 1  (symmetrised coefficient of z^mu z^nu in v^rho d_rho),
        W_mu f = 2 b_sigma Sigma^{mu sigma} f - 2 w b_mu f, on every field paired with its antifield with the signs of
        Def c7-vdot ((v.A,A*) - (v.lam,lam*) - (v.c,c*)), Sigma on A_V (vector), lam (spinor, (1/2) rho_mu rho_nu).
        Upsilon^0_K = -sum_{mu<=nu} c_{mu nu} d_{h_mu} d_{h_nu} D[V_{mu nu}] - i sum_mu d_{h_mu} D[W_mu]."""
import itertools
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ_I
import qb_core as B
import qb_job as QJ

V = B.V
I = sp.I
CJ = QJ.CJ


def register_K_terms(b, tag):
    """the vertex terms V_{mu nu} (mu <= nu, z^mu z^nu with multiplicity) and W_mu of S^(1)(K_{b,x_O})."""
    def spin(mu, nu):
        return sp.Rational(1, 2) * (B.Rbm[mu] * B.Rm[nu]) if mu != nu else sp.zeros(16, 16)
    names = {}
    for mu in range(4):
        for nu in range(mu, 4):
            # coefficient of z^mu z^nu (mu <= nu) in v^rho(z) d_rho: v^rho = z^2 b^rho - 2 (b.z) z^rho
            def dop(x, mu=mu, nu=nu):
                p = x['p']
                if mu == nu:
                    return sum((I * b[r] * p[r] for r in range(4)), 0) - 2 * b[mu] * I * p[mu]
                return -2 * (b[mu] * I * p[nu] + b[nu] * I * p[mu])
            def LA(s, dop=dop):
                a, As = s
                return dop(a) * sum((a['u'][M] * As['u'][M] for M in range(10)), 0)
            def Ll(s, dop=dop):
                l, ls = s
                return dop(l) * (sp.Matrix(l['u']).T * B.Mp * sp.Matrix(ls['u']))[0]
            def Lc(s, dop=dop):
                c, cs = s
                return dop(c) * c['u'][0] * cs['u'][0]
            nm = [f'KV{mu}{nu}A_{tag}', f'KV{mu}{nu}l_{tag}', f'KV{mu}{nu}c_{tag}']
            V.TERMS[nm[0]] = (sp.Integer(1), ('A', 'As'), LA, ('k', (0, 1)))
            V.TERMS[nm[1]] = (sp.Integer(-1), ('lam', 'lams'), Ll, ('k', (0, 1)))
            V.TERMS[nm[2]] = (sp.Integer(-1), ('c', 'cs'), Lc, ('k', (0, 1)))
            names[('V', mu, nu)] = nm
    for mu in range(4):
        # W_mu f = 2 b_sigma Sigma^{mu sigma} f - 2 w b_mu f
        def LA(s, mu=mu):
            a, As = s
            u = a['u']
            Su = [0] * 10
            for sg in range(4):
                if sg == mu:
                    continue
                # (Sigma^{mu sg} A)_m = delta_{m mu} A_sg - delta_{m sg} A_mu
                Su[mu] += 2 * b[sg] * u[sg]
                Su[sg] -= 2 * b[sg] * u[mu]
            tot = sum((Su[M] * As['u'][M] for M in range(10)), 0)
            tot += -2 * 1 * b[mu] * sum((u[M] * As['u'][M] for M in range(10)), 0)
            return tot
        def Ll(s, mu=mu):
            l, ls = s
            uv = sp.Matrix(l['u'])
            Sv = sp.zeros(16, 1)
            for sg in range(4):
                if sg != mu:
                    Sv += 2 * b[sg] * (spin(mu, sg) * uv)
            Sv += -2 * sp.Rational(3, 2) * b[mu] * uv
            return (Sv.T * B.Mp * sp.Matrix(ls['u']))[0]
        nm = [f'KW{mu}A_{tag}', f'KW{mu}l_{tag}']
        V.TERMS[nm[0]] = (sp.Integer(1), ('A', 'As'), LA, ('k', (0, 1)))
        V.TERMS[nm[1]] = (sp.Integer(-1), ('lam', 'lams'), Ll, ('k', (0, 1)))
        names[('W', mu)] = nm
    return names


class IJob(QJ.DJob):
    """DJob with a momentum h injected at the S vertex; gen 'xi' or 's0' as DJob, or 'terms' with an explicit list
    of quadratic vertex terms (even generator pieces, parity 1)."""

    def __init__(self, slot_a, slot_b, outs, gen, gvec, hmax=1, quad_terms=None, Spar=None, **kw):
        self.hmax = hmax
        if gen == 'terms':
            gsyms = sorted(set().union(*[sp.sympify(x).free_symbols for x in gvec]), key=lambda s: s.name)
            key = ('terms', tuple(quad_terms))
            if key not in QJ._STEMPL:
                t2 = CJ._block_template(quad_terms, [], [], extra=gsyms)
                QJ._STEMPL[key] = (t2, None, gsyms, Spar)
            self._key = key
            orig = QJ.s_templates
            QJ.s_templates = lambda g, v, t: QJ._STEMPL[key]
            try:
                super().__init__(slot_a, slot_b, outs, 'xi', gvec, **kw)
            finally:
                QJ.s_templates = orig
        else:
            super().__init__(slot_a, slot_b, outs, gen, gvec, **kw)
        # extend the ring by h0..h3
        names = self.names + [f'h{m}' for m in range(4)]
        Rn = ring(names, QQ_I)[0]
        old = self.R
        self.R = Rn
        g = dict(zip(names, Rn.gens))
        conv = lambda P: Rn.from_dict({tuple(list(m) + [0, 0, 0, 0]): c for m, c in P.items()}) if P else Rn.zero
        self.g, self.names = g, names
        self.q = [g[f'q{m}'] for m in range(4)]
        self.al = [g[f'al{e}'] for e in range(self.ntmax)]
        self.gimg = [g[nm] for nm in self.gnames]
        self.opol = [[conv(x) for x in v] for v in self.opol]
        self.omom = [[conv(x) for x in v] for v in self.omom]
        self.zav = [conv(x) for x in self.zav]; self.zbv = [conv(x) for x in self.zbv]
        self.hv = [g[f'h{m}'] for m in range(4)]
        self.h_idx = [names.index(f'h{m}') for m in range(4)]

    def trunc(self, P):
        if not P:
            return P
        hi = self.h_idx
        keep = {m: c for m, c in P.items() if sum(m[i] for i in hi) <= self.hmax}
        return self.R.from_dict(keep) if len(keep) != len(P) else P

    def config(self, seq, sigma):
        R = self.R
        q = self.q
        nv = len(seq)
        cur = [R.zero] * 4
        built = []
        vinfo = []
        spos = [i for i, (kd, bl) in enumerate(seq) if kd.startswith('S')][0]
        curs = []
        for bidx, (kind, blk) in enumerate(seq):
            curs.append(list(cur))
            built.append(('E', bidx, list(cur)))
            inj = [sum((self.omom[jj][m] for jj in blk), R.zero) for m in range(4)]
            if bidx == spos:
                inj = [inj[m] + self.hv[m] for m in range(4)]
            nxt = [cur[m] + inj[m] for m in range(4)]
            kk = [cur[m] - q[m] for m in range(4)]
            kq = [q[m] - nxt[m] for m in range(4)]
            if kind == 'I':
                if len(blk) == 1:
                    jj = blk[0]
                    templ = CJ.cubic_template(B.SLOT[self.outs[jj]][0])
                    ext = list(self.opol[jj]) + list(self.omom[jj])
                else:
                    j1, j2 = blk
                    templ = CJ.quartic_template()
                    ext = list(self.opol[j1]) + list(self.omom[j1]) + list(self.opol[j2]) + list(self.omom[j2])
                built.append(('H', self._hmat(templ, ext + kk + kq)))
            elif kind == 'S2':
                built.append(('H', self._hmat(self.t2, kk + kq + self.gimg)))
            else:
                jj = blk[0]
                ext = list(self.opol[jj]) + list(self.omom[jj])
                built.append(('H', self._hmat(self.t3[B.SLOT[self.outs[jj]][0]], ext + kk + kq + self.gimg)))
            vinfo.append((kind, blk))
            cur = nxt
        curs.append(list(cur))
        built.append(('E', nv, list(cur)))
        kk = [cur[m] - q[m] for m in range(4)]
        Mm = self._Mmat(kk, list(q))
        mats = []
        types = []
        for item in built:
            if item[0] == 'H':
                mats.append(item[1])
                continue
            _, j, cj = item
            mom = [q[m] - cj[m] for m in range(4)]
            if j == sigma:
                mats.append({'': self._emat('S', mom)}); types.append('S')
            else:
                mats.append({'': self._emat('P', mom)}); types.append('P')
        nodd = [0]
        for kind, blk in seq:
            nodd.append(nodd[-1] + sum(self.odd_out[jj] for jj in blk))
        ks = (-1) ** (nodd[sigma] + self.Spar * (nodd[spos] + (1 if spos < sigma else 0) + QJ.SIGN_EVEN_EXTRA))
        odd_seq = [jj for kind, blk in seq for jj in blk if self.odd_out[jj]]
        sgn = V.perm_sign(odd_seq) * ks
        fac = QQ_I(Fr(sgn, 2), 0)
        result = {}
        keylists = [list(m.keys()) for m in mats]
        allcols = sorted(set(j for Mrows in Mm.values() for rws in Mrows.values() for (j, v) in rws))
        isH = [x[0] == 'H' for x in built]
        for combo in itertools.product(*keylists):
            NumK = {mk: R.zero for mk in Mm}
            for i0 in allcols:
                vec = {i0: R.one}
                for km, key in zip(mats, combo):
                    rows = km[key]
                    new = {}
                    for i, vi in vec.items():
                        for (j, mij) in rows.get(i, ()):
                            new[j] = new.get(j, R.zero) + vi * mij
                    vec = {j: self.trunc(v) for j, v in new.items() if v}
                    vec = {j: v for j, v in vec.items() if v}
                    QJ.check_mem()
                    if not vec:
                        break
                if not vec:
                    continue
                QJ.check_mem()
                ps = -1 if B.PARITY[i0] else 1
                for mkey, Mrows in Mm.items():
                    for i, vi in vec.items():
                        for (j, mij) in Mrows.get(i, ()):
                            if j == i0:
                                NumK[mkey] += self.trunc(vi * (mij * ps))
            hkeys = [key for (km, key), h in zip(zip(mats, combo), isH) if h]
            for mkey in Mm:
                Num = NumK[mkey]
                if not Num:
                    continue
                words = {(): 1}
                for (kind, blk), key in zip(vinfo, hkeys):
                    if kind == 'S2':
                        W = {(): 1}
                    elif len(blk) == 1:
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
                    result[w2] = result.get(w2, R.zero) + Num * (QQ_I(c2, 0) * fac)
        return types, curs, {w: v for w, v in result.items() if v}

    def hcoeff(self, P, exps):
        """coefficient of prod h^exps (exps: tuple of 4) in P."""
        hi = self.h_idx
        d = {}
        for m, c in P.items():
            if tuple(m[i] for i in hi) == tuple(exps):
                mm = list(m)
                for i in hi:
                    mm[i] = 0
                d[tuple(mm)] = c
        return self.R.from_dict(d)
