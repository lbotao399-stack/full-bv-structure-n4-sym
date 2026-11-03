"""Route B, PART Q: the defect cycles D_x(X_a,X_b;eps) of qsk_theory Def qs-defect in momentum space.

Formula (derived in notes.txt): with the flow sign of the special-edge vertex ({f,g}_eps - {f,g} = minus f and g joined
by K_eps - K_0) and W_1 = (hbar/2) sum_n (1/n) STr[(P Pi H)^n] (Step-1 engine B convention, reproduced Ch5),
    D_x(X_a,X_b) = -(hbar/2) sum_{ordered vertex sequences, S^(1)(x) at position s} sum_{sigma}
                    s_x(s,sigma) STr[E_0 Pi h_1 E_1 ... Pi h_n E_n Pi H_M],
E_sigma = K_eps - K_0 = (e^{-eps k^2} - 1) k0 (Schwinger time in [0, eps], numerator -|k|^2), the other E_j = P(eps, oo)
= N(k) e^{-eps k^2}/k^2; h_v the Hessians of I at the outputs and, at position s, the Hessian of the quadratic part of
S^(1)(x) (no output) or of its cubic part T5 (one output A); H_M the Hessian of M = int omega X_a X_b (twoletter
convention, ODD_M_SIGN = +1);
    s_x(s,sigma) = (-1)^{n(<sigma) + |S^(1)(x)| (n(<s) + [s<sigma])}   (Koszul parameters of the odd special edge and,
for even x, of the odd vertex S^(1)(x), moved to the front; n(<j) = number of odd outputs at the vertices before j);
times the sign of the reordering of the odd outputs from chain order to label order (twoletter convention).
Exact coefficients of eps^{-m} (m >= 1) and eps^0 by the Feynman-parameter integrals of composite engine B (ce_int.J,
one special parameter in [0,1], ordinary parameters >= 1), units hbar/(16 pi^2), colour words as tl_engine.
Ring generators: q (loop), al (alphas), output polarisations o<j>_<c> and momenta r<j>_<m>, letter components
za_<c>, zb_<c>, letter-derivative variables (s,u), generator parameters (xi<j> or b<m>)."""
import os, sys, time, itertools, math
from fractions import Fraction as Fr
import qb_core as B
import ce_job as CJ
import ce_int as CI
import sympy as sp
from sympy.polys.rings import ring
from sympy.polys.domains import QQ_I
import resource

C = B.C
V = B.V
NT = B.NT
OFF = B.OFF
PARITY = B.PARITY
SLOT = B.SLOT
MEM_CAP_MB = 1000
SIGN_EVEN_EXTRA = 0      # printed Koszul order (special-edge parameter first); with the overall + below (G-pole, G-reg)
SIGN_VERSION = 2         # v1 (before 05:10, 2026-10-30): overall -(hbar/2) and t-first for even x; v2 = -(v1) for odd x, = v1 for even x
T5_LAMS_OUTPUT = True    # revision 2026-10-30 (finding Q4): T5 = 1/2 (rho(1/2[A^A]) xi, lam*) also placed with its lam* leg as the
                         # output (cycle legs A, A); False reproduces the omission of the first production pass
subst = CJ.subst
_K = CJ._K
_KP = CJ._KP


def _watchdog(limit_mb=1300, period=0.3):
    """hard memory guard: polls the current RSS of this process (ps) and exits with code 7 above limit_mb."""
    import threading, subprocess, os as _os
    pid = _os.getpid()

    def run():
        while True:
            try:
                out = subprocess.run(['ps', '-o', 'rss=', '-p', str(pid)], capture_output=True, text=True).stdout.strip()
                if out and int(out) / 1024 > limit_mb:
                    _os.write(2, f'WATCHDOG: RSS {int(out) / 1024:.0f} MB > {limit_mb} MB, exiting\n'.encode())
                    _os._exit(7)
            except Exception:
                pass
            time.sleep(period)
    th = threading.Thread(target=run, daemon=True)
    th.start()


_watchdog()


def check_mem():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20
    if rss > MEM_CAP_MB:
        raise MemoryError(f'peak RSS {rss:.0f} MB exceeds the cap {MEM_CAP_MB} MB')


def peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20


_STEMPL = {}


def s_templates(gen, gvec, tag):
    """(quadratic template, cubic template or None, list of generator symbols, parity of S^(1)(x))."""
    key = (gen, tuple(map(str, gvec)), tag)
    if key in _STEMPL:
        return _STEMPL[key]
    gsyms = sorted(set().union(*[sp.sympify(x).free_symbols for x in gvec]), key=lambda s: s.name)
    if gen == 'xi':
        quad, cub = B.register_xi_terms(gvec, tag=tag)
        par = 0
    elif gen == 's0':
        quad = B.register_s0_term(gvec, tag=tag); cub = []
        par = 0
    elif gen == 'Pb':
        quad = B.register_translation_term(gvec, tag=tag); cub = []
        par = 1
    else:
        raise ValueError(gen)
    t2 = CJ._block_template(quad, [], [], extra=gsyms)
    t3 = None
    if cub:
        # T5 with its output on each of its slot species: A (cycle legs A, lam*) and lam* (cycle legs A, A)
        t3 = {}
        for spc in ('A', 'lams'):
            e = sp.symbols(f'e0_0:{V.DIM[spc]}'); r = sp.symbols('r0_0:4')
            ext = [dict(sp=spc, u=list(e), p=list(r))]
            t3[spc] = CJ._block_template(cub, ext, list(e) + list(r), extra=gsyms)
    _STEMPL[key] = (t2, t3, gsyms, par)
    return _STEMPL[key]


class DJob:
    """Defect cycles of one pair (slot_a, slot_b) with one output leg list along one generator.
    gen: 'xi' (gvec: 16 Sigma components, sympy symbols or numbers), 'Pb' (gvec: 4 components of b), 's0' (gvec: 16
    Sigma^* components; only the intrinsic term s_0 of S^(1)(eta)).
    pol: 'sym' or list of dicts comp -> number; za, zb: None (symbolic) or dict comp -> number;
    deriv: None or ('gen', N) (point splitting, letters at x+y_a, x+y_b as in tl_engine).
    terms: which parts of S^(1)(x) to use: 'quad', 'cubic' or 'all'."""

    def __init__(self, slot_a, slot_b, outs, gen, gvec, pol='sym', za=None, zb=None, deriv=None, ntmax=None,
                 terms='all', tag=None, zero_mom=False):
        self.sa, self.sb, self.outs = slot_a, slot_b, list(outs)
        self.gen, self.deriv, self.terms = gen, deriv, terms
        tag = tag or (gen + '_' + str(abs(hash(tuple(map(str, gvec)))) % 10 ** 8))
        self.t2, self.t3, gsyms, self.Spar = s_templates(gen, gvec, tag)
        n = len(outs)
        self.ntmax = ntmax if ntmax is not None else n + 3
        names = [f'q{m}' for m in range(4)] + [f'al{e}' for e in range(self.ntmax)]
        for j, o in enumerate(outs):
            spc, comps = SLOT[o]
            if pol == 'sym':
                names += [f'o{j}_{c}' for c in comps]
            if not zero_mom:
                names += [f'r{j}_{m}' for m in range(4)]
        spa, ca = SLOT[slot_a]; spb, cb = SLOT[slot_b]
        if za is None:
            names += [f'za_{c}' for c in ca]
        if zb is None:
            names += [f'zb_{c}' for c in cb]
        if deriv is not None:
            names += [f's{m}' for m in range(4)] + [f'u{m}' for m in range(4)]
        self.gnames = [s.name for s in gsyms]
        names += self.gnames
        self.R = ring(names, QQ_I)[0]
        R = self.R
        g = dict(zip(names, R.gens))
        self.g, self.names = g, names
        self.q = [g[f'q{m}'] for m in range(4)]
        self.al = [g[f'al{e}'] for e in range(self.ntmax)]
        self.al_idx = [names.index(f'al{e}') for e in range(self.ntmax)]
        self.q_idx = [names.index(f'q{m}') for m in range(4)]
        self.gimg = [g[nm] for nm in self.gnames]
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
            self.omom.append([R.zero] * 4 if zero_mom else [g[f'r{j}_{m}'] for m in range(4)])

        def zvec(z, slot, nm):
            spc, comps = SLOT[slot]
            vec = [R.zero] * V.DIM[spc]
            for c in comps:
                vec[c] = g[f'{nm}_{c}'] if z is None else R(QQ_I.from_sympy(sp.nsimplify(z.get(c, 0))))
            return vec
        self.zav = zvec(za, slot_a, 'za'); self.zbv = zvec(zb, slot_b, 'zb')
        self.odd_out = [V.PAR[SLOT[o][0]] for o in outs]
        if deriv is not None:
            di = [names.index(x) for x in ('s0', 's1', 's2', 's3', 'u0', 'u1', 'u2', 'u3')]
            self.d_idx = di

    # ------------------------------------------------------------ letter derivative factors (as tl_engine 'gen')
    def dfac(self, p, which):
        R = self.R
        if self.deriv is None:
            return R.one
        N = self.deriv[1]
        sv = [self.g[(f's{m}' if which == 'a' else f'u{m}')] for m in range(4)]
        x = R(QQ_I(0, 1)) * sum((sv[m] * p[m] for m in range(4)), R.zero)
        out = R.one; t = R.one
        for kk in range(1, N + 1):
            t = t * x * R(QQ_I(Fr(1, kk), 0))
            out += t
        return out

    def trunc(self, P):
        if self.deriv is None or not P:
            return P
        N = self.deriv[1]
        di = self.d_idx
        keep = {m: c for m, c in P.items() if sum(m[i] for i in di) <= N}
        return self.R.from_dict(keep) if len(keep) != len(P) else P

    # ------------------------------------------------------------ matrices
    def _emat(self, name, mom):
        ent = CJ.mat_template(name)
        R = self.R
        rows = {}
        for (i, j, terms) in ent:
            v = subst(terms, mom, R)
            if v:
                rows.setdefault(i, []).append((j, v))
        return rows

    def _hmat(self, templ, images):
        syms, blocks = templ
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

    def _Mmat(self, k, kp):
        spa, _ = SLOT[self.sa]; spb, _ = SLOT[self.sb]
        syms, blocks = CJ.M_template(spa, spb)
        R = self.R
        out = {}
        for key, ent in blocks.items():
            if key == 'm01':
                fac = self.trunc(self.dfac(k, 'a') * self.dfac(kp, 'b'))
            else:
                fac = self.trunc(self.dfac(kp, 'a') * self.dfac(k, 'b'))
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

    # ------------------------------------------------------------ sequences with one S vertex
    def iblock_sequences(self, labels):
        isA = {j: SLOT[self.outs[j]][0] == 'A' for j in labels}
        res = []

        def rec(remaining, cur):
            if not remaining:
                res.append(list(cur))
                return
            for j in remaining:
                rec([x for x in remaining if x != j], cur + [('I', (j,))])
            for a, b in itertools.combinations(remaining, 2):
                if isA[a] and isA[b]:
                    rec([x for x in remaining if x not in (a, b)], cur + [('I', (a, b))])
        rec(list(labels), [])
        return res

    def sequences(self):
        n = len(self.outs)
        out = []
        if self.terms in ('all', 'quad'):
            for seq in self.iblock_sequences(range(n)):
                for s in range(len(seq) + 1):
                    out.append(seq[:s] + [('S2', ())] + seq[s:])
        if self.terms in ('all', 'cubic') and self.t3 is not None:
            for j in range(n):
                spc = SLOT[self.outs[j]][0]
                if spc not in self.t3 or (spc == 'lams' and not T5_LAMS_OUTPUT):
                    continue
                rest = [x for x in range(n) if x != j]
                for seq in self.iblock_sequences(rest):
                    for s in range(len(seq) + 1):
                        out.append(seq[:s] + [('S3', (j,))] + seq[s:])
        return out

    # ------------------------------------------------------------ one configuration
    def config(self, seq, sigma):
        """returns (types, curs, {word: Num}) with Num the STr numerator (polynomial in q), sign s_x(s,sigma), the
        reordering sign, the factor -(1/2) included; the special-edge factor -|k_sigma|^2 is added in integrate."""
        R = self.R
        q = self.q
        nv = len(seq)
        cur = [R.zero] * 4
        curs = []
        built = []
        vinfo = []
        spos = None
        for bidx, (kind, blk) in enumerate(seq):
            curs.append(list(cur))
            built.append(('E', bidx, list(cur)))
            inj = [sum((self.omom[jj][m] for jj in blk), R.zero) for m in range(4)]
            nxt = [cur[m] + inj[m] for m in range(4)]
            kk = [cur[m] - q[m] for m in range(4)]
            kq = [q[m] - nxt[m] for m in range(4)]
            if kind == 'I':
                if len(blk) == 1:
                    jj = blk[0]
                    templ = CJ.cubic_template(SLOT[self.outs[jj]][0])
                    ext = list(self.opol[jj]) + list(self.omom[jj])
                else:
                    j1, j2 = blk
                    templ = CJ.quartic_template()
                    ext = list(self.opol[j1]) + list(self.omom[j1]) + list(self.opol[j2]) + list(self.omom[j2])
                built.append(('H', self._hmat(templ, ext + kk + kq)))
            elif kind == 'S2':
                spos = bidx
                built.append(('H', self._hmat(self.t2, kk + kq + self.gimg)))
            else:
                spos = bidx
                jj = blk[0]
                ext = list(self.opol[jj]) + list(self.omom[jj])
                built.append(('H', self._hmat(self.t3[SLOT[self.outs[jj]][0]], ext + kk + kq + self.gimg)))
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
                mats.append({'': self._emat('S', mom)})
                types.append('S')
            else:
                mats.append({'': self._emat('P', mom)})
                types.append('P')
        # signs
        nodd = [0]
        for kind, blk in seq:
            nodd.append(nodd[-1] + sum(self.odd_out[jj] for jj in blk))
        # the Koszul factor of Def qs-defect (ii) without its leading minus (see fac below)
        ks = (-1) ** (nodd[sigma] + self.Spar * (nodd[spos] + (1 if spos < sigma else 0) + SIGN_EVEN_EXTRA))
        odd_seq = [jj for kind, blk in seq for jj in blk if self.odd_out[jj]]
        sgn = V.perm_sign(odd_seq) * ks
        fac = QQ_I(Fr(sgn, 2), 0)        # D_x = +(hbar/2) sum s_x STr (route-B STr normalisation; fixed by G-pole and G-reg)
        result = {}
        keylists = [list(m.keys()) for m in mats]
        allcols = sorted(set(j for Mrows in Mm.values() for rws in Mrows.values() for (j, v) in rws))
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
                    check_mem()
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
                hk = [key for km, key in zip(mats, combo) if key != '']
                hinfo = [vi for vi in vinfo]
                # keys of H matrices in order (E matrices carry key '')
                hkeys = []
                it = iter(combo)
                for km, key in zip(mats, combo):
                    pass
                hkeys = [key for (km, key), isH in zip(zip(mats, combo), [x[0] == 'H' for x in built]) if isH]
                for (kind, blk), key in zip(hinfo, hkeys):
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
        result = {w: v for w, v in result.items() if v}
        return types, curs, result

    # ------------------------------------------------------------ Gaussian + exact integrals
    def integrate(self, Num, curs, types, mlist=(0, 1)):
        R = self.R
        q = self.q
        nt = len(types)
        assert nt <= self.ntmax
        js = types.index('S')
        order = [js] + [j for j in range(nt) if j != js]
        al = self.al[:nt]
        ks = [q[m] - curs[js][m] for m in range(4)]
        Num = Num * (-sum((x * x for x in ks), R.zero))
        cent = [sum((al[p] * curs[order[p]][m] for p in range(nt)), R.zero) for m in range(4)]
        Phi = sum((al[p] * sum((curs[order[p]][m] ** 2 for m in range(4)), R.zero) for p in range(nt)), R.zero) - \
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
        cuts = tuple([1] * (nt - 1))
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
                check_mem()
                G = self._Japply(F, nt, m, cuts)
                for bkey, gv in G.items():
                    d = out.setdefault(m, {})
                    d[bkey] = d.get(bkey, R.zero) + self.trunc(gv * Cmu)
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

    # ------------------------------------------------------------ the kernel
    def kernel(self, mlist=(0, 1), verbose=False, by_chain=False):
        t0 = time.time()
        R = self.R
        total = {}
        chains = []
        nconf = 0
        for seq in self.sequences():
            nv = len(seq)
            for sigma in range(nv + 1):
                types, curs, words = self.config(seq, sigma)
                nconf += 1
                if not words:
                    continue
                cres = {}
                for w, Num in words.items():
                    co = self.integrate(Num, curs, types, mlist)
                    for m, d in co.items():
                        for b, P in d.items():
                            tw = total.setdefault(w, {}).setdefault(m, {})
                            tw[b] = tw.get(b, R.zero) + P
                            if by_chain:
                                cw = cres.setdefault(w, {}).setdefault(m, {})
                                cw[b] = cw.get(b, R.zero) + P
                if by_chain and cres:
                    chains.append((seq, sigma, cres))
            if verbose:
                print('   seq', seq, round(time.time() - t0, 1), 's', flush=True)
        total = {w: {m: {b: v for b, v in d.items() if v} for m, d in dm.items()} for w, dm in total.items()}
        total = {w: {m: d for m, d in dm.items() if d} for w, dm in total.items()}
        total = {w: dm for w, dm in total.items() if dm}
        return dict(total=total, chains=chains, seconds=time.time() - t0, configs=nconf)


# ------------------------------------------------------------------ formatting (tl_engine names)
def comp_name(slot, i):
    return CJ.T.comp_name(slot, i)


def poly_str(job, p):
    if not p:
        return '0'
    ren = {}
    for nm in job.names:
        if nm.startswith('o') and '_' in nm:
            j, c = nm[1:].split('_')
            ren[nm] = f'Y{j}[{comp_name(job.outs[int(j)], int(c))}]'
        elif nm.startswith('za_'):
            ren[nm] = f'Xa[{comp_name(job.sa, int(nm[3:]))}]'
        elif nm.startswith('zb_'):
            ren[nm] = f'Xb[{comp_name(job.sb, int(nm[3:]))}]'
        elif nm.startswith('r') and '_' in nm:
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


def table_str(job, tab):
    out = {}
    for w, dm in tab.items():
        for m, d in dm.items():
            for b, P in d.items():
                out[f'word {list(w)} | eps^{-m} | {b}'] = poly_str(job, P)
    return out
