"""Engine B, action side: the terms of S^(1)(x) as plane-wave vertex functions (conventions of vendor/vertices.py,
Def c1-vertex), written from the theory (action_theory.tex eq at-S1list, Def at-gens; manuscript Def c7-vdot):
  even x:  S^(1)(x) = int ((v.A, A*) - (v.lam, lam*) - (v.c, c*)),
           v.phi = v^mu d_mu phi + (w/4)(div v) phi + spin part;  weights w(A) = 1, w(lam) = 3/2, w(c) = 0;
           Delta: v = 2x, S(Delta) = sum_mu x^mu Vdir(e_mu) + V0, Vdir(a) = 2[(a.dA,A*) - (a.d lam,lam*) - (a.d c,c*)],
                  V0 = 2 (A,A*) - 3 (lam,lam*);  P_b: v = b, S(P_b) = Vdir(b)/2.
  xi in Sigma: T1+T2 = -(Gamma(xi,lam),A*), T3+T4 = (1/2)(rho(dA)xi, lam*), T5 = (1/2)(rho((1/2)[A^A])xi, lam*),
           rho(F) = (1/2) F_MN rho_M rho_N; T5 = (1/4) sum_{M != N} ([A_M,A_N], (rho_M rho_N xi, lam*)) as in qb_core.
Generator data enter linearly (template leading index)."""
import sympy as sp
import numpy as np
import eb_core as E
import ce_core as C

V = E.V
I = sp.I
Mp = C.Mp
Rm = C.Rm
Rbm = C.Rbm
MR = [Mp * Rm[M] for M in range(10)]


def _bd(b, x):
    return sum((I * b[m] * x['p'][m] for m in range(4)), 0)


def reg_translation(b, tag, scale=1):
    """scale * [(b.dA, A*) - (b.d lam, lam*) - (b.d c, c*)]."""
    def L_A(s):
        a, As = s
        return _bd(b, a) * sum((a['u'][M] * As['u'][M] for M in range(10)), 0)

    def L_l(s):
        l, ls = s
        return _bd(b, l) * (sp.Matrix(l['u']).T * Mp * sp.Matrix(ls['u']))[0]

    def L_c(s):
        c, cs = s
        return _bd(b, c) * c['u'][0] * cs['u'][0]
    V.TERMS['TA_' + tag] = (sp.Integer(scale), ('A', 'As'), L_A, ('k', (0, 1)))
    V.TERMS['Tl_' + tag] = (sp.Integer(-scale), ('lam', 'lams'), L_l, ('k', (0, 1)))
    V.TERMS['Tc_' + tag] = (sp.Integer(-scale), ('c', 'cs'), L_c, ('k', (0, 1)))
    return ['TA_' + tag, 'Tl_' + tag, 'Tc_' + tag]


def reg_weights(wA, wl, wc, tag, sel=None):
    """wA (A,A*) - wl (lam,lam*) - wc (c,c*); sel restricts A to 'V' or 'p' components."""
    def L_A(s):
        a, As = s
        rng = range(10) if sel is None else (range(4) if sel == 'V' else range(4, 10))
        return sum((a['u'][M] * As['u'][M] for M in rng), 0)

    def L_l(s):
        l, ls = s
        return (sp.Matrix(l['u']).T * Mp * sp.Matrix(ls['u']))[0]

    def L_c(s):
        c, cs = s
        return c['u'][0] * cs['u'][0]
    names = []
    if wA:
        V.TERMS['WA_' + tag] = (sp.nsimplify(wA), ('A', 'As'), L_A, ('k', (0, 1))); names.append('WA_' + tag)
    if wl:
        V.TERMS['Wl_' + tag] = (-sp.nsimplify(wl), ('lam', 'lams'), L_l, ('k', (0, 1))); names.append('Wl_' + tag)
    if wc:
        V.TERMS['Wc_' + tag] = (-sp.nsimplify(wc), ('c', 'cs'), L_c, ('k', (0, 1))); names.append('Wc_' + tag)
    return names


def rho_two(Fmn, xi):
    v = sp.zeros(16, 1)
    xv = sp.Matrix(xi)
    for M in range(10):
        for N in range(10):
            if Fmn[M][N] != 0:
                v += sp.Rational(1, 2) * Fmn[M][N] * (Rbm[M] * (Rm[N] * xv))
    return v


def _dA(a):
    F = [[0] * 10 for _ in range(10)]
    for M in range(10):
        for N in range(10):
            v = 0
            if M < 4:
                v += I * a['p'][M] * a['u'][N]
            if N < 4:
                v -= I * a['p'][N] * a['u'][M]
            F[M][N] = v
    return F


def reg_xi(xi, tag, parts=('T1', 'T2', 'T3', 'T4', 'T5')):
    """terms of S^(1)(xi): T1 (A*_V), T2 (A*_perp), T3 (dA_V), T4 (dA_perp), T5 cubic; xi a list of 16 symbols."""
    xv = sp.Matrix(xi)
    MRxi = [(xv.T * MR[M]) for M in range(10)]

    def mk_lam(rng):
        def L_lam(s):
            l, As = s
            u = sp.Matrix(l['u'])
            return sum(((MRxi[M] * u)[0] * As['u'][M] for M in rng), 0)
        return L_lam

    def mk_F2(rng):
        def L_F2(s):
            a, ls = s
            aa = dict(u=[a['u'][M] if M in rng else 0 for M in range(10)], p=a['p'])
            w = rho_two(_dA(aa), xi)
            return (w.T * Mp * sp.Matrix(ls['u']))[0]
        return L_F2
    RRxi = [[(Rbm[M] * (Rm[N] * xv)) for N in range(10)] for M in range(10)]

    def L_F3(s):
        a, b, ls = s
        lv = Mp * sp.Matrix(ls['u'])
        tot = 0
        for M in range(10):
            for N in range(10):
                if M == N:
                    continue
                tot += a['u'][M] * b['u'][N] * (RRxi[M][N].T * lv)[0]
        return tot
    quad, cub = [], []
    if 'T1' in parts:
        V.TERMS['X1_' + tag] = (sp.Integer(-1), ('lam', 'As'), mk_lam(range(4)), ('k', (0, 1))); quad.append('X1_' + tag)
    if 'T2' in parts:
        V.TERMS['X2_' + tag] = (sp.Integer(-1), ('lam', 'As'), mk_lam(range(4, 10)), ('k', (0, 1))); quad.append('X2_' + tag)
    if 'T3' in parts:
        V.TERMS['X3_' + tag] = (sp.Rational(1, 2), ('A', 'lams'), mk_F2(range(4)), ('k', (0, 1))); quad.append('X3_' + tag)
    if 'T4' in parts:
        V.TERMS['X4_' + tag] = (sp.Rational(1, 2), ('A', 'lams'), mk_F2(range(4, 10)), ('k', (0, 1))); quad.append('X4_' + tag)
    if 'T5' in parts:
        V.TERMS['X5_' + tag] = (sp.Rational(1, 4), ('A', 'A', 'lams'), L_F3, ('f', (0, 1, 2))); cub.append('X5_' + tag)
    return quad, cub


# ------------------------------------------------------------------ generator letters
_G = {}
BSYM = sp.symbols('gb0:4')
XISYM = sp.symbols('gxi0:16')


def tmpl_dir():
    """Vdir(b) = 2[(b.dA,A*) - (b.d lam,lam*) - (b.d c,c*)] with the odd CE coordinate in front; gen = b (4)."""
    if 'dir' not in _G:
        names = reg_translation(BSYM, 'dir', scale=2)
        _G['dir'] = E.build_template(names, [], extra_odd=1, gen_syms=BSYM, name='G_dir')
    return _G['dir']


def tmpl_weights(wA, wl, wc, tag, sel=None, codd=1):
    key = ('w', wA, wl, wc, sel, codd)
    if key not in _G:
        names = reg_weights(wA, wl, wc, tag, sel)
        _G[key] = E.build_template(names, [], extra_odd=codd, name=f'G_w_{tag}')
    return _G[key]


def tmpl_xi(parts=('T1', 'T2', 'T3', 'T4', 'T5')):
    key = ('xi', tuple(parts))
    if key not in _G:
        tag = 'xi' + ''.join(p[1] for p in parts)
        quad, cub = reg_xi(XISYM, tag, parts)
        tq = E.build_template(quad, [], extra_odd=0, gen_syms=XISYM, name=f'G_{tag}_q') if quad else None
        t3 = None
        if cub:
            t3 = {spc: E.build_template(cub, [spc], extra_odd=0, gen_syms=XISYM, name=f'G_{tag}_c_{spc}') for spc in ('A', 'lams')}
        _G[key] = (tq, t3)
    return _G[key]


# ------------------------------------------------------------------ rotations: so(V) spin part and so(W_R^perp)
def reg_rotation(rmat, rng, tag, scale=1):
    """scale * [ (r A, A*) - (r lam, lam*) ],  (r A)_P = sum_N r_PN A_N (P, N in rng),
    r lam = (1/4) sum_MN r_MN rho_M rho_N lam (spin lift fixed by the covariance of (lam, rho(A) lam), checked)."""
    idx = list(rng)
    RR = {}
    for M in idx:
        for N in idx:
            if M != N:
                RR[(M, N)] = Rbm[M] * Rm[N]

    def L_A(s):
        a, As = s
        return sum((rmat[P][N] * a['u'][N] * As['u'][P] for P in idx for N in idx if P != N), 0)

    def L_l(s):
        l, ls = s
        lv = sp.Matrix(l['u'])
        w = sp.zeros(16, 1)
        for (M, N), G in RR.items():
            if rmat[M][N] != 0:
                w += sp.Rational(1, 4) * rmat[M][N] * (G * lv)
        return (w.T * Mp * sp.Matrix(ls['u']))[0]
    V.TERMS['RA_' + tag] = (sp.Integer(scale), ('A', 'As'), L_A, ('k', (0, 1)))
    V.TERMS['Rl_' + tag] = (sp.Integer(-scale), ('lam', 'lams'), L_l, ('k', (0, 1)))
    return ['RA_' + tag, 'Rl_' + tag]


OMSYM = sp.symbols('gom0:6')          # omega_{mu nu}, mu < nu in V: (01,02,03,12,13,23)
RSYM = sp.symbols('gr0:15')           # r_{IJ}, 4 <= I < J <= 9
VPAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
WPAIRS = [(I, J) for I in range(4, 10) for J in range(I + 1, 10)]


def antisym(syms, pairs):
    m = [[0] * 10 for _ in range(10)]
    for s, (i, j) in zip(syms, pairs):
        m[i][j] = s
        m[j][i] = -s
    return m


def tmpl_Mspin():
    """V0 of S(M_omega): (L_v A)_mu spin part omega_{nu mu} A_nu and the spinor part -(1/4) omega rho rho, i.e. the
    rotation r = -omega; gen = omega (6 components, pairs VPAIRS).  CE coordinate odd."""
    if 'Mspin' not in _G:
        om = antisym(OMSYM, VPAIRS)
        r = [[-om[i][j] for j in range(10)] for i in range(10)]
        names = reg_rotation(r, range(4), 'Msp')
        _G['Mspin'] = E.build_template(names, [], extra_odd=1, gen_syms=OMSYM, name='G_Mspin')
    return _G['Mspin']


def tmpl_R():
    """S(R_r) = (r A_perp, A*_perp) - (r lam, lam*), gen = r (15 components, WPAIRS).  CE coordinate odd."""
    if 'R' not in _G:
        names = reg_rotation(antisym(RSYM, WPAIRS), range(4, 10), 'Rr')
        _G['R'] = E.build_template(names, [], extra_odd=1, gen_syms=RSYM, name='G_R')
    return _G['R']
