"""compare_delta_b / y1a_check.py

The engine conflict on the remainder Q(Y1a) (engine A: nonlocal; theory: the complete composite set is local;
engine B: K5 local).  Independent evaluation by this QA agent:

 (a) That_nu(p1,p2) := 16 pi^2 int d^4l/(2pi)^4 l_nu/(l^4 (l+p2)^2 (l-p1)^2)
       = - int_{x+y+z=1} x (y p2 - z p1)_nu / Delta^2,  Delta = y|p2|^2 + z|p1|^2 - |y p2 - z p1|^2
     (Feynman parameters re-derived here; own mpmath integrator, smooth parametrisation x = 1-u, y = u t, z = u(1-t)).
     Bhat(p) = 16 pi^2 int l/(l^4 (l+p)^2) = -p/|p|^2.
 (b) engine A's raw numbers beta(Y1a)/T_O at its three su(2) points against 8 p2/|p2|^2 + 8 |p1|^2 That(p1,p2)
     (p2 = -k - p1): checks A's closed form, independently of A's own fit.
 (c) exact identities: That(p2,p1) = -That(p1,p2) (l -> -l); div_q That(q,-k-q) = 2/(|q|^2 |q+k|^2)
     (div x/|x|^4 = 2 pi^2 delta^4 in R^4); numerically checked.
 (d) Q(Y1a): A* leg -> -(1/2) d*dA: alpha_Q = p2.beta, beta_Q = -(1/2)|p2|^2 beta (A's map; it is the map of
     two_point_sectors.py).  The output has two identical A|_Wperp legs, so only its symmetrisation is a functional:
       G_sym = (1/2)(G + G_swap),  alpha''(p) = alpha(p_d) - (div beta)(p_d), beta''(p) = beta(p_d), p_d = -k - p.
     For beta(Y1a) = T_O[l1 p2/|p2|^2 + l2 |p1|^2 That(p1,p2)] the identities give exactly
       plane-wave value of Q(Y1a) = 2 G_sym = (l2, (l1/2) k) T_O,
     i.e. -2 D_M O for l1 = l2 = 8 (D_M O <-> (-4, -2k) T_O), a LOCAL functional.  Checked numerically.
 (e) the decisive check of the theory note (delta_predictions.json): (hbar/2) d_Psi(Q O_tree) on (A_perp,A_perp)
     = Q(Y1a) + Q(Y1b) by [Q, d_Psi] = 0; its composite special-edge cycles through O are
       Q(Y1a) + {Pi_int|_perp, O}   (engine A)   vs   K5 - 4 x (log eps coefficient of d_delta I^CT)   (engine B).
     Engine B's three raw exact K5 samples are refitted here.
 (f) translations: Y1a_T (S^(1)(P_mu) in place of S_M) = (i/2) beta_mu(Y1a): Q(Y1a_T) = -2 D_P O.
Output: y1a_check.json
"""
import json, os, random
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 25


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def That(p1, p2):
    p1 = [mp.mpf(v) for v in p1]; p2 = [mp.mpf(v) for v in p2]
    a11 = dot(p1, p1); a22 = dot(p2, p2); a12 = dot(p1, p2)
    out = []
    for nu in range(4):
        def f(u, t):
            vv2 = t * t * a22 + (1 - t) ** 2 * a11 - 2 * t * (1 - t) * a12      # |t p2 - (1-t) p1|^2
            D = t * a22 + (1 - t) * a11 - u * vv2
            return -(1 - u) * (t * p2[nu] - (1 - t) * p1[nu]) / D ** 2
        out.append(mp.quad(f, [0, 1], [0, 1]))
    return out


res = {}
# ------------------------------------------------------------------ (b) engine A raw numbers
A = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_a/delta_results.json'))
fits = A['task2_component']['pieces']['Y1a fits']
rb = []
for alg in ['su2', 'gl3']:
    for pt in fits[alg]:
        p1 = pt['p1']; k = pt['k']; p2 = [-a - b for a, b in zip(k, p1)]
        raw = [complex(s.replace('(', '').replace(')', '')).real for s in pt['beta/T_O']]
        T = That(p1, p2)
        cf = [8 * mp.mpf(p2[m]) / dot(p2, p2) + 8 * dot(p1, p1) * T[m] for m in range(4)]
        rb.append({'alg': alg, 'p1': p1, 'k': k, 'engineA_raw_beta/T_O': raw,
                   'closed_form 8p2/|p2|^2 + 8|p1|^2 That': [mp.nstr(c, 15) for c in cf],
                   'max_abs_dev': mp.nstr(max(abs(cf[m] - raw[m]) for m in range(4)), 3)})
res['b_engineA_raw_vs_closed_form'] = rb

# ------------------------------------------------------------------ (c) identities
random.seed(20261028)
idc = []
for trial in range(2):
    q = [random.choice([-3, -2, -1, 1, 2, 3]) for _ in range(4)]
    k = [random.choice([-2, -1, 1, 2]) for _ in range(4)]
    pd = [-a - b for a, b in zip(k, q)]
    T12 = That(q, pd); T21 = That(pd, q)
    anti = max(abs(T12[m] + T21[m]) for m in range(4))
    h = mp.mpf('1e-6')
    dv = mp.mpf(0)
    for mu in range(4):
        qp = list(map(mp.mpf, q)); qm = list(map(mp.mpf, q))
        qp[mu] += h; qm[mu] -= h
        Tp = That(qp, [-a - b for a, b in zip(k, qp)]); Tm = That(qm, [-a - b for a, b in zip(k, qm)])
        dv += (Tp[mu] - Tm[mu]) / (2 * h)
    qk = [a + b for a, b in zip(q, k)]
    idc.append({'q': q, 'k': k, 'max |That(q,p_d) + That(p_d,q)|': mp.nstr(anti, 3),
                'div_q That(q,-k-q) (finite difference)': mp.nstr(dv, 15),
                '2/(|q|^2|q+k|^2)': mp.nstr(mp.mpf(2) / (dot(q, q) * dot(qk, qk)), 15)})
res['c_identities'] = idc


# ------------------------------------------------------------------ (d) Q(Y1a) symmetrised
def betaY(q, k, l1, l2):
    pd = [-a - b for a, b in zip(k, q)]
    T = That(q, pd)
    return [l1 * mp.mpf(pd[m]) / dot(pd, pd) + l2 * dot(q, q) * T[m] for m in range(4)]   # / T_O


def QY_sym(p, k, l1=8, l2=8, h=mp.mpf('1e-6')):
    pd = [-a - b for a, b in zip(k, p)]
    bp = betaY(p, k, l1, l2); bpd = betaY(pd, k, l1, l2)
    aQ = lambda q, b_: dot([-a - c for a, c in zip(k, q)], b_)                 # alpha_Q = p2.beta
    bQ = lambda q, b_: [-dot([-a - c for a, c in zip(k, q)], [-a - c for a, c in zip(k, q)]) / 2 * x for x in b_]
    # div_q beta_Q at q = p_d
    dv = mp.mpf(0)
    for mu in range(4):
        qp = list(map(mp.mpf, pd)); qm = list(map(mp.mpf, pd))
        qp[mu] += h; qm[mu] -= h
        dv += (bQ(qp, betaY(qp, k, l1, l2))[mu] - bQ(qm, betaY(qm, k, l1, l2))[mu]) / (2 * h)
    alpha_unsym = aQ(p, bp)
    alpha_sym = (aQ(p, bp) + aQ(pd, bpd) - dv) / 2
    beta_sym = [(x + y) / 2 for x, y in zip(bQ(p, bp), bQ(pd, bpd))]
    return alpha_unsym, bQ(p, bp), alpha_sym, beta_sym


pts = [(fits['su2'][i]['p1'], fits['su2'][i]['k']) for i in range(3)]
for trial in range(2):
    pts.append(([random.choice([-3, -2, -1, 1, 2, 3]) for _ in range(4)], [random.choice([-2, -1, 1, 2]) for _ in range(4)]))
rd = []
for (p1, k) in pts:
    au, bu, asym, bsym = QY_sym(p1, k)
    pw_alpha = 2 * asym; pw_beta = [2 * x for x in bsym]
    rd.append({'p1': p1, 'k': k,
               'unsymmetrised alpha_Q/T_O (engine A compared this)': mp.nstr(au, 12),
               'unsymmetrised beta_Q/(T_O) ': [mp.nstr(x, 10) for x in bu],
               'plane-wave Q(Y1a)/T_O alpha = 2 alpha_sym': mp.nstr(pw_alpha, 15),
               'plane-wave Q(Y1a)/T_O beta = 2 beta_sym': [mp.nstr(x, 15) for x in pw_beta],
               'target -2 D_M O: alpha': 8, 'target -2 D_M O: beta = 4k': [4 * x for x in k],
               'max_dev_from_target': mp.nstr(max([abs(pw_alpha - 8)] + [abs(pw_beta[m] - 4 * k[m]) for m in range(4)]), 3)})
res['d_QY1a_symmetrised'] = rd
# general (l1,l2) at one point: confirm plane wave = (l2, (l1/2) k)
p1, k = pts[3]
gen = {}
for (l1, l2) in [(8, 0), (0, 8), (3, 5)]:
    _, _, asym, bsym = QY_sym(p1, k, l1, l2)
    gen['l1=%d,l2=%d' % (l1, l2)] = {'pw alpha': mp.nstr(2 * asym, 12), 'expected l2': l2,
                                       'pw beta': [mp.nstr(2 * x, 12) for x in bsym], 'expected (l1/2) k': [l1 * x / 2 for x in k]}
res['d_general_coefficients'] = gen
res['d_exact_statement'] = ('For beta(Y1a) = T_O[l1 p2/|p2|^2 + l2 |p1|^2 That(p1,p2)], alpha(Y1a) = 0: plane-wave value of the '
                            'symmetrised Q(Y1a) = (alpha, beta) = (l2, (l1/2) k) T_O exactly (the That terms cancel by antisymmetry '
                            'in beta, by antisymmetry plus the contact identity in alpha). Engine A: l1 = l2 = 8 => (8, 4k) T_O = '
                            '-2 D_M O, local.  Engine A compared the unsymmetrised kernel of one ordering of the two identical '
                            'A|_Wperp output legs; that kernel carries That and is not a functional by itself.')

# ------------------------------------------------------------------ (e) decisive check against engine B's K5
Bk = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_b/delta_K5_ApAp_s31.json'))
from fractions import Fraction as Fr
refit = []
for smp in Bk['ApAp']:
    w1 = [Fr(x) for x in smp['w1']][4:]; w2 = [Fr(x) for x in smp['w2']][4:]
    u1 = [Fr(x) for x in smp['u1']][4:]; u2 = [Fr(x) for x in smp['u2']][4:]
    kk = [Fr(x) for x in smp['k']]
    T1 = dot(w1, u1) * dot(w2, u2) + dot(w1, u2) * dot(w2, u1)
    T2 = dot(w1, w2) * dot(u1, u2)
    a0 = Fr(smp['res']['0']['alpha']['1'])
    b0 = [Fr(smp['res']['0']['beta'][m]['1']) for m in range(4)]
    logs_zero = all(smp['res'][o]['alpha'][L] == '0' for o in ['-1', '0'] for L in ['L4', 'L3']) and \
        all(smp['res'][o]['beta'][m][L] == '0' for o in ['-1', '0'] for m in range(4) for L in ['L4', 'L3']) and \
        smp['res']['-1']['alpha']['1'] == '0' and all(smp['res']['-1']['beta'][m]['1'] == '0' for m in range(4))
    refit.append({'T1': str(T1), 'T2': str(T2), 'alpha': str(a0), '16T1-16T2': str(16 * T1 - 16 * T2),
                  'beta': [str(x) for x in b0], '(2/3)k T1': [str(Fr(2, 3) * x * T1) for x in kk],
                  'alpha_ok': a0 == 16 * T1 - 16 * T2, 'beta_ok': b0 == [Fr(2, 3) * x * T1 for x in kk],
                  'eps^-1 and log parts zero': logs_zero})
res['e_engineB_K5_refit'] = refit
# composite cycles through O, in plane-wave units of T_O (T1 = T_O for the traceless structure):
#   engine A: Q(Y1a)_sym + {Pi_int|_perp, O},  {Pi_int,O} = (5/3) D_M O - (2/3) {int(A,A*),O}, {int(A,A*),O} = 2 O
QY = (Fr(8), Fr(4))          # (alpha, beta/k)
DMO = (Fr(-4), Fr(-2))
PiO = (Fr(5, 3) * DMO[0] - Fr(2, 3) * 2, Fr(5, 3) * DMO[1])
A_cycles = (QY[0] + PiO[0], QY[1] + PiO[1])
B_cycles = (Fr(16) - Fr(16), Fr(2, 3))       # K5 alpha/T1 minus 4 x log coefficient (16 T1 - 16 T2 at T2 = 0); beta/(k T1)
res['e_decisive_check'] = {
    'theory': 'lim (hbar/2) d_Psi(Q O_tree[eps]) on (A_perp,A_perp) = Q(Y1a) + Q(Y1b), local, special-edge cycles through O plus the tree insertion of Q Pi_Psi',
    'engine A (Y1a symmetrised): Q(Y1a) + {Pi_int,O}': {'alpha/T_O': str(A_cycles[0]), 'beta/(k T_O)': str(A_cycles[1])},
    'engine B: K5 - 4 x (log eps coefficient of d_delta I^CT)': {'alpha/T1': str(B_cycles[0]), 'beta/(k T1)': str(B_cycles[1])},
    'agree': A_cycles == B_cycles,
    'equivalently': {'Q(Y1a)+Q(Y1b) (engine A)': '(4/3) O - (1/3) D_M O',
                     'K5 composite part + C_P(Q Pi_Psi, O) (engine B, theory)': '-(4/3) O - (1/3) D_M O + (8/3) O = (4/3) O - (1/3) D_M O'},
    'what it fixes': 'with l1, l2 free, agreement requires l2 - 8 = 0 and l1/2 - 10/3 = 2/3: both of engine A\'s Y1a coefficients are confirmed by engine B\'s K5',
}

# ------------------------------------------------------------------ (f) translations
tr = []
for (p1, k) in pts[:3]:
    pd = [-a - b for a, b in zip(k, p1)]
    bp = betaY(p1, k, 8, 8); bpd = betaY(pd, k, 8, 8)
    # alpha_T(q) = (i/2) beta(q) ; Q: alpha -> -(1/2)|p2|^2 alpha ; symmetrise (no beta): pw = alpha_Q(p) + alpha_Q(p_d)
    aQp = [-dot(pd, pd) / 2 * x / 2 for x in bp]       # divided by i
    aQpd = [-dot(p1, p1) / 2 * x / 2 for x in bpd]
    pw = [x + y for x, y in zip(aQp, aQpd)]           # plane wave / (i T_O)
    tr.append({'p1': p1, 'k': k, 'Q(Y1a_T) plane wave/(i T_O)': [mp.nstr(x, 12) for x in pw],
               'target -2 D_P O = 2 i k T_O -> 2k': [2 * x for x in k]})
res['f_translation'] = tr
res['f_translation_total'] = {
    'P_1(O) on (A_perp,A_perp)': 'Q(Y1a_T) + Q(Y1b_T) - {F(P),O} = -2 D_P O + (5/3) D_P O - f_P D_P O = -(1/3 + f_P) unit {S^(1)(P_b),O}',
    'engine A as reported': 'nonlocal Q-exact term (same artefact as for Delta)',
    'engine B': '-(1/3) unit {S^(1)(P),O} (f_P = 0)',
    'theory': 'R_{O,v} = c {S^(1)(P),O}, c = -1/3 unit (taken from engine B)'}
json.dump(res, open(os.path.join(HERE, 'y1a_check.json'), 'w'), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
