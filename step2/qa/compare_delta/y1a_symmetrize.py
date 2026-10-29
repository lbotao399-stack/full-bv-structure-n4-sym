"""QA compare_delta (Part D4, bilinear component): is engine A's Q(Y1a) nonlocal?

Engine A (d4_y1a.py / d4_y1a_ward.py): Y1a on legs (A_P at p1 [polarisation a], A^*_P at p2 [b]), p1 + p2 + k = 0,
alpha(Y1a) = 0, beta(Y1a) = T_O [ 8 p2/|p2|^2 + 8 |p1|^2 That(p1,p2) ],  T_O = (w1,a)(w2,b) + (w2,a)(w1,b) (vertex value),
That_nu(p1,p2) = (1/pi^2) int d^4l l_nu / (l^4 (l+p2)^2 (l-p1)^2) = - int_{x+y+z=1} x (y p2 - z p1)_nu / Delta^2.
Q (A^* leg -> -(1/2) d^*d A): alpha_Q = p2.beta, beta_Q = -(1/2)|p2|^2 beta  (one-ordering kernel K_Q).
Q(Y1a) is a quadratic functional of A_P alone: only K_Q + swap(K_Q) (its plane-wave value) is meaningful.
swap (bosonic, p1 fixed): alpha'(p1) = alpha(p2) - (div_q beta)(q = p2), beta'(p1) = beta(p2), with beta(q) the
first-leg function (second leg -q-k).  Claims tested here:
 (1) That(p2,p1) = -That(p1,p2);   (2) div_q That(q,-q-k) = 2/(|q|^2 |q+k|^2)  (from d.(l/l^4) = 2 pi^2 delta^4(l));
 (3) plane-wave value of Q(Y1a) = (alpha, beta) = (8, 4k) T_O = -2 x (plane-wave value of D_M O = (-4, -2k) T_O);
 (4) engine A's raw beta/T_O at its three kinematic points equals the closed form (fit check);
 (5) translation analogue: Y1a_T = (i/2) beta(Y1a) (s_mu = 2 x translation density), Q(Y1a_T) = -2 D_P O (D_P O <-> -i k T_O).
"""
import json
import mpmath as mp

mp.mp.dps = 20


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def That(p1, p2):
    a11 = dot(p1, p1); a22 = dot(p2, p2); a12 = dot(p1, p2)
    res = []
    for nu in range(4):
        def f(x, t):
            y = (1 - x) * t; z = (1 - x) * (1 - t)
            D = y * a22 + z * a11 - (y * y * a22 - 2 * y * z * a12 + z * z * a11)
            return -x * (y * p2[nu] - z * p1[nu]) / D ** 2 * (1 - x)
        res.append(mp.quad(f, [0, 1], [0, 1]))
    return res


def divThat_fd(q, k, h=mp.mpf('1e-4')):
    """finite-difference divergence of F(q) = That(q, -q-k), Richardson-extrapolated"""
    def F(qq, nu):
        p2 = [-qq[i] - k[i] for i in range(4)]
        a11 = dot(qq, qq); a22 = dot(p2, p2); a12 = dot(qq, p2)
        def f(x, t):
            y = (1 - x) * t; z = (1 - x) * (1 - t)
            D = y * a22 + z * a11 - (y * y * a22 - 2 * y * z * a12 + z * z * a11)
            return -x * (y * p2[nu] - z * qq[nu]) / D ** 2 * (1 - x)
        return mp.quad(f, [0, 1], [0, 1])
    tot = mp.mpf(0)
    for nu in range(4):
        def cd(hh):
            qp = list(q); qm = list(q); qp[nu] += hh; qm[nu] -= hh
            return (F(qp, nu) - F(qm, nu)) / (2 * hh)
        d1 = cd(h); d2 = cd(h / 2)
        tot += (4 * d2 - d1) / 3
    return tot


def beta_closed(q, k):
    p2 = [-q[i] - k[i] for i in range(4)]
    Th = That(q, p2)
    n2 = dot(p2, p2); n1 = dot(q, q)
    return [8 * p2[m] / n2 + 8 * n1 * Th[m] for m in range(4)], Th


pts = [((1, 2, 0, 1), (0, -1, 1, 2)), ((2, -1, 1, 0), (1, 1, -2, 1)), ((1, 1, 1, -1), (-2, 0, 1, 1))]
engA = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_a/delta_results.json'))
rawA = engA['task2_component']['pieces']['Y1a fits']['su2']
out = {'points': []}
for ip, (p1, k) in enumerate(pts):
    p1 = [mp.mpf(x) for x in p1]; k = [mp.mpf(x) for x in k]
    p2 = [-p1[i] - k[i] for i in range(4)]
    b1, Th12 = beta_closed(p1, k)
    Th21 = That(p2, p1)
    antisym = max(abs(Th12[m] + Th21[m]) for m in range(4))
    dv_fd = divThat_fd(p2, k)
    dv_th = 2 / (dot(p2, p2) * dot(p1, p1))
    # (4) fit check against engine A's raw numbers
    raw = [complex(s.strip('()').replace('-0j', '+0j')).real for s in rawA[ip]['beta/T_O']]
    fitdev = max(abs(mp.mpf(raw[m]) - b1[m]) for m in range(4))
    # (3) plane-wave value of Q(Y1a) in units of T_O
    n2 = dot(p2, p2); n1 = dot(p1, p1)
    # one-ordering kernel at first-leg momentum q: alpha_Q(q) = p2(q).beta(q), beta_Q(q) = -(1/2)|p2(q)|^2 beta(q)
    bq1 = b1
    aQ1 = dot(p2, bq1); bQ1 = [-n2 / 2 * x for x in bq1]
    bq2, _ = beta_closed(p2, k)          # first-leg function at q = p2 (second leg -p2-k = p1)
    aQ2 = dot(p1, bq2); bQ2 = [-n1 / 2 * x for x in bq2]
    # div_q of beta_Q(q) = -(1/2)|q+k|^2 beta(q) at q = p2, beta(q) = 8(-q-k)/|q+k|^2 + 8|q|^2 That(q,-q-k):
    # beta_Q(q) = 4(q+k) - 4 |q|^2 |q+k|^2 That(q,-q-k);  div = 16 - 4[(2q|q+k|^2 + 2(q+k)|q|^2).That + |q|^2|q+k|^2 div That]
    q = p2; qk = [q[i] + k[i] for i in range(4)]
    ThF = That(q, [-x for x in qk])
    grad = [2 * q[i] * dot(qk, qk) + 2 * qk[i] * dot(q, q) for i in range(4)]
    divbQ = 16 - 4 * (dot(grad, ThF) + dot(q, q) * dot(qk, qk) * dv_fd)
    alpha_pw = aQ1 + aQ2 - divbQ
    beta_pw = [bQ1[m] + bQ2[m] for m in range(4)]
    # targets: -2 D_M O = (8, 4k) T_O
    out['points'].append({
        'p1': [str(x) for x in p1], 'k': [str(x) for x in k],
        'That_antisymmetry_max_dev': mp.nstr(antisym, 5),
        'div_That_finite_difference': mp.nstr(dv_fd, 12), 'div_That_formula_2/(|p1|^2|p2|^2)': mp.nstr(dv_th, 12),
        'engineA_raw_beta_vs_closed_form_max_dev': mp.nstr(fitdev, 5),
        'Q(Y1a)_plane_wave_alpha/T_O': mp.nstr(alpha_pw, 12),
        'Q(Y1a)_plane_wave_beta/T_O': [mp.nstr(x, 12) for x in beta_pw],
        'target_-2D_MO_alpha/T_O': '8', 'target_-2D_MO_beta/T_O': [mp.nstr(4 * x, 6) for x in k],
        'translation: Q(Y1a_T)_plane_wave/(i T_O) = (1/2)beta_pw': [mp.nstr(x / 2, 12) for x in beta_pw],
        'translation target -2 D_P O = 2 i k T_O -> /(i T_O)': [mp.nstr(2 * x, 6) for x in k]})
    print(json.dumps(out['points'][-1]), flush=True)
out['conclusion'] = ('Q(Y1a) restricted to the bilinear (A_P,A_P) sector is local: plane-wave value (8, 4k) T_O = '
                     '-2 D_M O exactly (That-terms cancel by antisymmetry except the contact term 2/(|p1|^2|p2|^2)). '
                     'Engine A compared the unsymmetrised one-ordering kernel; its Ward split (-4O - D_M O) also used the '
                     'half-normalised average.  Translation analogue: Q(Y1a_T) = -2 D_P O.')
json.dump(out, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/compare_delta/y1a_symmetrize.json', 'w'), indent=1)
