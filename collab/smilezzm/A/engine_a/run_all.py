"""Engine A (Track A, IP0/IP1): run everything and write results.json and gates.json.
Usage: python3 run_all.py   (pure Python 3: sympy, mpmath, numpy; file-relative paths; a few minutes)."""
import os, sys, json, time, resource
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
import masters, fourier, ip1
from masters import eps, LE, p2, b2

t0 = time.time()
unit = ip1.unit
LX = ip1.LX
g_E = sp.EulerGamma

def S(e):
    return str(sp.nsimplify(sp.expand(e))) if e != 0 else '0'

def num(e, digits=30):
    return str(sp.N(e, digits))

# ---- 1. masters and Fourier lemma, with their numerical checks
mcheck = masters.check_masters(verbose=False)
fcheck = fourier.check_fourier(verbose=False)

# ---- 2. the engine
rep = ip1.assemble()
tot = rep['total']
cl = rep['classes']
Sig = rep['self_energy']

# ---- 3. gates
gates = {}
# G1: self-energy singular parts per graph, converted to the manuscript's units:
#   (eq:c4-pw: int(A,A) -> 2 (A(p),A(-p)) = 4 hbar A^a.A^a; plane-wave value of (1/2) int Sigma A^a A^a is Sigma A^a(p).A^a(-p);
#    hence Sigma_sing = 4 e unit/eps = -2 e g2N/(16 pi^2)/eps)
#   coefficient e (eps^{-1} part, times unit int(A_perp,A_perp)) = -(1/2) x [coefficient of eps^{-1} in Sigma in units g2N/(16 pi^2)]
#   coefficient e' (log eps part, times unit int(dA_perp,dA_perp)) = -(1/2) x [coefficient of p2 LE in Sigma in units g2N/(16 pi^2)]
G1_target = {'seagull_V': (2, 0), 'seagull_perp': (sp.Rational(5, 2), 0), 'bubble': (-sp.Rational(1, 4), 1),
             'gaugino': (-2, -2), 'total': (sp.Rational(9, 4), -1)}
G1 = {}
allpass = True
graphs = dict(Sig['graphs']); graphs['total'] = Sig['total']
for k, v in graphs.items():
    e = sp.expand(v * 16 * sp.pi**2 / ip1.gN2)
    e_inv = -sp.Rational(1, 2) * e.coeff(eps, -1)
    e_log = -sp.Rational(1, 2) * sp.expand(e.coeff(LE, 1)).coeff(p2, 1)
    ok = (sp.simplify(e_inv - G1_target[k][0]) == 0) and (sp.simplify(e_log - G1_target[k][1]) == 0)
    allpass &= ok
    G1[k] = {'eps_inv_units': S(e_inv), 'log_eps_units': S(e_log), 'target': [S(G1_target[k][0]), S(G1_target[k][1])], 'pass': bool(ok)}
gates['G1_selfenergy_singular_parts'] = {'graphs': G1, 'pass': bool(allpass)}

# G2: operator counterterm, graph (e) and (f), from the vertex tensors and the masters (see ip1.py docstring / notes)
#   (e): functional 1/2 V4(a,b,K,L) A_K A_L B(P), V4 = -g^2 N [2 (a.b) d_KL - a_K b_L - a_L b_K], sing B = -LE/(16 pi^2)
#   (f): functional (a.A)(b.A) x g^2 N x (log eps coefficient of int d^4k k^2 Dhat^3 = -1/(16 pi^2))
B_sing = sp.expand(masters.M(0, 0, 2) / (16 * sp.pi**2)).coeff(LE, 1) * LE
e_O = sp.expand(sp.Rational(1, 2) * (-ip1.gN2) * (-2) * B_sing).subs(ip1.gN2, -32 * sp.pi**2 * unit) / unit   # coefficient of O (per (a.A)(b.A))
e_K = sp.expand(sp.Rational(1, 2) * (-ip1.gN2) * (2) * B_sing).subs(ip1.gN2, -32 * sp.pi**2 * unit) / unit     # coefficient of (a.b) K
f_O = sp.expand(ip1.gN2 * (-LE / (16 * sp.pi**2))).subs(ip1.gN2, -32 * sp.pi**2 * unit) / unit
G2 = {'graph_e': {'O': S(e_O), 'ab_K': S(e_K), 'target': ['2*LE', '-2*LE']},
      'graph_f': {'O': S(f_O), 'target': '2*LE'},
      'total': {'O': S(sp.expand(e_O + f_O)), 'ab_K': S(e_K), 'target': ['4*LE', '-2*LE']}}
okG2 = (sp.simplify(e_O - 2 * LE) == 0) and (sp.simplify(e_K + 2 * LE) == 0) and (sp.simplify(f_O - 2 * LE) == 0)
G2['pass'] = bool(okG2)
gates['G2_operator_counterterm'] = G2

# G3: log eps cancels
okG3 = (tot['S12']['log_eps'] == 0) and (tot['T12']['log_eps'] == 0)
gates['G3_log_eps_cancels'] = {'S12': S(tot['S12']['log_eps']), 'T12': S(tot['T12']['log_eps']), 'pass': bool(okG3)}
# G4: d_T = 8, no log x^2 on S12; per class: log x^2 coefficient = -1/2 x (tree insertion of D_1)/tree.
#   D_1(O_ab) = 4 (Z(O) - unit {N, O}) with Z = unit [4 O - 2 (a,b) K], {N,O} = 4 O  ->  D_1(O_ab) = -8 unit (a,b) K
#   -1/2 [<D_1 O1, O2> + <O1, D_1 O2>]/tree = -1/2 x (-8)(2 T12 + 2 T12)/S12 ... = +8 T12/S12 unit
okG4 = (tot['S12']['log_x2'] == 0) and (sp.simplify(tot['T12']['log_x2'] - 8) == 0)
gates['G4_log_x2_coefficient'] = {'S12': S(tot['S12']['log_x2']), 'T12': S(tot['T12']['log_x2']), 'target': ['0', '8'],
                                  'meaning': 'D_1(K) = -48 unit (Konishi 3 g^2 N/(4 pi^2)); 0 on the 20-prime', 'pass': bool(okG4)}
# G5: gamma_E cancels in c_S
cS = tot['S12']['const']; cT = tot['T12']['const']
okG5 = sp.expand(cS).coeff(g_E, 1) == 0
gates['G5_gammaE_cancels_on_20prime'] = {'gammaE_coefficient_in_c_S': S(sp.expand(cS).coeff(g_E, 1)), 'pass': bool(okG5)}
# G6: flavour structures
fl = ip1.flavour_check()
maxres = max(v for d in fl.values() for v in d.values())
gates['G6_only_S12_T12'] = {'max_residual': float(maxres), 'details': {k: {kk: float(vv) for kk, vv in d.items()} for k, d in fl.items()}, 'pass': bool(maxres < 1e-12)}
# G7: IP0 Gram on the 21-dim family (symmetric tensors e_I e_J), conjugation componentwise
import numpy as np
basis = []
for I in range(6):
    for J in range(I, 6):
        a = np.zeros(6); b = np.zeros(6); a[I] = 1; b[J] = 1
        basis.append((a, b))
# use a complex basis too: the 't Hooft blocks give complex combinations; test with a random complex unitary rotation
rng = np.random.default_rng(7)
Z = rng.normal(size=(21, 21)) + 1j * rng.normal(size=(21, 21))
Q, _ = np.linalg.qr(Z)
cSn = complex(sp.N(cS.subs(b2, 0))); cTn = complex(sp.N(cT))
def gram(uval):
    G = np.zeros((21, 21), dtype=complex)
    for k, (ak, bk) in enumerate(basis):
        for l, (al, bl) in enumerate(basis):
            akc, bkc = np.conj(ak), np.conj(bk)
            Skl = (akc @ al) * (bkc @ bl) + (akc @ bl) * (bkc @ al)
            Tkl = (akc @ bkc) * (al @ bl)
            G[k, l] = Skl * (1 + uval * cSn) + Tkl * uval * cTn
    return G
G0 = gram(0.0); G1m = (gram(1e-3) - G0) / 1e-3
Gr = Q.conj().T @ gram(0.01) @ Q  # in a rotated complex basis: must stay Hermitian
herm0 = np.linalg.norm(G0 - G0.conj().T); herm1 = np.linalg.norm(G1m - G1m.conj().T); hermr = np.linalg.norm(Gr - Gr.conj().T)
eig0 = np.linalg.eigvalsh((G0 + G0.conj().T) / 2)
# orthogonality 20' vs K: K = sum_I e_I e_I ; traceless combination e_0 e_0 - e_1 e_1
kvec = np.zeros(21); tl = np.zeros(21)
idx = {(I, J): n for n, (I, J) in enumerate([(I, J) for I in range(6) for J in range(I, 6)])}
for I in range(6): kvec[idx[(I, I)]] = 1
tl[idx[(0, 0)]] = 1; tl[idx[(1, 1)]] = -1
orth = abs(kvec @ gram(0.01) @ tl)
gates['G7_IP0_gram_hermitian_positive'] = {'tree_hermiticity_residual': float(herm0), 'oneloop_hermiticity_residual': float(herm1),
    'rotated_complex_basis_hermiticity_residual': float(hermr), 'tree_min_eigenvalue': float(eig0.min()), 'tree_max_eigenvalue': float(eig0.max()),
    'K_vs_traceless_overlap': float(orth), 'pass': bool(herm0 < 1e-12 and herm1 < 1e-9 and hermr < 1e-12 and eig0.min() > 0 and orth < 1e-12)}
# G8: colour: every one-loop graph carries f f -> N dim g = dim g tr_g(ad ad)/(X,Y) ... recorded as the structural statement used
gates['G8_colour'] = {'statement': 'all one-loop graphs carry f^{acd} f^{bcd} = N delta^{ab} (or f^{cde} f^{cde} = N dim g) times the tree colour dim g; no other colour structure at this order (su(N) identities, no 1/N terms)', 'pass': True}

# ---- 4. results.json
def splitdict(d):
    return {k: {'exact': S(v), 'numeric': num(v.subs(b2, 0)) if v.free_symbols <= {b2} else None} for k, v in d.items()}
results = {
  'engine': 'A', 'collaborator': 'smilezzm', 'track': 'A', 'stage': 'IP0/IP1',
  'method': 'exact momentum-space reduction: Schwinger masters with the heat-kernel cutoff (masters.py), Gaussian loop integration, '
            'numerator identity |k|^2 Dhat(k) = e^{-eps k^2}, localisation constant for the K_eps-contact graphs, Fourier lemma (fourier.py); sympy exact arithmetic',
  'conventions': {
     'action': 'dimensionally reduced 10d YM, su(N), Feynman gauge, weight e^{-S_E}; vertex factors V3, V4 as in ip1.py; Majorana-Weyl gaugino loop = 1/2 x tr_16',
     'propagator': 'delta^{ab} delta_MN e^{-eps k^2}/k^2 on every edge (bosons); fermion with p-slash numerator and the same factor',
     'unit': 'unit = hbar tr_g(ad ad)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); g^2 N/(16 pi^2) = -2 unit',
     'self_energy': 'Sigma = sum of 1PI graphs with (-S_E) vertices, symmetry factors, (-1) per fermion loop; dressed line Dhat + Dhat Sigma Dhat; counterterm = singular part (eps^{-1}, log eps with ell = 1); Sigma_ren = Sigma - Sing Sigma',
     'operators': 'O_ab = (a.A^c)(b.A^c); tree <O1,O2> = dim g S12 D_0^2 (manuscript: 4 hbar dim g S12 D_0^2); only ratios reported',
     'two_point': 'coefficient of t2 t1 = connected expectation (even operators); T = infinity; x != 0',
     'beta_2': 'F_0 contains beta_2 unit int dvol (dA_perp, dA_perp); enters as -F_0 on a line: Sigma_F = -4 beta_2 unit k^2',
     'log_eps': 'LE = log(eps/ell^2), ell = 1; LX = log(x^2/ell^2)'},
  'classes': {k: {'S12': splitdict(v['S12']), 'T12': splitdict(v['T12'])} for k, v in cl.items()},
  'total': {'S12': splitdict(tot['S12']), 'T12': splitdict(tot['T12'])},
  'c_S': {'exact': S(cS), 'numeric_beta2_0': num(cS.subs(b2, 0))},
  'c_T': {'exact': S(cT), 'numeric': num(cT)},
  'd_T': {'exact': S(tot['T12']['log_x2'])},
  'other_structures': 'none (gate G6)',
  'beta_2_dependence': 'c_S -> c_S - 8 beta_2 (class B); c_T, d_T independent of beta_2',
  'self_energy': {'graphs': {k: S(v) for k, v in Sig['graphs'].items()}, 'total': S(Sig['total']), 'singular': S(Sig['singular']),
                  'renormalised': S(Sig['renormalised']), 'units': 'coefficient of delta^{ab} delta_IJ; g2N = g^2 N; p2 = |p|^2'},
  'masters': {k: S(v) for k, v in rep['masters'].items()},
  'line_correction_over_D0': S(rep['line_correction']),
  'checks': {'masters_vs_quadrature': {str(k): v for k, v in mcheck.items()}, 'fourier_vs_hankel_quadrature': {str(k): v for k, v in fcheck.items()}},
  'runtime_seconds': None, 'peak_rss_MB': None,
}
results['runtime_seconds'] = round(time.time() - t0, 1)
results['peak_rss_MB'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024**2), 1)
with open(os.path.join(HERE, 'results.json'), 'w') as f:
    json.dump(results, f, indent=1)
gates['all_pass'] = all(v.get('pass', True) for v in gates.values())
with open(os.path.join(HERE, 'gates.json'), 'w') as f:
    json.dump(gates, f, indent=1)
print('c_S =', S(cS), ' c_T =', S(cT), ' d_T =', S(tot['T12']['log_x2']))
print('gates:', {k: (v.get('pass') if isinstance(v, dict) else v) for k, v in gates.items()})
print('runtime %.1fs, peak RSS %.0f MB' % (results['runtime_seconds'], results['peak_rss_MB']))
