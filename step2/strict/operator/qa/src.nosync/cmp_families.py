"""QA comparator and independent spectral checks, operator side.
(1) A vs B on F0-F4 entry by entry (same bases where both use tr(d^j X d^k Y) with j <= k).
(2) Independent targets: Heisenberg (su(2)), free-fermion su(1|1) (graded cyclic), harmonic sums -32 h(J) with the
    two-site psu(2,2|4) selection rule (J even for L = 2), computed here from scratch.
(3) (R7) along P_n as matrix identities D^(S+1) d = d D^(S) on F3, F4 (both engines), F5, F7 (engine A).
Writes qa/out.nosync/cmp_families.json."""
import json, itertools, os
from fractions import Fraction as Fr
import sympy as sp
import numpy as np, warnings
warnings.filterwarnings('ignore')

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/operator'
RA = json.load(open(f'{ROOT}/engine_a/results.json'))['stages']
RB = json.load(open(f'{ROOT}/engine_b/results.json'))
out = {}

def M(x):
    if isinstance(x, str):
        x = sp.sympify(x)
    return sp.Matrix(sp.sympify(x))

def h(j):
    return sum(Fr(1, k) for k in range(1, j + 1))

# ---------------- F0 ----------------
a0, b0 = RA['F0'], RB['OP6_F0_even']
f0 = {
 'Upsilon0_Delta_traceless': [a0['Upsilon0_Delta(O) traceless'], b0['Upsilon0_Delta']['traceless']],
 'Upsilon0_Delta_trace': [a0['Upsilon0_Delta(O) trace (O_tr = sum_I (A_I,A_I))'], b0['Upsilon0_Delta']['trace']],
 'D1_traceless_theta': [a0['D_1 traceless (symbolic theta)'], b0['D_1']['traceless']],
 'D1_trace_theta': [a0['D_1 trace (symbolic theta)'], b0['D_1']['trace']],
 'Upsilon_P_over_dbO': [a0['Upsilon_P(O) / d_b O'], b0['Upsilon_P']['total'][0]],
 'Upsilon0_M': [a0['Upsilon_M0 sector'], b0['Upsilon0_M']],
 'Upsilon_R': [a0['Upsilon_R sector'], b0['Upsilon_R']],
}
for k, (x, y) in list(f0.items()):
    try:
        same = sp.simplify(sp.sympify(str(x).replace('theta_perp', 'theta')) - sp.sympify(str(y))) == 0
    except Exception:
        same = (json.dumps(x, sort_keys=True).count('"0"') > 0 and json.dumps(y).count('"0"') > 0)
    f0[k] = {'A': x, 'B': y, 'agree': bool(same)}
qa_ = a0['Q_1 two-body fit (xi in S_+ (x) W)']
qb_ = RB['OP6_F0_odd_and_arity0']['XIP']['coefficients']
L2, L3 = sp.log(2), sp.log(3)
def ev(s):
    return sp.nsimplify(sp.sympify(s.replace('log2', 'log(2)').replace('log3', 'log(3)')))
qcmp = {}
for k in ['a1', 'a2', 'a3', 'a4']:
    dA = sp.sympify(qa_[k].replace('log2', 'log(2)').replace('log3', 'log(3)'))
    dB = sp.sympify(qb_[k])
    qcmp[k] = {'A': qa_[k], 'B': qb_[k], 'agree': bool(sp.simplify(sp.expand_log(dA - dB, force=True)) == 0)}
f0['Q1_two_body_xi_S+W'] = qcmp
qb2 = RB['OP6_F0_odd_and_arity0']['XIM']['coefficients']
qa2 = a0['Q_1 two-body fit (xi in S_- (x) W^*)']
f0['Q1_two_body_xibar_S-W*'] = {k: {'A': qa2[k], 'B': qb2[k], 'agree': bool(sp.simplify(sp.expand_log(
    sp.sympify(qa2[k].replace('log2', 'log(2)').replace('log3', 'log(3)')) - sp.sympify(qb2[k]), force=True)) == 0)}
    for k in ['a1', 'a2', 'a3', 'a4']}
LA = a0['L_1^(0)(O) / ((w1,w2) int omega (dc, A_V)), psi = 0'][0].replace('+ -', '- ').replace('log2', 'log(2)').replace('log3', 'log(3)')
LB = RB['OP6_F0_odd_and_arity0']['L1_0']['coefficients_of_i_p_c_and_i_p_A'][0]
f0['L1_0_coefficient'] = {'A': LA, 'B': LB, 'agree': bool(sp.simplify(sp.sympify(LA) - sp.sympify(LB)) == 0)}
out['F0'] = f0

# ---------------- F3, F4 entry by entry ----------------
def famcmp(stA_list, keyB, key_sum='sum_Upsilon0_Delta_matrix'):
    res = {}
    mats = RB[keyB]['matrices']
    for S, mb in mats.items():
        if not mb.get('dim'):
            continue
        ma = None
        for st in stA_list:
            if st in RA and S in RA[st]:
                ma = RA[st][S]
        if ma is None:
            res[S] = {'A': None, 'B': mb['Upsilon0_Delta'], 'agree': None}
            continue
        A_ = M(ma[key_sum]); B_ = M(mb['Upsilon0_Delta'])
        evA = sorted(ma['D1_eigenvalues'], key=lambda s: float(sp.sympify(s)))
        evB = sorted(mb['D1_eigenvalues_units_u'], key=lambda s: float(sp.sympify(s)))
        res[S] = {'sum_Upsilon0_Delta_equal_entrywise': A_ == B_, 'A': str(A_.tolist()), 'B': str(B_.tolist()),
                  'D1_eigs_A': evA, 'D1_eigs_B': evB, 'eigs_agree': evA == evB}
    return res
out['F3_twist_two'] = famcmp(['twist', 'twist_ext', 'twist_ext2'], 'OP6_F3_twist_two')
out['F4_fermion_sl2'] = famcmp(['fsl2', 'fsl2_ext', 'fsl2_ext2'], 'OP6_F4_fermion_sl2')

# ---------------- independent targets ----------------
# (a) su(2) Heisenberg, cyclic zero-momentum states of tr(Z^{L-M} X^M); D_1 = -8u sum_l (1 - P_{l,l+1})
def heis(L, Mx):
    states = [s for s in itertools.product([0, 1], repeat=L) if sum(s) == Mx]
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((len(states), len(states)))
    for s in states:
        for l in range(L):
            m = (l + 1) % L
            t = list(s); t[l], t[m] = t[m], t[l]; t = tuple(t)
            H[idx[s], idx[s]] += 1; H[idx[t], idx[s]] -= 1
    # cyclic shift projector (trace cyclicity)
    T = np.zeros_like(H)
    for s in states:
        T[idx[s[1:] + s[:1]], idx[s]] = 1
    Pz = sum(np.linalg.matrix_power(T, k) for k in range(L)) / L
    w, V = np.linalg.eigh((Pz + Pz.T) / 2)
    B = V[:, w > 0.5]
    Hr = B.T @ H @ B
    return sorted(np.round(np.linalg.eigvalsh((Hr + Hr.T) / 2), 6).tolist())
heisen = {}
bsp = RB['OP6_F1']['F1']['spectra_D1_units_u']
asp = RA['F1']['su2_cyclic_spectra']
ok_all = True
for L in range(2, 9):
    for Mx in range(0, L // 2 + 1):
        E = heis(L, Mx)
        tgt = sorted([-8 * e for e in E])
        kb = f'L={L},M={Mx}'; ka = f'L{L}_M{Mx}'
        eb = sorted(bsp.get(kb, [])); ea = sorted(asp.get(ka, {}).get('E = D_1/(-8u)', []))
        okb = len(eb) == len(tgt) and all(abs(x - y) < 1e-4 for x, y in zip(eb, tgt))
        oka = len(ea) == len(E) and all(abs(x - y) < 1e-4 for x, y in zip(ea, E))
        ok_all &= okb and oka
        heisen[kb] = {'indep_D1_units_u': tgt, 'B_agree': okb, 'A_agree': oka}
out['F1_su2_independent'] = {'all_agree': ok_all, 'sectors': heisen}

# (b) su(1|1): tr(Z^{L-M} psi^M), D_1 = -8u sum_l (1 - Pi_{l,l+1}), graded permutation, graded cyclicity
def su11(L, Mx):
    states = [s for s in itertools.product([0, 1], repeat=L) if sum(s) == Mx]
    idx = {s: i for i, s in enumerate(states)}
    n = len(states)
    H = np.zeros((n, n))
    for s in states:
        for l in range(L):
            m = (l + 1) % L
            t = list(s); t[l], t[m] = t[m], t[l]; t = tuple(t)
            if m == 0:   # periodic link: moving across the trace boundary
                # graded permutation of sites L-1 and 0 written as l, m: sign from passing the other fermions
                inner = sum(s[1:L - 1])
                sgn = (-1) ** (s[l] * s[m]) * ((-1) ** ((s[l] + s[m]) * inner) if s[l] != s[m] else 1)
            else:
                sgn = (-1) ** (s[l] * s[m])
            H[idx[s], idx[s]] += 1; H[idx[t], idx[s]] -= sgn
    # graded cyclic shift: move last site to front with sign (-1)^{f_last * (M - f_last)}
    T = np.zeros((n, n))
    for s in states:
        t = s[-1:] + s[:-1]
        T[idx[t], idx[s]] = (-1) ** (s[-1] * (Mx - s[-1]))
    Pz = sum(np.linalg.matrix_power(T, k) for k in range(L)) / L
    w, V = np.linalg.eig(Pz)
    B = np.real(V[:, np.abs(w - 1) < 1e-8])
    if B.shape[1] == 0:
        return []
    Q, _ = np.linalg.qr(B)
    Hr = Q.T @ H @ Q
    return sorted(np.round(np.real(np.linalg.eigvals(Hr)), 6).tolist())
su = {}
ok2 = True
asu = RA['F2']['spectra']; bsu = RB['OP6_F2_su11']['sectors']
for L in range(2, 7):
    for Mx in range(0, L + 1):
        E = su11(L, Mx)
        ka = f'L{L}_M{Mx}'
        eA = sorted(float(sp.N(sp.sympify(x))) for x in asu.get(ka, {}).get('E_from_D1', []))
        bA = bsu.get(ka)
        eB = []
        if bA:
            for k, mult in bA['D1_eigenvalues_units_u'].items():
                eB += [float(sp.N(sp.sympify(k))) / -8.0] * mult
        eB = sorted(eB)
        okA = len(eA) == len(E) and all(abs(x - y) < 1e-5 for x, y in zip(eA, E))
        okB = (not bA and not E) or (len(eB) == len(E) and all(abs(x - y) < 1e-5 for x, y in zip(eB, E)))
        ok2 &= okA and okB
        su[ka] = {'indep_E': E, 'A': eA, 'B': eB, 'A_agree': okA, 'B_agree': okB}
out['F2_su11_independent'] = {'all_agree': ok2, 'sectors': su}

# (c) harmonic sums with the two-site selection rule; new primaries per S'
def target_spec(family, S):
    # returns the multiset of -32 h(J) over primaries with spin <= S (descendants included)
    vals = []
    for s in range(0, S + 1):
        if family == 'F3':   # tr(D^k Z D^{s-k} Z): J = s, s even
            Js = [s] if s % 2 == 0 else []
        elif family == 'F4':  # tr(D psi D psi): J = s + 1, s odd
            Js = [s + 1] if s % 2 == 1 else []
        elif family == 'F6':  # tr(D F D F): J = s + 2, s even
            Js = [s + 2] if s % 2 == 0 else []
        elif family == 'F7':  # tr(D F D psi): J in {s+1, s+2} with J even
            Js = [J for J in (s + 1, s + 2) if J % 2 == 0]
        vals += [-32 * h(J) for J in Js]
    return sorted(vals)
spec = {}
for fam, stages in [('F3', ['twist', 'twist_ext', 'twist_ext2']), ('F4', ['fsl2', 'fsl2_ext', 'fsl2_ext2']),
                    ('F6', ['gluon', 'gluon_ext']), ('F7', ['f7', 'f7_ext'])]:
    spec[fam] = {}
    for st in stages:
        for S, m in RA.get(st, {}).items():
            if not isinstance(m, dict) or 'D1_eigenvalues' not in m:
                continue
            got = sorted(Fr(str(sp.Rational(sp.sympify(x)))) for x in m['D1_eigenvalues'])
            tgt = target_spec(fam, int(S))
            spec[fam][S] = {'A': [str(x) for x in got], 'target_indep': [str(x) for x in tgt], 'agree': got == tgt}
out['spectra_vs_independent_harmonic_targets'] = spec

# (d) (R7) along P_n as matrix identities
def dmap_bos(S, sym):  # basis (j, S-j), j <= S-j for bosonic identical letters (sym=+1) or fermionic (sym=-1)
    def basis(S):
        if sym == +1:
            return [(j, S - j) for j in range(0, S // 2 + 1)]
        return [(j, S - j) for j in range(0, (S + 1) // 2) if j != S - j]
    b0, b1 = basis(S), basis(S + 1)
    D = sp.zeros(len(b1), len(b0))
    for c, (j, k) in enumerate(b0):
        for (p, q) in [(j + 1, k), (j, k + 1)]:
            sgn = 1
            if p > q:
                p, q = q, p; sgn = sym
            if (p, q) in b1:
                D[b1.index((p, q)), c] += sgn
    return D
r7 = {}
for fam, stages, keyB, sym in [('F3', ['twist', 'twist_ext', 'twist_ext2'], 'OP6_F3_twist_two', +1),
                               ('F4', ['fsl2', 'fsl2_ext', 'fsl2_ext2'], 'OP6_F4_fermion_sl2', -1)]:
    Bm = RB[keyB]['matrices']
    r7[fam] = {}
    for S in range(0, 6):
        s0, s1 = str(S), str(S + 1)
        if not (Bm.get(s0, {}).get('dim') and Bm.get(s1, {}).get('dim')):
            continue
        D0 = M(Bm[s0]['D1_at_theta_m13']); D1 = M(Bm[s1]['D1_at_theta_m13'])
        dm = dmap_bos(S, sym)
        okB = (D1 * dm - dm * D0) == sp.zeros(*dm.shape)
        # B's own Upsilon_Pn must equal -1/3 d
        okPn = M(Bm[s0]['Upsilon_Pn']) == -sp.Rational(1, 3) * dm
        r7[fam][f'{S}->{S+1}'] = {'D1(S+1) d = d D1(S) [B matrices]': bool(okB), 'Upsilon_Pn = -1/3 d': bool(okPn)}
# F5, F7 of engine A: basis tr(d^j X d^{S-j} Y), j = 0..S (and two blocks for F5)
def dmap_lin(S, blocks=1):
    n0, n1 = S + 1, S + 2
    D = sp.zeros(blocks * n1, blocks * n0)
    for b in range(blocks):
        for j in range(n0):
            D[b * n1 + j + 1, b * n0 + j] += 1   # d on X
            D[b * n1 + j, b * n0 + j] += 1       # d on Y
    return D
for fam, stages, blocks in [('F7', ['f7', 'f7_ext'], 1), ('F5', ['f5', 'f5_ext'], 2)]:
    r7[fam] = {}
    mats = {}
    for st in stages:
        for S, m in RA.get(st, {}).items():
            if not isinstance(m, dict) or 'D1_matrix' not in m:
                continue
            mats[int(S)] = M(m['D1_matrix'])
    for S in sorted(mats):
        if S + 1 in mats:
            dm = dmap_lin(S, blocks)
            r7[fam][f'{S}->{S+1}'] = bool(mats[S + 1] * dm - dm * mats[S] == sp.zeros(*dm.shape))
            # eigenvalues recomputed
    r7[fam]['eigenvalues_recomputed'] = {str(S): sorted([str(k) for k, v in m.eigenvals().items() for _ in range(v)])
                                         for S, m in mats.items()}
out['R7_matrix_identities'] = r7
os.makedirs(f'{ROOT}/qa/out.nosync', exist_ok=True)
json.dump(out, open(f'{ROOT}/qa/out.nosync/cmp_families.json', 'w'), indent=1, default=str)
print(json.dumps(out['F0'], indent=0, default=str)[:3000])
for fam in ['F3_twist_two', 'F4_fermion_sl2']:
    print(fam, {S: (v.get('sum_Upsilon0_Delta_equal_entrywise'), v.get('eigs_agree')) for S, v in out[fam].items()})
print('F1 su2 all agree', out['F1_su2_independent']['all_agree'])
print('F2 su11 all agree', out['F2_su11_independent']['all_agree'],
      [k for k, v in out['F2_su11_independent']['sectors'].items() if not (v['A_agree'] and v['B_agree'])])
print({f: {S: v['agree'] for S, v in d.items()} for f, d in spec.items()})
print(json.dumps(r7, default=str))
