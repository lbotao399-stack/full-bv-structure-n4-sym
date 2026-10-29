"""Builds predictions.json for the two-letter engines (twoletter_theory.tex, Section sec:tl-pred).
All numbers exact (fractions / sympy).  Units: u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)); for su(N) with
(X,Y) = -(4 hbar/g^2) tr(XY): u = -g^2 N/(32 pi^2), so g^2 N/(8 pi^2) = -4 u.
Spectrum dictionary (Beisert hep-th/0307015 = B03c, step2/literature/n4_untwisted.json):
  D(g) = D_0 + (g^2 N/8 pi^2) H,  H = sum_l H_{l,l+1},  H_12 = 2 h(J_12)        [B03c (1.1)-(1.3)]
  anomalous dimension gamma = (g^2 N/8pi^2) E_H = -4 E_H u = 2 eig(lc^natural)  =>  eig(lc^natural) = -2 E_H u,
  D_1 eigenvalue (Delta = 2E normalisation, Thm c6-dil) = 2 gamma = -8 E_H u  (at theta = -1/3, prediction P7)."""
import json, os, itertools
from fractions import Fraction as Fr
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))


def h(j):
    return sum((Fr(1, k) for k in range(1, j + 1)), Fr(0))


def fs(x):
    return str(x)


# ---------------- su(2) Heisenberg chain on cyclic single traces ----------------
def heisenberg_cyclic(L, M):
    """H = sum_l (1 - P_{l,l+1}) (periodic) on the cyclic (trace) subspace of words with L-M letters Z and M letters X.
    Basis: necklaces (orbits under rotation); the trace identifies rotations (bosonic letters)."""
    words = [w for w in itertools.product('ZX', repeat=L) if w.count('X') == M]
    def canon(w):
        return min(tuple(w[i:] + w[:i]) for i in range(L))
    neck = sorted(set(canon(w) for w in words))
    idx = {n: i for i, n in enumerate(neck)}
    # state |n> = sum of the distinct rotations of n (trace); H acts on words, project back
    Hm = sp.zeros(len(neck), len(neck))
    for j, n in enumerate(neck):
        rots = set(tuple(n[i:] + n[:i]) for i in range(L))
        for w in rots:
            for l in range(L):
                a, b = l, (l + 1) % L
                # (1 - P) w
                w2 = list(w); w2[a], w2[b] = w2[b], w2[a]; w2 = tuple(w2)
                # contribution to coefficient of rotation-class: count words, divide by orbit size later
                Hm[idx[canon(w)], j] += 1
                Hm[idx[canon(w2)], j] -= 1
    # convert word counts to coefficients of the orbit sums: coefficient of |n'> = (count of words in orbit n')/|orbit n'|
    for i, n in enumerate(neck):
        orb = len(set(tuple(n[k:] + n[:k]) for k in range(L)))
        for j in range(len(neck)):
            Hm[i, j] = sp.Rational(Hm[i, j], orb)
    x = sp.Symbol('E')
    cp = sp.factor_list(Hm.charpoly(x).as_expr())
    factors = []
    roots = []
    for fac, mult in cp[1]:
        pol = sp.Poly(fac, x)
        rts = [complex(r) for r in pol.nroots(n=30)]
        for r in rts:
            assert abs(r.imag) < 1e-20
        exact = [str(sp.nsimplify(r)) for r in sp.roots(pol, x).keys()] if pol.degree() <= 2 else []
        factors.append(dict(factor=str(fac), multiplicity=int(mult), exact_roots=exact))
        roots += sorted([r.real for r in rts]) * int(mult)
    return dict(charpoly_factors=factors, roots=sorted(roots))


def bethe_two_magnon(L):
    """two-magnon cyclic states, Bethe momenta p = 2 pi n/(L-1): E = 8 sin^2(p/2), n = 0..floor((L-2)/2) (p = pi excluded)"""
    return sorted([sp.nsimplify(sp.simplify(8 * sp.sin(sp.pi * n / (L - 1)) ** 2)) for n in range(0, (L - 2) // 2 + 1)],
                  key=lambda x: float(x))


def sl2_matrix(n):
    """B03c (3.14): H'_12 on (a1^dag)^k (a2^dag)^{n-k}|00>, k = 0..n: diag h(k)+h(n-k), off-diagonal -1/|k-k'|"""
    M = sp.zeros(n + 1, n + 1)
    for k in range(n + 1):
        for kp in range(n + 1):
            M[kp, k] = (sp.Rational(h(k) + h(n - k))) if k == kp else -sp.Rational(1, abs(k - kp))
    return M


def main():
    P = {}
    P['description'] = __doc__
    P['units'] = dict(unit='hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y))', su_N='-g^2 N/(32 pi^2)',
                      gYM2N_over_8pi2_in_units='-4', gamma_from_lc_natural='gamma = 2 eig(lc^natural)',
                      D1_from_lc_natural='D_1 = 4 lc^natural + Xi_Delta (Thm thm:tl-D1); Xi_Delta = 0 on bilinears at theta=-1/3')
    # P1: Ch5 graph-by-graph and totals (scheme specific eps^-1, universal log)
    P['P1_scalar_bilinears_Ch5'] = dict(
        operator='O = int omega (iota_{w1}A, iota_{w2}A)',
        graphs_log_units={'(e) quartic bubble, outputs A_perp A_perp': '2 O - 2 (w1,w2) int omega (A_perp,A_perp)',
                          '(f) triangle, outputs A_perp A_perp': '2 O',
                          '(c) quartic bubble, outputs A_V A_V': '-2 (w1,w2) int omega (A_V,A_V)',
                          '(d) triangle, outputs A_V A_V': '+2 (w1,w2) int omega (A_V,A_V)',
                          '(b) cubic, output A_V': '0 (colour)'},
        lc_units=dict(traceless='4 O', trace='-8 O'),
        eps_minus_1='2 hbar dim g (w1,w2)/(16 pi^2) int omega  (k = 0 self-contraction, scheme specific)',
        source='Prop c5-table, Thm c5-Oeps')
    # P2: Callan-Symanzik operator (anomalous dimensions)
    P['P2_lc_natural'] = dict(
        rule='lc^natural(O) = lc(O) - u {F, O},  {F,-} multiplies c,A,lambda,A*,lambda*,c* by 1,2,4,-2,-4,-1',
        values_units={'traceless bilinear (20prime)': '0', 'trace bilinear (Konishi)': '-12 O',
                      'letter A_perp': '-2', 'letter A_V (mod terms from sO)': '-2', 'letter lambda': '-4', 'letter c': '-1'},
        anomalous_dimension_units={'20prime': '0', 'Konishi': '-24 = 3 g^2 N/(4 pi^2)', 'A_perp': '-4', 'A_V': '-4',
                                   'lambda': '-8', 'c': '-2'},
        source='Thm c6-dil, Thm thm:tl-brst')
    # P3: Ch9 dilatation and translation components
    P['P3_D1_P1_Ch9'] = dict(
        D1_traceless='-(1/3+theta) u {S1(Delta),O}', D1_trace='-48 u O - (1/3+theta) u {S1(Delta),O}',
        theta_fixed='-1/3 (Thm c9-kappa, arity two along (Delta,P_b))',
        P1='-(1/3+theta_P) u {S1(P_b),O}',
        composite_triangle='C2_Delta(pair A_perp A_perp, no derivatives) = -2 u {S1(Delta),O}',
        composite_two_point_perp='Pi_Delta|perp = u [5/3 int(v.A_perp, A*_perp) - 2/3 int(A_perp, A*_perp)]',
        one_body_D1_on_A_perp_letter='u[-8 X + (5/3 - theta) v.X]',
        F1_Delta_quadratic='u[22/3 int(A,A*) - 47/3 int(lambda*,lambda) - 8/3 int(c*,c)] + theta-direction',
        Ob1_Delta_quadratic='u[-11/3 int(dA,dA) + 47/3 int(lambda,dslash lambda) - 14/3 int(dc,A*)]',
        source='Thm c9-quadratic, Thm c9-D1, Prop c9-triangle, Thm c9-P1, Cor thm:tl-D1bil')
    # P4: one-body Ward identities (scheme independent, exact)
    P['P4_W1_gates'] = dict(
        statement='lc(s N) = {L, N} for one-letter N = int omega <zeta, X>; L = u s(F)',
        multipliers={'-1/4(dA,dA)': 4, '1/2(lambda,dslash lambda)': 8, '(dc,A*)': -1, 'cubic gauge': 6,
                     '1/2(lambda,(dslash_A - dslash)lambda)': 10, '([A,c],A*)': 1, '([c,lambda],lambda*)': 1,
                     '1/2([c,c],c*)': 1, 'quartic': 8},
        cases={'A_perp': 'lc({I,N}) = u {I,N}  (pair (A_perp,c), output A_perp c, W=1)',
               'lambda': 'lc({I,N}) = u {I,N}  (pair (lambda,c))',
               'c': 'lc({I,N}) = u {I,N}  (pair (c,c), outputs c c, W=0)',
               'A_V': 'lc({I,N}) = u ({I,N} - QN), QN = -int omega <zeta, dc>  (pair (A_V,c): one-field output dc)',
               'A*_perp': 'lc({I,N}) = u (4 QN + 6 {I_AAA,N} + 8 {I_AAAA,N} + 10 {I_llA,N} + {I_AcA*,N})'},
        two_body='K^log(s1 X_a, X_b) + (-1)^{|X_a|} K^log(X_a, s1 X_b) = s K^log(X_a,X_b)',
        source='Thm thm:tl-brst, Cor cor:tl-ward, Prop prop:tl-W1')
    # P5: spectrum, Beisert
    hs = {str(j): fs(h(j)) for j in range(0, 11)}
    two_site = {str(j): fs(2 * h(j)) for j in range(0, 11)}
    tw = {}
    for S in range(0, 12, 2):
        EH = 4 * h(S)
        tw[str(S)] = dict(E_H=fs(EH), lc_natural_units=fs(-2 * EH), gamma_units=fs(-4 * EH),
                          gamma_gYM2N_over_2pi2=fs(h(S)), D1_units=fs(-8 * EH))
    P['P5_spectrum_Beisert'] = dict(
        dictionary='eig(lc^natural) = -2 E_H units, gamma = -4 E_H units = (g^2 N/8 pi^2) E_H, D_1 = -8 E_H units',
        harmonic_numbers=hs, two_site_energies_2h_j=two_site,
        twist_two=dict(statement='length 2, spin S even: E_H = 4 h(S), deltaD = g^2 N h(S)/(2 pi^2) [B03c (2.31)]', values=tw),
        konishi=dict(E_H='6', lc_natural_units='-12', gamma_units='-24', D1_units='-48'),
        half_BPS=dict(E_H='0'),
        scheme='spectrum: independent of scheme, theta, finite parts (Thm c8-invariants (iii)); B03c regulator: point splitting',
        verified_sectors='lc^natural spectra verified on: twist two N <= 8 (both engines, referee), su(2) L <= 8 (both engines, '
                         'referee), Konishi, su(1|1) two-site density (engine A gate G8f, referee); fermion families beyond the '
                         'su(1|1) two-site density and field-strength families: not tested',
        D1_status='eig(D_1) = -8 E_H u: derived on the planar scalar sector at theta_perp = -1/3 (Prop prop:tl-scalarsector); '
                  'elsewhere conditional on Conjecture P7 (conj:tl-P7)',
        source='B03c (1.1)-(1.3), (2.31), (4.27)-(4.29); step2/literature/n4_untwisted.json')
    # su(2) Heisenberg
    heis = {}
    for L in range(2, 9):
        for M in range(0, L // 2 + 1):
            ev = heisenberg_cyclic(L, M)
            heis['L=%d,M=%d' % (L, M)] = ev
    bethe = {str(L): [str(x) for x in bethe_two_magnon(L)] for L in range(4, 9)}
    # consistency: two-magnon Bethe energies appear in the M = 2 cyclic spectra
    cons = {}
    for L in range(4, 9):
        spec = heis['L=%d,M=2' % L]['roots']
        cons[str(L)] = all(any(abs(float(b) - s_) < 1e-12 for s_ in spec) for b in bethe_two_magnon(L))
    P['P6_su2_Heisenberg'] = dict(
        statement='su(2) sector Tr(Z^{L-M} X^M): H = sum_l (1 - P_{l,l+1}) (H_12 = 2h(1) on the antisymmetric pair); '
                  'eig(lc^natural) = -2 E units',
        cyclic_spectra_E=heis, two_magnon_bethe_E=bethe, bethe_in_spectra=cons,
        konishi_descendant='L=4,M=2 contains E=6 (Tr[Z,X][Z,X])')
    # sl(2) two-site matrix (B03c (3.14))
    sl2 = {}
    for n in range(0, 6):
        Mx = sl2_matrix(n)
        ev = Mx.eigenvals()
        evl = sorted([sp.nsimplify(x) for x, m in ev.items() for _ in range(m)], key=float)
        expect = sorted([sp.Rational(2 * h(j)) for j in range(0, n + 1)], key=float)
        sl2[str(n)] = dict(matrix=[[str(Mx[i, j]) for j in range(n + 1)] for i in range(n + 1)],
                           eigenvalues=[str(x) for x in evl], equals_2h_j=bool(evl == expect))
    P['P6b_sl2_two_site'] = dict(
        statement="B03c (3.14): H'_12 (a1^dag)^k (a2^dag)^{n-k}|00> = (h(k)+h(n-k)) - sum_{k'!=k} 1/|k-k'|; "
                  "matrix-level prediction (conditional): on planar single traces, modulo s-exact terms, "
                  "lc^natural = -2 sum_l H_{l,l+1} units in the letter basis of B03c",
        levels=sl2)
    # P7 prediction on composite part
    P['P7_composite_part'] = dict(
        statement='Conjecture conj:tl-P7: the parameter-independent part of Xi_Delta = -{F_Delta prime,-} + C_Delta induces '
                  '-(1/3) u {S1(Delta),-} on cohomology; Xi_Delta = -u[(1/3){S1(Delta),-} + sum_{s=perp,lambda,V} theta_s {S_s,-}], '
                  's = perp, lambda, V; at theta_s = -1/3: Spec D_1 = Spec 4 lc^natural = 2 x anomalous dimensions',
        status='proved on the scalar bilinears and on the planar scalar sector (letters A_perp without derivatives; '
               'Prop prop:tl-scalarsector, inputs Thm c9-D1, Prop c9-triangle, Thm c9-P1); conjecture elsewhere',
        planar_scalar_sector=dict(C2_Delta_pair='-u v.(X_a X_b) mod cyclically trivial', C2_Pb_pair='-u d_b(X_a X_b)',
                                  Xi_Delta='-(1/3+theta_perp) u v.', P1='-(1/3+theta_P) u {S1(P_b),-}',
                                  D1_pair_density='-8u(1 - P + K/2) - (1/2)(1/3+theta_perp) u v.(X_a X_b)'),
        open_question='whether the cubic part of eq c9-arity1 forces theta_perp = theta_lambda = theta_V')
    P['P9_W2'] = dict(statement='two-body Ward identity (Cor cor:tl-ward): K(s1 X_a, X_b) + (-1)^{|X_a|} K(X_a, s1 X_b) = s K(X_a,X_b)',
                      gate='engine A G3-W2: full identity lc(sM) = s lc(M) + {L,M} on M = tr(zeta X_a X_b), 741/741, 232 nontrivial')
    # P8 scheme facts
    P['P8_scheme'] = dict(
        log_universal='log eps coefficients independent of the cutoff shape h (Prop prop:tl-shape)',
        examples={'bubble N=|l|^4': dict(per_edge='|p|^2/8 eps^-1 + 1/4 eps^-2', total_time='|p|^2 eps^-1 + 3 eps^-2', log='0'),
                  'triangle N=|l|^4': dict(per_edge_eps_m1='1/3', total_time_eps_m1='3', log='equal (-31/25 at the sample)')},
        special_edge='no log eps; finite part -int_{R_s} F_0 log(alpha_s/min_{f!=s} alpha_f) (Lemma lem:tl-special); '
                     'toy: -1/(2 eps) - 3|p|^2/8',
        k0='self-contraction: pure eps^{-(w_a+w_b)/2}, (d^gamma S_eps)(0) = (-1)^{m+1}/(16 pi^2 m 2^{m-1} eps^m) x pairings')
    with open(os.path.join(HERE, 'predictions.json'), 'w') as f:
        json.dump(P, f, indent=1)
    print('predictions.json written;', 'sl2 levels equal 2h(j):', all(v['equals_2h_j'] for v in sl2.values()),
          '; Bethe in spectra:', cons)
    print('L=4,M=2', heis['L=4,M=2'], '\n L=6,M=2', heis['L=6,M=2'])


if __name__ == '__main__':
    main()
