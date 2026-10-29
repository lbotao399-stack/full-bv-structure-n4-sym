"""Audit script a3: outer derivation B, non-exactness argument for G_B, inner-Delta freedom,
and the diagonal derivations of the conformal algebra.

Part 1  B as ad of a diagonal element of gl(4|4) on the blocks of [SW,(17)] (block table of ch02):
        design def:o-B / targets F5 say ad diag(1_4,0_4); ch08 def:c8-B says ad diag(0_4,1_4).
Part 2  non-exactness of G_B (ch08 prop:c8-B, design prop:o-GB (ii)): linear system for
        phi^(0) = c S^(1)(Delta) + sum_f a_f N_f (both invariant under the compact Spin(4) x SU(4)),
        using {S^(1)(xi), S^(1)(Delta)} = S^(1)([xi,Delta]) = S^(1)(xi) for xi in Sigma (lem:c2-components,
        prop:c2-special) and the four field-reweighting factors of the proof.
Part 3  derivations D of the conformal algebra with D P = kP P, D K = kK K, D Delta = z Delta, D M = mu M.
Part 4  shifts of the components under the field redefinition phi = kappa S^(1)(Delta) (type (a),
        x-dependent phi^(0), F_0 unchanged), from (s phi)^(1)(x) = -kappa e_x S^(1)(x), e_x = ad_Delta eigenvalue.
"""
import json, os
import sympy as sp

res = {}
# Part 1 -------------------------------------------------------------------------------
blocks = {  # ch02 table: subspace -> (row block, column block) of [SW,(17)]
    'translations': (2, 1), 'so(V) on S+': (1, 1), 'so(V) on S-': (2, 2), 'special conformal': (1, 2),
    'so(W)': (3, 3), 'S+ (x) W (Sigma)': (3, 1), 'S- (x) W* (Sigma)': (2, 3),
    'S- (x) W (Sigma*)': (3, 2), 'S+ (x) W* (Sigma*)': (1, 3)}
B_values = {'translations': 0, 'so(V) on S+': 0, 'so(V) on S-': 0, 'special conformal': 0, 'so(W)': 0,
            'S+ (x) W (Sigma)': 1, 'S- (x) W* (Sigma)': -1, 'S- (x) W (Sigma*)': 1, 'S+ (x) W* (Sigma*)': -1}
delta_eig = {'translations': -2, 'so(V) on S+': 0, 'so(V) on S-': 0, 'special conformal': 2, 'so(W)': 0,
             'S+ (x) W (Sigma)': -1, 'S- (x) W* (Sigma)': -1, 'S- (x) W (Sigma*)': 1, 'S+ (x) W* (Sigma*)': 1}


def ad_diag(d, rc):
    r, c = rc
    return d[r] - d[c]


d10 = {1: 1, 2: 1, 3: 0}   # diag(1_4, 0_4)
d01 = {1: 0, 2: 0, 3: 1}   # diag(0_4, 1_4)
dDelta = {1: 1, 2: -1, 3: 0}  # Delta = diag(1,-1|0)
p1 = {}
for k, rc in blocks.items():
    p1[k] = {'B (design/ch08 values)': B_values[k], 'ad diag(1_4,0_4)': ad_diag(d10, rc),
             'ad diag(0_4,1_4)': ad_diag(d01, rc), 'ad Delta from blocks': ad_diag(dDelta, rc),
             'ad Delta (ch02 table)': delta_eig[k]}
res['part1_block_rule'] = p1
res['part1_verdict'] = {
    'B == ad diag(0_4,1_4) (ch08 def:c8-B)': all(B_values[k] == ad_diag(d01, rc) for k, rc in blocks.items()),
    'B == ad diag(1_4,0_4) (design def:o-B, targets F5)': all(B_values[k] == ad_diag(d10, rc) for k, rc in blocks.items()),
    'B == -ad diag(1_4,0_4)': all(B_values[k] == -ad_diag(d10, rc) for k, rc in blocks.items()),
    'block table reproduces ad_Delta eigenvalues': all(delta_eig[k] == ad_diag(dDelta, rc) for k, rc in blocks.items())}

# Part 2 -------------------------------------------------------------------------------
c, aV, aW, ap, am, sig = sp.symbols('c a_V a_W a_lp a_lm sigma')
# four linear terms of S^(1)(xi) for xi in S+ (x) W; factors from design prop:o-GB (ii) / ch08 prop:c8-B
fac_plus = [am - aV, ap - aW, aV - ap, aW - am]
# for xibar in S- (x) W*: (Gamma(xibar,lam)|_V,A*_V) pairs with lambda_+, (Gamma|_W, A*_W) with lambda_-,
# (rho(F_V) xibar, lam*) with lam*_- (antifield of lambda_-), (rho(dA_W) xibar, lam*) with lam*_+
fac_minus = [ap - aV, am - aW, aV - am, aW - ap]
# {S^(1)(x), S^(1)(Delta)} = S^(1)([x,Delta]) = -e_x S^(1)(x); for x in Sigma, e_x = -1: coefficient +c
# equation along a block with B-value b: c + sigma*factor_i = b for every term (sigma = +-1 convention)
sols = {}
for sigma in (1, -1):
    eq_plus = [sp.Eq(c + sigma * f, 1) for f in fac_plus]
    eq_both = eq_plus + [sp.Eq(c + sigma * f, -1) for f in fac_minus]
    s_plus = sp.solve(eq_plus, [c, aV, aW, ap, am], dict=True)
    s_both = sp.solve(eq_both, [c, aV, aW, ap, am], dict=True)
    # same with c forced to 0 (translation-invariant phi^(0): only field countings)
    s_plus0 = sp.solve([e.subs(c, 0) for e in eq_plus], [aV, aW, ap, am], dict=True)
    sols[f'sigma={sigma}'] = {'S+W block only, phi0 = c S1(Delta) + sum a_f N_f': [str(s) for s in s_plus],
                              'S+W block only, c = 0 (field countings only)': [str(s) for s in s_plus0],
                              'both blocks of Sigma': [str(s) for s in s_both]}
res['part2_GB_nonexactness_system'] = sols
res['part2_translation_component'] = ('(s phi)^(1)(P_b) = -{S^(1)(P_b), c S^(1)(Delta)} = c S^(1)([P_b,Delta]) = 2c S^(1)(P_b); '
                                      'G_B^(1)(P_b) = 0 forces c = 0 since S^(1)(P_b) is not (Q+{I,-})-exact')

# Part 3 -------------------------------------------------------------------------------
kP, kK, z, mu = sp.symbols('kP kK z mu')
# [Delta,P] = -2P ; [Delta,K] = 2K ; [K,P] = Delta + R (R rotation) ; [M,M'] = M''
eqs = []
# D[Delta,P] = [D Delta,P] + [Delta,D P]:  -2 kP = -2 z - 2 kP
eqs.append(sp.Eq(-2 * kP, -2 * z - 2 * kP))
# D[Delta,K]: 2 kK = 2 z + 2 kK
eqs.append(sp.Eq(2 * kK, 2 * z + 2 * kK))
# D[K,P] = z Delta + mu R  vs (kK + kP)(Delta + R)
eqs.append(sp.Eq(z, kK + kP)); eqs.append(sp.Eq(mu, kK + kP))
# D[M,M'] = mu M'' vs 2 mu M''
eqs.append(sp.Eq(mu, 2 * mu))
res['part3_diagonal_derivations'] = str(sp.solve(eqs, [kP, kK, z, mu], dict=True))
res['part3_reading'] = ('only kappa*(1,-1,0,0) = -(kappa/2) ad_Delta on (P,K,Delta,M): a universal P_1 = kappa delta_P is '
                        'consistent with the homotopy representation only with K_1 containing -kappa delta_K and no rescaling '
                        'of the classical Delta-action (z = 0); rotations cannot be rescaled (mu = 0)')

# Part 4 -------------------------------------------------------------------------------
kap = sp.Symbol('kappa')
shift = {}
for k in blocks:
    e = delta_eig[k]
    # component change: -{(s phi)^(1)(x), O} = +kappa e_x {S^(1)(x), O}
    shift[k] = str(kap * e) + ' x (classical transformation)'
shift['Delta'] = '0'
res['part4_inner_Delta_field_redefinition_shifts'] = shift
res['part4_notes'] = [
    '(s phi)^(0) = {S, S^(1)(Delta)} = 0: arity-zero finite counterterm unchanged',
    'P_1 -> P_1 - 2 kappa delta_P, K_1 -> K_1 + 2 kappa delta_K, Q_1 -> Q_1 - kappa delta_xi (both chiralities), '
    'S_1 -> S_1 + kappa delta_eta (both chiralities), D_1, M_1, R_1 unchanged',
    'twisted: r_xi(s phi) = Q_0^xi(kappa S^(1)(Delta)): [Q_1^xi] unchanged',
    'reflection: ad_Delta commutes with dagger_n (dagger_n preserves Sigma and Sigma*), so real kappa is reflection '
    'compatible; it is absent from the freedom tables of the design (field redefinitions listed only for the five '
    'field-counting phi and the arity-one phi^(1))']

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a3_outer_B_freedoms.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
