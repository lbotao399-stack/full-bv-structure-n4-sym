"""Engine A, part Q: Weinberg-Srednicki two-component forms of the reported formulas (Def c8-twocomp, Rem qs-twocomp and
step2/notation/dictionary.json: (rho(w) xi)^alpha_j = w_ji xi^{alpha i}, (rho(v) xi)^{alphadot i} = -i v_m sigmabar^{m alphadot alpha}
xi_alpha^i, (xi, eta) = xi_alpha^i eta^alpha_i - xibar_{alphadot i} etabar^{alphadot i}, iota_w A = 1/4 w^{ij} A_ij,
(w1, w2) = 1/4 w1^{ij} w2_ij, xi^{alpha i} = eps^{alpha beta} xi_beta^i, eps^{12} = +1).

Orientation (revision 2026-10-30, finding Q13): under the identity on vector indices the engine's lam_+ (+1 eigenspace of
rho(e0 e1 e2 e3) in cliff.py) is the dictionary's S_- (x) W^* (dotted), i.e. route B's lam_-, and the engine's
"xi in S_+ (x) W" is the dictionary's xibar_{alphadot i}.  The theory is invariant under the reflection of one axis of V
combined with W <-> W^* (exact in route B's code, qa/out.nosync/mirror_B.json), and engine A finds the same numbers for
both chiralities (TQ4); the forms below are therefore written in the dictionary orientation, for xi_alpha^i (S_+ (x) W)."""


def forms(a1, a2, a3, a4, cQ, cocycle=None):
    out = {
        'orientation': 'dictionary orientation (S_+ = undotted for dx^1 dx^2 dx^3 dx^4); engine A labels mirrored, numbers identical by the reflection symmetry',
        'O': 'int omega (iota_w1 A, iota_w2 A) = 1/16 int omega (w1^{ij} A_ij, w2^{kl} A_kl)',
        'F1': 'int omega [(iota_w1 A, (rho(w2) xi, lam)) + (1<->2)] = 1/4 int omega [(w1^{kl} A_kl, w2_{ji} xi^{alpha i} lam_alpha^j) + (1<->2)]',
        'F2': '(w1,w2) int omega (rho(A|_perp) xi, lam) = 1/4 w1^{ij} w2_ij int omega (A_{ji}, xi^{alpha i} lam_alpha^j)',
        'F3': '(w1,w2) int omega (rho(A|_V) xi, lam) = 1/4 w1^{ij} w2_ij int omega i (lambdabar_{alphadot i}, A_m sigmabar^{m alphadot alpha} xi_alpha^i)',
        'F4': '(w1,w2) int omega (c, (xi, lam^*)) = 1/4 w1^{ij} w2_ij int omega (c, xi_alpha^i lam^{*alpha}_i)',
        'classical': '{S^(1)(xi), O} = +F1  (A_ij -> +Gamma(xi,lam)_ij, Gamma_kl(xi, lam) = -2 eps_{jikl} (xi^j lam^i), dictionary entry Gamma; Cor c7-Q)',
        'Q1_two_body_psi0': 'Q_1(O)|_two-body = (%s) F1 + (%s) F2 + (%s) F3 + (%s) F4  units u = -g^2 N/(32 pi^2)  (psi = 0)' % (a1, a2, a3, a4),
        'c_Q': 'Q_1(O)|_traceless = c_Q {S^(1)(xi), O},  c_Q = %s - theta_perp - p_perp' % cQ,
        'trace_psi0': 'Q_1(O)|_trace, two-body (psi = 0) = (%s) (1/4 w1^{ij} w2_ij) int omega [(A_ji, xi^{alpha i} lam_alpha^j) + i (lambdabar_{alphadot i}, A_m sigmabar^{m alphadot alpha} xi_alpha^i)]' % a2,
        'Konishi': 'Q_1( int omega (u'"'"'^alpha lam_alpha^i, A_ij z^j) ) contains 0 x int omega (lambdabar_{alphadot i} z^i, lambdabar^{alphadot}_j z^j): c_K = 0 (psi = 0)',
        'bc_lamstar': 'Q_1( int omega (c, u_alpha lam^{*alpha}_4) ) on lambdabar_{alphadot 4} lambdabar^{alphadot}_4: 0',
        'betabeta': 'Q_1( int omega (u'"'"'^alpha lam_alpha^I, u'"'"'^beta lam_beta^J) ) on (lambdabar_{alphadot 4}, d^{alphadot} A_{K4}), K the complement of I, J: '
                    'nonzero, holomorphic derivatives d_alphadot = u'"'"'^alpha sigma^m_{alpha alphadot} d_m only, the same coefficient for both colour orders '
                    '(per unit weight-vector letters 2/3 on d^{alphadot} gamma^K lambdabar_{alphadot 4}, modulo the Q^xi_0-exact Dirac component)',
        'bgamma': 'Q_1( int omega (u'"'"'^alpha F_mn (sigma^m sigmabar^n)_alpha^beta u'"'"'_beta, A_{I4}) ) on (lambdabar_{alphadot 4}, d^{alphadot} A_{I4}): '
                  'nonzero, holomorphic derivatives only; modulo Im Q^xi_0 the two colour orders have unequal coefficients containing log(3/2) '
                  '(engine weight-vector normalisation 64/3 - 128 log(3/2) and 32 + 64 log(3/2)); no reduction modulo total derivatives',
        'eta_killing': 'eps_alpha^i(x) = -i x_m sigma^m_{alpha alphadot} etabar^{alphadot i},  epsbar_{alphadot i}(x) = i x_m eta^alpha_i sigma^m_{alpha alphadot}; '
                       'intrinsic vertex s0 = (A^{kj} eps_{alpha beta} eta^beta_j, lam^{*alpha}_k) - (A_kj eps_{betadot alphadot} etabar^{alphadot j}, lambdabar^{*betadot k})',
        'K_b': 'v^m = x^2 b^m - 2 (b.x) x^m; intrinsic z-polynomial vertex of Rem qs-twocomp with z = x - x_O',
    }
    if cocycle is not None:
        out['cocycle_scheme'] = ('psi(O) = (%s) u 1/4 w1^{ij} w2_ij int omega (A_m, A_m), psi(F1) = (%s) u F3; then Q_1(O)|_two-body = (%s) F1 + (%s) F2: '
                                 'the term i (lambdabar_{alphadot i}, A_m sigmabar^{m alphadot alpha} xi_alpha^i) drops out' %
                                 (cocycle['kappa'], cocycle['kappa_prime'], a1, a2))
    return out
