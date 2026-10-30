"""Collect the referee computations into findings.json (run after the other scripts).

Inputs (this directory): xspace_J_dps20.json, xspace_R_dps20.json, xspace_R_dps26.json, localisation_dps20.json,
localisation_ii_dps20.json, selfenergy_dps22.json, fourier_opct_fourier.json, fourier_opct_fourier3.json,
fourier_opct_algebra.json.
Also assembles c_S, c_T, d_T of thm:result from the referee's OWN ingredient values (J, R constants from the
x-space fits, c_loc numeric, Sigma_ren constant from the direct momentum integrals, Fourier entries), using the
class formulas whose combinatorial factors were re-derived in fourier_opct.py (EX, QU) and selfenergy.py.
"""
import json
import os

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 20
ld = lambda n: json.load(open(os.path.join(HERE, n)))
J = ld('xspace_J_dps20.json'); R = ld('xspace_R_dps20.json'); R26 = ld('xspace_R_dps26.json')
LOC = ld('localisation_dps20.json'); LOC2 = ld('localisation_ii_dps20.json')
SE = ld('selfenergy_dps22.json'); FO = ld('fourier_opct_fourier.json'); FO3 = ld('fourier_opct_fourier3.json')
AL = ld('fourier_opct_algebra.json')
g, L2, L3, pi = mp.euler, mp.log(2), mp.log(3), mp.pi
M = mp.mpf

bJ = M(J['fits']['6-param']['b_const']); aJ = M(J['fits']['6-param']['a_logeps'])
bR = M(R['fits']['6-param']['b_const']); aR = M(R['fits']['6-param']['a_logeps'])
cl = M(LOC['c_loc']['momentum_radial'])
cSig = M(SE['direct_momentum']['Sigma_ren_const_minus_2g_2log2'])     # constant of Sigma_ren/(g^2N/16pi^2)/p^2 at p=1 minus 2g+2log2

# refits with richer bases (no claimed values used)
def refit(rows, xkey, ykey, basis):
    X = [[f(M(r[xkey])) for f in basis] for r in rows]
    Y = [M(r[ykey]) for r in rows]
    A = mp.matrix(X); y = mp.matrix(Y)
    return mp.lu_solve(A.T * A, A.T * y)
solB = refit(SE['B_hankel']['rows'], 'eps', 'b_hankel',
             [lambda e: mp.log(e), lambda e: 1, lambda e: e, lambda e: e * mp.log(e), lambda e: e * e, lambda e: e * e * mp.log(e)])
solL2 = refit(LOC2['rows'], 'eps', 'LHSii_minus_J_over_D0sq',
              [lambda e: 1, lambda e: e * mp.log(e)**2, lambda e: e * mp.log(e), lambda e: e])

# --- assembly at |x| = 1 (L = 0), in units, relative to the tree, beta_2 kept apart ---
SEc = 8 * g - 24 * L2 - 4 * cSig                     # two lines: 2 * (-2) * [-2(2g - 2log2) + 2g + 2log2 + cSig]
EXc = 4 * bJ - 8 * bR - 64 * pi**2 * cl              # EX = 4 log eps + const ; 32pi^2[J/D0^2 - 4R/D0 - 2 c_loc]
QUc = 4 * bJ                                         # QU = (2 tau - 1)(-4 log eps + QUc)
cS = SEc + EXc - QUc                                 # tau = 0 part (log eps cancels: 4 + 4 - 8 = 0)
cT = 2 * QUc
assembled = dict(SE=mp.nstr(SEc, 12), EX_const=mp.nstr(EXc, 12), QU_const=mp.nstr(QUc, 12),
                 c_S_without_beta2=mp.nstr(cS, 12), c_S_claimed=mp.nstr(mp.mpf(27) / 2 + 4 * mp.log(mp.mpf(3) / 2), 12),
                 c_T=mp.nstr(cT, 12), c_T_claimed=mp.nstr(8 * (g - 3 * L2), 12),
                 SE_claimed=mp.nstr(mp.mpf(43) / 2 + 8 * g - 24 * L2, 12),
                 EX_claimed=mp.nstr(-8 - 4 * g + 8 * L2 + 4 * L3, 12),
                 QU_claimed=mp.nstr(4 * g - 12 * L2, 12),
                 log_eps_coefficients=dict(J_fit_a=mp.nstr(aJ, 10), R_fit_a=mp.nstr(aR, 10),
                                           EX=mp.nstr(4 * aJ - 8 * aR, 10), QU_times_minus_2tau_plus_1=mp.nstr(4 * aJ, 10),
                                           note='EX log eps coefficient = 4 a_J - 8 a_R = 4 ; QU: (2tau-1)*4 a_J = -(2tau-1)*4 ; OC: -8(1-tau) (G2) ; sum 0'))
print(json.dumps(assembled, indent=1))

def agree(x, y, tol):
    return abs(M(x) - M(y)) <= tol

F = []
def add(i, claim, method, found, claimed, ok, sev, comment):
    F.append(dict(id=i, claim=claim, method=method, value_found=found, value_claimed=claimed,
                  agreement=ok, severity=sev, comment=comment))

add('R1-J', 'J(x)=D0^2/(8pi^2)(log(x^2/eps)+gamma_E-3log2)+o(1)',
    'direct x-space quadrature (two-centre (r,s) reduction, mpmath dps 20), |x|=1, eps=1e-3..1e-5 (7 values); free fit a log eps+b+c eps log eps+d eps(+eps^2 terms)',
    dict(a=mp.nstr(aJ, 12), b=mp.nstr(bJ, 12), diff_at_eps_1em5=J['rows'][-1]['diff'], diff_over_eps_at_eps_1em5=mp.nstr(M(J['rows'][-1]['diff']) / M('1e-5'), 8)),
    dict(a=-1, b=J['b_claimed']), agree(bJ, J['b_claimed'], 1e-7), 'none',
    'difference to the claimed asymptotics is -6 eps + O(eps^2 log eps) (smaller than the allowed O(eps log eps)); 6-parameter fit reproduces the constant to 2e-8.')
add('R2-R', 'R(x)=D0/(16pi^2)(log(x^2/eps)+1+gamma_E-3log2)+o(1)', 'as R1, also at dps 26 (identical to 18 digits)',
    dict(a=mp.nstr(aR, 12), b=mp.nstr(bR, 12), diff_at_eps_1em5=R['rows'][-1]['diff'], dps26_b=R26['fits']['6-param']['b_const']),
    dict(a=-1, b=R['b_claimed']), agree(bR, R['b_claimed'], 1e-7), 'none',
    'difference is -4 eps + O(eps^2 log eps).')
add('R3-cloc', 'c_loc=int d^4k/(2pi)^4 e^{-2k^2}(e^{-k^2}-1)/k^4=log(2/3)/(16pi^2)',
    'radial momentum quadrature and, independently, the x-space form int d^4w D_1(w)[D_2(w)-D_1(w)]',
    dict(momentum=LOC['c_loc']['momentum_radial'], position=LOC['c_loc']['position_space']), LOC['c_loc']['claimed'], True, 'none', 'agreement to 21 digits.')
add('R4-loc-i', 'Lemma lem:local (i): LHS -> D0[R + D0 c_loc]',
    'exact reduction LHS-D0R = Term1 + Term2 (heat semigroup K*D_eps = D_2eps; harmonic multipoles l=1,2 of D0 around x, Hobson identity), 2d quadratures, eps=1e-2..1e-5',
    dict(rows=[(r['eps'], r['LHS_minus_D0R_over_D0sq'], r['difference']) for r in LOC['lemma_i']['rows']], Term1_over_D0sq_at_eps_1em5=LOC['lemma_i']['rows'][-1]['Term1_over_D0sq']),
    LOC['lemma_i']['claimed'], True, 'none',
    'Term1 equals D0^2 c_loc up to exponentially small terms (mean-value property of the harmonic D0); the remainder Term2 ~ -0.02 eps log eps -> 0 (2.4e-6 relative to 2.6e-3 at eps=1e-5). Terms l>=3 of the multipole series are neglected (O(eps^3 log eps)).')
add('R5-loc-ii', 'Lemma lem:local (ii): LHS -> J + 2 D0^2 c_loc',
    'LHS-J = <P,(K*P)-P>, P=D(x-z)D(z), symmetric half-space reduction + multipoles l<=2 + heat semigroup; eps=1e-2..1e-5',
    dict(rows=[(r['eps'], r['LHSii_minus_J_over_D0sq'], r['difference']) for r in LOC2['rows']], fit_limit_with_eps_log2eps=mp.nstr(solL2[0], 10)),
    LOC2['claimed'], True, 'none', 'difference behaves like eps log^2 eps (-4.4e-6 at eps=1e-5) -> 0.')
pA = SE['partA']
add('R6-G1', 'G1: singular parts of the scalar self-energy graphs = ch04 table (2, 5/2, -1/4, -2 for eps^-1; 0, 0, 1, -2 for log eps)',
    'own Feynman rules: vertex functions by polarisation of -S_E (su(2), su(3)), seagull contractions, 10d chiral trace; singular parts from the direct momentum integrals (free fit); conversion Sigma=-2e g^2N/16pi^2 derived from eq:c4-pw (factor 2) and (X,Y)=-(4hbar/g^2)tr',
    dict(seagull_factors=[pA['su3_seagull_V_factor_over_minus_N'], pA['su3_seagull_Wperp_factor_over_minus_N']],
         bubble=SE['direct_momentum']['fits']['bubble'], gaugino=SE['direct_momentum']['fits']['gaugino'],
         Sigma_sing_total=dict(inv_eps=SE['direct_momentum']['total_inv_eps'], log_eps=SE['direct_momentum']['total_log_eps'])),
    dict(Sigma_units_gsqN_16pi2=dict(seagull_V=-4, seagull_perp=-5, bubble=[0.5, -2], gaugino=[4, 4], total=[-4.5, 2])), True, 'none',
    'all eight singular coefficients agree; total -9/(2 eps) + 2 p^2 log eps = table 9/4, -1.')
add('R7-fermion-sign', 'sign and normalisation of the gaugino loop (Majorana-Weyl factor 1/2 x tr_16)',
    'tr_16(G^I G^K G^J G^L)=-16 d^IJ d^KL with explicit 32x32 gammas; supersymmetry test: with an identical cutoff on every loop momentum the quadratic divergences must cancel',
    dict(chiral_trace_maxdev=pA['chiral_trace_tr16(GI GK GJ GL)+16dIJdKL_maxdev'],
         quadratic_divergences=pA['quadratic_divergence_identical_cutoff_seagulls_bubble_gaugino']),
    'gaugino numerator 8 g^2N l.(l-p)', True, 'none', 'seagulls -9, bubble +1, gaugino +8: sum 0 -> sign and factor 1/2 confirmed independently of Euclidean Majorana conventions.')
add('R8-B', 'M[0,0,2]=1-gamma_E-log2-log(eps p^2) (bubble master)',
    'one-dimensional position-space Hankel transform of D_eps^2, mpmath quadosc, eps=1e-2..1e-6',
    dict(diff_rows=[(r['eps'], r['diff']) for r in SE['B_hankel']['rows']], refit_log=mp.nstr(solB[0], 12), refit_const=mp.nstr(solB[1], 12)),
    SE['B_hankel']['claimed_const'], True, 'none', 'difference = 7.5 eps + O(eps^2) exactly; the constant is confirmed (refit const within 1e-9).')
add('R9-SigmaRen', 'Sigma_ren(p) = (g^2N/16pi^2) p^2 (2 log p^2 + 2gamma_E + 2log2 - 43/8)',
    'referee reduction Sigma/(g^2N) = -9 T - 2 p^2 B + 9 C with C exact, plus direct numerical momentum integrals of bubble and gaugino at 8 eps (1e-4..2e-2) and free 7-parameter fit',
    SE['direct_momentum']['Sigma_ren_const_minus_2g_2log2'], '-5.375 (= -43/8)', True, 'none',
    'analytic: 9C = 9/(32pi^2 eps) - 27 p^2/(128 pi^2) + O(eps), -2 p^2 b = p^2(2 log(eps p^2)+2g+2log2-2) -> -2-27/8 = -43/8. Numeric agreement 2e-8.')
add('R10-pergraph', 'Prop prop:selfenergy per-graph finite constants: bubble +167/72, gaugino -13/9 (with the -2/+4 (log(eps p^2)+g+log2) terms)',
    'direct numerical momentum integrals (R9)',
    dict(bubble_const=SE['direct_momentum']['fits']['bubble']['const'], gaugino_const=SE['direct_momentum']['fits']['gaugino']['const'],
         referee_exact='bubble 13/8, gaugino -7'),
    dict(bubble=SE['direct_momentum']['predictions']['bubble_const_note(167/72 - 2g - 2log2)'], gaugino=SE['direct_momentum']['predictions']['gaugino_const_note(-13/9 + 4g + 4log2)']),
    False, 'minor',
    'The per-graph constants printed in the proposition are wrong: numerically the bubble constant is 13/8 (=117/72) and the gaugino constant -7, matching the referee reduction (bubble = 2p^2B + C, gaugino = 4(2C - p^2B)). The printed ones sum to 7/8, not to the stated total -43/8 (the -13/9 is only the -8*13/72 part from M[1,1,4]; the 16 M[0,0,3] part -50/9 is missing; 167/72 exceeds 13/8 by 25/36). The total -43/8, which is all the theorem uses, is correct.')
fo = {}
fo.update(FO['fourier']); fo.update(FO3['fourier'])
add('R11-fourier', 'Lemma lem:fourier: FT of log q^2, log^2 q^2, 1/q^2, log q^2/q^2 at x != 0',
    '4d Hankel transform with exponential (Abel) damping e^{-eta q}, Neville extrapolation eta->0, |x|=1,2,3',
    {k: v['extrapolated'] for k, v in fo.items()}, {k: v['claimed'] for k, v in fo.items()}, True, 'none',
    'agreement to 1e-9..1e-11 at |x|=3, 7e-7..2e-10 at |x|=2; at |x|=1 the extrapolation is less converged (rel. 2e-5..2e-7, consistent with the scatter of the one-point-less extrapolant). The |x| dependence (log x^2 terms) is confirmed by the three radii.')
G2 = AL.get('G2', {})
add('R12-G2', 'ch05 thm:c5-Oeps / prop:c5-table: (e) unit log eps [2O - 2(a,b)K], (f) unit log eps [2O]',
    'own vertex functions (polarisation of -S_E), explicit su(2), su(3), random complex a,b; least-squares projection on O_ab and (a,b)K over many external legs; log coefficient -(1/16pi^2)',
    G2, '(e): O 2, (a,b)K -2; (f): O 2, (a,b)K 0', True, 'none', 'agreement to machine precision; colour factor scales as N (no 1/N).')
CB = AL.get('combinatorics', {})
add('R13-EX-QU', 'EX and QU vertex, colour, flavour and symmetry factors: EX/tree = -g^2N x [(k+kp).(2q-k-kp) integral]/D0^2 ; QU/tree = -g^2N (2T-S)/S x J/D0^2',
    'symmetrised vertex functions of O and -S_E contracted with explicit |Aut| (tree 2, EX 2, QU 4), random momenta and flavour vectors, su(2), su(3)',
    CB, dict(EX='-N', QU='-N(2T-S)/S'), True, 'none',
    'confirms the two labelled exchange graphs give S12 with colour N dim g and relative factor -g^2N; the quartic factor 2T-S and symmetry factor.')
add('R14-assembly', 'thm:result: c_S = 27/2 + 4 log(3/2) - 8 beta_2, c_T = 8(gamma_E - 3log2), d_T = 8; per class values',
    'assembled from the referee-computed J, R constants, c_loc, Sigma_ren constant (numeric) and the Fourier entries, with the class formulas checked in R13/R6/R12',
    assembled, dict(c_S='27/2+4log(3/2) = ' + mp.nstr(mp.mpf(27) / 2 + 4 * mp.log(1.5), 12), c_T=mp.nstr(8 * (g - 3 * L2), 12), d_T=8),
    agree(cS, mp.mpf(27) / 2 + 4 * mp.log(mp.mpf(3) / 2), 1e-6) and agree(cT, 8 * (g - 3 * L2), 1e-6), 'none',
    'per-class SE, EX, QU constants and log eps cancellation reproduced; beta_2 enters as -8 beta_2 (two lines x -4 beta_2 unit, plane-wave value of int(dA,dA) = 2|p|^2 (A(p),A(-p)) per eq:c4-pw); d_T = 8 from QU log x^2 (8 tau) consistent with the Konishi value 3g^2N/(4pi^2).')
gt = AL.get('gram_tree', {})
add('R15-gram-eigs', 'Prop prop:ip0: Gram_tree has eigenvalues (4hbar dim g/16pi^4) x {1 (twenty times), 2 (once, on K/sqrt6)} in the basis e_I e_J',
    'explicit 21x21 matrix S12 in the basis O_{e_I e_J}, I<=J',
    gt.get('eigenvalues_in_units_of_4hbar_dimg_over_16pi4'), '1 x20, 2 x1', False, 'minor',
    'In the basis e_I e_J the matrix is diagonal with 2 on the six O_{e_I e_I} and 1 on the fifteen O_{e_I e_J}, I<J (the proof itself says so); the eigenvalues are {1 x15, 2 x6}. K/sqrt6 has norm 2 but so do the five orthonormalised traceless diagonal combinations. Positivity, Hermiticity and the orthogonality of K to the 20-prime (also checked) are unaffected.')
add('R16-contact', '2q^2 x (five-line graph) is a contact term at x != 0',
    'analytic review: no 2-edge cycles, triangle subgraphs and the whole graph are UV convergent (degree -2), so F5_eps(q) = q^-2 f(eps q^2), f(y) = 6 zeta(3)/(16pi^2)^2 + O(y log^k y); lemma fourier (last statement) then gives O(eps log^2 eps) at x != 0',
    'justified', 'contact term', True, 'none', 'Not computed numerically; the argument is sound and only needs the remainder bound of [tp, lem:t-ft], which I could not read.')
add('R17-herm', 'Prop prop:herm (Hermiticity of the reflection form)',
    'line-by-line check against ch08 def:c8-Theta, thm:c8-fields (ii),(iii), def:c8-2pt',
    'valid', 'valid', True, 'minor',
    'Each step follows from thm:c8-fields (ii) (W^Theta = W(J^Theta) for even J-prime = J + t1 a^Theta + t2 b), Theta^2 = 1, F^Theta(0) = F(0)^c and the order reversal of c. Gap: the two-point function uses the renormalised O(eps) = O - d_delta I^CT(O); reflection compatibility of the OPERATOR counterterms (d_delta I^CT(O^Theta) = d_delta I^CT(O)^Theta) is needed and is not among the stated hypotheses (it follows by the uniqueness/reality argument of thm:c8-fields (iii), but should be stated). The graded symmetry used is [tp, lem:t-elem (iii)], not readable by me.')
add('R18-hermrad', 'Thm thm:hermrad (Hermiticity of the one-loop radial form)',
    'line-by-line check against ch08 def:c8-radial, prop:c8-C, def:c8-inner',
    'conclusion follows under its hypotheses; one justification is incorrect as written', 'valid', True, 'minor',
    'The claim that the formal exponential e^{s delta^hbar} "terminates on each component of fixed weight" is not correct: the classical part of P_n - K_n acts on operators at x != 0 as a first-order differential operator mixing all derivative orders (finite conformal flow), so the series in s does not terminate. The needed identity <e^{-s delta}F, G> = <F, e^{s delta}G> follows instead from d/ds <e^{-s delta}F, e^{s delta}G> = 0 using the infinitesimal Ward identity on the flowed observables, PROVIDED the flow exists on renormalised observables and keeps the supports disjoint for s in [0, pi/4] (true classically by prop:c8-C). The Ward identity itself is [tp, thm:t-ward] (H1-H4), and the vanishing of Lambda^(0)-exact two-point functions is [tp, cor:t-exact]; neither could be read by me.')
add('R19-ip0-frame', 'Prop prop:ip0: radial Gram matrix = minimal-scheme two-point coefficient at unit separation',
    'consistency with ch08 def:c8-inner and def:c8-radial',
    'consistent only if the frame P_1 = 0 is reached through the arity-one finite counterterm theta_P = -1/3 (which does not enter the arity-zero two-point function), not through a finite part psi of the operators', 'unconditional', None, 'minor',
    'def:c8-radial gives P_1(O) = -(1/3 + theta_P) units on the scalar bilinears; if instead P_1 = 0 were imposed by the finite part psi\' = (a_P/2){S(Delta),-}, the two-point coefficient would change (remark: a finite part m adds 2m), and in a scheme with P_1 = a_P{S(P_b),-} def:c8-inner multiplies by 1 - 2 a_P w. The proposition defers to [ip, prop:ip-frames (a)] and [ip, thm:ip-origin], which I could not read; this step is unverified.')
add('R20-spec-pw', 'spec section 1: "int dvol (A,A) -> (A(p),A(-p)) (ch04 eq:c4-pw)"',
    'comparison with ch04 eq:c4-pw', 'eq:c4-pw has 2 (A(p),A(-p))', '(A(p),A(-p))', False, 'minor',
    'The spec misquotes eq:c4-pw by a factor 2; the theory note uses a compensating kernel convention, and its final conversions (Sigma_sing = -2e g^2N/16pi^2, Sigma_F = -4 beta_2 unit p^2) agree with the referee derivation from eq:c4-pw. No effect on results, but the convention should be stated once.')
add('R21-unverified', 'steps resting on documents outside the referee scope',
    'listing', '[tp, lem:t-ft] remainder estimate; [tp, lem:t-elem (iii)]; [tp, thm:t-ward], [tp, cor:t-exact]; [ip, lem:ip-W], [ip, thm:ip-origin], [ip, prop:ip-frames (a)]; ch06 thm:c6-dil values of D_1 (consistent with the remark c6-konishi and with d_T = 8)',
    None, None, 'none', 'Not verified (not in the allowed reading list). The numerical consequences that could be tested (Fourier remainder via B, J, R convergence) behave as claimed.')

allok = all(f['agreement'] in (True, None) for f in F if f['severity'] != 'minor') and all(f['severity'] in ('none', 'minor') for f in F)
verdict = dict(verdict='partial',
               summary='Theorem thm:result (c_S = 27/2 + 4 log(3/2) - 8 beta_2, c_T = 8(gamma_E - 3 log 2), d_T = 8, per-class values, no log eps) is CONFIRMED by an independent third method on every recomputed ingredient: J, R, c_loc, Lemma local (i) and (ii), B master, Sigma_ren constant -43/8, G1, G2, Fourier entries, EX/QU combinatorics, fermion-loop sign. No fatal or major error found. Minor errors: (1) the per-graph finite constants of the bubble (167/72) and gaugino (-13/9) in prop:selfenergy are wrong (correct: 13/8 and -7; the printed ones do not sum to the stated total); (2) the eigenvalue statement for Gram_tree in prop:ip0 is wrong ({1 x15, 2 x6} in the basis e_I e_J); (3) thm:hermrad justifies the finite flow by a false termination claim (repairable by an ODE argument); (4) prop:herm should list the reflection compatibility of the operator counterterms; (5) steps resting on [tp], [ip] and the P_1 = 0 frame were not verifiable.',
               all_non_minor_agree=allok)
json.dump(dict(findings=F, verdict=verdict, assembled=assembled), open(os.path.join(HERE, 'findings.json'), 'w'), indent=1, default=str)
print(json.dumps(verdict, indent=1))
