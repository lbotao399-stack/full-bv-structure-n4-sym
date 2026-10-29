# Builds step2/literature/sweep.json (literature sweep on the one-loop correction Q_1 of the
# holomorphically twisted differential and its untwisted N=4 counterparts).
# Every statement below was read from the arXiv text (pdftotext of the arXiv PDF of the stated
# version, or the arXiv abs page where marked "abstract only"). Equation numbers are those of the
# stated version. Formulas are transcriptions (plain text / LaTeX-like); they are data, not inputs.
import json, os

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sweep.json")

meta = {
    "task": "Step 2 literature sweep: primary papers defining/computing loop corrections to the differential "
            "of the holomorphically twisted 4d gauge theory on local operators (Q = Q_0 + Q_1 + ...), and the "
            "untwisted N=4 papers they are compared with (dynamic supercharges, 1/16-BPS Q, Konishi anomaly).",
    "date": "2026-10-28",
    "method": "WebFetch of arXiv abs/HTML/PDF; the PDF text was read with pdftotext (no further files written). "
              "Equation numbers refer to the version listed in 'version_read'.",
    "name_clash_warning": [
        "Manuscript (this project): Q_1 in Sigma, Q_2 in Sigma^* are spinor parameters (Prop prop:c7-susy). Not used here.",
        "Holomorphic-twist literature (2306.01039, 2512.07771, 2606.27955, 2610.22332, 2610.20340): Q_n = n-loop part "
        "of the twisted differential (Q_0 = tree-level interacting differential, Q_1 = one-loop correction).",
        "Beisert hep-th/0310252 and Zwiebel hep-th/0511109: J_k = coefficient of g^k (g^2 = g_YM^2 N/(8 pi^2)); "
        "Beisert's (Q_1)_a^alpha is the O(g) classical length-changing supercharge, i.e. the analogue of the twisted Q_0, "
        "NOT of the twisted Q_1.",
        "2606.27955 (2.13) and 2610.20340 footnote 2 print Q = Q_free + g_YM Q_0 + g_YM^2 Q_1 + O(g_YM^3): the twisted "
        "one-loop Q_1 sits at order g_YM^2 in the unrescaled normalisation."
    ],
    "equation_numbering_note": "2306.01039: numbers are those of the v2 PDF (JHEP version), e.g. (3.27)-(3.45), (C.1)-(C.9). "
                               "The arXiv HTML rendering numbers the same equations consecutively ((58)-(76) for (3.27)-(3.45)).",
}

leads = [
    {"lead": "Budzik-Gaiotto-Kulp-Williams-Wu-Yu, Semi-chiral operators in 4d N=1 gauge theories, 2306.01039",
     "verdict": "confirmed (v1 2023-06-01, v2 2024-09-03; JHEP05(2024)245)."},
    {"lead": "Gaiotto-Kulp-Wu, Higher operations in perturbation theory, 2403.13049",
     "verdict": "confirmed (v1 2024-03-19; cited as JHEP05(2025)230 in 2512.07771 ref [17])."},
    {"lead": "2510.26608 (Gaiotto/Williams/Wu?)",
     "verdict": "authorship in the lead is wrong: 2510.26608 is Gui, Wang, Williams, 'Higher-dimensional Chiral Algebras in "
                "the Jouanolou Model and free-field realization' (v1 2025-10-30, v2 2026-06-02). It contains a recursion for "
                "Laman type I' graph weights in complex dimension 2 (Theorem 4.13) and explicit one-, two-, three-loop weights; "
                "it does not compute Q_1 of any gauge theory."},
    {"lead": "Budzik-Kulp, 2512.07771",
     "verdict": "confirmed: 'Loop Corrected Supercharges from Holomorphic Anomalies' (v1 2025-12-08, v2 2026-03-19, "
                "v3 2026-07-20; JHEP08(2026)089 per search result). This is the paper with the complete one-loop Q_1 for N=4."},
    {"lead": "Beisert hep-th/0310252",
     "verdict": "confirmed: 'The su(2|3) Dynamic Spin Chain', Nucl. Phys. B682 (2004) 487."},
    {"lead": "Chang-Yin arXiv 2209.06728",
     "verdict": "wrong attribution: 2209.06728 is Chang, Lin, 'Words to describe a black hole' (JHEP02(2023)109). "
                "Chang-Yin '1/16 BPS States in N=4 SYM' is arXiv 1305.6314 (PRD 88 (2013) 106005). Both are in the sweep."},
]

papers = []

# ---------------------------------------------------------------- twisted side: framework + integrals
papers.append({
    "id": "2207.14321",
    "title": "Feynman Diagrams in Four-Dimensional Holomorphic Theories and the Operatope",
    "authors": ["Kasia Budzik", "Davide Gaiotto", "Justin Kulp", "Jingxiang Wu", "Matthew Yu"],
    "version_read": "v1 (2022-07-28); published JHEP 07 (2023) 127 (per 2306.01039 ref [32])",
    "category": "twisted: universal one-loop (triangle) integral",
    "read_depth": "full text: Sec. 1, 2.2-2.7, 4.3, 5, App. B",
    "conventions": {
        "setting": "4d holomorphic theories (holomorphic twist of N=1); superspace propagators = (0,1)-forms (Bochner-Martinelli); "
                   "external holomorphic momenta lambda_v; shifts z_e on edges.",
        "coupling_hbar": "no coupling or hbar in the integrals; theory data enter only through the pairing eta_ab (Sec. 5, (5.4),(5.7)).",
        "twisting_supercharge": "not fixed in this paper (companion papers).",
    },
    "defined": [
        "Generating-function Feynman integrals I_Gamma[lambda_v; z_e] for graphs Gamma (Sec. 2); operations "
        "{O_n,...,O_0}_{lambda} = Pi[Q_free int prod O_k e^{lambda_k x_k}/(2 pi i)^2 O_0(0)] (5.3), contribution of a graph (5.7).",
        "Regularisation: Schwinger time with UV cutoff eps on each propagator and IR cutoff L (Sec. 2.2 around (2.13)); "
        "integration region identified with the positive part of RP^{|Gamma_1|-1} as L/eps -> infinity (Fig. 3, Sec. 2.6).",
    ],
    "computed": [
        "Triangle (one-loop) master integral by direct integration (Sec. 4.3, App. A): "
        "I[lambda;z] = -e^{lambda_2.z_12+lambda_3.z_13} (lambda_2 ^ lambda_3) [ e^{-lambda_2.Z}/((lambda_1.Z)(lambda_2.Z)) "
        "+ e^{lambda_3.Z}/((lambda_1.Z)(lambda_3.Z)) + 1/((lambda_2.Z)(lambda_3.Z)) ]  (4.19), Z = z_12 + z_23 - z_13; "
        "power series (4.20): (lambda_2^lambda_3) e^{...} sum_{n,m>=0} (lambda_3.Z)^n (-lambda_2.Z)^m/(n+m+2)!.",
        "Quadratic relations (1.1)/(Sec. 3) used to bootstrap two- and three-loop integrals (Sec. 4.4-4.8).",
    ],
    "scheme": {
        "claimed_scheme_independent": [
            "Final integrals finite, no counterterms needed although UV/IR regulators are used at intermediate steps (Sec. 1, bullet list p.3).",
        ],
        "scheme_dependent": [
            "App. B: replacing the combined propagator by a more general F(y,x) d^2y (B.1) with constraint (B.5) leaves the tree "
            "diagram unchanged; shifting P -> P + d_y H changes I_Gamma by boundary terms (B.6) that are 'essentially an "
            "infinitesimal operator redefinition' relating the operations before and after the deformation.",
        ],
    },
    "explicit_Q1": ["No Q_1 of a specific theory; the triangle (4.19) is the universal coefficient used by 2306.01039 App. C and 2512.07771 (2.45)."],
    "other_generators": [],
    "account": "Defines the universal Feynman integrals of 4d holomorphic (twisted N=1) theories with Schwinger-time UV/IR "
               "regularisation, proves quadratic (associativity) identities among them and evaluates the one-loop triangle in "
               "closed form (4.19)-(4.20). This triangle is the only integral entering Q_1. Appendix B is the only place in the "
               "swept twisted literature where a change of regularisation (propagator gauge) is shown to act on the operations "
               "by an operator redefinition (B.6).",
    "caveats": ["Gauge here means the choice of propagator form/Schwinger projective gauge (Fig. 3 'Feynman gauge t1+t2+t3=1'), not a field gauge."],
})

papers.append({
    "id": "2306.01039",
    "title": "Semi-Chiral Operators in 4d N=1 Gauge Theories",
    "authors": ["Kasia Budzik", "Davide Gaiotto", "Justin Kulp", "Brian R. Williams", "Jingxiang Wu", "Matthew Yu"],
    "version_read": "v2 (2024-09-03), JHEP 05 (2024) 245",
    "category": "twisted: definition of Q = Q_0 + Q_1 + ... and first one-loop computation (pure N=1 SYM)",
    "read_depth": "full text: Sec. 2 (start), 3 (3.1)-(3.45), 4.2.1, 5.2, App. C, App. E",
    "conventions": {
        "spinors": "Euclidean; supercharges Q_{alphadot}, Qbar_alpha with {Q_{alphadot}, Qbar_beta} = P_{alphadot beta} (2.1); "
                   "dotted/undotted roles switched relative to standard notation (footnote 7).",
        "twisting_supercharge": "Q = Q_{-dot} (Sec. 2); z^alpha = x^{+dot alpha} holomorphic, {Q, Qbar_alpha} = d_{zbar^alpha} (2.3).",
        "fields": "b in Omega^{0,*}(C^2, g^vee), c in Omega^{0,*}(C^2, g)[1] (3.1); action int d^2z Tr b(dbar c - 1/2[c,c]) (3.2); "
                  "dictionary: F_{+dot+dot} <-> d^2z b^(0); gaugino lambda-bar_alpha <-> d_{z^alpha} c^(0) (Sec. 3, bullets after (3.3)).",
        "coupling_hbar": "BF action has no coupling; holomorphic coupling tau only in int d^2z tau Tr d_alpha c d^alpha c (3.4), a total "
                         "derivative (affects only instantons). Loop counting by hbar in (3.27).",
        "su2_normalisation": "Q b^A = eps_ABC c^B b^C, Q c^A = 1/2 eps_ABC c^B c^C (3.22); b^2 = b^A b^A = -2 Tr b^2 (3.23); "
                             "psi^alpha chi_alpha = psi_2 chi_1 - psi_1 chi_2 = -psi_alpha chi^alpha (footnote 54).",
    },
    "defined": [
        "Full differential dbar + hbar Delta_BV + {I, .}_BV (3.27); expansion Q = Q_0 + Q_1 + ... with Q_n from n-loop diagrams (3.28).",
        "Q_0 O = Q_free[int_{|z|>eps} I d^2z O(0)] = int_{|z|=eps} I d^2z O(0) (3.30) = {I, O} (3.31): tree level, excision regulator.",
        "Q_n O = 1/(n+1)! {I, I, ..., I, O} (3.33) (higher bracket with n+2 slots); Maurer-Cartan condition 1/2{I,I} + 1/6{I,I,I} + ... = 0 (3.34).",
        "Five-step evaluation procedure (after (3.34)): regularised propagators, Q_free integrated by parts onto propagators, project to "
        "semi-chiral fields at the origin, integrals 'finite and independent of the regulator' (step 5). Propagator (3.36), dbar P = d^2xbar delta^4 (3.37).",
        "Nilpotency order by order: Q^2 = 0 (3.41); sum_{n1+n2=n} Q_{n1} Q_{n2} = 0 (3.42); first condition {Q_0, Q_1} = 0.",
        "Selection rules: Q_n V(N_b,N_c,J_1,J_2) = 0 if N_b + N_c < n+1 (3.43); Q_n: V(n+1,0,J1,J2) -> V(1,1,J1+n,J2+n) (3.44); "
        "Q_n: V(n,1,J1,J2) -> V(0,2,J1+n,J2+n) (3.45); Q_n acts as a derivation on longer monomials; U(1)_b filtration.",
    ],
    "computed": [
        "One loop: only triangle diagrams (Sec. 3.2.3); Q_1 O = 1/2 {I, I, O} (C.1); action on longer words = sum over pairs (C.2); "
        "triangle taken from 2207.14321 ([32]).",
        "Bootstrap of Q_1 from {Q_0,Q_1}=0 up to two constants C_1, C_2 (C.5)-(C.7); direct evaluation gives C_2 = 0, C_1 = -2 (text after (C.5)).",
        "Q_1 is not nilpotent; two-loop Q_2 bootstrapped: Q_2 b^A b^A b^B b^B = (20 C_1^2/3) eps_ABC b^A d_alpha d_beta b^B d^alpha d^beta c^C (C.8), (C.9).",
        "Planar pure SYM: Q_1 A_n ~ n Tr b^{n-2} d_alpha b d^alpha c (5.20), Q_1 A_n = (2n/(n-1)) d_alpha B^alpha_{n-1} in tree-level cohomology (5.21); "
        "only the C_n tower survives (Sec. 5.2).",
        "Holomorphic confinement: Q_1 Tr b^2 proportional to d_alpha(Tr b d^alpha c) (3.40), so the derivative of the stress tensor becomes Q-exact; "
        "no higher corrections act on Tr b^2 by (3.43).",
    ],
    "scheme": {
        "claimed_scheme_independent": [
            "Cohomology expected scheme independent: Q obtained in different schemes should be related by quasi-isomorphisms (Sec. 3.2.1, text before (3.28)).",
            "Integrals in step 5 'finite and independent of the regulator' (Sec. 3.2.1).",
        ],
        "scheme_dependent": [
            "'The resulting Q operator will be scheme-dependent' (Sec. 3.2.1, paraphrase).",
            "App. E: parameterisation of deformation space scheme dependent beyond first order; BRST anomaly as odd vector field Q = eta_{1;i}^j g_j d_{g_i} + ... (E.3).",
            "Footnote 28: Leibniz rule holds only for the loop-corrected (regularised) product, which is not associative.",
        ],
    },
    "explicit_Q1": [
        "SU(2), eps_ABC normalisation (3.22): Q_1(b^A b^B) = d_alpha b^A d^alpha c^B + d_alpha b^B d^alpha c^A - 2 delta^AB d_alpha b^C d^alpha c^C (C.3).",
        "Q_1(b^A c^B) = (-1)[ d_alpha c^A d^alpha c^B - delta^AB d_alpha c^D d^alpha c^D ] (C.4).",
        "General derivative insertions: (C.6)-(C.7) with functions F, H, G, K (App. C).",
        "Gauge-invariant: Q_1 Tr b^2 proportional to d_alpha(Tr b d^alpha c) (3.40); planar towers (5.20)-(5.21).",
        "N=4 (Sec. 4.2.1): only the tree-level Q_0 is given, (4.14); no one-loop formula for N=4 in this paper.",
    ],
    "other_generators": [
        "Holomorphic stress tensor S_alpha not closed at one loop: (Q + dbar) S_alpha ~ d_alpha Tr(d_beta c d^beta c) (3.24), "
        "described as the direct translation of the Konishi anomaly; a holomorphic vector field shifts the action by "
        "int d^2z v^alpha d_alpha Tr(d_beta c d^beta c) (3.25) and tau by a multiple of d_alpha v^alpha (holomorphic scale anomaly); "
        "divergence-free (symplectic) vector fields stay symmetries.",
        "Footnote 20: on a curved complex surface Y the one-loop quantisation has an anomaly proportional to c_1(TY) times the adjoint quadratic Casimir.",
        "Footnote 27: even translations can receive quantum corrections (P_mu J^mu for an anomalous current).",
    ],
    "account": "Introduces the loop-corrected twisted differential Q = Q_0 + Q_1 + ... acting on the free semi-chiral cohomology, "
               "with Q_n = (n+1)!^{-1}{I,...,I,O} computed by n-loop Laman diagrams. For pure N=1 SYM it computes Q_1 on pairs of "
               "letters (C.3)-(C.4) from the triangle of 2207.14321, finds Q_1 Tr b^2 proportional to the divergence of the stress "
               "tensor (holomorphic confinement), shows Q_1^2 != 0 so that Q_2 is required (C.8)-(C.9), and computes the planar "
               "one-loop cohomology. It states that Q is scheme dependent while its cohomology should be scheme independent up to "
               "quasi-isomorphism. For N=4 it gives only the tree-level differential (4.14).",
    "caveats": [
        "Relative normalisation of (C.4) versus 2512.07771 (3.35) [Q_1(b^A c^B) = 1/2 f_ACD f_BCE d c^D d c^E with Killing-normalised generators (3.12)] "
        "not reconciled here; 2512.07771 footnote 23 records a corrected factor 1/2 in its z=w=0 formulas.",
        "Twisting supercharge is antichiral-labelled Q_{-dot} with switched dotted/undotted conventions; 2512.07771 uses Q = Q^1_- with dotted holomorphic indices.",
    ],
})

papers.append({
    "id": "2403.13049",
    "title": "Higher Operations in Perturbation Theory",
    "authors": ["Davide Gaiotto", "Justin Kulp", "Jingxiang Wu"],
    "version_read": "v1 (2024-03-19); cited as JHEP 05 (2025) 230 in 2512.07771",
    "category": "twisted: general framework (BRST anomaly brackets, schemes)",
    "read_depth": "full text: Sec. 2.2, 4.1-4.5 (incl. 4.4.3), 6.1-6.2 (start)",
    "conventions": {
        "setting": "holomorphic-topological theories on C^H x R^T, action int[(Phi, d Phi) + I(Phi)] d^H x (1.1), d = dx_R d/dx_R + dxbar d/dxbar (1.2).",
        "coupling_hbar": "perturbative couplings g^i of interactions; no specific normalisation.",
    },
    "defined": [
        "Action of BRST on local operators in the deformed theory: odd map Q with L_eta^2 = {eta, Q} + Q^2 = 0 (2.26); "
        "Q O = {O} + {I, O} + 1/2 {I, I, O} + ... (2.27) (right L_infinity-module); {I,O}(x_2) = int_{|x_12|=eps} J(x_1) O(x_2) (2.28).",
        "BRST anomaly brackets from Wick contractions with regularised propagators P_eps, d P_eps = K_eps (Sec. 4.1, (4.9)-(4.12)).",
        "4d holomorphic lambda-bracket {f_lambda g} = oint_{S^3} omega_BM(w) e^{lambda.w} e^eta(f(w) (x) g(0)) (4.40), omega_BM (4.41); "
        "associativity relations (4.42), (4.43).",
        "Interacting brackets {O_1,...,O_n}_I = {O_1..O_n} + {I,O_1..O_n} + 1/2{I,I,O_1..O_n} + ... (4.44).",
        "Point-splitting schemes labelled by configurations of little disks p in D_d(n) (Sec. 6.1, (6.2)-(6.5)); integrated composites over chains (6.18)-(6.21).",
    ],
    "computed": [
        "Quadratic axioms for the brackets proved geometrically (Sec. 3.5, 4.2); non-renormalisation for >= 2 topological directions (Sec. 3.8).",
    ],
    "scheme": {
        "claimed_scheme_independent": ["Operations from closed chains in H_*(Conf_n) are scheme-independent in cohomology (end of Sec. 6.1, paraphrase)."],
        "scheme_dependent": ["An actual choice of contour gamma_n in D_d(n) is needed to prescribe a renormalisation scheme (Sec. 6.1)."],
    },
    "explicit_Q1": ["None for a 4d gauge theory; supplies the definition used by 2306.01039 (3.33) and 2512.07771 (2.22)-(2.27)."],
    "other_generators": ["Section 2.3: translations P_mu expected to receive no quantum corrections in the space of position-dependent couplings ((2.31) and footnote 13)."],
    "account": "General formalism: the loop-corrected BRST action on local operators is the L_infinity-module structure (2.27) whose "
               "coefficients are BRST-anomaly brackets computed by regularised Feynman diagrams; gives the 4d holomorphic lambda-brackets "
               "and their quadratic axioms (Wess-Zumino consistency) and a point-splitting description of schemes. It is the framework "
               "cited for the scheme statements of 2306.01039 and 2512.07771.",
    "caveats": [],
})

papers.append({
    "id": "2512.07771",
    "title": "Loop Corrected Supercharges from Holomorphic Anomalies",
    "authors": ["Kasia Budzik", "Justin Kulp"],
    "version_read": "v3 (2026-07-20); v1 2025-12-08, v2 2026-03-19; JHEP 08 (2026) 089 (search result)",
    "category": "twisted: complete one-loop Q_1 for Lagrangian 4d gauge theories incl. N=4",
    "read_depth": "full text: Sec. 1-3, App. A, App. B (B.1)-(B.6), bibliography",
    "conventions": {
        "spinors": "Euclidean; {Q^i_alpha, Qbar_{j betadot}} = delta^i_j P_{alpha betadot} (2.1), {Q,Q} = {Qbar,Qbar} = 0 (2.2); "
                   "z^{alphadot} := x^{+ alphadot}, zbar^{alphadot} := x^{- alphadot} (2.4); '-dot, +dot = 1, 2' (footnote 20).",
        "twisting_supercharge": "Q := Q^1_- for N=1 (1.2)/(2.3); Q := Q^4_- for N=4 (Sec. 3.4, before (3.51)). Anti-holomorphic translations Q-exact: {Q, Qbar^1_{alphadot}} = P_{- alphadot} (2.5).",
        "superfields": "O = e^{i dzbar^{alphadot} Qbar^1_{alphadot}} O^(0) (2.7); D_- = Q + dbar (2.8); descent (2.9)-(2.10).",
        "fields_and_dictionary": "BV action int_{C^2} d^2z [Tr b(dbar c + 1/2[c,c]) + beta_I(dbar gamma^I + c gamma^I) + W[gamma]] (3.1) (bracket placement reconstructed from the extracted text; note the sign +1/2[c,c] versus -1/2[c,c] in 2306.01039 (3.2)); "
                                 "N=1 dictionary gamma ~ phibar, beta ~ psi_+, d_{alphadot} c ~ lambdabar_{alphadot}, b ~ F_{++}, d_{alphadot} ~ D_{+ alphadot} (3.6); "
                                 "N=4 dictionary gamma^I ~ Phi^{4I}, beta_I ~ Psi^I_+, d_{alphadot} c ~ Psibar_{4 alphadot}, b ~ F_{++}, d_{alphadot} ~ D_{+ alphadot} (3.51).",
        "N4_superpotential": "W[gamma] = Tr gamma^1[gamma^2, gamma^3] = 1/3! eps_IJK f_ABC (gamma^I)^A (gamma^J)^B (gamma^K)^C (3.48).",
        "lie_algebra": "[t_A, t_B] = f_AB^C t_C, indices raised/lowered with delta^AB, generators normalised with respect to the Killing form (3.12).",
        "coupling_hbar": "Q = Q_free + Q_tree + hbar Q_{1-loop} + hbar^2 Q_{2-loop} + ... (1.3); Q O = {I,O}_0 + hbar/2! {I,I,O}_0 + ... (2.22); "
                         "no relation of hbar to g_YM is given; propagator P = (zbar^2 dzbar^1 - zbar^1 dzbar^2)/|z|^4 (2.31) with dbar P = pi^2 delta^4 d^2zbar (2.32); "
                         "bracket measure d^2z/(2 pi i)^2 (2.27).",
    },
    "defined": [
        "n-loop correction Q_n O = 1/(n+1)! {I,...,I,O} (1.6); L_infinity Maurer-Cartan: sum_{n1+n2=n} Q_{n1}Q_{n2} = 0 (2.25), "
        "{Q_0,Q_1} = 0, Q_1^2 + {Q_0,Q_2} = 0 (2.26)/(1.12).",
        "Higher lambda-brackets (2.27); graph contribution as differential operator I_Gamma[lambda + d_w; w] (2.41); only Laman graphs contribute; "
        "at one loop only the triangle (Sec. 2.3, Fig. 1).",
        "Regulator: Schwinger parametrisation with UV cutoff t in [eps, infinity) on every propagator, P -> P_eps(z+w, zbar, t) (2.42); combined "
        "propagator on spacetime x Schwinger space (A.2) with heat kernel K_t(x) = t^{-2} e^{-|x|^2/t}; integral over RP^2_> (A.3)-(A.9).",
        "Generating operator D_{w,z}(f,g) (3.9), D_{0,0}(f,g) = 1/2 d_{alphadot} f d^{alphadot} g (3.10).",
    ],
    "computed": [
        "One-loop master (triangle) integral (2.45), series (2.46)-(2.47), I[lambda;0] = 1/2 lambda_1 ^ lambda_2 (2.48); App. A (A.1)-(A.12).",
        "Q_1 acts on pairs of letters: Q_1(f_1...f_n) = sum_{i<j} (-1)^{eps_ij} Q_1(f_i f_j) f_1 ... (3.7)-(3.8); vanishing pairs "
        "Q_1(cc) = Q_1(gamma gamma) = Q_1(gamma c) = Q_1(beta c) = 0 (3.11).",
        "General gauge group + chirals in R + W: (3.13), (3.14), (3.16)-(3.20); adjoint matter (3.21)-(3.26); N=4 with shifts (3.27)-(3.31).",
        "N=4 at zero shifts (3.64)-(3.69) (also obtained by acting with tree-level psl(3|3) generators on Q_1(bc), (3.55)-(3.63)); "
        "superfield form (3.70); derivative insertions App. B (B.1)-(B.6), coefficients C_mn (B.3)-(B.4).",
        "N=1 SYM (3.35)-(3.37), N=2 SYM (3.41)-(3.47).",
        "Planar N=4 (Sec. 3.4.1): Q_1 splits/joins traces (3.71)-(3.76); Q_1(gamma gamma) = 0 implies Q_1 annihilates Tr gamma^{(I_1}...gamma^{I_n)} (3.77) "
        "and, through footnote 24, their psl(3|3) descendants; the single-trace-to-single-trace part of Q_1 does not lift the infinite-N "
        "single-trace tree-level cohomology; the trace-splitting/joining part is left open.",
    ],
    "scheme": {
        "claimed_scheme_independent": [
            "Lifting of cohomology classes from the final perturbative cohomology is scheme independent (Sec. 1, paragraph after (1.12); "
            "footnote 9: under the assumption that non-perturbative corrections can only further reduce the cohomology).",
            "Integrands manifestly UV and IR finite; holomorphic theories finite, no counterterms (footnote 22, Sec. 2.2.2 end).",
        ],
        "scheme_dependent": [
            "Higher lambda-brackets, hence Q_n, depend on the scheme (Sec. 1, after (1.12)); 'generically scheme dependent beyond first-order "
            "in perturbation theory'; a scheme in which Q_n vanishes at some order may exist, depending on the theory; equivalence between schemes = "
            "quasi-isomorphism of the L_infinity conformal algebras (end of Sec. 2.2.2).",
            "Footnote 21: generically no composite-operator scheme exists that is associative and in which the interacting Q obeys Leibniz; "
            "its existence would be a one-loop non-renormalisation theorem.",
            "Footnote 8: corrections on an operator of length n truncate at most at Q_{n+1} (graph combinatorics), expected Q_{n-1} from explicit integrals.",
        ],
    },
    "explicit_Q1": [
        "Pure gauge sector (any g): Q_1(b^A(z) c^B(w)) = f_ACD f_BCE D_{w,z}(c^D, c^E) (3.13); Q_1(b^A(z) b^B(w)) = f_ACD f_BCE [D_{w,z}(c^D,b^E) - D_{w,z}(b^D,c^E)] (3.14).",
        "N=4, zero shifts: Q_1(b^A c^B) = 1/2 f_ACD f_BCE d_{alphadot} c^D d^{alphadot} c^E (3.64).",
        "Q_1(b^A b^B) = 1/2 f_ACD f_BCE [ d c^D d b^E - d b^D d c^E + d (beta_I)^D d (gamma^I)^E - d (gamma^I)^D d (beta_I)^E ] (3.65) "
        "(each 'd X d Y' = d_{alphadot} X d^{alphadot} Y).",
        "Q_1((beta_I)^A (gamma^J)^B) = 1/2 delta_I^J f_ACD f_BCE d c^D d c^E (3.66) (also computed in 2511.09519).",
        "Q_1((beta_I)^A (beta_J)^B) = 1/2 eps_IJK f_ACD f_BCE [ d c^D d (gamma^K)^E - d (gamma^K)^D d c^E ] (3.67).",
        "Q_1(b^A (gamma^I)^B) = 1/2 f_ACD f_BCE [ d c^D d (gamma^I)^E - d (gamma^I)^D d c^E ] (3.68).",
        "Q_1((beta_I)^A b^B) = 1/2 f_ACD f_BCE [ d (beta_I)^D d c^E + d c^D d (beta_I)^E + eps_IJK d (gamma^J)^D d (gamma^K)^E ] (3.69).",
        "Superfield C = c + theta_I gamma^I + 1/2 eps_IJK theta^I theta^J beta_K + theta^1 theta^2 theta^3 b (3.52); Q_0 C = 1/2[C,C] (3.54); "
        "Q_1(C^A(theta) C^B(theta')) = -1/2 f_ACD f_BCE (theta_1-theta'_1)(theta_2-theta'_2)(theta_3-theta'_3) d_{alphadot} C^D(theta) d^{alphadot} C^E(theta') (3.70).",
        "Derivatives: Q_1(d_{alphadot} f g) = d_{alphadot} Q_1(fg) - Q_1(f d_{alphadot} g) (B.1); "
        "Q_1(b^A d_1^m d_2^n c^B) = f_ACD f_BCE 1/(m+n+2) sum_{k=0}^m sum_{l=0}^n 1/(k+l+1) binom(m,k) binom(n,l) d_1^k d_2^l d_{alphadot} c^D d_1^{m-k} d_2^{n-l} d^{alphadot} c^E (B.5); (B.6) and analogues.",
        "Tree level N=4: Q_0 c = 1/2[c,c], Q_0 b = [c,b] - [gamma^I, beta_I] (3.49); Q_0 gamma_I = [c,gamma_I], Q_0 beta_I = [c,beta_I] + 1/2 eps_IJK [gamma^J, gamma^K] (3.50).",
    ],
    "other_generators": [
        "Residual psl(3|3) (commutant of Q^4_- in psl(4|4)) used at tree level, super vector fields Q^I_+ ~ d_{theta_I}, S^I_+ ~ theta^I(z.d + theta.d_theta), "
        "Qtilde_{I alphadot} ~ theta^I d_{alphadot}, Stilde^{I alphadot} ~ z^{alphadot} d_{theta_I} (3.55)-(3.56); component actions (3.57)-(3.58).",
        "Footnote 24: the tree-level psl(3|3) action may itself receive loop corrections, but Q_1 (anti)commutes with the tree-level generators R_0 on "
        "Q_0-cohomology: {Q_1, R_0} + {Q_0, R_1} = 0. The corrections R_1 (one-loop S, Stilde, ... of the commutant) are not computed.",
        "Sec. 2.2.1: corrections to Q identified with a 'twice-generalized Konishi anomaly'; if classical scale/U(1)_R symmetry is broken at one loop, "
        "one-loop Q lifts the semi-chiral S_{alphadot} (2.18)-(2.19); for N=4 no such breaking is discussed.",
        "Intro (1.4)-(1.5): analogy with D acting on :phi1 phi2: and with P_mu J_5^mu = g^2 hbar 2C(R)/(16 pi^2) tr F Ftilde.",
        "Footnote 6 / (1.11): one-loop correction to the BRST operator of the chiral-algebra subsector of N=4 computed in 2310.20086: "
        "Q_1^{2d}(C^A(theta) C^B(theta')) = -2 f_ABC (theta^1-theta'^1)(theta^2-theta'^2) d C^C(theta).",
    ],
    "account": "Gives the complete one-loop correction Q_1 to the twisted supercharge for 4d N=1 Lagrangian gauge theories with chirals and "
               "superpotential, specialised to N=1, N=2 and N=4 SYM, with the Schwinger-time UV cutoff regularisation of 2207.14321. Q_1 is a "
               "second-order operator acting on pairs of letters, with coefficients from the one-loop triangle; on N=4 it is (3.64)-(3.69), fixed "
               "from Q_1(bc) by tree-level psl(3|3), and packaged in C^{2|3} superspace as (3.70). It shows Q_1 annihilates the half-BPS-tower "
               "representatives Tr gamma^{(I1}...gamma^{In)} and leaves the infinite-N single-trace cohomology unlifted by its single-trace part. "
               "Q_n are scheme dependent; liftings are claimed scheme independent. Other generators (psl(3|3)) are used only at tree level; their "
               "one-loop corrections are not computed. No comparison with Beisert/Zwiebel spin-chain supercharges appears (no such citation).",
    "caveats": [
        "Footnote 23: a missing factor 1/2 in the z=w=0 formulas of Secs. 3.2-3.4 was corrected (thanks to J. Park); earlier versions differ by this factor.",
        "Twisting supercharge differs between sections: Q^1_- (general N=1) vs Q^4_- (N=4, matching Chang-Yin 1305.6314).",
        "Absolute normalisation of Q_1 relative to g_YM is not fixed in this paper (hbar loop counting only).",
    ],
})

papers.append({
    "id": "2509.16737",
    "title": "Unravelling the Holomorphic Twist II: Anomalies and Extended Supersymmetry",
    "authors": ["Pieter Bomans", "Niklas Garner", "Brian R. Williams", "Jingxiang Wu"],
    "version_read": "v1 (2025-09-20)",
    "category": "twisted: one-loop anomalies of symmetries as higher brackets; equivariant QME",
    "read_depth": "full text: Sec. 3.2 end, 3.3.1-3.3.2, 4.3.1",
    "conventions": {
        "setting": "holomorphic twist of 4d N>=1; background fields mu (Beltrami), A (flavour); L = symmetry dg Lie algebra.",
        "N4_superfield": "A = C + theta^3 B = c + theta^I phi_I + 1/2 eps_IJK theta^I theta^J lambda_K + theta^1 theta^2 theta^3 b (4.81); "
                         "S_BV = int_{C^{2|3}} d^2z d^3theta Tr(1/2 A dbar A + 1/3 A^3) (4.80).",
    },
    "defined": [
        "L-equivariant classical master equation d_CE S^L + 1/2{S^L, S^L}_BV = 0 (3.46) for S^L[mu,A] = S_0 + int d^2z(A^a J_a - mu^i S_i) (3.45).",
        "Anomaly = failure of the L-equivariant QME (d_CE + Q) I^L + 1/2{I^L, I^L}_BV = 0 (3.47) ('inner quantum L background' if satisfied).",
        "Background-coupled quantum differential Q^L = d_CE + Q + {I^L_0 -} + hbar/2 {I^L_0 I^L_0 -} + hbar^2/6 {I^L_0 I^L_0 I^L_0 -} + ... (3.49).",
    ],
    "computed": [
        "Flavour triangle: (d_CE + Q) I^L + 1/2{I^L,I^L}_BV = -k_FFF/(6 pi^2) int tr(A_F dA_F dA_F) (3.48).",
        "One-loop non-semi-chirality of currents from ternary brackets: Q^L J_Fa (3.50); with ABJ anomaly Q J_f = dbar J_f - hbar k_FFf/(2 pi^2) eps^{ij} tr_F(d_i c_F d_j c_F) (3.51); "
        "Q S_i = dbar S_i - hbar k_RFF/(4 pi^2) eps^{jk} d_i tr_F(d_j c_F d_k c_F) (3.52) for the stress tensor when the twisting U(1)_R has an ABJ anomaly.",
    ],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": ["No Q_1 on general operators; (3.51)-(3.52) are Q_1 on the currents."],
    "other_generators": [
        "N=4 (Sec. 4.3.1): classical currents S_i, R^I_J (4.82), L^I_i (4.83), G^I_i, Gtilde^I (4.84), M_IJ, N (4.85) of the twisted superconformal symmetry; "
        "no loop corrections computed for them.",
    ],
    "account": "Formulates anomalies of (infinite-dimensional) symmetries of holomorphically twisted theories as failures of the "
               "L-equivariant quantum master equation (3.47) and as ternary-bracket (one-loop) terms in the action of Q on the symmetry "
               "currents (3.49)-(3.52). This is the twisted-side analogue of an equivariant one-loop QME along symmetry ghosts; for N=4 "
               "it lists only the classical currents.",
    "caveats": [],
})

papers.append({
    "id": "2510.26608",
    "title": "Higher-dimensional Chiral Algebras in the Jouanolou Model and free-field realization",
    "authors": ["Zhengping Gui", "Minghao Wang", "Brian R. Williams"],
    "version_read": "v2 (2026-06-02)",
    "category": "twisted: Feynman graph weights (recursion), no gauge-theory Q_1",
    "read_depth": "intro (Theorem statement before (0.5)), Sec. 4.4 ((4.1)-(4.5)), App. C start",
    "conventions": {"setting": "complex dimension d (mostly 2), Jouanolou model of configuration-space cohomology; propagators P, residue."},
    "defined": ["Higher chiral operations; Laman graphs of type I' and Henneberg I' moves (Sec. 4); recursion Theorem 4.13."],
    "computed": [
        "One-loop triangle weight (4.1), agreeing with 2207.14321 Sec. 4.3 (4.22)-(4.23) (stated after (4.1)); leading derivative terms (4.2)-(4.5), e.g. "
        "mu(P_1o P_12 P_2o dz) = 1/2 Lambda dz, Lambda = lambda_1 ^ lambda_2 (4.2).",
        "Two-loop Theta-diagram weight to lowest order (0.5): 1/24 (lambda_1 ^ (lambda_3 + 2 lambda_2)).(lambda_3 ^ (lambda_1 + 2 lambda_2)) dz.",
        "Higher Kac-Moody/Virasoro central extensions via one-loop coefficients (Prop. 3.13, Thm 3.15, App. C).",
    ],
    "scheme": {"claimed_scheme_independent": ["Weights of type I' Laman graphs completely fixed by the L_infinity relations (Thm 4.13, paraphrase of intro)."], "scheme_dependent": []},
    "explicit_Q1": ["None (no gauge theory)."],
    "other_generators": [],
    "account": "Mathematical construction of higher chiral operations; its relevance to Q_1 is that it rederives the one-loop triangle "
               "weight used by 2512.07771 and gives two- and three-loop weights by an L_infinity recursion (relevant for Q_2).",
    "caveats": ["Not by Gaiotto/Williams/Wu as the lead suggested."],
})

papers.append({
    "id": "1809.02661",
    "title": "Renormalization for holomorphic field theories",
    "authors": ["Brian R. Williams"],
    "version_read": "v2 (2020-02-24)",
    "category": "Costello-formalism foundation for holomorphic theories (one-loop finiteness, one-loop QME obstruction)",
    "read_depth": "Sec. 3 (Thm 3.1, Thm 3.4), Sec. 4.2 (Lemma 4.3)",
    "conventions": {"formalism": "Costello's homotopy RG flow W(P_{eps<L}, I) with heat-kernel propagators; holomorphic theories on C^d."},
    "defined": ["One-loop obstruction hbar Theta[L] = Q^hol I[L] + 1/2 lim_{eps->0} e^{-I/hbar} e^{hbar d_{P_{eps<L}}} {I,I}_eps e^{I/hbar} mod hbar^2 (23) (Lemma 4.3)."],
    "computed": ["Theorem 3.1 / 3.4: a holomorphic theory on C^d admits a one-loop prequantisation with no counterterms (the eps -> 0 limit of W(P_{eps<L}, I) mod hbar^2 exists)."],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": ["None; supplies the counterterm-free one-loop setting in which twisted Q_1 is finite."],
    "other_generators": ["Classification of holomorphic one-loop anomalies (Sec. 4)."],
    "account": "In Costello's BV renormalisation framework, holomorphic theories on C^d are one-loop finite (no counterterms, Thm 3.4) and "
               "the one-loop QME obstruction reduces to (23). This is the twisted-side counterpart of the project's UV-boundary analysis: "
               "on the twisted side no singular counterterms exist at one loop, so all scheme freedom of Q_1 lies in the regularisation/"
               "propagator choice and finite redefinitions.",
    "caveats": [],
})

papers.append({
    "id": "2310.20086",
    "title": "On 1/8-BPS black holes and the chiral algebra of N=4 SYM",
    "authors": ["Chi-Ming Chang", "Ying-Hsuan Lin", "Jingxiang Wu"],
    "version_read": "v2 (2024-03-06), Adv. Theor. Math. Phys. 28 (2024) 2431 (per 2512.07771 ref [54])",
    "category": "twisted (holomorphic-topological / Omega-deformed): loop corrections to the BRST differential",
    "read_depth": "Sec. 3 ((3.1)-(3.23))",
    "conventions": {
        "twisting_supercharge": "Q_HT = Q^4_- + Q_{-dot}^3 (3.1) (index placement as extracted: 'Q4- + Q-dot3'); one-form supercharge (3.2).",
        "fields": "C in Omega^{0,*}_C (x) Omega^*_{R^2} (x) g[theta^1, theta^2] (3.4), action int dz dtheta^1 dtheta^2 Tr(1/2 C dC + 1/6 C[C,C]) (3.5).",
        "coupling_hbar": "Q_int = Q_0 + hbar Q_1 + hbar^2 Q_2 + ... (Sec. 3.1, text after (3.11)).",
    },
    "defined": ["Q_0 C = 1/2[C,C] (3.9); components (3.10)."],
    "computed": [
        "4d HT twist: all loop diagrams contributing to the BRST differential vanish (H=1, T=2: omega_Gamma vanishes on 3-Laman graphs), so Q_int = Q_0 "
        "to all orders and obeys Leibniz (Sec. 3.1, text after (3.11)).",
        "2d Omega-deformed theory: Q^{2d}_1(q_k^A q_l^B) = 2 eps_kl f^ABC dc^C, Q^{2d}_1(b^A b^B) = -2 f^ABC db^C (3.17); Q^{2d}_1(q_k^A b^B) = 2 f^ABC dq_k^C, "
        "Q^{2d}_1(b^A c^B) = -2 f^ABC dc^C (3.18); Q^{2d}_1(C^A(z,theta) C^B(z,theta')) = -2(theta^1-theta'^1)(theta^2-theta'^2) f^ABC dC^C (3.19); "
        "Q^{2d}_0 + Q^{2d}_1 = (1/hbar) oint dz/(2 pi i) Tr(b[c,c] + eps^{ij} c[q_i,q_j]) (3.20); Q^{2d}_2 trivial.",
    ],
    "scheme": {
        "claimed_scheme_independent": [],
        "scheme_dependent": ["Composite operators: normal-ordered (4d) vs nested regularised (2d) products give differentials Q^{4d}_0 vs Q^{2d}_BRST with (expected) isomorphic cohomologies, left open (end of Sec. 3.2)."],
    },
    "explicit_Q1": ["4d HT twist: Q_1 = 0. 2d chiral algebra: (3.17)-(3.19)."],
    "other_generators": [],
    "account": "Shows that for the holomorphic-topological twist of N=4 (Schur sector) the 4d BRST differential receives no loop corrections, "
               "while the 2d Omega-deformed chiral algebra has a one-loop Q_1^{2d} acting on pairs (3.17)-(3.19), reproducing the BRST current of "
               "the VOA. Footnote 5 reads Zwiebel hep-th/0511109 as showing the leading loop correction 'takes the form [Q_0, -]'.",
    "caveats": ["Footnote 5 identification with Zwiebel's O(g^3) result is interpretive (different sector and generator)."],
})

papers.append({
    "id": "2002.10517",
    "title": "A Taxonomy of Twists of Supersymmetric Yang-Mills Theory",
    "authors": ["Chris Elliott", "Pavel Safronov", "Brian R. Williams"],
    "version_read": "project source [ESW] (local copy of the arXiv PDF)",
    "category": "tree-level twist dictionary (project primary source)",
    "read_depth": "Sec. 10.3, 10.3.1 (Theorem 10.9)",
    "conventions": {"twisting_supercharge": "rank (1,0) Q in S_+ (x) W spanning C w_1 inside the scalar summand (Sec. 10.3.1)."},
    "defined": [],
    "computed": ["Theorem 10.9: the holomorphic twist of 4d N=4 SYM on R^4 is perturbatively equivalent to BF theory with fields T^*[-1] Map(C^2_Dol, Bg); "
                 "field collections listed after the proof (c + A_{0,1} + B, phi + psi + sigma^*, ...)."],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": ["None (classical)."],
    "other_generators": [],
    "account": "Classical (perturbative, tree-level) identification of the N=4 holomorphic twist used by the twisted-side literature; "
               "no loop corrections. Relevant to R2: a twisted statement may appear only as a consequence of the untwisted computation.",
    "caveats": [],
})

# ---------------------------------------------------------------- twisted Q_1 applied to N=4 lifting
papers.append({
    "id": "2511.09519",
    "title": "Konishi lifts a black hole",
    "authors": ["Jaehyeok Choi", "Eunwoo Lee"],
    "version_read": "v1 (2025-11-12)",
    "category": "untwisted N=4: one-loop Q inferred from the generalized Konishi anomaly; lifting in SO(7)",
    "read_depth": "full text",
    "conventions": {
        "twisting_supercharge": "Q = Q^4_- of N=4 (N=1 language Q = Q_-, psi_+ = Q_+ Phi) (Sec. 2 end).",
        "loop_counting": "'n-loop' = O(g_YM^{2n}) corrections to the dilatation operator; classical Q-cohomology called tree level (footnote 2).",
        "coupling": "tau = theta/(2 pi) + 4 pi i/g_YM^2 (Sec. 1); Tr T_R^A T_R^B = C(R) delta^AB (Sec. 2).",
    },
    "defined": [
        "Q = Q_0 + Q_1 + ...; cohomology computed iteratively; Q-cohomology = cohomology of Q~_1 = pi_0 Q_1 iota_0 on Q_0-cohomology (Sec. 1, via 2306.01039 App. D).",
    ],
    "computed": [
        "Konishi anomaly (adjoint chiral): Dbar^2 Tr Phibar e^{[V,.]} Phi |_anomaly = g_YM^2 C(adj)/(16 pi^2) Tr W_alpha W^alpha (2.1); generalised, adjoint indices: "
        "Dbar^2 Phibar^A (e^{[V,.]} Phi)^B |_anomaly = g_YM^2/(16 pi^2) (W^adj_alpha W^{adj alpha})^AB (2.2); SU(N) (2.3), SO(N) (2.4).",
        "One-loop correction inferred: Q_K = Abar^{ij,kl} g_YM^2/(32 pi^2) (d_psi)_ij (d_Phibar)_kl (2.5), A^{ij,kl} (2.6) (SO(N)), (d_psi)_ij = d/d(psi_+)^{ij}, (d_Phibar)_kl = d/dPhibar^{kl}; "
        "example (unnumbered, after (2.6)): Q Tr Phibar psi_+ = g_YM Tr Phibar dWbar/dPhibar + g_YM^2 (N-2)/(16 pi^2) Tr lambdabar_{alphadot} lambdabar^{alphadot}.",
        "SO(7): Q_K O_f = O_m + Q_0 Lambda with O_m proportional to M_1 - M_2 - 4 M_3 (Sec. 3, unnumbered), M_i (3.2), Q_0-exact combination M_1 + 2M_2 - 10M_3 (3.3); "
        "M_1 - M_2 - 4M_3 independent of the 741 Q_0-exact operators in the sector (J_a = 0, R_i = 3, 8 letters).",
    ],
    "scheme": {
        "claimed_scheme_independent": [],
        "scheme_dependent": ["Discussion (Sec. 4): only the Konishi-anomaly part of Q_1 is used; 'other quantum corrections' for more complicated composites not computed; derivative letters would need further terms (Sec. 2 end)."],
    },
    "explicit_Q1": [
        "Second-order operator Q_K (2.5)-(2.6) on (psi_+, Phibar) pairs producing lambdabar lambdabar; coefficient g_YM^2/(32 pi^2) fixed by the Konishi anomaly. "
        "This is the (beta_I, gamma^J) entry (3.66) of 2512.07771 in untwisted letters.",
    ],
    "other_generators": [],
    "account": "Uses the generalized Konishi anomaly to infer the one-loop correction to Q^4_- on (psi_+, Phibar) pairs (the only entry needed), "
               "with an explicit g_YM^2/(32 pi^2) normalisation, and shows that in SO(7) N=4 SYM it maps the classical fortuitous class O_f onto "
               "the non-Coulomb monotone O_m, lifting both and resolving the SO(7)/Sp(3) mismatch of 2510.24008. It is the only swept paper that "
               "ties the normalisation of Q_1 to g_YM through an untwisted anomaly coefficient.",
    "caveats": ["Partial Q_1 (Konishi entry only, no derivative letters); C(R) is called 'dual Coxeter number of R' in the text (it is the index defined by Tr T_R^A T_R^B = C(R) delta^AB)."],
})

papers.append({
    "id": "2610.22332",
    "title": "Non-Coulomb Q-cohomology and quantum lifting",
    "authors": ["Tanay Kibe", "Sam van Leuven", "Kayleigh Mathieson"],
    "version_read": "v2 (2026-10-20)",
    "category": "application of twisted Q_1 to N=4 lifting (general theorems)",
    "read_depth": "Sec. 1, 2 (start), 5 ((5.1)-(5.5))",
    "conventions": {"twisting_supercharge": "Q = Q^4_- with Chang-Yin letters phi_n = Phi^{4n}, psi_n = -i Psi^n_+, lambda_{alphadot} = Psibar^4_{alphadot}, f = -i F_{++}, D_{alphadot} = D_{+ alphadot} (2.1)."},
    "defined": ["Q = Q_0 + hbar Q_1 + hbar^2 Q_2 + O(hbar^3) (5.1); Q_0^2 = 0, {Q_0,Q_1} = 0, Q_1^2 = -{Q_0,Q_2} (5.2); [Q_1] nilpotent on Q_0-cohomology."],
    "computed": [
        "Pair action used: Q_1 (psi_m)_a (phi^n)_b = -kappa_Q delta^n_m (f f)_{ab}^{de} (lambda_{alphadot})_d (lambda^{alphadot})_e (5.5) (normalisation kappa_Q left free; index placement garbled in extraction).",
        "[Q_1][F_4] = c[S^4], c != 0, S = Tr(lambda_+ lambda_-), SO(7) (1.2)/(5.3).",
        "General (Sec. 6): [Q_1] annihilates every monotone class (multi-trace, finite rank, incl. derivatives); every nonzero [Q_1] image has vanishing Cartan restriction; "
        "key identity Delta^k F = Q_0 H for Q_0-closed F with vanishing Cartan restriction (1.3).",
    ],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": ["Uses Q_1 of 2512.07771; normalisation kappa_Q not needed (5.5)."]},
    "explicit_Q1": ["(5.5), from 2512.07771 (3.66) in Chang-Yin letters."],
    "other_generators": ["Classical commutant supercharges Qbar_{m alphadot} anticommute with Q_1 in Q_0-cohomology (Sec. 5, citing 2512.07771); used to transport pairings (5.4)."],
    "account": "Applies the Budzik-Kulp Q_1 to N=4 with general gauge group: monotones are [Q_1]-closed, [Q_1]-images are non-Coulomb, so "
               "every one-loop pairing has a fortuitous source and a non-Coulomb target; constructs the SO(7) partner of S^4.",
    "caveats": [],
})

papers.append({
    "id": "2606.27955",
    "title": "Quantum black hole cohomologies",
    "authors": ["Seok Kim", "Seongmin Kim", "Siyul Lee", "Jiyoo Park"],
    "version_read": "v2 (2026-06-29)",
    "category": "application of twisted Q_1 to N=4 SU(2) lifting",
    "read_depth": "Sec. 2 ((2.1)-(2.20))",
    "conventions": {
        "coupling": "Q = Q_0 + g_YM^2 Q_1 + O(g_YM^4) after Phi -> g_YM^{-1} Phi (2.2); originally Q = Q_free + g_YM Q_0 + g_YM^2 Q_1 + O(g_YM^3) (2.13).",
        "lie_algebra": "tr(T_A T_B) = kappa delta_AB, tr(T_A [T_B, T_C]) = kappa f_ABC, kappa = dual Coxeter number = N = 2 (2.4).",
        "tree_level": "[Q_0, phi_m] = 0, {Q_0, psi_m} = -i eps_mnl [phi_n, phi_l], {Q_0, lambda_alpha} = 0, [Q_0, f] = i[phi_m, psi_m], [Q_0, D_alpha] = -i[lambda_alpha, } (2.3).",
    },
    "defined": ["Lifting mechanisms: Q_1 O_0 not Q_0-exact, or O_0 = Q_1 Lambda_0 + Q_0 Lambda_1 (2.17)-(2.20); grading d_BMN increased by Q_1."],
    "computed": [
        "Q_1 pair actions (2.5)-(2.10), 'see (B.5-10) of [14]' (= 2512.07771) with overall factors adapted, e.g. Q_1((psi_I)^A d_1^m d_2^n (phi^J)^B) = 1/4 kappa^2 delta_I^J f_ACD f_BCE "
        "[1/(m+n+2)] sum_k sum_l [1/(k+l+1)] binom(m,k) binom(n,l) d_1^k d_2^l d_alpha lambda^D d_1^{m-k} d_2^{n-l} d^alpha lambda^E (2.7); derivative rule (2.11); pairwise (2.12).",
        "SU(2): lightest fortuitous cohomology and its 'hairy' versions unlifted; many heavier core fortuitous ones lifted (abstract, Sec. 3).",
    ],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": ["(2.5)-(2.10) in Chang-Yin letters with kappa normalisation; 'lambda without index' defined by lambda_alpha = D_alpha lambda."],
    "other_generators": ["Bonus U(1) (number of letters): Q_0 charge 1, Q_1 charge 0, Q_n charge 1-n (Sec. 2 after (2.12))."],
    "account": "Translates the Budzik-Kulp Q_1 into Chang-Yin letters with an explicit kappa normalisation and states the coupling order "
               "Q_1 ~ g_YM^2 (2.13); finds lifted and unlifted SU(2) fortuitous classes.",
    "caveats": ["The coupling order g_YM^2 of Q_1 (letter-number preserving) should be compared with Beisert's O(g^2) Q_2 (2 -> 2 sites) rather than his O(g^3) Q_3 (length-changing); not checked in any swept paper."],
})

papers.append({
    "id": "2610.20340",
    "title": "Black hole cohomologies in N=4 SYM",
    "authors": ["Seunggyu Kim", "Jehyun Lee", "Siyul Lee", "Jaewon Song"],
    "version_read": "v1 (2026-10-14)",
    "category": "application of twisted Q_1 to N=4 SU(N), N=2,3,4 lifting",
    "read_depth": "Sec. 2 ((2.1)-(2.5)), grep of lifting statements",
    "conventions": {"coupling": "Q = Q + g_YM^2 Q_1 + O(g_YM^4) (2.1); footnote 2: Q = Q_free + g_YM Q + g_YM^2 Q_1 + O(g_YM^3); kappa = 2 in (2.5)."},
    "defined": [],
    "computed": ["Q_1((psi_n)^A (phi^m)^B) = (kappa_1^2/16) delta_n^m f_ACD f_BCE lambda^D_alpha lambda^{alpha,E} (2.4) with kappa_1 left free; SU(3) brane-fused family not lifted at first quantum order; several black-hole cohomologies lifted (abstract, Sec. 2-3)."],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": ["Overall normalisation of Q_1 'does not play any role' (after (2.5))."]},
    "explicit_Q1": ["(2.4) only; full rules from 2512.07771 and 2606.27955."],
    "other_generators": [],
    "account": "Uses Q_1 to test lifting of BMN-sector black-hole cohomologies; conjectures parity rule.",
    "caveats": [],
})

papers.append({
    "id": "2510.24008",
    "title": "Violation of S-duality in classical Q-cohomology",
    "authors": ["Chi-Ming Chang", "Ying-Hsuan Lin"],
    "version_read": "v1 (2025-10-28)",
    "category": "untwisted N=4 classical cohomology; motivates loop lifting",
    "read_depth": "Sec. 2 ((2.1)-(2.9)), Discussion",
    "conventions": {"tree_level": "Q(Psi) = Psi^2 (2.2), Leibniz (2.3); H^n(g[A], g), A = C[z^+, z^-] (x) Lambda(theta^1,theta^2,theta^3) (2.4)."},
    "defined": [],
    "computed": ["Explicit pair of classes (fortuitous (1/2,1/2,5/2,5/2,5/2) and monotone chiral-ring (0,0,3,3,3)) present for SO(7) and absent for USp(6) at tree level (Table 1, Sec. 3-4)."],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": ["None; Discussion: either S-duality fails or the tree-level exactness conjecture fails, the extra pair acquiring an anomalous dimension from higher-loop or non-perturbative effects."],
    "other_generators": [],
    "account": "Classical cohomology computation whose S-duality mismatch is resolved by the one-loop Q_1 in 2511.09519.",
    "caveats": [],
})

papers.append({
    "id": "2306.04693",
    "title": "Following Black Hole States",
    "authors": ["Kasia Budzik", "Harish Murali", "Pedro Vieira"],
    "version_read": "v2 (2023-09-14)",
    "category": "untwisted N=4 one-loop dilatation; check against twisted Q_1",
    "read_depth": "abstract, footnotes 5 and 10",
    "conventions": {},
    "defined": [],
    "computed": ["Footnote 5: by explicit computation, the first higher-loop correction of 2306.01039 preserves the SU(2) BH cohomology O_0."],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": [],
    "other_generators": ["One-loop dilatation via {Q,Q^dagger} with classical Q; Konishi-like states followed in N (abstract)."],
    "account": "Follows 1/16-BPS states in N with the one-loop dilatation; records a check that the first twisted loop correction preserves the threshold black-hole cohomology.",
    "caveats": [],
})

# ---------------------------------------------------------------- untwisted N=4 dynamic generators
papers.append({
    "id": "hep-th/0310252",
    "title": "The su(2|3) Dynamic Spin Chain",
    "authors": ["Niklas Beisert"],
    "version_read": "v2 (2003-12-09), Nucl. Phys. B682 (2004) 487",
    "category": "untwisted N=4: interacting su(2|3) generators H, Q, S to O(g^5)/O(g^6) by symmetry",
    "read_depth": "Sec. 2.1, 2.4, 3, 4.1-4.4, 5 (start), 6.1",
    "conventions": {
        "fields": "phi_a (a=1,2,3), psi_alpha (alpha=1,2) (2.1); states |A_1...A_L> = Tr W_{A_1}...W_{A_L} (2.3).",
        "expansion": "J(g) = sum_k g^k J_k (3.1); J_k has n+m = k+2 legs (3.7); g^2 = g_YM^2 N/(8 pi^2) (6.7), alpha_1 = 1 (6.8).",
        "subsector": "Delta_0 = p + q/2 + 3q_2/2 (2.27); H = D, deltaH = D - D_0 (2.30).",
    },
    "defined": ["Interaction symbols {A_1..A_n / B_1..B_m} (3.2)-(3.4); gauge identity (3.6); similarity transformations J -> T J T^{-1} (3.13), (3.14), (3.39), (4.4); coupling redefinitions g -> alpha_1 g (3.40), g -> f(g) (4.9), (6.1)."],
    "computed": [
        "Tree level (3.10)-(3.12); O(g): (Q_1)_a^alpha = (1/sqrt2) alpha_1 e^{i beta_1 + i beta_2} eps^{alpha beta} eps_{abc} {beta / bc} (replaces a fermion beta by two bosons bc), "
        "(S_1)_alpha^a = (1/sqrt2) alpha_1 e^{-i beta_1 - i beta_2'} eps^{abc} eps_{alpha beta} {bc / beta} (3.37), symbol {A/B} replaces A by B (3.2)-(3.4); "
        "H_2, H_3 (3.37); constraints (3.35) incl. c_1 = -c_4 fixing H_2 at O(g^3).",
        "O(g^2): (Q_2)_a^alpha, (S_2)_alpha^a with free parameters gamma_1..gamma_4 (Table 1); H_4 (Table 1).",
        "O(g^3): Q_3, S_3 exist with 6 parity-even structures each (Sec. 4.1; Table 2 row Q_{k-2}, k=5); fixed by [Q, deltaH] = [S, deltaH] = O(g^6) (4.2), (5.1); "
        "'the resulting generators are too lengthy to be displayed' (Sec. 5): not printed.",
        "Three-loop H_6 in the bosonic subsector (Table 3); BMN scaling fixes alpha_3 = sigma_1 = sigma_2 = 0 (6.6).",
    ],
    "scheme": {
        "claimed_scheme_independent": ["Energies (eigenvalues) up to coupling redefinition; parameters gamma_{1..4} 'rotate only the eigenstates' (Sec. 4.4)."],
        "scheme_dependent": [
            "Similarity transformations T(g) (4.4)-(4.6): J_2 -> J_2 + [T_2, J_0]; at each order the free coefficients equal the number of H_{k-2} structures (Table 2).",
            "Coupling redefinition g -> f(g): H_4 -> alpha_1^4 H_4 + 2 alpha_1 alpha_3 H_2 (4.10); changes energies (6.2), fixed by BMN scaling (6.5)-(6.6).",
            "Gauge transformations (interactions vanishing on cyclic states, (3.6)) (Sec. 4.2).",
        ],
    },
    "explicit_Q1": [
        "No twisted Q_1. Leg counting: twisted Q_n maps n+1 letters to 2 (2306.01039 (3.44)-(3.45)), i.e. n+3 legs, which is Beisert's order k = n+1 by (3.7): "
        "twisted Q_0 <-> Beisert O(g) (Q_1), twisted Q_1 <-> Beisert O(g^2) (Q_2, 2 -> 2 sites, Table 1), twisted Q_2 <-> Beisert O(g^3) (Q_3, 3 <-> 2). "
        "This dictionary is an inference of this sweep, not a statement of either paper.",
    ],
    "other_generators": [
        "S_k listed alongside Q_k at every order (3.37), Table 1, Table 2: one-loop-order corrections of special supercharges exist in this sector.",
        "{S,Q} = L + R + delta(1/3 H + 1/6 deltaH) + O(g^4) (4.3) satisfied automatically; (3.41)-(3.42).",
    ],
    "account": "Fixes the planar su(2|3) generators (Hamiltonian, supercharges Q and special supercharges S) order by order in g from the "
               "algebra plus Feynman-diagram structure, up to similarity transformations and coupling redefinitions, which are the "
               "spin-chain form of the scheme freedom. O(g) and O(g^2) supercharges are printed; O(g^3) ones are only counted. Its O(g^2) "
               "Q_2 is the leg-count analogue of the twisted Q_1.",
    "caveats": ["The su(2|3) generators Q_a^alpha are not identical to Q^4_-; mapping to the 1/16-BPS supercharge needs a sector dictionary (not given)."],
})

papers.append({
    "id": "hep-th/0307015",
    "title": "The Complete One-Loop Dilatation Operator of N=4 Super Yang-Mills Theory",
    "authors": ["Niklas Beisert"],
    "version_read": "v5 (2011-02-27)",
    "category": "untwisted N=4: one-loop dilatation on all local operators",
    "read_depth": "Sec. 2 ((2.7)-(2.24)), Sec. 5 ((5.15)), Sec. 6",
    "conventions": {"expansion": "J(g) = sum g^k J_k (2.7); D(g) = D + g_YM^2 N/(8 pi^2) H + O(g^3) (2.11)."},
    "defined": ["[J_0, D_2] = 0 (2.10): D_2 invariant under classical psl(4|4); chirality B conserved at one loop, broken at higher loops by the Konishi anomaly (Sec. 2 after (2.10))."],
    "computed": ["D(g) = D - g_YM^2/(8 pi^2) sum_j h(j) (P_j)^{AB}_{CD} :Tr[W_A, Wcheck^C][W_B, Wcheck^D]: + O(g^3) (2.22), h(j) harmonic numbers (2.21); planar H_12 = sum_j 2 h(j) P_{12,j} (2.24); "
                 "Konishi submultiplet energies E = 4C_2 = 6C_1 = 6C_1 - 2C_0 (6.1)."],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": [],
    "other_generators": [
        "Sec. 6: the generators receiving corrections are P, K, Q, Qdot, S, Sdot; closure should constrain them (not computed).",
        "(5.15): conjecture deltaS(g) ~ g Tr Psi Q, deltaSdot(g) ~ g Tr Psibar Q for the first correction to the superboosts.",
    ],
    "account": "The one-loop dilatation generator D_2 on all (planar and non-planar) local operators of N=4, fixed by classical psl(4|4) "
               "invariance and harmonic-number coefficients (2.22). It is the standard of comparison for the project's D_1 on gauge-invariant "
               "operators; the paper states that P, K, Q, S also receive corrections but does not compute them.",
    "caveats": [],
})

papers.append({
    "id": "hep-th/0511109",
    "title": "N=4 SYM to Two Loops: Compact Expressions for the Non-Compact Symmetry Algebra of the su(1,1|2) Sector",
    "authors": ["Benjamin I. Zwiebel"],
    "version_read": "v3 (2006-03-15), JHEP 0602:055 (per 2310.20086 ref [32])",
    "category": "untwisted N=4: O(g^2) psu(1,1|2) generators incl. P, K, Q, S and O(g^3) dynamic supercharges",
    "read_depth": "Sec. 1, 2, 3, 4",
    "conventions": {
        "generators": "psu(1,1|2): J_0 = -L + 2D_0 + deltaD, J^{++} = P_22, J^{--} = K_22, eight Q/S (2.3); psu(1|1)^2: T^{->+} = Qdot_14, T^{<--} = Sdot_14, T^{->-} = S_11, T^{<-+} = Q_11 (2.5) (upper/lower index placement not recoverable from the extracted text); "
                      "[L, T^{+-}] = +-T^{+-} (2.6); {T^{->+}, T^{<--}} = 1/2 deltaD (2.7).",
        "coupling": "g^2 = g_YM^2 N/(8 pi^2) (2.15); T's expand in odd powers of g (Sec. 1).",
    },
    "defined": ["Harmonic-number generator h (3.3); two-site generator x (3.4); similarity J -> U J U^{-1}, U = 1 + g^2 y + ..., J_2 -> J_2 + [J_0, y] (3.7)."],
    "computed": [
        "O(g^2): X^{+-}_2 = +-[X^{+-}_0, x] + [X^{+-}_0, y] for X^{+-} in {J^{++}, J^{--}, eight Q's} (3.6): the one-loop corrections of P, K and the psu(1,1|2) supercharges are commutators.",
        "O(g^3): T^{+-}_3 = +-[T^{+-}_1, x] + [T^{+-}_1, y] + alpha T^{+-}_1 (4.1), alpha <-> g -> g + alpha g^3 (4.2); deltaD_4 from (4.4)-(4.5); deltaD_4 -> deltaD_4 + 2 alpha deltaD_2 (4.6), alpha = 0 to match field theory.",
    ],
    "scheme": {
        "claimed_scheme_independent": [],
        "scheme_dependent": ["y (similarity), alpha (coupling), gauge transformations (Sec. 3.2, 4.2); homogeneous solutions neither found nor excluded (3.8)-(3.9), (4.3)."],
    },
    "explicit_Q1": ["Spin-chain O(g^3) dynamic supercharge (4.1); O(g^2) psu(1,1|2) supercharges (3.6). Not the twisted Q_1; T's have no O(g^2) term because of (2.6)."],
    "other_generators": ["K_22 = J^{--} and P_22 = J^{++} corrected at O(g^2) by (3.6): an explicit one-loop special-conformal correction in this sector."],
    "account": "In the su(1,1|2) sector the one-loop (O(g^2)) corrections of all psu(1,1|2) generators, including the special conformal K and "
               "special supercharges, are commutators with a universal two-site generator x plus a similarity transformation, and the "
               "O(g^3) dynamic supercharges are T_3 = +-[T_1,x] + [T_1,y] + alpha T_1. 2310.20086 cites this as 'Q_1 takes the form [Q_0,-]'.",
    "caveats": ["Restricted sector; non-planar lift conjectured (Sec. 4.4). Correspondence of T^{<-+} = Q^1_1 with the 1/16-BPS Q^4_- up to R-relabelling is plausible but unchecked; if so, twisted Q_1 (2 -> 2, O(g_YM^2)) must vanish on this sector by (2.6)."],
})

papers.append({
    "id": "0806.1786",
    "title": "Iterative Structure of the N=4 SYM Spin Chain",
    "authors": ["Benjamin I. Zwiebel"],
    "version_read": "v2 (2008-06-19), JHEP 0807:114",
    "category": "untwisted N=4: all-order ansatz for coupling dependence of psu(1,1|2) generators",
    "read_depth": "Sec. 3.1-3.2 ((3.1)-(3.8))",
    "conventions": {},
    "defined": ["Ansatz: d J^{+-}/d lambda = +-[J^{+-}(lambda), X(lambda)] (3.1); X_ab = {Qhat_b, [Shat_a, h]}, X = eps^{ab} X_ab + 1/2 [H, h] (3.2); equations for h (3.3); J^{+-}(lambda) = U(-+X,lambda) J_0^{+-} U^dagger(+-X,lambda), U = P exp int X (3.4)."],
    "computed": ["Planar three-loop dilatation generator in the psu(1,1|2) sector (Sec. 4.5) passing spectral tests (Sec. 5)."],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": [],
    "other_generators": ["All raising/lowering generators (P, K, Q, S of psu(1,1|2)) deformed by X to all orders (3.1)-(3.4)."],
    "account": "Extends hep-th/0511109 to an all-order iterative ansatz in which the coupling dependence of every psu(1,1|2) generator is a flow generated by X.",
    "caveats": [],
})

papers.append({
    "id": "1305.6314",
    "title": "1/16 BPS States in N=4 SYM",
    "authors": ["Chi-Ming Chang", "Xi Yin"],
    "version_read": "v2 (2013-09-05), Phys. Rev. D 88 (2013) 106005",
    "category": "untwisted N=4: classical (O(g)) Q^4_- on BPS letters, relative Lie algebra cohomology",
    "read_depth": "Sec. 1-3, Sec. 4 (table)",
    "conventions": {
        "supercharge": "Q = Q^4_-, S = S_4^- = Q^dagger, 2{Q,Q^dagger} = Delta = E - 2J - H_1 - H_2 - H_3 (2.1).",
        "letters": "phi_n = Phi^{4n}, psi_n = -i Psi^n_+, lambda_{alphadot} = Psibar_{alphadot 4}, f = -i F_{++}, D_{alphadot} = D_{+ alphadot} (2.2)-(2.3); "
                   "[D_alphadot, D_betadot] = eps f, D_alphadot lambda^alphadot = [phi_n, psi_n] (2.4).",
    },
    "defined": ["Q on letters: [Q,phi_n] = 0, {Q,psi_n} = -i eps_nmp [phi_m, phi_p], {Q,lambda} = 0, [Q,f] = i[phi_n,psi_n] (2.5); [Q, D zeta] = -i[lambda, zeta] + D Q zeta (2.6); "
                "superfield Psi(z,theta) = -i[lambda(z) + 2 theta^n phi_n(z) + eps_mnp theta^m theta^n psi_p(z) + 4 theta^1 theta^2 theta^3 f(z)] (2.9); {Q, Psi} = Psi^2 (2.10); Psi|_{z=theta=0} = 0 (2.11)."],
    "computed": ["Q-cohomology = H^*(G_N, sl_N; C), G_N = C[z_+, z_-] (x) Lambda[theta^1,theta^2,theta^3] (x) sl_N (3.6); infinite-N cohomology = multi-gravitons (Sec. 4)."],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": ["None (tree level = twisted Q_0)."],
    "other_generators": [],
    "account": "The standard formulation of the classical 1/16-BPS cohomology: (2.10) is the untwisted counterpart of the twisted Q_0 C = 1/2[C,C] "
               "(2512.07771 (3.54)); its superfield (2.9) is the counterpart of (3.52) there, with different normalisations.",
    "caveats": ["Phases (-i) differ from GGKM (2.6) and from 2512.07771 (3.49)-(3.50); a sign/phase dictionary is needed before comparing Q_1 normalisations."],
})

papers.append({
    "id": "0803.4183",
    "title": "Comments on 1/16 BPS Quantum States and Classical Configurations",
    "authors": ["Lars Grant", "Pietro A. Grassi", "Seok Kim", "Shiraz Minwalla"],
    "version_read": "v1 (2008-03-28), JHEP 05 (2008) 049",
    "category": "untwisted N=4: classical Q^4_- cohomology; non-renormalisation conjecture",
    "read_depth": "Sec. 1, 2.1-2.2",
    "conventions": {"letters": "Phi^{4m} = phibar^m, Phi^{4p} = phi_p, Psi^4_alpha = lambda_alpha, Psi^m_alpha = psi_{m alpha} (2.3); Delta = 2{Q_-, S^-} = E - 2J - H_1 - H_2 - H_3 (2.4)."},
    "defined": ["[Q_-, phibar_n] = 0, {Q_-, psi_{n+}} = eps_nmp [phibar_m, phibar_p], {Q_-, lambdabar_betadot} = 0, [Q_-, f_{++}] = i[phibar_m, psi_{m+}] (2.6)."],
    "computed": ["Cohomology below energy N = Fock space of single-trace cohomology (Sec. 2.1)."],
    "scheme": {"claimed_scheme_independent": ["Conjecture (Sec. 1): the supersymmetric spectrum equals the 1/16-BPS states of the finite-N one-loop Beisert Hamiltonian."], "scheme_dependent": []},
    "explicit_Q1": [],
    "other_generators": [],
    "account": "Formulates 1/16-BPS counting as cohomology of the classical Q and states the one-loop-exactness conjecture later disproved (2510.24008, 2511.09519).",
    "caveats": [],
})

papers.append({
    "id": "2209.06728",
    "title": "Words to describe a black hole",
    "authors": ["Chi-Ming Chang", "Ying-Hsuan Lin"],
    "version_read": "v4 (2023-02-16), JHEP 02 (2023) 109",
    "category": "untwisted N=4: non-renormalisation argument (assuming Leibniz)",
    "read_depth": "Sec. 3 ((3.1)-(3.5), footnotes 8-12)",
    "conventions": {"expansion": "D = D_0 + g_YM^2 N/(8 pi^2) H + O(g_YM^3) (3.1); H proportional to {Q, Q^dagger} on BPS words (3.3); leading Q = Tr Psi Psi d/dPsi, Q^dagger = Tr Psi d/dPsi d/dPsi (3.2)."},
    "defined": [],
    "computed": ["H^*(G_N, sl_N; C[[g_YM]]) isomorphic to H^*(G_N, sl_N; C) (x) C[[g_YM]] (3.4)-(3.5), hence perturbative invariance of the 1/16-BPS spectrum, under the assumption that Q obeys Leibniz in perturbation theory."],
    "scheme": {
        "claimed_scheme_independent": [],
        "scheme_dependent": ["Footnote 10: a regularisation scheme can be chosen such that Q on the BPS superfield is undeformed (DRED, or their ref. [62]); footnote 11: the Konishi anomaly is the only known obstruction to Leibniz, stated absent in N=4 with their ref. [64]."],
    },
    "explicit_Q1": ["None; the argument is invalidated by Q_1 not obeying Leibniz (2512.07771 footnote 21; 2511.09519)."],
    "other_generators": [],
    "account": "Gives the perturbative non-renormalisation argument for the classical Q-cohomology, explicitly conditional on a Leibniz-rule Q; "
               "the twisted Q_1 violates exactly this assumption.",
    "caveats": ["Reference [64] (support for 'Konishi anomaly absent in N=4') is printed as A. Kostyuk, H. Stoecker, nucl-th/0501077 (an unrelated heavy-ion paper); "
                "plausibly a misprint for hep-th/0501077 (Eden-Jarczak-Sokatchev-Stanev), which however argues that in N=4 the classical/quantum split of the Konishi anomaly is scheme-meaningless, not that it is absent."],
})

# ---------------------------------------------------------------- Konishi anomaly
papers.append({
    "id": "hep-th/0211170",
    "title": "Chiral Rings and Anomalies in Supersymmetric Gauge Theory",
    "authors": ["Freddy Cachazo", "Michael R. Douglas", "Nathan Seiberg", "Edward Witten"],
    "version_read": "v2 (2002-12-17), JHEP 12 (2002) 071",
    "category": "Konishi anomaly: generalized (chiral-ring) form",
    "read_depth": "Sec. 3 ((3.1)-(3.19))",
    "conventions": {"setting": "N=1 U(N) with adjoint Phi; J = Tr Phibar e^{ad V} Phi (3.1); identities written with Dbar^{alphadot} (footnote 6)."},
    "defined": ["Generalized currents J_f = Tr Phibar e^{ad V} f(Phi, W_alpha) (3.11)."],
    "computed": [
        "Dbar^2 J = 1/(32 pi^2) tr_adj (ad W_alpha)(ad W^alpha) (3.3); with superpotential Dbar^2 J = Tr Phi dW/dPhi + N/(16 pi^2) Tr W_alpha W^alpha - 1/(16 pi^2) Tr W_alpha Tr W^alpha (3.4).",
        "Dbar^2 J_f = Tr f(Phi,W_alpha) dW/dPhi + 1/(32 pi^2) sum_{i,j} [W_alpha, [W^alpha, df/dPhi_{ij}]]_{ji} (3.14); Ward identity (3.15).",
    ],
    "scheme": {
        "claimed_scheme_independent": ["(3.4)/(3.14) exact for chiral-ring purposes: no perturbative or non-perturbative chiral corrections (argument after (3.14) and (3.4))."],
        "scheme_dependent": ["Higher-loop contributions to the anomaly are renormalisation-scheme dependent but non-chiral (text after (3.5))."],
    },
    "explicit_Q1": ["Untwisted chiral-ring analogue of Q_1: the second-order (anomaly) term of (3.14)."],
    "other_generators": [],
    "account": "The generalized Konishi anomaly: anomalous Qbar-variation of Tr Phibar e^V f(Phi,W) with a one-loop coefficient 1/(32 pi^2), exact in the chiral ring. "
               "2512.07771 calls its semi-chiral extension the 'twice-generalized Konishi anomaly'.",
    "caveats": [],
})

papers.append({
    "id": "hep-th/0104016",
    "title": "Properties of the Konishi multiplet in N=4 SYM theory",
    "authors": ["Massimo Bianchi", "Stefano Kovacs", "Giancarlo Rossi", "Yassen S. Stanev"],
    "version_read": "v1 (2001-04-02)",
    "category": "Konishi anomaly in N=4; anomalous term in a Konishi descendant",
    "read_depth": "Sec. 2-3 ((1)-(11))",
    "conventions": {"coupling": "standard g (footnote 3: g of their earlier papers = 2g here); N=1 superfield action (1); Fermi-Feynman gauge."},
    "defined": ["K_1 = 1/2 sum_i :tr(phi_i phi_i): naive (4)."],
    "computed": [
        "1/4 Dbar^2 :tr(e^{-2gV} Phibar_I e^{2gV} Phi^I): = tr(Phi^I dW/dPhi^I) + 6 g^2 N/(32 pi^2) tr(W^alpha W_alpha) (7).",
        "Level-two descendant (K_10*)_AB = 3 sqrt2 g t^{ijk}_AB tr(phi_i phi_j phi_k) + 6 g^2 N/(32 pi^2) tr(psibar_{alphadot A} psibar_B^{alphadot}) (9); "
        "protected (E_10*)_AB = -tr(psibar psibar) + sqrt2 g t^{ijk}_AB tr(phi phi phi) (10), descendant of Q_20' (11).",
    ],
    "scheme": {"claimed_scheme_independent": [], "scheme_dependent": []},
    "explicit_Q1": ["Untwisted: the g^2 N tr(psibar psibar) term of (9) is the one-loop (anomalous) part of the double Q-variation of K_1; it corresponds to the lambdabar lambdabar output of twisted Q_1."],
    "other_generators": [],
    "account": "Constructs the N=4 Konishi multiplet including the anomaly-induced bilinear fermion term in its 10-bar descendant (9), "
               "distinguishing it from the protected descendant (10) of the stress-tensor multiplet.",
    "caveats": [],
})

papers.append({
    "id": "hep-th/0501077",
    "title": "Operator mixing in N=4 SYM: The Konishi anomaly revisited",
    "authors": ["B. Eden", "C. Jarczak", "E. Sokatchev", "Ya. S. Stanev"],
    "version_read": "v2 (2005-02-03)",
    "category": "Konishi anomaly in N=4 as mixing; scheme dependence",
    "read_depth": "abstract, Sec. 1-2.1 ((1)-(4))",
    "conventions": {"setting": "N=1 superspace, W = (g/3) eps_IJK Tr(Phi^I Phi^J Phi^K); regulator: supersymmetric dimensional reduction (SSDR) for the g^4 computation (not usable in this project)."},
    "defined": ["Dbar^2 K = Tr Phi dW/dPhi + g^2/(32 pi^2) Tr W^alpha W_alpha = B + g^2/(32 pi^2) F (1); renormalised form [Dbar^2 K]_R = a(g)[B]_R + b(g)[F]_R (2)."],
    "computed": ["Two-point diagonalisation conditions (4); mixing matrix to O(g^4): protected O_10 unrenormalised; Konishi descendant K_10 changes form at every loop; at one loop the F coefficient of (1) is recovered, at two loops it becomes a divergent factor (Sec. 1)."],
    "scheme": {
        "claimed_scheme_independent": ["Only K_10 (with the Konishi anomalous dimension) can appear on the right of the anomaly equation (abstract)."],
        "scheme_dependent": ["a(g), b(g) depend on operator normalisation, i.e. the subtraction scheme (after (2)); in a superconformal scheme the split into 'classical' and 'quantum' anomaly terms is not possible and one-loop exactness of the Konishi anomaly is 'out of context' (abstract)."],
    },
    "explicit_Q1": [],
    "other_generators": [],
    "account": "Recasts the N=4 Konishi anomaly as superconformal operator mixing of the trilinear B and bilinear F operators and argues that "
               "the separation into classical and one-loop quantum parts is scheme dependent; the scheme-independent content is the "
               "conformal (anomalous-dimension) structure. Directly relevant to R3's 'what is scheme independent'.",
    "caveats": ["Uses DRED/SSDR; cite only for statements, not for regulator choices."],
})

# ---------------------------------------------------------------- abstract-level / peripheral
peripheral = [
    {"id": "2311.04304", "title": "Unravelling the Holomorphic Twist: Central Charges", "authors": ["Pieter Bomans", "Jingxiang Wu"],
     "read_depth": "abstract only",
     "note": "Binary and ternary lambda-brackets of the twisted N=1 SCFT; central extensions of higher Kac-Moody/Virasoro encode a, c and flavour central charges. Relevant to anomalies of the twisted conformal symmetry (K-type), not to Q_1 on operators."},
    {"id": "2401.08113", "title": "Feynman Graph Integrals on C^d", "authors": ["Minghao Wang"],
     "read_depth": "abstract only",
     "note": "UV finiteness of holomorphic graph integrals via compactification of metric-graph moduli; cited by 2512.07771 for finiteness."},
    {"id": "1906.04221", "title": "Twisted characters and holomorphic symmetries", "authors": ["Ingmar Saberi", "Brian R. Williams"],
     "read_depth": "abstract only",
     "note": "Classical twist of chiral/vector multiplets (betagamma with superpotential); cited by 2512.07771 for the tree-level dictionary (3.6)."},
    {"id": "2606.05314", "title": "Mass-Flow Invariance of Q-Cohomology in BMN Matrix Quantum Mechanics",
     "authors": ["Chi-Ming Chang", "Zhengyuan Du", "Sarthak Duary", "Kangning Liu", "Yi-Xiao Tao"], "read_depth": "abstract only",
     "note": "mu-derivative of the BMN supercharges is an adjoint action (similarity) => Q-cohomology mu-independent; algebraic, not a loop computation."},
    {"id": "0905.3738", "title": "Exacting N=4 Superconformal Symmetry", "authors": ["Till Bargheer", "Niklas Beisert", "Wellington Galleas", "Florian Loebbert", "Tristan McLoughlin"],
     "read_depth": "abstract only",
     "note": "Leg-changing corrections to the free superconformal generators acting on tree amplitudes (holomorphic anomaly); amplitudes, not local operators."},
    {"id": "1111.0083", "title": "From Scattering Amplitudes to the Dilatation Generator in N=4 SYM", "authors": ["Benjamin I. Zwiebel"],
     "read_depth": "abstract only",
     "note": "Leading-order interactions acting on L sites from tree amplitudes; leading-order generators."},
    {"id": "Konishi 1984 / Konishi-Shizuya 1985", "title": "Anomalous Supersymmetry Transformation of Some Composite Operators in SQCD (PLB 135 (1984) 439); Functional Integral Approach to Chiral Anomalies in Supersymmetric Gauge Theories (Nuovo Cim. A 90 (1985) 111)",
     "authors": ["K. Konishi", "K.-i. Shizuya"], "read_depth": "not read (not on arXiv); bibliographic data from 2306.01039 refs [30],[31]",
     "note": "Original Konishi anomaly."},
]

cross = {
    "what_Q1_is_in_the_twisted_literature": [
        "Definition: Q_n O = (n+1)!^{-1}{I,...,I,O} (2306.01039 (3.33); 2512.07771 (1.6)); one loop: only the triangle; Q_1 is a second-order operator on pairs of letters (2512.07771 (3.7)).",
        "Regulator: Schwinger-time UV cutoff on every propagator, integrals finite without counterterms (2207.14321 Sec. 2; 2512.07771 (2.42), App. A; 1809.02661 Thm 3.4).",
        "N=4 result: 2512.07771 (3.64)-(3.70) (+ App. B for derivatives); first entry (3.66) also in 2511.09519 (2.5) with physical normalisation g_YM^2/(32 pi^2).",
        "Conclusions to reproduce (R3): (i) Q_1 on pairs as (3.64)-(3.69)/(3.70); (ii) Q_1(gamma gamma) = 0, so Q_1 kills Tr gamma^{(I1}...gamma^{In)} and (via psl(3|3)) their descendants; "
        "(iii) {Q_0,Q_1} = 0, Q_1^2 = -{Q_0,Q_2} != 0; (iv) SO(7): Q_1 O_f = O_m mod Q_0 (2511.09519); (v) [Q_1] kills monotones, images non-Coulomb (2610.22332); "
        "(vi) N=1 SYM: Q_1 Tr b^2 proportional to d_alpha(Tr b d^alpha c) (2306.01039 (3.40)).",
    ],
    "scheme_statements": [
        "Q_n scheme dependent beyond first order; schemes related by L_infinity quasi-isomorphisms (2306.01039 Sec. 3.2.1; 2512.07771 Sec. 2.2.2).",
        "Change of propagator => operator redefinition (2207.14321 (B.6)).",
        "Lifting pattern scheme independent (2512.07771 Sec. 1, footnote 9).",
        "Spin-chain side: similarity transformations and coupling redefinitions (hep-th/0310252 (3.39), (4.4)-(4.10); hep-th/0511109 (3.7), (4.1)-(4.2)).",
        "Konishi anomaly: classical/quantum split scheme dependent in N=4 (hep-th/0501077 abstract, (2)); higher-loop anomaly terms scheme dependent but non-chiral (hep-th/0211170 after (3.5)).",
        "No swept paper displays an explicit family 'Q_1 -> Q_1 + [Q_0, h] + ...' parametrised by finite counterterms of the untwisted theory; the twisted papers state the quasi-isomorphism principle only.",
    ],
    "normalisation_map": [
        "Twisted: hbar loop counting only (2512.07771 (1.3), (2.22)); propagator P with dbar P = pi^2 delta^4 (2.32), bracket measure (2 pi i)^{-2} (2.27); Killing-normalised f (3.12).",
        "Untwisted Q^4_- conventions: Chang-Yin (2.5) (phases -i), GGKM (2.6), 2606.27955 (2.3) (-i), 2610.20340 (2.2) (no i, different field phases, footnote after (2.5)).",
        "Coupling order: twisted Q_1 at g_YM^2 relative to the unrescaled expansion Q = Q_free + g_YM Q_0 + g_YM^2 Q_1 + O(g_YM^3) (2606.27955 (2.13); 2610.20340 footnote 2).",
        "Physical coefficient: only 2511.09519 (2.5) fixes Q_1 (Konishi entry) in g_YM units: g_YM^2/(32 pi^2) with A^{ij,kl} (2.6).",
        "Relative factor between 2306.01039 (C.4) and 2512.07771 (3.35)/(3.64) not reconciled (generator normalisation + footnote 23 factor 1/2).",
    ],
    "leg_counting_dictionary_inferred": "Twisted Q_n maps n+1 letters to 2 (2306.01039 (3.44)-(3.45)), i.e. n+3 legs; Beisert's J_k has k+2 legs (hep-th/0310252 (3.7)). "
                                        "Hence twisted Q_n <-> spin-chain order g^{n+1}: Q_0 <-> O(g) dynamic supercharge, Q_1 <-> O(g^2) 2->2 terms, Q_2 <-> O(g^3). "
                                        "Inferred here; no swept paper compares the twisted Q_1 with Beisert's Q_2 or Zwiebel's generators (no citation of either in 2512.07771, 2306.01039, 2606.27955, 2610.22332).",
    "other_generators_at_one_loop_known": [
        "D: complete one-loop D_2 on all local operators (hep-th/0307015 (2.22)); su(2|3) to three loops (hep-th/0310252); psu(1,1|2) to two/three loops (hep-th/0511109 (4.5), 0806.1786).",
        "S (special supercharges): su(2|3) S_1 (O(g)) and S_2 (O(g^2)) printed (hep-th/0310252 (3.37), Table 1), S_3 counted only; psu(1,1|2) S's at O(g^2) and Sdot^4_1, S^1_1 at O(g^3) (hep-th/0511109 (3.6), (4.1)); conjecture (5.15) of hep-th/0307015.",
        "K (special conformal): psu(1,1|2) K_22 = J^{--} at O(g^2) (hep-th/0511109 (3.6)); hep-th/0307015 Sec. 6 states K is corrected, not computed.",
        "P, M, R: hep-th/0307015 Sec. 6 lists P among corrected generators in its basis; 2403.13049 footnote 13/(2.31) expects no corrections to translations acting on couplings.",
        "Twisted side: tree-level psl(3|3) only; one-loop R_1 not computed (2512.07771 footnote 24); current non-closure (2509.16737 (3.51)-(3.52)); holomorphic stress tensor (2306.01039 (3.24)).",
        "Not found in the sweep: one-loop S, K (or D on letters/BV components) as components of an untwisted equivariant BV cochain along psl(4|4) ghosts, or any full-theory (all-sector, non-planar) one-loop S or K on local operators.",
    ],
    "searches_run": [
        "WebSearch (extended): one-loop correction supercharge holomorphic twist N=4 SYM Q-cohomology arXiv 2026",
        "WebSearch (extended): semi-chiral / holomorphic twist loop corrections supercharge two-loop Q_2 N=4 fortuitous lifting 2026",
        "WebSearch: semi-chiral ring one-loop Q_1 holomorphic BF SQCD Konishi 2025",
        "INSPIRE citation lists: refersto:recid:3090379 (2512.07771; 12 hits), refersto:recid:2665272 or 3082243 (2306.01039, 2511.09519; 63 hits); scanned titles, opened the relevant ones.",
    ],
}

doc = {"meta": meta, "leads_verification": leads, "papers": papers, "peripheral_or_abstract_only": peripheral, "cross_paper_findings": cross}

with open(OUT, "w") as fh:
    json.dump(doc, fh, indent=1, ensure_ascii=False)
print(OUT, len(papers), "papers,", len(peripheral), "peripheral")
