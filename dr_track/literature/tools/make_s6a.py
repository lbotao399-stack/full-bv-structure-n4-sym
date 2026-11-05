import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rec_io import write_record

OUT = "/Users/libotao/N4_BV_DR/literature/parts/S6"
PROG = "/Users/libotao/N4_BV_DR/literature/parts/S6_PROGRESS.txt"
SEC = "6a BV master equation / quantum action principle in DR and regularization-free BV"

R = []

R.append(dict(
 id="S6.01_Barnich1998", section=SEC, authors="G. Barnich",
 title="On the quantum Batalin-Vilkovisky formalism and the renormalization of non linear symmetries",
 year=1998, arxiv="hep-th/9807210", journal="JHEP 9903 (1999) 010",
 claim="Builds an extended antifield formalism containing all observables, with an extended master equation. If an invariant regularization exists, subtracted and effective actions satisfy the extended master equation to all orders independently of power counting. For dimensional regularization the evanescent breaking of the regularized master equation is coupled to a global fermionic source rho*, and the operator Delta_d = (tau/(-i hbar)) d^L/d rho* takes the role of the BV Laplacian; the formal BV Delta term vanishes in DR because delta(0)=0.",
 loop_order="all orders (structural); explicit bookkeeping of one- and two-loop anomalies",
 theory="general Lagrangian gauge theories with non-linear (local or global) symmetries",
 scheme="Dimensional regularization with the conventions of the paper's ref. [32]; regulator parameter tau; breaking theta_tau = (1/2tau)(S_tau,S_tau) + (1/tau) Delta_c S_tau coupled to global source rho* (Sec. 2.1); renormalized quantum action principle taken from the cited literature.",
 symmetry_argument="Quantum action principle applied to the regularized extended master equation: (1/2)(W,W) + Delta_c W - i hbar Delta_d W = A (eq. 131), consistency condition (W,A) + Delta_c A - i hbar Delta_d A = 0 (eq. 132); anomalies are BRST cocycles at ghost number 1.",
 symmetry_status=dict(
  proved=["With an invariant regularization (theta_tau = 0) the subtracted and effective actions obey the extended master equation (Sec. 2.2, eqs. 72-73).",
          "In DR the would-be BV Laplacian term drops out (delta(0)=0); Delta_d is nilpotent, anticommutes with Delta_c and is a graded derivation of the antibracket (Sec. 2.8, eqs. 129-132).",
          "One-loop divergences are cocycles at ghost number 0, the one-loop anomaly a cocycle at ghost number 1 (Theorem 6, Sec. 2.4, eqs. 85-89); two-loop anomaly is local and obeys a consistency condition (Sec. 2.6, eqs. 109-115)."],
  assumed=["Renormalized quantum action principle in dimensional renormalization (taken from refs. [30-34] of the paper).",
           "Knowledge of the relevant local BRST cohomology."],
  checked_partially=[]),
 evidence=["Sec. 2.1: theta_tau definition; eqs. (68)-(71) regularized master equation with source rho* and its QAP form.",
           "Sec. 2.2: invariant case, eqs. (72)-(73).",
           "Sec. 2.4-2.6: one- and two-loop anomalies, eqs. (85)-(115).",
           "Sec. 2.8: quantum BV Delta operator in DR, eqs. (129)-(132); Delta_d defined after eq. (130).",
           "Footnote 3 (Introduction): evanescent breaking terms of the regularized master equation are what make Delta non-trivial in dimensional renormalization."],
 red_flags=["Delta_d exists only before the regulator is removed and with rho* kept non-zero; the formal BV Delta has no place in the final renormalized theory (Sec. 2.8).",
            "No supersymmetric application."],
 cross_checks=["S3.02_BreitenlohnerMaison1977 (quantum action principle in DR)", "S6.02_DPT1996 (BPHZ version)", "S6.04_FredenhagenRejzner2011 (regularization-free version)"],
 source_access="full text via ar5iv", confidence="high",
 notes="Gives the precise dictionary needed for a DR-based QME: in DR the entire 'quantum' part of the QME is the evanescent breaking insertion; this is the object that a BV computation of DRED or BMHV breaking has to evaluate."))

R.append(dict(
 id="S6.02_DPT1996", section=SEC, authors="F. De Jonghe, J. Paris, W. Troost",
 title="The BPHZ renormalised BV master equation and two-loop anomalies in chiral gravities",
 year=1996, arxiv="hep-th/9603012", journal="Nucl. Phys. B476 (1996) 559",
 claim="Using Zimmermann normal products and BPHZ, the BV master equation is turned into a local quantum operator equation valid to all orders; anomalies are computed with it at one and two loops for chiral W3 gravity; the anomaly-computation method is usable outside BV.",
 loop_order="all orders (structure); explicit one and two loops",
 theory="chiral W3 gravity (2d)",
 scheme="BPHZ subtraction with Zimmermann normal products (no regulator); quantum action principle of Lowenstein-Zimmermann type.",
 symmetry_argument="Renormalized master equation as a normal-product operator identity; anomaly = local insertion obeying the Wess-Zumino consistency condition, whatever loop order it first appears at.",
 symmetry_status=dict(proved=["Local all-order form of the renormalized BV master equation in BPHZ (abstract).", "Explicit two-loop anomalies in chiral W3 gravity (abstract)."], assumed=[], checked_partially=[]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["Not a supersymmetric or four-dimensional gauge theory example."],
 cross_checks=["S6.01_Barnich1998"],
 source_access="abstract only", confidence="medium",
 notes="Shows that a regulator-free (BPHZ) renormalized QME exists as a local operator statement; relevant as the scheme-independent benchmark against which DR-induced breakings are measured."))

R.append(dict(
 id="S6.03_Hollands2007", section=SEC, authors="S. Hollands",
 title="Renormalized Quantum Yang-Mills Fields in Curved Spacetime",
 year=2007, arxiv="0705.3340", journal="Rev. Math. Phys. 20 (2008) 1033",
 claim="Constructs perturbative quantum Yang-Mills theory on globally hyperbolic spacetimes: the algebra of gauge-invariant observables is the cohomology of the interacting BRST charge of an auxiliary gauge-fixed theory with ghosts and antifields; the key step is a hierarchy of renormalized Ward identities ensuring conservation of the interacting BRST current and nilpotency of the interacting BRST charge; OPE closes on gauge-invariant operators and RG flow stays in that space.",
 loop_order="all orders (formal power series)",
 theory="pure Yang-Mills (no matter) on curved spacetime",
 scheme="Renormalization of local covariant time-ordered products (abstract does not name the scheme; the construction is of Epstein-Glaser type, not dimensional regularization).",
 symmetry_argument="Renormalized (anomalous) Ward identities for the BRST current; anomalies controlled by BRST cohomology.",
 symmetry_status=dict(proved=["Ward identities at the renormalized level guaranteeing conserved BRST current and nilpotent BRST charge (abstract)."], assumed=[], checked_partially=[]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["No matter, no supersymmetry; details of how potential anomalies are removed were not read here."],
 cross_checks=["S6.04_FredenhagenRejzner2011"],
 source_access="abstract only", confidence="medium",
 notes="Rigorous regularization-independent existence result for a renormalized BRST/BV structure in 4d gauge theory; a supersymmetric analogue (constant SUSY ghosts) was not found in this survey."))

R.append(dict(
 id="S6.04_FredenhagenRejzner2011", section=SEC, authors="K. Fredenhagen, K. Rejzner",
 title="Batalin-Vilkovisky formalism in perturbative algebraic quantum field theory",
 year=2011, arxiv="1110.5232", journal="Commun. Math. Phys. (DOI 10.1007/s00220-012-1601-1; volume/page not verified)",
 claim="Constructs the BV complex of perturbatively renormalized QFT without path integral or intermediate regularization: the renormalized time-ordered product is shown equivalent to the classical pointwise product, the renormalized BV algebra is the classical one rewritten with the time-ordered product, and an operator replaces the ill-defined graded Laplacian; this operator is identified with the anomaly term of the anomalous Master Ward Identity of Brennecke and Duetsch.",
 loop_order="all orders (formal)",
 theory="general perturbative gauge theories in the pAQFT framework",
 scheme="Epstein-Glaser renormalization of time-ordered products (no regulator).",
 symmetry_argument="Renormalized QME expressed through the quantum BV operator; failure = anomaly of the anomalous Master Ward Identity.",
 symmetry_status=dict(proved=["Equivalence of renormalized time-ordered product with classical product; definition of the quantum BV operator; identification with the AMWI anomaly (abstract)."], assumed=[], checked_partially=[]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["No supersymmetric example; explicit loop computations not part of this work."],
 cross_checks=["S6.01_Barnich1998 (DR analogue: Delta_d)", "S6.03_Hollands2007"],
 source_access="abstract only", confidence="medium",
 notes="Supplies a scheme-free definition of 'the anomaly of the QME' against which the DR evanescent insertion (S6.01) and the BMHV/DRED breaking (S2.10, S3.08) can be compared."))

R.append(dict(
 id="S6.05_Anselmi2014AB", section=SEC, authors="D. Anselmi",
 title="Adler-Bardeen theorem and manifest anomaly cancellation to all orders in gauge theories",
 year=2014, arxiv="1402.6453", journal="Eur. Phys. J. C74 (2014) 3083",
 claim="In renormalizable, perturbatively unitary gauge theories with matter there is a subtraction scheme in which gauge anomalies cancel manifestly to all orders when they vanish at one loop (the scheme with no finite local subtractions from two loops on). The proof uses the BV formalism with a combined dimensional + higher-derivative regularization and no renormalization-group argument, hence covers conformal and finite theories.",
 loop_order="all orders",
 theory="non-abelian gauge theory with left-handed fermions (worked model); right-handed fermions, scalars, photons in Sec. 8 (not read)",
 scheme="'DHD' regularization: spacetime R^4 x R^(-eps), hatted evanescent components (Sec. 2); right-handed partners psi_R coupled through a matrix varsigma_IJ to regulate chiral propagators in all directions (eqs. 2.5-2.7); gauge-invariant higher-derivative factors (D^2/Lambda^2)^8 on F^2 and (D^2/Lambda^2)^3 on fermion kinetic terms, modified gauge fixing (eqs. 3.1-3.3); limits: eps -> 0 at fixed Lambda, then Lambda -> infinity (Sec. 3.1); gamma5 = i gamma0 gamma1 gamma2 gamma3 with continued Dirac algebra, no named gamma5 scheme.",
 symmetry_argument="BV antibracket, master equation with gauge fermion (eqs. 2.2-2.4); anomaly functional A = (Gamma,Gamma) = <(S,S)>_S (eq. 2.11); regularized violation (S_r0,S_r0) = O(eps) (eq. 2.9) and DHD version (eqs. 3.5, 3.10, 4.2); one-loop anomaly split into non-evanescent, divergent-evanescent and evanescent parts (eqs. 4.3-4.4, 5.1); higher-derivative theory has only one-loop divergences (eqs. 3.11, 4.1).",
 symmetry_status=dict(
  proved=["Counterterms confined to one loop in the HD-regularized theory; the subtracted anomaly functional is convergent (Secs. 3-5)."],
  assumed=["No analytically evanescent coefficients in the regularized action (Sec. 2.2).",
           "No products of two or more fermion bilinears in the relevant functionals, to avoid Fierz ambiguities (Secs. 2.2, 2.6).",
           "Perturbative unitarity and renormalizability."],
  checked_partially=["The higher-loop cancellation step (Sec. 6) and the Lambda -> infinity limit (Sec. 7) were not read here."]),
 evidence=["Sec. 2: eqs. (2.2)-(2.19) BV setup, chiral regularization, evanescent terms S_ev.",
           "Sec. 3: eqs. (3.1)-(3.11) DHD regularization and limit prescription.",
           "Secs. 4-5: eqs. (4.1)-(5.1) one-loop anomaly decomposition."],
 red_flags=["The Fierz restriction excludes exactly the structures on which supersymmetry invariance depends (Yukawa and gaugino Fierz identities); supersymmetry is not treated.",
            "Existence statement: the scheme is characterized implicitly (no finite subtractions beyond one loop), not constructed explicitly for practical computations."],
 cross_checks=["S6.06_Anselmi2014chiralDR", "S3.08_BIMS2020 and later BMHV works (explicit finite restoring counterterms instead of an existence proof)"],
 source_access="full text via ar5iv (Secs. 1-5 only)", confidence="medium",
 notes="Closest existing template for an all-order BV/QME argument combined with DR; a supersymmetric version would have to lift the Fierz restriction, which is where DRED breaks (S2.10)."))

R.append(dict(
 id="S6.06_Anselmi2014chiralDR", section=SEC, authors="D. Anselmi",
 title="Weighted power counting and chiral dimensional regularization",
 year=2014, arxiv="1405.3110", journal="Phys. Rev. D89 (2014) 125024",
 claim="Proposes a dimensional regularization for chiral and parity-violating theories in which spinors, vectors and tensors keep their physical-dimension components (standard gamma matrices), while propagators are regularized by evanescent higher-derivative kinetic terms (Majorana-type for chiral fermions); with weighted power counting the scheme is claimed consistent to all orders.",
 loop_order="all orders (claim)", theory="chiral gauge theories (examples: renormalizability, Adler-Bardeen theorem)",
 scheme="'Chiral dimensional regularization': physical-dimension index algebra, evanescent higher-derivative kinetic terms in the -eps directions.",
 symmetry_argument="Weighted power counting controls the evanescent sector; gauge symmetry handled in BV (per the companion papers).",
 symmetry_status=dict(proved=["All-order consistency of the scheme (claimed in abstract)."], assumed=[], checked_partially=[]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["Explicit computations become harder (stated in abstract); no supersymmetric discussion."],
 cross_checks=["S6.05_Anselmi2014AB", "S2.10_Stoeckinger2005 (DRED also keeps 4d index algebra, but formally)"],
 source_access="abstract only", confidence="medium",
 notes="Structurally intermediate between DRED (4d Clifford algebra kept) and BMHV (d-dim algebra); untested for supersymmetry."))

R.append(dict(
 id="S6.07_Costello2007", section=SEC, authors="K. J. Costello",
 title="Renormalisation and the Batalin-Vilkovisky formalism",
 year=2007, arxiv="0706.1533", journal="arXiv math.QA (later expanded into the book 'Renormalization and Effective Field Theory', AMS 2011; book reference from general knowledge)",
 claim="A method to renormalize certain QFTs on compact manifolds inside the BV formalism; examples Yang-Mills (4d), Chern-Simons and holomorphic Chern-Simons; for Chern-Simons the renormalization respects all symmetries up to homotopy, yielding a manifold invariant.",
 loop_order="all orders (formal)", theory="Yang-Mills (4d), Chern-Simons, holomorphic Chern-Simons",
 scheme="Not stated in the abstract. (General knowledge, not verified here: length-scale effective actions with heat-kernel propagators and a scale-dependent QME; counterterms; obstruction to the QME in a local cohomology group.)",
 symmetry_argument="QME at each scale; obstruction theory in the local deformation complex (general knowledge, not verified in this record).",
 symmetry_status=dict(proved=["Renormalization respecting symmetries up to homotopy for Chern-Simons (abstract)."], assumed=[], checked_partially=[]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["No dimensional regularization at all; comparison with DR-based results requires an explicit scheme map."],
 cross_checks=["S6.08_EWY2017"],
 source_access="abstract only", confidence="medium",
 notes="Included as the reference framework for a regulator other than DR in which the QME and its one-loop obstruction are well defined."))

R.append(dict(
 id="S6.08_EWY2017", section=SEC, authors="C. Elliott, B. Williams, P. Yoo",
 title="Asymptotic Freedom in the BV Formalism",
 year=2017, arxiv="1702.05973", journal="J. Geom. Phys. (DOI 10.1016/j.geomphys.2017.08.009; volume not verified)",
 claim="Defines the beta function in Costello's BV renormalization framework as a cohomology class in the obstruction-deformation complex; for translation-invariant classically scale-invariant theories the one-loop beta function is a well-defined element of the local deformation complex, computable either as a rescaling anomaly or as the logarithmic one-loop counterterm; computes it for first-order Yang-Mills (asymptotic freedom).",
 loop_order="1", theory="first-order Yang-Mills",
 scheme="Costello's framework (regularization not named in the abstract).",
 symmetry_argument="Scale (rescaling) anomaly as a class in the local deformation complex.",
 symmetry_status=dict(proved=["One-loop beta function is a well-defined class, equal to the rescaling anomaly / log counterterm; value for YM (abstract)."], assumed=[], checked_partially=[]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["One loop only; no fermions or supersymmetry in the abstract."],
 cross_checks=["S1d.12_DHKS2007 (DR analogue: dilatation variation of the action proportional to eps times the Lagrangian)"],
 source_access="abstract only", confidence="medium",
 notes="BV analogue of the DR statement that scale breaking is the eps-proportional Lagrangian insertion; useful for the conformal part of the synthesis."))

R.append(dict(
 id="S6.09_BGKWWY2023", section=SEC, authors="K. Budzik, D. Gaiotto, J. Kulp, B. R. Williams, J. Wu, M. Yu",
 title="Semi-Chiral Operators in 4d N=1 Gauge Theories",
 year=2023, arxiv="2306.01039", journal="JHEP 05 (2024) 245",
 claim="Studies quarter-BPS local operators of 4d N=1 SYM through the holomorphic twist, including loop corrections to the space of operators and to the algebraic operations of the twisted theory; classifies planar single-trace quarter-BPS operators for SU(N) and perturbative ones for SU(2), SU(3) at small charges; proposes a holographic dual of the twisted theory.",
 loop_order="perturbative (loop corrections to the twisted differential; order not stated in abstract)",
 theory="4d N=1 SYM (and outline with matter)",
 scheme="Not stated in the abstract.",
 symmetry_argument="Holomorphic twist: the supercharge becomes a BV/BRST-type differential whose loop corrections act on composite operators.",
 symmetry_status=dict(proved=[], assumed=[], checked_partially=["Classification of quarter-BPS operators at small charges (abstract)."]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["Regularization of the loop corrections not checked in this record."],
 cross_checks=["S1a.25_Beisert2003su23 (O(g) dynamic supercharge in the untwisted theory)"],
 source_access="abstract only", confidence="medium",
 notes="Example of loop corrections to a supercharge acting on composite operators computed in a BV-type framework without DR; relevant to formulating one-loop Q-tests for N=4 composite operators."))

for r in R:
    write_record(r, OUT, PROG)
print(len(R), "records written")
