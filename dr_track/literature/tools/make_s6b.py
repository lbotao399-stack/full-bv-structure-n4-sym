import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rec_io import write_record

OUT = "/Users/libotao/N4_BV_DR/literature/parts/S6"
PROG = "/Users/libotao/N4_BV_DR/literature/parts/S6_PROGRESS.txt"
SECb = "6b conformal symmetry in d != 4 (breaking and restoration in DR)"
SECc = "6c all-order protection results with manifestly supersymmetric or algebraic schemes; BPS protection checks in DR"

R = []

R.append(dict(
 id="S6.10_BraunManashov2013", section=SECb, authors="V. M. Braun, A. N. Manashov",
 title="Evolution equations beyond one loop from conformal symmetry",
 year=2013, arxiv="1306.5644", journal="Eur. Phys. J. C73 (2013) 2544",
 claim="At the critical coupling of a theory in non-integer dimension (O(n) phi^4 in d=4-2eps, su(n) phi^3 in d=6-2eps) the theory is taken to be exactly conformal; the collinear conformal generators acquire coupling-dependent corrections, and commutation relations with the evolution kernel determine the non-conformal part of the l-loop kernel from the (l-1)-loop correction to the special conformal generator; MS kernels are eps-independent, so results carry over to integer dimension.",
 loop_order="method to all orders; explicit checks at two/three loops in scalar theories",
 theory="O(n) phi^4 (d=4-2eps), su(n) phi^3 (d=6-2eps); gauge theories announced for later work",
 scheme="Dimensional regularization at the Wilson-Fisher critical point u*, MS(-bar); kernels eps-independent in minimal subtraction (Sec. 2.2).",
 symmetry_argument="Exact conformal invariance of the critical theory in d dimensions; special conformal Ward identity derived from the finite renormalized generating functional.",
 symmetry_status=dict(
  proved=["Generator corrections: S_- uncorrected (eqs. 44-46); Delta S_0 = -eps + (1/2) H(u*) (eqs. 39, 51); Delta S_+ = (z1+z2)(-eps + (1/2) u* H^(1)) + O(eps^2) at one loop (eqs. 39, 52-62, 80).",
          "Consistency relations [S_+^(0), H(u*)] = [H(u*), Delta S_+] and [Delta S_+, S_-] = 2 Delta S_0 (eqs. 88-91).",
          "A naive all-order guess for S_+ (eq. 81) fails; explicit extra term (1/4) u*^2 [H^(2), z1+z2] (eq. 83)."],
  assumed=["The critical theory in non-integer d is exactly conformally invariant (Secs. 1, 2.3).",
           "eps-singular terms of the Ward identity cancel because the renormalized generating functional is finite (Sec. 3).",
           "Transfer to integer d via eps-independence of MS kernels (Secs. 1, 2.2)."],
  checked_partially=["Only scalar theories in this paper."]),
 evidence=["Sec. 2.1 eqs. (2)-(4), (12)-(13): models and critical points.", "Sec. 2.3 eq. (39); Sec. 3 eqs. (44)-(62), (80)-(83); Sec. 4 eqs. (88)-(91); eq. (103) eps-u* relation."],
 red_flags=["Non-integer-d theories are non-unitary, so 'scale implies conformal' is not available; conformality is assumed (motivates S6.11).",
            "Does not apply directly to N=4 SYM: with beta_4 = 0 the d-dimensional beta function -2eps a has only the trivial zero, so no interacting critical point exists in d = 4-2eps."],
 cross_checks=["S6.11_BMMS2018", "S6.13_BelitskyMueller1998", "S6.14_BKM2005"],
 source_access="full text via ar5iv (first 100k characters; Secs. 1-4)", confidence="high",
 notes="Provides the precise form of 'conformal symmetry broken by eps and restored at the critical point': generator deformations of order eps and of order the kernel."))

R.append(dict(
 id="S6.11_BMMS2018", section=SECb, authors="V. M. Braun, A. N. Manashov, S. Moch, M. Strohmaier",
 title="Conformal symmetry of QCD in d dimensions",
 year=2018, arxiv="1810.04993", journal="journal not verified (DESY 18-175)",
 claim="Gives a perturbative argument that QCD in d = 4-2eps at the critical coupling a* is conformal in the gauge-invariant sector: the inhomogeneous term in the conformal variation of a renormalized gauge-invariant operator is itself gauge-invariant up to BRST-exact and EOM terms that drop out of correlators.",
 loop_order="all orders (perturbative argument)", theory="QCD in d = 4-2eps at the critical point",
 scheme="MS in d = 4-2eps (eq. 5); multiplicative renormalizability assumed (Sec. 2); Landau gauge used for scale invariance of all Green functions (Sec. 5).",
 symmetry_argument="Scale and conformal Ward identities: the bulk term carries beta(a) and vanishes at a* (eq. 24); products of the Lagrangian with renormalized operators decompose into gauge-invariant, BRST and EOM pieces (eqs. 26-28, 35; Joglekar-Lee); renormalized transformation law delta_K O_j = (delta_jj' K_Dj + 2 x^mu gamma_jj') O_j' + O^mu_j (eq. 44).",
 symmetry_status=dict(
  proved=["Vanishing of bulk contributions at beta(a*)=0 (eq. 24); dilatation identity F_j = - sum gamma_jj' O_j' (eq. 42); gauge invariance of the inhomogeneous term up to BRST/EOM terms (Secs. 4-5)."],
  assumed=["Multiplicative renormalizability; existence of a*; diagonalizability of the anomalous-dimension matrix (fn. 2 warns it can fail with fermions).",
           "Joglekar-Lee structure of gauge-invariant operator mixing."],
  checked_partially=[]),
 evidence=["Sec. 1 eqs. (1)-(3); Sec. 2 eqs. (4), (9)-(14), (21)-(23); Sec. 3 eqs. (16)-(25); Sec. 4 eqs. (26)-(35); Sec. 5 eqs. (38)-(44)."],
 red_flags=["Gauge fixing and ghosts break conformal symmetry already in d=4; only gauge-invariant correlators transform covariantly.",
            "The authors flag that F Ftilde = dK shows inhomogeneous terms that are not gauge invariant can arise.",
            "Evanescent operators, the d-dimensional energy-momentum tensor and gamma5 are not discussed."],
 cross_checks=["S6.10_BraunManashov2013", "S6.12_BMMS2021 (gamma5 in Larin's scheme)"],
 source_access="full text via ar5iv", confidence="high",
 notes="The BRST-exact/EOM bookkeeping here is exactly what a BV formulation organizes (antifield-dependent terms); a BV proof would replace the Joglekar-Lee input."))

R.append(dict(
 id="S6.12_BMMS2021", section=SECb, authors="V. M. Braun, A. N. Manashov, S. Moch, M. Strohmaier",
 title="Three-loop off-forward evolution kernel for axial-vector operators in Larin's scheme",
 year=2021, arxiv="2101.01471", journal="Phys. Rev. D103 (2021) 094018",
 claim="Extends the critical-point conformal technique to axial-vector leading-twist operators with gamma5 in Larin's scheme; three-loop off-forward kernels are reconstructed from anomalous dimensions and the two-loop special conformal anomaly; an explicit finite renormalization kernel, encoding the axial-minus-vector difference, is needed to return to the conventional MS-bar scheme.",
 loop_order="3", theory="QCD, flavour-nonsinglet axial-vector twist-two operators",
 scheme="DR at the Wilson-Fisher critical point in d = 4-2eps; Larin's gamma5 prescription plus finite renormalization.",
 symmetry_argument="Conformal symmetry at the critical point; gamma5 prescription-induced violations of the vector/axial relation absorbed in a finite renormalization kernel.",
 symmetry_status=dict(proved=["Three-loop axial kernel and the finite renormalization kernel (abstract)."], assumed=["Conformal symmetry at the critical point (as in S6.10/S6.11)."], checked_partially=[]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["Shows that a non-symmetric gamma5 prescription interferes with symmetry relations used to bootstrap higher loops; the needed finite kernel is scheme data, not physics."],
 cross_checks=["S5.04_Larin1993 (finite Z5)", "S6.11_BMMS2018"],
 source_access="abstract only", confidence="medium",
 notes="Template for a 'symmetry-restoring finite counterterm' in a conformal context."))

R.append(dict(
 id="S6.13_BelitskyMueller1998", section=SECb, authors="A. V. Belitsky, D. Mueller",
 title="Broken conformal invariance and spectrum of anomalous dimensions in QCD",
 year=1998, arxiv="hep-ph/9804379", journal="Nucl. Phys. B537 (1999) 397",
 claim="Uses the conformal algebra and conformal Ward identities to constrain the dilatation and special conformal anomalies of twist-two operators; the anomalies are computed at leading order in MS; the conformal consistency relation then gives the off-diagonal part of the two-loop anomalous-dimension matrix, and the NLO off-diagonal ERBL kernels in the flavour-singlet parity-even, parity-odd and chiral-odd sectors.",
 loop_order="special conformal anomaly at one loop; off-diagonal anomalous dimensions at two loops",
 theory="QCD twist-two operators",
 scheme="Minimal subtraction (MS) in dimensional regularization (per abstract; details not verified).",
 symmetry_argument="Broken conformal Ward identities: the special conformal anomaly at order l-1 fixes the non-diagonal (conformal-symmetry-breaking) part of the anomalous dimensions at order l.",
 symmetry_status=dict(proved=["LO special conformal anomaly; two-loop off-diagonal anomalous dimensions; NLO ERBL off-diagonal kernels (abstract)."], assumed=[], checked_partially=[]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["Parity-odd sector involves gamma5; its treatment was not checked in this record."],
 cross_checks=["S6.14_BKM2005", "S6.10_BraunManashov2013"],
 source_access="abstract only", confidence="medium",
 notes="Original form of the conformal-anomaly method in DR."))

R.append(dict(
 id="S6.14_BKM2005", section=SECb, authors="A. V. Belitsky, G. P. Korchemsky, D. Mueller",
 title="Integrability of two-loop dilatation operator in gauge theories",
 year=2005, arxiv="hep-th/0509121", journal="Nucl. Phys. B735 (2006) 17",
 claim="Computes one- and two-loop dilatation operators in the non-compact SL(2) sector (three quark or gaugino fields of equal helicity with light-cone derivatives) of QCD and of N=1,2,4 SYM; two-loop spectra are integrable in planar SYM, and with extended supersymmetry the N-dependence factorizes into a c-number; conformal symmetry is violated from two loops on.",
 loop_order="2", theory="QCD (SU(3)) and N=1,2,4 SYM (SU(Nc)), SL(2) sector",
 scheme="Computed in parallel in DREG (MS-bar) and DRED (DR-bar) (Sec. 4.1), light-like axial gauge with principal-value prescription (eq. 4.3); conversion eqs. (5.48)-(5.53).",
 symmetry_argument="In d = 4-2eps the coupling runs with beta_eps(g^2) = -2eps + beta(g^2) (eq. 4.6), so conformal symmetry is broken for eps != 0 even for N=4; the breaking enters the dilatation operator from two loops; at L=3, V_DR^(1) = V^(1) - (1/6) V^(0) - (3/4)(4-N) 1l (eq. 5.54), the c-number term vanishing for N=4.",
 symmetry_status=dict(
  proved=["Two-loop dilatation operators in both schemes; DR-bar and MS-bar operators differ by a c-number/normalization at two loops (Sec. 5.4.1, eqs. 5.48-5.54)."],
  assumed=["DRED preserves supersymmetry at two loops (stated as the reason DRED is mandatory, not verified; Sec. 4.1).",
           "Principal-value prescription for light-like gauge poles (eq. 4.3)."],
  checked_partially=["Sec. 6.1.3 (conformal Ward identity discussion for N=4) only partly read."]),
 evidence=["Sec. 2 eqs. (2.3)-(2.9): field content, beta = 0 for N=4.", "Sec. 4.1, eq. (4.6).", "Sec. 5.4.1, eqs. (5.48)-(5.54)."],
 red_flags=["SUSY preservation of DRED is an assumption here.", "Light-cone gauge with PV prescription is a non-covariant ingredient."],
 cross_checks=["S6.13_BelitskyMueller1998", "S1a.25_Beisert2003su23 (closure assumptions in the compact sectors)"],
 source_access="full text via ar5iv (up to part of Sec. 6.1.3)", confidence="medium",
 notes="For N=4 the only conformal-breaking parameter in d dimensions is the explicit -2eps (no critical point), exactly the term a d-dimensional BV conformal Ward identity has to carry."))

R.append(dict(
 id="S6.15_BHMP2010", section=SECb, authors="N. Beisert, J. Henn, T. McLoughlin, J. Plefka",
 title="One-Loop Superconformal and Yangian Symmetries of Scattering Amplitudes in N=4 Super Yang-Mills",
 year=2010, arxiv="1002.1733", journal="JHEP 1004 (2010) 085",
 claim="In dimensional regularization, one-loop amplitudes violate tree-level superconformal and Yangian symmetry (IR divergences, holomorphic-anomaly collinear terms); the anomalies are computed by acting with tree generators on unitarity cuts and are cancelled by a universal one-loop deformation of the generators built from tree data; checked for MHV and the 6-point NMHV amplitude.",
 loop_order="1", theory="planar N=4 SYM amplitudes",
 scheme="Dimensional regularization of one-loop amplitudes (variant not verified in abstract).",
 symmetry_argument="Deformed generators J = J^(0) + g^2 J^(1) restoring the symmetry on regulated amplitudes.",
 symmetry_status=dict(proved=["Universal one-loop generator deformation (abstract)."], assumed=[], checked_partially=["Verified for MHV amplitudes and 6-point NMHV (abstract)."]),
 evidence=["Abstract only; section/eq. no. not verified."],
 red_flags=["On-shell, IR-regulated quantities; not a statement about the UV regularization of composite operators."],
 cross_checks=["S1d.12_DHKS2007 (dual conformal anomaly in DR)", "S1a.25_Beisert2003su23 (deformed generators on composite operators)"],
 source_access="abstract only", confidence="medium",
 notes="On-shell analogue of 'symmetry broken by the regulator, restored by O(g^2) generator corrections'."))

R.append(dict(
 id="S6.16_BPS2015", section=SECc, authors="I. L. Buchbinder, N. G. Pletnev, K. V. Stepanyantz",
 title="Manifestly N=2 supersymmetric regularization for N=2 supersymmetric field theories",
 year=2015, arxiv="1509.08055", journal="Phys. Lett. B (2015) (DOI 10.1016/j.physletb.2015.10.071; volume/page not verified)",
 claim="Constructs a higher covariant derivative regularization in N=2 harmonic superspace plus N=2 Pauli-Villars superfields, preserving N=2 supersymmetry and background gauge invariance at every step; proves the N=2 non-renormalization theorem: divergences only at one loop, hypermultiplets not renormalized, beta function one-loop exact; N=4 SYM finite to all orders.",
 loop_order="all orders", theory="N=2 gauge theories with hypermultiplets; N=4 SYM as special case",
 scheme="N=2 higher-derivative term S_Lambda (eq. 13), background field split and gauge fixing (eqs. 14-25), Nielsen-Kallosh sector (eqs. 29-31), commuting analytic PV superfields with masses proportional to Lambda (eqs. 34-35); generating functional eq. (36).",
 symmetry_argument="Manifest N=2 supersymmetry and background gauge invariance of the regularized theory; power counting omega_Lambda = omega - 2(P-V) (eq. 33); analyticity gives gamma_phi = 0 (eq. 40); beta one-loop exact (eqs. 45-46); second proof via NSVZ (eqs. 37-38).",
 symmetry_status=dict(
  proved=["Finiteness beyond one loop and one-loop exactness of beta in the HD+PV scheme (Secs. 3-5)."],
  assumed=["N_D > 0 beyond one loop and manifest N=2 SUSY of the regularized effective action taken from ref. [5] (fn. 6).",
           "NSVZ relation in the second proof supported by evidence, not derived (Sec. 4).",
           "One-loop computation done at m0 = 0, xi = 1 (Sec. 5)."],
  checked_partially=[]),
 evidence=["Sec. 3 eqs. (13)-(36); Sec. 3.3 eqs. (32)-(33); Sec. 4 eqs. (37)-(38), (40); Sec. 5 eqs. (45)-(46); Introduction: earlier proofs (Howe-Stelle-Townsend, Howe-Stelle-West) implicitly assumed an N=2 supersymmetric gauge-invariant regularization."],
 red_flags=["Anomalies are not discussed.",
            "The Introduction cites DRED as breaking SUSY at three loops (Avdeev 1982); that claim is disputed by 2603.02892 (S2.13)."],
 cross_checks=["S1e.11_Stepanyantz2020 (N=1 HD all-loop NSVZ)", "S2.13_ChakrabortyMoch2026"],
 source_access="full text via ar5iv", confidence="high",
 notes="Establishes all-order N=4 finiteness with a regulator that is supersymmetric by construction; this is the cleanest 'genuinely established' protection result found, and it is not a DR-type scheme."))

R.append(dict(
 id="S6.17_LPS1988", section=SECc, authors="C. Lucchesi, O. Piguet, K. Sibold",
 title="Vanishing beta functions in N=1 supersymmetric gauge theories",
 year=1988, arxiv="none (pre-arXiv)", journal="Helv. Phys. Acta 61 (1988) 321",
 claim="For N=1 SYM with simple gauge group: if the gauge anomaly vanishes, the one-loop gauge beta function and one-loop matter anomalous dimensions vanish, and the one-loop Yukawa solution lambda = rho g is isolated and non-degenerate, then the gauge and Yukawa beta functions vanish to all orders (one-coupling theory via reduction of couplings).",
 loop_order="all orders", theory="N=1 SYM with chiral matter, simple gauge group",
 scheme="Algebraic renormalization with renormalized Ward identities for R, SUSY, BRS (Slavnov) and chiral symmetries (eqs. 3.5-3.10 of hep-ph/9604216), assuming a supersymmetric superspace renormalization scheme; not DR.",
 symmetry_argument="Supercurrent anomaly expansion relates the R-axial anomaly coefficient r to beta functions and anomalous dimensions (eq. 3.31 of hep-ph/9604216); non-renormalization of the axial anomaly coefficients (eq. 3.32) then forces beta_g = 0 to all orders (eq. 3.46).",
 symmetry_status=dict(
  proved=["All-order vanishing of beta functions under one-loop hypotheses, within the stated scheme (as summarized in hep-ph/9604216 Sec. 3)."],
  assumed=["A renormalization scheme preserving supersymmetry (superspace scheme of the cited ref. [13]).",
           "Non-renormalization of the R-axial/chiral anomaly coefficients."],
  checked_partially=[]),
 evidence=["Secondary: Lucchesi-Zoupanos hep-ph/9604216, Sec. 3, criterion after eq. (3.32); eqs. (3.5)-(3.10), (3.14)-(3.17), (3.31)-(3.37), (3.43)-(3.46)."],
 red_flags=["Statement is about vanishing beta functions in a specific scheme; scheme requirements stated explicitly for phases (hep-ph/9604216 Sec. 3).",
            "Relies on a SUSY-preserving scheme that is not DRED."],
 cross_checks=["S4.03_MaggiorePiguetWolf1996 (constant-ghost algebraic renormalization)", "S6.16_BPS2015"],
 source_access="secondary: cited in hep-ph/9604216 (full text via ar5iv)", confidence="medium",
 notes="Example of an all-order protection theorem whose regularization input is an assumed supersymmetric scheme."))

R.append(dict(
 id="S6.18_PSZ1999", section=SECc, authors="S. Penati, A. Santambrogio, D. Zanon",
 title="Two-point functions of chiral operators in N=4 SYM at order g^4",
 year=1999, arxiv="hep-th/9910197", journal="JHEP 9912 (1999) 006",
 claim="The two-point function of the chiral primary Tr Phi^3 in N=4 SU(N) SYM receives no correction through order g^4 for any N; at order g^2 this follows from stress-tensor-multiplet non-renormalization, at g^4 the cancellation is between zeta(3) and zeta(5) pole terms of different diagrams.",
 loop_order="2 (order g^4)", theory="N=4 SU(N) SYM in N=1 superspace",
 scheme="N=1 superspace; dimensional regularization n = 4-2eps for momentum integrals (App. B); minimal subtraction keeping only 1/eps poles; whether the D-algebra is done in 4d is not stated explicitly (Sec. 2 'standard techniques'); uniqueness method for integrals (eqs. B.4-B.10).",
 symmetry_argument="Conformal invariance fixes the x-dependence (eq. 2.4); vanishing shown for the 1/eps poles only; finite momentum-space parts claimed to be contact terms (Secs. 2, 5).",
 symmetry_status=dict(
  proved=["Cancellation of all 1/eps poles at O(g^4): 6a -(3/2)zeta3, 6b +3zeta3, 6i -5zeta5, 6j +5zeta5-(3/2)zeta3 (eqs. 4.6, 4.8, 4.9, 4.21); non-planar colour factors vanish (App. A)."],
  assumed=["The regularization preserves the relevant (super)symmetry; not discussed.", "Finite parts irrelevant (contact terms)."],
  checked_partially=["Only one operator (Tr Phi^3) and only pole terms."]),
 evidence=["Sec. 2 eq. (2.4); Sec. 4 eqs. (4.3)-(4.21); App. A eqs. (A.18)-(A.21); App. B eqs. (B.4)-(B.10)."],
 red_flags=["Protection is checked on 1/eps poles only; finite scheme-dependent terms are not tested.",
            "N=4 SUSY only non-linearly realized in N=1 superspace."],
 cross_checks=["S1a.07_Sieg2010 (Sec. 4 links BPS protection and conformal invariance)"],
 source_access="full text via ar5iv", confidence="high",
 notes="Representative two-loop DR check of BPS protection for a composite operator."))

for r in R:
    write_record(r, OUT, PROG)
print(len(R), "records written")
