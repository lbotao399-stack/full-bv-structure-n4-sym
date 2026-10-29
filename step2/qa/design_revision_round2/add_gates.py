"""Add explicit acceptance gates to the K, P and S items of step2/theory/engine_plan.json (round 2)."""
import json, os
P = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'theory', 'engine_plan.json')
d = json.load(open(P))
C = {c['id']: c for c in d['computations']}
C['E11']['gates'] = [
    'G11.1 complete composite sums: the subset with S^(1)(x) on the composite line reproduces the full translation (rotation) part; complement 0 exactly (as engine B part2.translation_part_check)',
    'G11.2 locality: the composite remainder is a local functional (homogeneous polynomial in the leg momenta, no log|p|) after the complete sum',
    'G11.3 bilinears: P_1 = -(1/3 + f_P) unit ad_{S1(P_b)} O reproduced by both engines (engine A with the corrected composite term, engine B); letters -f_P',
    'G11.4 R_1 = 0 exactly; M_1 rescaling coefficient m reported (prediction m = 0 on cohomology)',
    'G11.5 the scheme used downstream is recorded in the output JSON (translation covariant, or P_1 = 0 with the family finite parts and their p_F)']
C['E8']['gates'] = [
    'G8.1 S_1(eta) - Q_1(rho(x)eta) + p ad_{S1(eta)} is translation invariant and linear in 4 eta (U_eta), zero on the scalar bilinears',
    'G8.2 the same -p for both chiralities (S^alpha_i and Sbar^{alphadot i})',
    'G8.3 two engines agree exactly']
C['E10']['gates'] = [
    'G10.1 eps -> 0 at fixed T: the limit is identical (exact arithmetic, or 40 digits) for at least three values of T; the O(eps) approach eps(a + b log eps) has no eps^{-1} and no log eps term after the counterterm vertices',
    'G10.2 locality of the limit: homogeneous polynomial in the leg momenta (p -> 2p scaling exact), no log|p|',
    'G10.3 complete composite sums (R_{O,K}); single composite graphs are never reduced by Q-exactness separately',
    'G10.4 read-off by the corrected eq:o-L1graphs (letters: trees of hbar L_1^(0) with S^(1)(K_b) subtracted)',
    'G10.5 comparison with design thm:o-KD: K_1 + p ad_{S1(K_b)} == (1/4) div v gamma O with the p of E11 and the gamma of Part D, in both schemes',
    'G10.6 relation at (eta1,eta2) with S_1 of E8 (E10d) and K = P^dagger (E10e) as cross-checks; disagreement is reported, not fitted',
    'G10.7 two engines agree; every job deterministic, one process, peak RSS recorded (< 1 GB)']
json.dump(d, open(P, 'w'), indent=1, ensure_ascii=False)
for c in d['computations']:
    print(c['id'], '|', c.get('resources'))
