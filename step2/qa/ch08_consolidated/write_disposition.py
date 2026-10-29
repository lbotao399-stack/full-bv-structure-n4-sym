#!/usr/bin/env python3
import json, sys
from pathlib import Path
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from edits import E
ROOT = HERE.parents[2]
ED = {e['id']: e for e in E}
D = 'step2/qa/ch08_corrections_from_delta.json'
G = 'step2/qa/ch08_corrections_from_design.json'
P = 'step2/qa/ch08_corrections_from_ip.json'
C = 'step2/qa/ip_completeness2.json'
OUT7 = 'outside Chapter 8: the item edits ch07 rem:c7-internal; only manuscript/ch08_definitions.tex is edited. The Chapter 8 statement that relied on it (def:c8-radial, P_1=0 attributed to prop:c7-cov) is corrected in R21/R22, which describe P_1 through thm:c9-P1 and prop:c9-scheme.'
X = {}
def put(src, i, status, reason, edits=(), merged_with=None):
    X[(src, i)] = dict(status=status, merged_with=merged_with, reason=reason, edits=list(edits))

# ---- delta stream
put(D,0,'rejected','Its order |y|^{-4} along Delta (convergence only by oddness) is wrong: the Euler and rotation operators preserve homogeneity, (2(y,d_y)+2)|y-a|^{-2} = -2|y|^{-2}+O(|y|^{-3}), so the Delta insertion is O(|y|^{-5}) and absolutely convergent (delta#10, #16, #21, #27; audit_delta_2b/insertion_decay.json). Its correct parts (A x A edge of order |y|^{-2}; composite kernel of order |x-y|^{-1}) enter through delta#21 (R15, R16).',['R15','R16'])
put(D,1,'rejected','The T = infinity case is the first identity of eq:c8-transfer at (T_0, infinity) with K_infinity = 0, Delta_infinity = 0 (Definitions def:c3-K, def:c3-delta), an exact identity of kernels; it needs only the convergence of the flow with P(T_0,infinity), and no limit of the K_T-bracket terms enters (delta#27). The proof of thm:c8-L (ii) is made explicit on this point (R09).',['R09'])
put(D,2,'rejected','Its claim that the extra terms leave the spectrum unchanged is false: they shift the eigenvalue at Weyl weight 2 by -4(1/3+theta) units (delta#11); the terms z, G_Delta, R_O are superseded by ch09 thm:c9-D1 (exact formula) and thm:c9-kappa (theta fixed). Replaced by delta#12 (R10).')
put(D,3,'merged-with','Core point (P_1 = 0 is a scheme condition; prop:c7-cov gives only regulator covariance) merged into the def:c8-radial text of design#0 (R21, R22). Its finite part "c times the Weyl weight" is rejected: it differs from (c/2){S^(1)(Delta),-} = c(E+w) by the position scaling, is not a map of local functionals and changes the Delta and K components by non-multiples (ch09 prop:c9-scheme; design_revision_round2/scheme_shift_round2.json); R22 uses (a_P/2){S^(1)(Delta),-}. Its requested covariance statement in cor:c8-gram is already present ("Under a finite part a -> (1+hbar psi\')a, with psi\' reflection compatible, ... the relation is preserved").',['R21','R22'],G+'#0')
put(D,4,'rejected','Vacuous: by prop:c3-H the gauge-fixing operators with the arrows of def:c3-qgf and H the componentwise Laplacian are unique, so a smooth family with those properties is constant, and under its hypothesis (difference by (a), (b)) the conclusion is item (iii); the comparison gauge with (Q^GF)^2 = 0 is not available in the minimal complex (delta#26 evidence; delta#27).')
put(D,5,'rejected',OUT7)
put(D,6,'merged-with','Same formula D_1(O) = 2cO - (1/3+f) unit {S^(1)(Delta),O}; R10 uses delta#12, in which theta is fixed by the (Delta,P_b) component (thm:c9-kappa) rather than "by the whole equation since it changes the spectrum".',['R10'],D+'#12')
put(D,7,'merged-with','Same correction of the parenthetical; R27 uses delta#13 (cites the ch09 identification).',['R27'],D+'#13')
put(D,8,'merged-with','def:c8-radial part merged into design#0 (R21, R22); thm:c8-adjoint parenthetical superseded by design#1 (hypothesis P_1=M_1=R_1=0 removed, R30); cor:c8-gram wording "In the scheme ... with P_1=0" merged into R32. Its finite part "-1/3 units times the Weyl weight" is replaced by (a_P/2){S^(1)(Delta),-} (prop:c9-scheme; see delta#3).',['R21','R22','R30','R32'],G+'#0')
put(D,9,'rejected',OUT7)
put(D,10,'merged-with','Same correction of the Tsafe proof (A x A edge of order |y|^{-2}, derivative of the insertion on that edge); R15 uses delta#21. Its odd insertion lacked the minus sign of -(Gamma(epsilon,lambda),A^*) (eq:c2-S1Q); R15 has it.',['R15'],D+'#21')
put(D,11,'applied','No text of its own. Its recommendation not to adopt the thD replacement of the Remark is followed (delta#2 rejected); its neutral items (thm:c8-L (ii), thm:c8-invariants (iv)) are decided by delta#27.')
put(D,12,'applied','Applied as the text of Remark (What is computed), consistent with ch09 thm:c9-D1, thm:c9-kappa, cor:c9-ident. Modified: the hypothesis of thm:c9-kappa is stated in full ("the one-loop quantisation and the local representatives exist through arity two"), theta called a scheme parameter, cor:c9-ident cited for eq:c6-step2.',['R10'])
put(D,13,'applied','Applied; hypothesis completed as in delta#12 and cor:c9-ident cited.',['R27'])
put(D,14,'merged-with','Its def:c8-radial text (thm:c9-P1 value -(1/3+theta_P), finite part (a/2){S^(1)(Delta),-} with a = -(1/3+theta_P) units, changes a, -a/2, a/2, reflection compatible) is the content of R22 (a written a_P). Its thm:c8-adjoint parenthetical is superseded by design#1 (hypothesis removed). Its cor:c8-gram wording is merged into R32.',['R21','R22','R32'],G+'#0')
put(D,15,'rejected',OUT7)
put(D,16,'merged-with','Same correction as delta#21, which is the fuller version (adds the K_b origin and the kernel orders); R15 uses delta#21 with the odd term -(Gamma(epsilon(y),lambda),A^*) as here.',['R15'],D+'#21')
put(D,17,'applied','Applied (R01); the citation in the proof of lem:c8-transfer is made consistent with it (R02).',['R01','R02'])
put(D,18,'merged-with','Same convention as delta#17 (R01 uses the delta#17 wording, which cites Def c3-P and Lemma c9-ops).',['R01','R02'],D+'#17')
put(D,19,'merged-with','Kernel normalisation adopted where the translation composite kernel is named: R22 writes 2(d_b (P(eps,infty)*)^2 (x) 1)K_0 (ch09 Translations).',['R22'],G+'#0')
put(D,20,'rejected','It corrects R_O inside the thD replacement text of the Remark, which is not adopted (delta#2); R_O does not appear in Chapter 8.')
put(D,21,'applied','Applied (R15 for the edge decay, R16 for the K_b origin and the kernel orders). Modified: the index expression v^m d_m written index-free as "the derivative along v" (ESW v(-)); v.A expanded by Definition def:c7-vdot; P(eps,T)* defined at first use; Prop c9-composite (iv) cited for the small-momentum orders.',['R15','R16'])
put(D,22,'rejected','Same as delta#1: the T = infinity transfer uses eq:c8-transfer at (T, infinity) with Delta_infinity = 0 and needs only convergence of the flow (delta#27); the same holds for prop:c8-Tindep (ii), whose proof already uses that identity.',['R09'])
put(D,23,'merged-with','Its point (the identification with Chapter 6 is conditional) is retained: R10 and R27 state the exact formula of thm:c9-D1 and the conditional fixing of theta by thm:c9-kappa; its "not established in Feynman gauge" is superseded by those ch09 theorems.',['R10','R27'],D+'#12')
put(D,24,'merged-with','Merged into design#0 (R21, R22). Its finite part "c times the Weyl weight" is rejected (see delta#3); its reflection-compatibility requirement is in R22 (real a_P) and already in cor:c8-gram.',['R21','R22'],G+'#0')
put(D,25,'rejected',OUT7)
put(D,26,'rejected','Same as delta#4 (vacuous by prop:c3-H; the (Q^GF)^2 = 0 comparison gauge is outside the minimal complex, as the item itself records).')
put(D,27,'applied','Its three withdrawals are followed (delta#0, #1, #4, #22, #26 rejected on its grounds); its argument for thm:c8-L (ii) is written into the proof (R09).',['R09'])
put(D,28,'merged-with','Same content as delta#12 (the multiple vanishes under the (Delta,P_b) component); R10 uses delta#12.',['R10'],D+'#12')
put(D,29,'merged-with','Same parenthetical correction; R27 uses delta#13.',['R27'],D+'#13')
put(D,30,'rejected','Superseded by ch09 thm:c9-kappa, which states and proves the (Delta,P_b) consequence with its hypotheses (arity-two quantisation and representatives); R10 and R27 cite it. The proposed item is also incomplete: thm:c9-kappa needs, besides densities of weight 3, the position-linear densities of weight 4 in the arity-two output.')
put(D,31,'merged-with','Same finite part (c/2){S^(1)(Delta),-}, real, as in R22 (c written a_P); the "-1/3 units on the scalar bilinears" is stated with theta_P as in thm:c9-P1.',['R21','R22'],G+'#0')
put(D,32,'rejected',OUT7)
put(D,33,'merged-with','Followed for thm:c8-L (ii) and thm:c8-invariants (iii) (statements unchanged). For prop:c8-Tsafe the exponents -6/-5/-4 are kept, but the proof wording is corrected because the A x A component carries no derivative (delta#21).',['R09'],D+'#27')
# ---- design stream
put(G,0,'applied','Applied (R21 removes P_1 = 0 from the assumption; R22 adds the dependence on P_1 and the renormalisation-condition statement). Modified: no new symbols (the family F and projection pi_F replaced by "a family stable under Lambda" and "the finite part of the operators of the family", the hypothesis of prop:c9-scheme); p written a_P (name of ch09 thm:c9-kappa); value cited from thm:c9-P1 with theta_P; the letters value omitted (not stated in the manuscript); kernel written 2(d_b(P*)^2 (x) 1)K_0 (delta#19); rem:c7-internal not cited (it states the superseded "no one-loop correction"); the field-redefinition sentence condensed; {S^(1)(Delta),-} = 2(E+w) cited from lem:c9-SDelta (i).',['R21','R22'])
put(G,1,'applied','Applied: the proof uses only reflection compatibility, the Ward identity and cor:c8-homotopy for the full Lambda^(1). Combined in one replacement with the Ward-identity wording of ip#3.',['R30'])
put(G,2,'applied','Applied (wording "the first in a scheme with P_1=0"); ip#5 merged.',['R23'])
put(G,3,'merged-with','Same formula as ipc2#3 for real a_P (pbar = p); R32 uses ipc2#3, which also names the Gram matrix of that scheme, and keeps this item\'s reason (ad_{S(P_b)})^adj = -ad_{S(K_bbar)}.',['R32'],C+'#3')
put(G,4,'applied','Applied (R05); the parenthetical "this step is needed: ..." stated as a mathematical clause.',['R05'])
put(G,5,'applied','Applied verbatim (R04); same position-degree convention as design#4.',['R04'])
put(G,6,'applied','Partly applied: the product of functional-valued cochains with total-parity Koszul signs, e^{J/hbar}, associativity and graded commutativity (R07). Rejected part: "|O| below are total parities". With the total parity, on an odd O of degree 0 and even x, y one gets [Lambda^(1)(x),Lambda^(1)(y)] = +Lambda^(1)([x,y]), against cor:c8-homotopy (ii); |O| is the ghost-number degree (ip#2, ipc2#6, ip_completeness2_checks.json K4; user decision), R08.',['R07','R08'])
put(G,7,'applied','Applied verbatim.',['R12'])
put(G,8,'applied','Applied (three replacements).',['R03','R06','R11'])
put(G,9,'applied','Applied verbatim; consistent with the proof of thm:c8-twist (N S^xi = -S^(1)(xi)+S^(2)(xi,xi)).',['R13'])
put(G,10,'applied','Statement text applied verbatim (R14); its fibre-factor derivation of 3/(8 pi^2) placed in the proof (R16), together with the orientation cancellation of even-even chains.',['R14','R16'])
put(G,11,'applied','Applied; last clause stated affirmatively ("so the defect graphs are evaluated at finite T").',['R17'])
put(G,12,'applied','Applied verbatim.',['R18'])
put(G,13,'applied','Applied with dedupe and one correction: the per-graph log T coefficient refers to the cut tree of prop:c8-Tsafe (stated there by design#10); "the limit eps->0 is therefore local and the same at every T" replaced by "the same at every T, and the existence of the local representative at fixed T is the locality of that limit" (T d/dT = O(eps) gives T-independence; locality is the separate short-distance statement, as in the original remark); the composite-line sentence shortened to refer to prop:c8-Tsafe; v named as the vector field of K_b.',['R19'])
# ---- ip stream
put(P,0,'applied','Applied verbatim.',['R24'])
put(P,1,'applied','Applied verbatim, in one replacement with ipc2#0.',['R25'])
put(P,2,'merged-with','Same content as ipc2#6; R08 uses the ipc2#6 wording with "in the sense of Definition def:c2-dce" from this item.',['R08'],C+'#6')
put(P,3,'merged-with','Ward-identity part (two-point function of Definition def:c8-2pt; |F|, |Y| total parities) applied verbatim in R30. The part keeping "P_1=M_1=R_1=0 (frame ...)" is rejected: these hypotheses are unused by the proof and removed (design#1); P_1 = 0 is described as a renormalisation condition in def:c8-radial.',['R30'],G+'#1')
put(P,4,'merged-with','Merged into design#0. Its "frame reached by the field redefinition generated by a real multiple of S^(1)(Delta)" is replaced by the family-wise finite part: a field redefinition changes the components of every operator by the same multiples, so it cannot set P_1 = 0 on families with different translation components (R22, last sentence).',['R21','R22'],G+'#0')
put(P,5,'merged-with','Same content as design#2.',['R23'],G+'#2')
put(P,6,'applied','Applied verbatim.',['R31'])
put(P,7,'applied','Applied; the expansion was re-derived from Y = (1+sA+s^2B)^{-1} Z (1+sA+s^2B), A = Gram_tree^{-1}Gram_{1/2}, B = Gram_tree^{-1}Gram_1. Modified: "the tree norms are of relative order (unit)^{L-2}" replaced by "the tree norms of the two lengths differ by one power of the unit" (each further propagator carries hbar(-,-)^{-1}, proportional to the unit).',['R34'])
put(P,8,'applied','Applied verbatim.',['R33'])
# ---- ip_completeness2 ch08_issues
put(C,0,'applied','Applied, in one replacement with ip#1; "(at one loop the arity-zero part of a solution ...)" kept.',['R25'])
put(C,1,'applied','Applied verbatim.',['R26'])
put(C,2,'applied','Applied; p written a_P.',['R28'])
put(C,3,'applied','Applied; p written a_P; opening phrased "In a scheme with P_1=0 on the family (Definition def:c8-radial)" (delta#8, #14); the reason (ad_{S(P_b)})^adj = -ad_{S(K_bbar)} of design#3 added.',['R32'])
put(C,4,'applied','Applied; "it does not satisfy (ii) and (iii) there, since ... are not inversion covariant" stated affirmatively ("(ii) and (iii) there hold for Theta_n and fail for Theta_I, since ... change under the inversion"). The two-component formulas were re-derived from -(x,x)^{-2} rho(x) conj(lambda(Ix)), (x,x)^{-3} rho(x) conj(lambda*(Ix)) with Lemma c8-cc and eq:c8-pairing.',['R20'])
put(C,5,'applied','Applied verbatim.',['R29'])
put(C,6,'applied','Applied (R08); its last clause on products of functional-valued cochains realised in R07 (with design#6).',['R08','R07'])
put(C,7,'merged-with','Same removal of P_1 = 0 from the assumption and the same finite part (p/2){S^(1)(Delta),-}, reflection compatible exactly when p is real; merged into R21, R22.',['R21','R22'],G+'#0')

src_items = {}
d = json.loads((ROOT/D).read_text()); src_items[D] = [x['text'][:160] for x in d]
g = json.loads((ROOT/G).read_text()); src_items[G] = [x['label'] for x in g['entries']]
p = json.loads((ROOT/P).read_text()); src_items[P] = [x['label'] for x in p['corrections']]
c = json.loads((ROOT/C).read_text()); src_items[C] = [x['id']+' '+x['label'] for x in c['ch08_issues']]
entries = []
for src in (D, G, P, C):
    for i, lab in enumerate(src_items[src]):
        v = X[(src, i)]
        entries.append(dict(source=src, index=i, item=lab, status=v['status'], merged_with=v['merged_with'], reason=v['reason'],
                            edits=[dict(id=k, label=ED[k]['label'], old=ED[k]['old'], new=ED[k]['new']) for k in v['edits']]))
assert len(entries) == len(X) == 65
counts = {}
for e in entries: counts[e['status']] = counts.get(e['status'], 0) + 1
out = dict(file='step2/qa/ch08_consolidated/disposition.json', target='manuscript/ch08_definitions.tex',
           backup='step2/qa/ch08_consolidated/ch08_before_consolidation.tex', date='2026-10-29',
           counts=counts, edits_applied=[dict(id=e['id'], label=e['label'], sources=e['sources']) for e in E],
           entries=entries)
(HERE/'disposition.json').write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(counts)
