"""Build step2/qa/ch09_12_review/findings.json. Every `replace` is cut from the source file (exact text),
checked to occur exactly once; every \\label in a replaced span must survive in `by`."""
import json, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent.parent
ROOT = HERE.parents[2]
MS = ROOT / 'manuscript'
F9 = 'manuscript/ch09_cubic_and_dilatation.tex'
F10 = 'manuscript/ch10_two_letter_kernels.tex'
F11 = 'manuscript/ch11_composite_and_odd_kernels.tex'
F12 = 'manuscript/ch12_strict_construction.tex'
S12 = 'manuscript/supplements/ch12.tex'
_cache = {}


def text(f):
    if f not in _cache:
        _cache[f] = (ROOT / f).read_text(encoding='utf-8')
    return _cache[f]


def lit(f, s):
    n = text(f).count(s)
    if n != 1:
        raise SystemExit(f'count {n} for literal in {f}: {s[:120]!r}')
    return s


def span(f, start, end):
    t = text(f)
    if t.count(start) != 1:
        raise SystemExit(f'start count {t.count(start)} in {f}: {start[:100]!r}')
    i = t.index(start)
    j = t.index(end, i) + len(end)
    s = t[i:j]
    return lit(f, s)


def proof_of(f, label):
    """The proof environment following the statement that carries \\label{label}."""
    t = text(f)
    i = t.index('\\label{' + label + '}')
    a = t.index('\\begin{proof}', i)
    b = t.index('\\end{proof}', a) + len('\\end{proof}')
    return lit(f, t[a:b])


findings = []


def add(file, location, severity, claim, evidence, replace, by):
    labs_old = set(re.findall(r'\\label\{([^}]*)\}', replace))
    labs_new = set(re.findall(r'\\label\{([^}]*)\}', by))
    if labs_old - labs_new:
        raise SystemExit(f'label lost in {location}: {labs_old - labs_new}')
    assert replace != by
    findings.append(dict(file=file, location=location, severity=severity, claim=claim, evidence=evidence,
                         replace=replace, by=by))


# ---------------------------------------------------------------- R10 leftovers
add(F9, 'chapter introduction, first paragraph, last sentence', 'major',
    'R10: the introduction names Theorem c9-kappa (operator-level L o L = 0 on the bilinears, under existence hypotheses) as the place where the multiple of the classical transformation is fixed. In the strict route it is fixed by the action-level arity-two obstructions along (Delta,P_b) and (Delta,xi) on two legs, Theorem c9-theta2; Theorem c9-kappa is now labelled a cross-check in its own statement.',
    'Theorem c9-kappa statement: "a cross-check of that computation"; Theorem c9-theta2 (i),(ii); results/summary.txt sec. 2 and 6 ("theta = -1/3 -> computed from the construction twice ... Thm c9-kappa = cross-check").',
    lit(F9, r'and the arity-two equation along $(\Delta,P_b)$ fixes when the one-loop quantisation and the local representatives exist through arity two (Theorem \ref{thm:c9-kappa}).'),
    r'and the action-level arity-two equations along $(\Delta,P_b)$ and $(\Delta,\xi)$ fix on two legs (Theorem \ref{thm:c9-theta2}).')

add(F10, 'Theorem c10-brst, proof, first item, after eq. c10-finitesp', 'major',
    'Finite-special-vertex double count survives here: the finite vertex (eq c10-finitesp) is said to "enter L_1(M)". In the corrected reading it is the epsilon^0 coefficient of the counterterm partners of the bridges and cancels against the bare bridges, so it enters no component (Theorem c11-classes (ii), Lemma c12-bridge, Remark c10-scheme1 (d) of the same chapter).',
    'Theorem c11-classes (ii): "they are used up here"; Lemma c12-bridge eq c12-bridge; Remark c10-scheme1 (d); results/operator_final.json master_formula.revision ("F_x and the vacuum-pole term ... are cancelled by the bridges"); disposition.json QA-F1.',
    lit(F10, r'with $\{-,-\}_{H^m}$ the bracket built from $(H^m\otimes1)K_0$ in place of $K_0$; it depends on the regularisation and enters $L_1(M)$.'),
    r'with $\{-,-\}_{H^m}$ the bracket built from $(H^m\otimes1)K_0$ in place of $K_0$; it is the coefficient of $\varepsilon^0$ of the counterterm partners of the graphs with the special edge on a bridge, and with these graphs it tends to $0$ (Lemma \ref{lem:c12-bridge}).')

add(F10, 'Proposition c10-arity1, statement after the displayed formula', 'major',
    'R10: the two-body kernel is described as the limit of cycles with "exactly one defect insertion ... together with their images under Q", including "the arity-one counterterm vertex of O", and the one-body kernel as collecting "trees of the arity-zero finite counterterm with S(x)" and "defect or composite two-point functionals", with the scaling Ward identity cited for x = Delta. This is the pre-strict description. In the construction the two-body kernel is [eps^0] of the defect cycles of the pair (Definition c11-defect (iv)) and the one-body term is -{hbar F_1^(1)(x),X} alone (master formula, Theorem c11-classes (iii), Theorem c12-master (iii)); Theorem c10-D1 is a second evaluation. Chapter 11, Section 11.6 opens by citing this proposition as the definition of the two-body kernel.',
    'Theorem c11-classes (iii) eq c11-decompx; Theorem c12-master eq c12-master; operator_final.json master_formula.statement; summary.txt sec. 6.',
    span(F10, r'The two-body kernel is the $\varepsilon\to0$ limit of the local parts of the cycles through the pair', r'for $x=P_b$ it is Proposition \ref{prop:c10-P1}.'),
    r'The two-body kernel of $(X_a,X_b)$ along $x$ is the coefficient of $\varepsilon^0$ of the defect cycles of the pair, the cycles through the two letters with exactly one special edge $K_\varepsilon-K_0$ and exactly one vertex $S^{(1)}(x)$, placed on a line (for even $x$ the composite line of Definition \ref{def:c9-composite} with the vector field of $x$ in place of $v$), at the special edge, or as a cubic vertex (Definition \ref{def:c11-defect}); the one-body kernel of a letter $X$ is $-\{\hbar F_1^{(1)}(x),X\}$. This is the master formula of Theorem \ref{thm:c11-classes} (iii). For $x=\Delta$ and $x=P_b$, Theorem \ref{thm:c10-D1} and Proposition \ref{prop:c10-P1} evaluate the same flow a second time, through the scaling of the regulator and the composite kernel.')

add(F10, 'Remark c10-scheme1, item (b)', 'minor',
    'In the master formula hbar F_1^(0) enters the arity-one component only through the arity-one equation that fixes hbar F_1^(1)(x); the graphs of hbar F_1^(0) with S^(1)(x) belong to class (Ac), which vanishes by that equation. "enters through trees, one-body" restates the pre-strict one-body kernel.',
    'Theorem c11-classes (ii) (Ac); Theorem c12-master (ii); Proposition c12-xi2 (vi).',
    lit(F10, r'(b) $\hbar F_1^{(0)}$ enters through trees, one-body.'),
    r'(b) $\hbar F_1^{(0)}$ enters through the arity-one equation, which contains $\{\hbar F_1^{(0)},S^{(1)}(x)\}$, and therefore through $\hbar F_1^{(1)}(x)$ in the one-body kernel.')

# ---------------------------------------------------------------- finite special vertex remnants
add(F12, 'Theorem c12-master (iv)', 'minor',
    'Eq c11-pole no longer contains shift terms: in the rewritten Theorem c11-classes (iv) the special-edge parts of the counterterm brackets are grouped with the bridges at every power of epsilon. "with the shift terms j>=1 removed" refers to the former form of c11-pole.',
    'Theorem c11-classes (iv), text after eq c11-pole: "the special-edge parts of the counterterm brackets, with all their powers of epsilon, having been grouped with the bridges in (ii)".',
    lit(F12, r'and the coefficients of $\varepsilon^{-m}$ by \eqref{eq:c11-pole} with the shift terms $j\ge1$ removed, these being the poles of the bridge partners grouped in (ii);'),
    r'and the coefficients of $\varepsilon^{-m}$ by \eqref{eq:c11-pole}, the poles of the bridge partners being grouped with the bridges in (ii);')

add(F12, 'Remark c12-finitesp, third sentence', 'minor',
    'The remark compares the kernels of Definition c11-defect (iv) and of eq c12-kernel as if they could differ (they "coincide" on D<2). After the rewrite Definition c11-defect (iv) is Upsilon_x = [eps^0] D_x, identical to eq c12-kernel on every leg set; the comparison is a leftover of the double-count reading.',
    'Definition c11-defect (iii)-(iv) (current text); eq c12-kernel.',
    lit(F12, r'On the leg sets with superficial degree $D<2$ the defect cycles have no poles and the finite special vertex vanishes (Proposition \ref{prop:c11-oddbound}), so there the kernels of Definition \ref{def:c11-defect} (iv) and of \eqref{eq:c12-kernel} coincide; the composite pair kernels of Definition \ref{def:c11-sp} are coefficients of $\varepsilon^0$ of cycles and are the same in both readings.'),
    r'The kernel \eqref{eq:c12-kernel} is that of Definition \ref{def:c11-defect} (iv); on the leg sets with superficial degree $D<2$ the defect cycles have no poles and the finite special vertex vanishes (Proposition \ref{prop:c11-oddbound}), and the composite pair kernels of Definition \ref{def:c11-sp} are coefficients of $\varepsilon^0$ of cycles.')

add(F11, 'Proposition c11-freedom (v)', 'minor',
    'The cutoff shape is said to enter "the finite special vertices". These cancel against the bridges and enter no component (Theorem c11-classes (ii)); listing them among the scheme-dependent inputs of the components is a remnant of the double count.',
    'Theorem c11-classes (ii); Lemma c12-bridge.',
    lit(F11, r'of the special-edge cycles of $L_1^{(0)}$ and of the finite special vertices;'),
    r'and of the special-edge cycles of $L_1^{(0)}$;')

add(S12, 'supplement after_label thm:c12-master, last sentence', 'minor',
    'Audit narrative ("In Chapters 10 and 11 a finite vertex ... was added to the density"): describes a former error of the text, which FORMAT.md excludes, and is false once Chapters 10 and 11 are corrected.',
    'FORMAT.md: "No audit notes"; Definition c11-defect (iii), Theorem c11-classes (ii).',
    lit(S12, r'In Chapters 10 and 11 a finite vertex built from the regulator insertion on an external line and the counterterm was added to the density; the graph with the regulator insertion on that external line cancels it exactly (Lemma \ref{lem:c12-bridge}), so the density is the finite part of the graphs through the two fields alone.'),
    r'The finite vertex built from the regulator insertion on an external line and the counterterm cancels exactly against the graph with the regulator insertion on that external line (Lemma \ref{lem:c12-bridge}), so the density is the finite part of the graphs through the two fields alone.')

# ---------------------------------------------------------------- claims stronger than the evidence
add(F12, 'Remark c12-method, first sentence', 'major',
    '"Every number is computed by two independent evaluations" contradicts the caption of Table c12-spectra (last three rows by the first evaluation only), Remark c12-next (second evaluation on F5-F7 listed as open) and the results (F5-F7, P_1 on F5-F7: one engine).',
    'operator_final.json families_spectra F5/F6/F7 and "P_1 on F3-F7": "computed by one engine + referee on F5-F7"; Table c12-spectra caption.',
    lit(F12, r'Every number is computed by two independent evaluations. The first evaluates'),
    r'Every number is computed by two independent evaluations, except where a statement names one of them. The first evaluates')

add(F12, 'Remark c12-method, penultimate sentence', 'major',
    'Single-pair tables: "at most three letter derivatives in the first evaluation". The first route covers n <= 2; the n = 3 tables are open. Proposition c11-masterfam (v) of the same text says "at most two letter derivatives in the simplex evaluation".',
    'operator_final.json F8_tables.evidence ("engine A 4634/4634 rows D<2 and 646/646 rows D>=2 (n <= 2); engine B every complete row n <= 1") and open list ("engine A n=3 tables ... resumable"); Proposition c11-masterfam (v).',
    lit(F12, r'the single-pair tables cover the $45$ pairs of slots with defect cycles with at most three letter derivatives in the first evaluation and one in the second.'),
    r'the single-pair tables cover the $45$ pairs of slots with defect cycles with at most two letter derivatives in the first evaluation and one in the second, on the rows on which the second is complete.')

add(F12, 'Remark c12-method, last sentence', 'minor',
    '"fixed ... by a regression to the values of Chapters 4 and 9" is engine language (FORMAT.md: no audit notes); state the normalisation as such. The a priori global sign of the defect cycles is an open item, so the sentence must stay a normalisation statement.',
    'operator_final.json open: "a priori overall sign of defect cycles (gate-fixed)"; FORMAT.md.',
    lit(F12, r'The overall sign of the defect cycles is fixed in each evaluation by a regression to the values of Chapters 4 and 9.'),
    r'The overall sign of the defect cycles is normalised in each evaluation to the arity-zero obstruction of Theorem \ref{thm:c4-obstruction} and to the cubic obstruction of Section 9.1.')

add(F12, 'Proposition c12-spectra, second sentence', 'major',
    'The translation sums on the three families with gaugino or field-strength letters (F5-F7) come from one evaluation; under the convention of the chapter introduction ("computed by the two independent evaluations ... unless a statement says otherwise") the sentence claims two.',
    'operator_final.json "P_1 on F3-F7": "computed by two engines + referee on F3, F4; computed by one engine + referee on F5-F7"; Table c12-spectra caption.',
    lit(F12, r'On twist two, the fermion $\mathfrak{sl}(2)$ sector and the three families with gaugino or field-strength letters, $\sum_{a<b}\pm\Upsilon_{P_n}=-\tfrac13\partial_n$ units, so $P_1=-(\tfrac13+t_P)\partial_n$ units;'),
    r'On twist two and the fermion $\mathfrak{sl}(2)$ sector, and in the first evaluation of Remark \ref{rem:c12-method} on the three families with gaugino or field-strength letters, $\sum_{a<b}\pm\Upsilon_{P_n}=-\tfrac13\partial_n$ units, so $P_1=-(\tfrac13+t_P)\partial_n$ units;')

add(F12, 'Proposition c12-spectra, proof', 'major',
    'R10 / unverified locality: the proof derives the spectra from Theorem c12-D1 (4 lc + composite pair kernels), whose sum rule rests on Theorem c12-local; the identities G-Cid and G-typeI behind Lemma c12-Cid and Proposition c12-typeI are untested, and Remark c12-next lists them as next computations. The spectra were computed from the master formula (defect cycles + one-body term), which is how Proposition c11-masterfam states them, and Remark c12-specref says "the numbers come from eq c12-master". The proof should rest on the master formula, with Theorem c12-D1 as second evaluation.',
    'operator_final.json open: "G-Cid, G-typeI (locality theorem untested)"; families_spectra statuses; Proposition c11-masterfam and its proof; Remark c12-specref.',
    proof_of(F12, 'prop:c12-spectra'),
    '\\begin{proof}\n' + r'Theorem \ref{thm:c12-master} along $\Delta$ and $P_n$: the two-body kernels are the coefficients of $\varepsilon^0$ of the defect cycles of the pairs of letters of the families, evaluated through the two-field parts (Proposition \ref{prop:c10-twist} and the remark before Proposition \ref{prop:c11-targets}), and the one-body term is the bracket with $\hbar F_1^{(1)}(\Delta)$ at $t_\Delta=-\tfrac13$ (Proposition \ref{prop:c11-masterfam}). Theorem \ref{thm:c12-D1} with the composite pair kernels of Theorem \ref{thm:c11-families} gives the same spectra, a second evaluation. Scheme independence is Theorem \ref{thm:c8-invariants} (iii), $t_\Delta$ being fixed by Theorem \ref{thm:c12-ar2}. The translation identity is (R7) of Corollary \ref{cor:c12-rel} along $P_n$.' + '\n\\end{proof}')

add(F12, 'Remark c12-pair, last two sentences', 'major',
    'The pair form (Conj op-pair) is conjectural beyond (A_perp,A_perp); the remark states that it "holds entry by entry in both evaluations" on every pair of the tables and that "a proof pair by pair follows" from Theorem c12-conj. The second route has 11 incomplete rows and covers n <= 1 only, and the pair-by-pair argument needs the agreement of the s int omega X_a X_b terms of the two evaluations, which is not established.',
    'operator_final.json F8_tables: "conjectured (proved on (A_perp,A_perp))", "11 engine-B rows incomplete (no verdict)", "engine B every complete row n <= 1"; summary.txt sec. 4 "Conj op-pair [Cj]".',
    span(F12, r'On every pair of the tables of Remark \ref{rem:c12-method} the pair form holds entry by entry in both evaluations.', r'by the classes (F), (C0) in the first.'),
    r'On every complete row of the tables of Remark \ref{rem:c12-method} the pair form holds entry by entry, in both evaluations at one letter derivative and in the first at two. For a pair whose functional $\int\omega\,X_aX_b$ is not $s$-closed, Theorem \ref{thm:c12-conj} gives the pair form once the terms in $s\int\omega\,X_aX_b$ of the two evaluations, Theorem \ref{thm:c12-local} (vi) in the second and the classes (F), (C0) in the first, are shown to agree.')

add(F11, 'Proposition c11-masterfam (iv)', 'minor',
    'The matrix identity D_1(d_n O) = d_n D_1(O) is claimed on "the families of (ii), (iii)", which include the gluon sl(2) sector; (R7) was computed on twist two, the fermion sector, the two-gaugino family and the mixed family, not on the gluon family (Corollary c12-rel and Proposition c12-spectra list exactly these).',
    'operator_final.json "(R7) along P_n": "exact matrix identity on F3 (0..6), F4 (1..5), F5 (0..3), F7 (0..4)"; Corollary c12-rel, last sentence.',
    lit(F11, r"and $D_1(\partial_nO)=\partial_nD_1(O)$ on twist two and the families of (ii), (iii), as an identity of matrices from $S'$ to $S'+1$;"),
    r"and $D_1(\partial_nO)=\partial_nD_1(O)$ on twist two, the fermion $\mathfrak{sl}(2)$ sector, the family of two gaugino components with $Z'F_{11}$ and the family $(F_{11},\psi)$, as an identity of matrices from $S'$ to $S'+1$;")

add(F12, 'Remark c12-freetable, Table c12-free, rows Sigma* and K_b', 'minor',
    'Row Sigma*: the difference is fixed by (xi_+,eta_-) (computed); (xi_-,eta_+) is open, so "(xi_pm,eta_mp)" claims an uncomputed pair. Row K_b: "up to the inner direction" names no equation; t_K is fixed only by (K_b,P_b\') and (eta_+,eta_-), both open (Theorem c12-fix table, Remark c12-next).',
    'summary.txt sec. 2 and 8 (open: (K_b,P_b\'), (eta,eta\'), (xi_-,eta_pm)); Table c12-ar2 rows.',
    span(F12, r'$\Sigma^*$ & $S^V(\eta)$, $S^\perp(\eta)$; difference fixed by $(\xi_\pm,\eta_\mp)$ & up to the inner and outer directions\\', r'$K_b$ & $S^{(1)}(K_b)_d$, $G^K_1,\dots,G^K_4$ & up to the inner direction\\'),
    r'$\Sigma^*$ & $S^V(\eta)$, $S^\perp(\eta)$; difference fixed by $(\xi_+,\eta_-)$ & up to the inner and outer directions\\' + '\n' + r"$K_b$ & $S^{(1)}(K_b)_d$, $G^K_1,\dots,G^K_4$ & $(K_b,P_{b'})$, $(\eta_+,\eta_-)$, up to the inner direction\\")

add(F12, 'Theorem c12-algebra, last sentence', 'minor',
    '{Q,S} is computed on two legs along (xi_+,eta_+) and (xi_+,eta_-) only; (xi_-,eta_pm) is open.',
    'summary.txt sec. 2 (computed pairs) and sec. 8 (open (xi_-,eta_pm)); Table c12-ar2.',
    lit(F12, r'Computed on two legs: $[D,P]$, $[D,Q]$, $[Q,P]$, $[Q,Q]$ and $\{Q,S\}$, by Theorem \ref{thm:c12-ar2}.'),
    r'Computed on two legs: $[D,P]$, $[D,Q]$, $[Q,P]$, $[Q,Q]$, and $\{Q,S\}$ along $(\xi_+,\eta_\pm)$, by Theorem \ref{thm:c12-ar2}.')

add(F12, 'Theorem c12-fix (iv)', 'minor',
    'Two problems in (iv). (a) Theorem c12-fix is structural (rho left unspecified), but t_{eta_pm} = t_{xi_mp} follows only from the computed classes: t_{eta_+} = t_{xi_-} from (xi_+,eta_+) and (xi_+,xi_-) of Theorem c12-ar2; t_{eta_-} = t_{xi_+} needs (xi_-,eta_-), which is open. (b) Sign: (ii) says phi^(0) = c S^(1)(Delta) units generates t_z -> t_z + c lambda_z, i.e. t_P -> t_P - 2c, while (iv) says the same phi shifts hbar F_1^(1)(P_b) by +2c S^(1)(P_b). Lemma c9-SDelta (i) gives {S^(1)(Delta),S^(1)(x)} = v.S^(1)(x) = -lambda_x S^(1)(x) (weight 4 - lambda_x/2), so the proof\'s "c{S(Delta),S(x)} = c lambda_x S(x)" also has the opposite sign; the same convention is used in Theorem c9-theta2 (iii) and Theorem c11-oddar2 (ii). The direction is unaffected (c is free); the replacement makes (iv) consistent with (ii) and the sign of the generating phi should be fixed once from Lemma c9-SDelta (i).',
    'Lemma c9-SDelta (i); proof of Theorem c9-theta2 (i) ({S(Delta)_d,S(P_b)_d} = 2 S(P_b)_d); Corollary c12-inv; summary.txt sec. 8 open (xi_-,eta_pm).',
    lit(F12, r'(iv) in a translation-covariant scheme $t_P=t_\Delta$ and $t_{\eta_\pm}=t_{\xi_\mp}$ up to the outer direction, the inner direction being excluded since $\phi^{(0)}=c\,S^{(1)}(\Delta)$ units shifts $\hbar F_1^{(1)}(P_b)$ by $2c\,S^{(1)}(P_b)$ units and $\hbar F_1^{(1)}(\Delta)$ by $0$;'),
    r'(iv) in a translation-covariant scheme $t_P=t_\Delta$, the inner direction being excluded since it shifts $t_P$ by $c\lambda_{P_b}=-2c$ and $t_\Delta$ by $0$; with the computed classes of Theorem \ref{thm:c12-ar2}, $t_{\eta_+}=t_{\xi_-}$, a relation invariant under the outer direction;')

add(S12, 'supplement after_label thm:c12-sol, last sentence', 'minor',
    'The coefficients of the four second-order free-field symmetries G^K_j are open; the three-leg classes are not computed, so "fixed at three legs" asserts an unperformed computation (Theorem c12-cubic gives a criterion only).',
    'summary.txt sec. 8 open: "G^K_j and SU(4)-15 classes", "three/four-leg classes"; Remark c12-next.',
    lit(S12, r'whose coefficients are fixed at three legs.'),
    r'whose coefficients enter the parts of the equations with three legs.')

add(F9, 'Remark c9-methodD, fourth sentence', 'minor',
    'Remark c9-methodD counts Remark c9-twoeval among the independent second evaluations that "agree with the defect cycles", but Proposition c9-triangle and Remark c9-twoeval state that the two coefficients 8 of the composite triangle are fixed by that agreement; on the bilinears the second evaluation is therefore a determination, not a check.',
    'Proof of Proposition c9-triangle: "its two coefficients are 8 and 8, the values fixed by Remark c9-twoeval"; Remark c9-twoeval last sentence; Remark c10-inputs.',
    lit(F9, r'The second evaluation of Proposition \ref{prop:c9-structure} (iii), Remark \ref{rem:c9-Pi} and Remark \ref{rem:c9-twoeval} evaluates the composite cycles and the graphs with $S^{(1)}(\Delta)$ on a leg directly as $\varepsilon\to0$ limits in a Feynman-parameter representation, with $E$ acting on Fourier components as $-(4+(p,\partial_p))$; it agrees with the defect cycles in every sector and placement class by placement class.'),
    r'The second evaluation of Proposition \ref{prop:c9-structure} (iii) and Remark \ref{rem:c9-Pi} evaluates the composite cycles and the graphs with $S^{(1)}(\Delta)$ on a leg directly as $\varepsilon\to0$ limits in a Feynman-parameter representation, with $E$ acting on Fourier components as $-(4+(p,\partial_p))$; it agrees with the defect cycles in every sector and placement class by placement class. On the scalar bilinears the agreement of Remark \ref{rem:c9-twoeval} fixes the two coefficients of the composite triangle.')

add(F11, 'Section 11.6, end of the opening paragraph', 'minor',
    'Notation clash: Section 11.6 sets L := log(3/2), while Theorem c11-classes (iv) of the same section (eq c11-logcov) and Lemma c11-Cid use L^(1)(x) for the coefficient of log eps in the arity-one counterterm cochain (Chapter 12 uses L, L^(0), L^(1) for the logarithmic coefficient cochain and L(O), L_1 for the local representative). The minimal fix makes the two uses explicit in Section 11.6; a rename of the cochain throughout Chapters 11-12 (e.g. to the lc of the operator side) is the clean fix.',
    'Theorem c11-classes (iv), eq c11-logcov; Lemma c11-Cid; Definition c12-comp; Prop c11-oddsol (L = log 3/2).',
    lit(F11, r'In this section $\psi$ denotes the finite part of the operators of Definition \ref{def:c8-renorm}, and $L:=\log\tfrac32$.'),
    r'In this section $\psi$ denotes the finite part of the operators of Definition \ref{def:c8-renorm}, and $L:=\log\tfrac32$; with its superscript, $L^{(1)}(x)$ is the coefficient of $\log\varepsilon$ in $I^{CT}_{\mathfrak{psl}(4|4)}[\varepsilon]^{(1)}(x)$ (Lemma \ref{lem:c11-Cid}).')

# ---------------------------------------------------------------- duplication: Chapters 9-11 prove, Chapter 12 cites
add(F12, 'Lemma c12-wz, proof', 'minor',
    'Duplication: Lemma c12-wz (i)-(iii) is Lemma c9-sform with the arity-zero components added, and its proof repeats that proof line by line. Proposed division: Lemma c9-sform proves (Chapter 9 uses it first, Chapter 11 cites it); Chapter 12 states and cites. The Wess-Zumino signs agree in both places and with results/summary.txt sec. 0.',
    'Lemma c9-sform and proof; summary.txt sec. 0 (corrected WZ).',
    proof_of(F12, 'lem:c12-wz'),
    '\\begin{proof}\n' + r'Lemma \ref{lem:c9-sform} and its proof, with $\Lambda^{(0)}=Q+\{I,-\}$ and $S(x)=S^{(1)}(x)$; the components of arity zero are $(sF)^{(0)}=(Q+\{I,-\})F^{(0)}$ and $(s\,\mathrm{Ob})^{(0)}=0$, and the closedness of $\mathcal R^{(1)}(x)$ is the arity-one component of $s\big(\hbar\,\mathrm{Ob}-s(\hbar F_1)\big)=0$ when the arity-zero component of $\hbar\,\mathrm{Ob}-s(\hbar F_1)$ vanishes.' + '\n\\end{proof}')

add(F12, 'Theorem c12-theta, proof', 'minor',
    'Duplication: Theorem c12-theta and its proof repeat Theorem c9-theta2 (i) and its proof (same bracket count 2 theta_d S(P_b)_d, same independence argument). Proposed division: Theorem c9-theta2 proves; Chapter 12 cites.',
    'Theorem c9-theta2 (i) and proof.',
    proof_of(F12, 'thm:c12-theta'),
    '\\begin{proof}\n' + r'Theorem \ref{thm:c9-theta2} (i) and its proof, in the notation $\theta^\Delta_d$, $\theta^{P_b}_d$ of Definition \ref{def:c12-gens}.' + '\n\\end{proof}')

add(F12, 'Theorem c12-ar2, proof', 'minor',
    'Duplication: Table c12-ar2 and eq c12-t repeat Table c11-oddar2 and eq c11-trel (Theorem c11-oddar2), with the same proof. Proposed division: Theorem c9-theta2 and Theorem c11-oddar2 prove; Chapter 12 keeps the statement with the explicit obstructions eq c12-ob2 and cites. All entries agree with results/action_final.json and summary.txt sec. 2.',
    'Theorem c11-oddar2, Table c11-oddar2, eq c11-trel; summary.txt sec. 2.',
    proof_of(F12, 'thm:c12-ar2'),
    '\\begin{proof}\n' + r'Theorem \ref{thm:c9-theta2} for the rows along $\Delta$, with Theorem \ref{thm:c12-theta} and $o_d=-\tfrac23$ units, and Theorem \ref{thm:c11-oddar2} for the other rows and for \eqref{eq:c12-t}; the obstructions \eqref{eq:c12-ob2} are the sums \eqref{eq:c12-local2} over the cycles of Proposition \ref{prop:c12-classes2}.' + '\n\\end{proof}')

add(F12, 'Theorem c12-ob1, proof', 'minor',
    'Duplication: the two-leg arity-one obstructions of Theorem c12-ob1 are those of Theorem c4-obstruction, Theorem c9-quadratic (ii), Theorem c9-P1 (i) and Proposition c11-oddsol (i)-(iv), recomputed in the proof. Proposed division: Chapters 9 and 11 prove; Chapter 12 collects and cites.',
    'Theorem c9-quadratic (ii) eq c9-ObDelta; Theorem c9-P1 (i); Proposition c11-oddsol; summary.txt sec. 1 (all entries agree).',
    proof_of(F12, 'thm:c12-ob1'),
    '\\begin{proof}\n' + r'Theorem \ref{thm:c4-obstruction} in arity zero, Theorem \ref{thm:c9-quadratic} (ii) along $\Delta$, Theorem \ref{thm:c9-P1} (i) along $P_b$, and Proposition \ref{prop:c11-oddsol} (i)--(iv) along $\xi$, $\eta$, $M_\omega$ and $R_r$, each obstruction being the sum \eqref{eq:c12-local1} over the cycles of Theorem \ref{thm:c12-classes} on the bases of Propositions \ref{prop:c12-xi2} and \ref{prop:c12-two}; the position part along $\eta$ is Proposition \ref{prop:c12-pos}. $L^{(1)}(\xi)|_2$ is the coefficient of $\log\varepsilon$ of the composite and cubic cycles, and $\{\mathcal N,-\}$ multiplies the arrows of $T_1$, $T_2$ ($\lambda\to A^*$) by $-2+4=2$ and those of $T_3$, $T_4$ ($A\to\lambda^*$) by $2-4=-2$.' + '\n\\end{proof}')

add(F12, 'Theorem c12-sol, proof', 'minor',
    'Duplication: the two-leg solutions of Theorem c12-sol are those of Theorem c9-quadratic (iii), Theorem c9-P1 (i) and Proposition c11-oddsol. Proposed division: Chapters 9 and 11 prove; Chapter 12 collects and cites.',
    'Theorem c9-quadratic (iii) eq c9-F1Delta; Theorem c9-P1 (i); Proposition c11-oddsol (i)-(v); summary.txt sec. 1.',
    proof_of(F12, 'thm:c12-sol'),
    '\\begin{proof}\n' + r'Theorem \ref{thm:c9-quadratic} (iii) along $\Delta$, Theorem \ref{thm:c9-P1} (i) along $P_b$, and Proposition \ref{prop:c11-oddsol} along $M_\omega$, $R_r$, $\xi$, $\eta$ and $K_b$, with the free directions of Propositions \ref{prop:c12-xi2} and \ref{prop:c12-two}. Along $\xi$, \eqref{eq:c12-xisol} with $r_0=0$, the coefficient $\tfrac{27}8$ of $g_2$ in the obstruction being the contribution $\tfrac{27}{16}\cdot2$ of the finite counterterm, and $-2r_1=-2r_3=\log\tfrac32-\tfrac{79}{12}$.' + '\n\\end{proof}')

add(F12, 'Section 12.6, opening paragraph, last sentence', 'minor',
    'Duplication of the locality theorem: Definition c12-Xi, Lemma c12-Cid, Lemma c12-BD, Proposition c12-typeI and Theorem c12-local repeat Definition c11-Xi, Lemma c11-Cid, Lemma c11-trees, Proposition c11-typeI and Theorem c11-local, which Chapter 11 now proves in full (not an outline). Proposed division: Chapter 11 proves (for Delta, P_b; it is used there by Corollary c11-D1local and Sections 11.3-11.5); Chapter 12 states the theorem for x in {Delta,P_b,M_omega}, with the rotations added by Remark c12-localM, and cites. The next four findings replace the Chapter 12 proofs accordingly.',
    'Chapter 11 Section 11.2 (Lemma c11-Cid, Lemma c11-trees, Prop c11-typeI, Thm c11-local and proofs); rewrite_log.json ch11 def:c11-Xi insert and proof:thm:c11-local replace.',
    lit(F12, r'The theorem of this section completes the outline of Theorem \ref{thm:c11-local}.'),
    r'Lemma \ref{lem:c12-Cid}, Lemma \ref{lem:c12-BD}, Proposition \ref{prop:c12-typeI} and Theorem \ref{thm:c12-local} are Lemma \ref{lem:c11-Cid}, Lemma \ref{lem:c11-trees}, Proposition \ref{prop:c11-typeI} and Theorem \ref{thm:c11-local}, proved there for $\Delta$ and $P_b$, together with the rotations of Remark \ref{rem:c12-localM}.')

add(F12, 'Lemma c12-Cid, proof', 'minor',
    'Duplication (see the finding on the opening of Section 12.6): the proof repeats that of Lemma c11-Cid.',
    'Lemma c11-Cid and proof.',
    proof_of(F12, 'lem:c12-Cid'),
    '\\begin{proof}\n' + r'Lemma \ref{lem:c11-Cid} and its proof, which use of the composite line only $[Q,P\star x.P\star]=0$, true along $P_b$, $\Delta$ and $M_\omega$ (Proposition \ref{prop:c12-composite} (i), (ii)); the absence of poles of $I^{CT,(1)}(x)$ is Proposition \ref{prop:c12-poles}.' + '\n\\end{proof}')

add(F12, 'Lemma c12-BD, proof', 'minor',
    'Duplication: the proof repeats that of Lemma c11-trees.',
    'Lemma c11-trees and proof.',
    proof_of(F12, 'lem:c12-BD'),
    '\\begin{proof}\n' + r'Lemma \ref{lem:c11-trees} and its proof; (B) is Lemma \ref{lem:c12-bridge} with $G$ the composite cycle through the pair with its trees.' + '\n\\end{proof}')

add(F12, 'Proposition c12-typeI, proof', 'minor',
    'Duplication: the proof repeats that of Proposition c11-typeI.',
    'Proposition c11-typeI and proof.',
    proof_of(F12, 'prop:c12-typeI'),
    '\\begin{proof}\n' + r'Proposition \ref{prop:c11-typeI} and its proof, with Lemma \ref{lem:c12-Cid} in place of Lemma \ref{lem:c11-Cid}.' + '\n\\end{proof}')

add(F12, 'Theorem c12-local, proof', 'minor',
    'Duplication: the proof repeats that of Theorem c11-local (items (i)-(vi) correspond one to one).',
    'Theorem c11-local and proof.',
    proof_of(F12, 'thm:c12-local'),
    '\\begin{proof}\n' + r'Theorem \ref{thm:c11-local} (i)--(vi) and its proof, with Lemmas \ref{lem:c12-Cid}, \ref{lem:c12-BD} and Proposition \ref{prop:c12-typeI} in place of Lemmas \ref{lem:c11-Cid}, \ref{lem:c11-trees} and Proposition \ref{prop:c11-typeI}; along $M_\omega$ every step applies by Remark \ref{rem:c12-localM}.' + '\n\\end{proof}')

out = HERE / 'findings.json'
out.write_text(json.dumps(findings, indent=1, ensure_ascii=False) + '\n', encoding='utf-8')
# final verification: every replace occurs exactly once in its file
for k, f in enumerate(findings):
    n = text(f['file']).count(f['replace'])
    assert n == 1, (k, f['location'], n)
# no two findings overlap in one file
for a in range(len(findings)):
    for b in range(a + 1, len(findings)):
        fa, fb = findings[a], findings[b]
        if fa['file'] == fb['file']:
            t = text(fa['file'])
            ia, ib = t.index(fa['replace']), t.index(fb['replace'])
            if not (ia + len(fa['replace']) <= ib or ib + len(fb['replace']) <= ia):
                raise SystemExit(f'overlap {a} {b}')
from collections import Counter
print(len(findings), 'findings', Counter(f['severity'] for f in findings))
