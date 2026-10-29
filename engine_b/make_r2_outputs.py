"""Assemble engine B round-2 results into results_r2.json (T1, T2, T3, formula_chains, checks)."""
import glob, json
import sympy as sp

U = 'u = hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2)'
SETFILES = [('(A|_V, A|_V, A|_V)', 'r2_T1_AV_AV_AV.json'), ('(A|_V, A|_Wperp, A|_Wperp)', 'r2_T1_AV_Ap_Ap.json'),
            ('(lambda, lambda, A|_V)', 'r2_T1_lam_lam_AV.json'), ('(lambda, lambda, A|_Wperp)', 'r2_T1_lam_lam_Ap.json'),
            ('(c, A|_V, A*|_V)', 'r2_T1_c_AV_As_V.json'), ('(c, A|_Wperp, A*|_Wperp)', 'r2_T1_c_Ap_As_p.json'),
            ('(c, lambda, lambda*)', 'r2_T1_c_lam_lams.json'), ('(c, c, c*)', 'r2_T1_c_c_cs.json'),
            ('(A|_V)^4', 'r2_T1_AV4.json'), ('(A|_V)^2 (A|_Wperp)^2', 'r2_T1_AV2Ap2.json'), ('(A|_Wperp)^4', 'r2_T1_Ap4.json')]
TREE = {
    '(A|_V, A|_V, A|_V)': 'the (A|_V)^3 part of int dvol(-1/2 (dA, 1/2[A^A]))',
    '(A|_V, A|_Wperp, A|_Wperp)': 'the A|_V (A|_Wperp)^2 part of int dvol(-1/2 (dA, 1/2[A^A]))',
    '(lambda, lambda, A|_V)': 'the A|_V part of int dvol 1/2 (lambda, (dslash_A - dslash) lambda)',
    '(lambda, lambda, A|_Wperp)': 'the A|_Wperp part of int dvol 1/2 (lambda, (dslash_A - dslash) lambda)',
    '(c, A|_V, A*|_V)': 'the A|_V part of int dvol ([A,c],A*)',
    '(c, A|_Wperp, A*|_Wperp)': 'the A|_Wperp part of int dvol ([A,c],A*)',
    '(c, lambda, lambda*)': 'int dvol ([c,lambda],lambda*)',
    '(c, c, c*)': 'int dvol 1/2 ([c,c],c*)',
    '(A|_V)^4': 'the (A|_V)^4 part of int dvol(-1/4 (1/2[A^A], 1/2[A^A]))',
    '(A|_V)^2 (A|_Wperp)^2': 'the (A|_V)^2 (A|_Wperp)^2 part of int dvol(-1/4 (1/2[A^A], 1/2[A^A]))',
    '(A|_Wperp)^4': 'the (A|_Wperp)^4 part of int dvol(-1/4 (1/2[A^A], 1/2[A^A]))',
}
WEIGHT_TEXT = {
    'AAA': '{I,Psi} on a vertex with n_A fields A, n_A* fields A*, n_lambda fields lambda, n_lambda* fields lambda*, n_c fields c, n_c* fields c*: '
           'weight 2(n_A - n_A*) + 4(n_lambda - n_lambda*) + (n_c - n_c*); AAA: 2*3 = 6',
    'lamlamA': '2*1 + 4*2 = 10', 'AcAs': '2*(1-1) + 1 = 1', 'clamlams': '4*(1-1) + 1 = 1', 'cccs': '(2-1) = 1', 'AAAA': '2*4 = 8'}
LEGEND = ('Plane-wave evaluation variables (code level, used only to state multilinear forms): leg j carries polarisation a_j (A|_V), b_j (A|_Wperp), '
          'lambda_j (Sigma), incoming momentum p_j; (x.y) is the metric pairing; colour structures: ([X_0,X_1],X_2) for three legs, '
          '([X_0,X_1],[X_2,X_3]) and ([X_0,X_3],[X_1,X_2]) for four legs (([X_0,X_2],[X_1,X_3]) = their sum by Jacobi), and '
          '(1/24) sum over orderings of tr_g(ad_{X_s1} ad_{X_s2} ad_{X_s3} ad_{X_s4}) (symmetrised adjoint trace, units hbar/(16 pi^2), '
          'not proportional to tr_g(ad_X ad_Y)/(X,Y)). All other colour structures are in units u.')

def struct_text(k):
    k = k.replace('cadj K12:', '([X_0,X_1],[X_2,X_3]) x ').replace('cadj K14:', '([X_0,X_3],[X_1,X_2]) x ')
    k = k.replace('str4:', '(1/24) sum_sigma tr_g(ad ad ad ad) x ')
    for j in range(4):
        k = k.replace(f'a{j}', f'a_{j}').replace(f'b{j}', f'b_{j}').replace(f'p{j}', f'p_{j}')
    return k

def t1():
    out = {}; chains = {}
    for name, fn in SETFILES:
        d = json.load(open(fn))
        graphs = []
        for g in d['graphs']:
            mult = g['log_eps_multiple_of_tree']
            if mult != 'not proportional':
                le = f'{mult} x tree' if mult not in ('0',) else '0'
            else:
                le = {struct_text(k): v for k, v in g['log_eps_structure'].items()}
            gid = g['id'].replace('[A:V]', '(A|_V)').replace('[A:perp]', '(A|_Wperp)').replace('[A:V,A:V]', '(A|_V,A|_V)') \
                .replace('[A:V,A:perp]', '(A|_V,A|_Wperp)').replace('[A:perp,A:perp]', '(A|_Wperp,A|_Wperp)') \
                .replace('[lam:all]', '(lambda)').replace('[c:all]', '(c)').replace('[As:V]', '(A*|_V)').replace('[As:perp]', '(A*|_Wperp)') \
                .replace('[lams:all]', '(lambda*)').replace('[cs:all]', '(c*)')
            gid = g['topology'] + ': ' + gid
            graphs.append(dict(id=gid, vertices=g['vertices'], internal_lines=g['internal_lines'], symmetry_factor=g['symmetry_factor'],
                               log_eps=le, eps_inv='0', eps_inv_half='0'))
            n = {'bubble': 2, 'triangle': 3, 'box': 4}[g['topology']]
            w = d['omega']
            chains[f'T1 {name} | {gid}'] = dict(
                amplitude=f'sum over the words of this graph class (all assignments of the external legs to its vertices, both orientations): '
                          f'(hbar/{2 * n}) sign(odd reordering) tr_g(M_{n}...M_1) int d^4q/(2pi)^4 STr[P(q-P_0) Pi d^2I_(v1) ... P(q-P_{n - 1}) Pi d^2I_(v{n})], '
                          f'vertices {g["vertices"]}, lines {g["internal_lines"]}',
                numerator=f'STr numerator homogeneous of degree 2E-4+w = {2 * n - 4 + w} in (q, external momenta); P numerators: 2 (-,-)^-1 (A), '
                          f'i rho(q)|_(Sigma*->Sigma) (lambda), -i q (c-A*); colour tr_g(M_n...M_1) with M = ad_X (cubic vertex), '
                          f'ad_[X,Y], -ad_Y ad_X, -ad_X ad_Y (quartic vertex, three Hodge contractions); tr_g(ad_X ad_Y ad_Z) = (1/2)(tr_g(ad ad)/( , ))([X,Y],Z)',
                schwinger='P(eps,infty)(q - P_e) = N(q - P_e) int_eps^infty dt_e e^{-t_e (q-P_e)^2}, one t_e per edge; Gaussian q-integral with T = sum_e t_e',
                eps_expansion=('eps d/deps X = -sum_e X|_{t_e = eps}; Taylor in the external momenta (lambda P): order lambda^j scales as eps^{(j - w)/2}; '
                               'radial int_0^infty q^{3+2m-2b} e^{-E eps q^2} dq and angular average exact (24 short roots of F4, a spherical 5-design): '
                               f'eps d/deps X = -(16 pi^2)^-1 < [lambda^{w}] Num(qhat, lambda P) prod_f |qhat - lambda P_f|^-2 > '
                               + ('(the lambda^0 term, an odd polynomial on the sphere, gives eps^{-1/2} x 0)' if w == 1 else '')),
                local_part=f'log eps coefficient (units u): {le if isinstance(le, str) else json.dumps(le)}; eps^-1 and eps^-1/2: 0')
        out[name] = dict(external_set=name, superficial_degree=d['omega'], tree_vertex=TREE[name], graphs=graphs,
                         total=f'{d["total_multiple_of_tree"]} x tree  (log eps, units u; eps^-1 = eps^-1/2 = 0)',
                         total_structure_fit={struct_text(k): v for k, v in d['total_fit'].items()},
                         predicted_from_A=f'{d["predicted_from_A"]} x tree  [{WEIGHT_TEXT[d["tree_term"]]}]',
                         verdict='agrees' if str(d['total_multiple_of_tree']) == str(d['predicted_from_A']) else 'DISAGREES',
                         n_random_exact_configurations=d['n_configs'])
    return out, chains

def coupling_direction():
    aV, aP, aL, aC, k = sp.symbols('alpha_AV alpha_AWperp alpha_lambda alpha_c kappa')
    data = {  # computed multiples of the free/tree terms (log eps, units u)
        'kinetic A|_V (-1/4 int (dA|_V,dA|_V))': (2 * aV + k, 4), 'kinetic A|_Wperp': (2 * aP + k, 4),
        'kinetic lambda (1/2 int (lambda, dslash lambda))': (-2 * aL + k, 8), 'kinetic ghost (int (dc, A*))': (-aV - aC + k, -1),
        '(A|_V)^3': (3 * aV + k, 6), 'A|_V (A|_Wperp)^2': (aV + 2 * aP + k, 6), 'lambda lambda A|_V': (-2 * aL + aV + k, 10),
        'lambda lambda A|_Wperp': (-2 * aL + aP + k, 10), 'A|_V c A*|_V': (-aC + k, 1), 'A|_Wperp c A*|_Wperp': (-aC + k, 1),
        'c lambda lambda*': (-aC + k, 1), 'c c c*': (-aC + k, 1), '(A|_V)^4': (4 * aV + k, 8), '(A|_V)^2 (A|_Wperp)^2': (2 * aV + 2 * aP + k, 8),
        '(A|_Wperp)^4': (4 * aP + k, 8)}
    eqs = [e - v for (e, v) in data.values()]
    sol = sp.solve(eqs, [aV, aP, aL, aC, k], dict=True)
    quad_only = sp.solve([data[x][0] - data[x][1] for x in list(data)[:4]], [aV, aP, aL, aC], dict=True)
    return dict(
        ansatz='log eps part of I^CT[eps] = log eps u [ (Q + {I,-})( alpha_AV int (A|_V,A*|_V) + alpha_AWperp int (A|_Wperp,A*|_Wperp) '
               '+ alpha_lambda int (lambda*,lambda) + alpha_c int (c*,c) ) + kappa S ],  S = the full action, kappa S = d/dt S[(1+t)(-,-)] at t=0 times kappa',
        weights='(Q+{I,-}) int alpha (A,A*) multiplies a monomial by alpha(n_A - n_A*); int alpha(lambda*,lambda) by -alpha(n_lambda - n_lambda*); '
                'int alpha(c*,c) by -alpha(n_c - n_c*)  [bracket of classical_spec.md: {f,g} = sum_i (d_xi_i f)(d_x_i g) + (-1)^|f| (d_x_i f)(d_xi_i g), '
                'pairs (x,xi) = (A,A*), (lambda*,lambda), (c*,c)]',
        equations={kk: f'{e} = {v}' for kk, (e, v) in data.items()},
        quadratic_sector_alone=f'{quad_only} (kappa free: the quadratic sector does not fix the coupling direction)',
        solution_all_sectors=str(sol),
        conclusion='unique solution alpha_AV = alpha_AWperp = 2, alpha_lambda = -4, alpha_c = -1, kappa = 0: the 15 equations (4 quadratic, 8 cubic, '
                   '3 quartic) are consistent and force kappa = 0; the coefficient along the coupling direction vanishes at one loop.')

def t2():
    raw = json.load(open('r2_T2_raw.json')); chk = json.load(open('r2_T2_check.json'))
    ids = {
        'v1=AcA*[c] v2=AAA[A|_V] special edge A-A* (K_eps - K_0), P edge A|_V': 'QME2-a',
        'v1=AcA*[c] v2=AAA[A|_V] special edge A-A* (K_eps - K_0), P edge A|_Wperp': 'QME2-b',
        'v1=AcA*[c] v2=AcA*[A|_V] special edge A-A* (K_eps - K_0), P edge c-A*': 'QME2-c',
        'v1=c lambda lambda*[c] v2=lambda lambda A[A|_V] special edge lambda-lambda* (K_eps - K_0), P edge lambda': 'QME2-d',
        'v1=ccc*[c] v2=AcA*[A|_V] special edge c-c* (K_eps - K_0), P edge c-A*': 'QME2-e'}
    graphs = []; chains = {}
    for cl, v in raw['classes_over_i_pa'].items():
        fin = v.get('0', '0').replace('(p0**2 + p1**2 + p2**2 + p3**2)', 'p^2')
        finc = sp.nsimplify(fin.replace('*p^2', '').replace('p^2', '1')) if 'p^2' in fin else sp.Integer(0)
        g = dict(id=ids[cl], description=cl, vertices=[cl.split()[0][3:], cl.split()[1][3:]],
                 internal_lines=['special edge ' + cl.split('special edge ')[1].split(' (')[0] + ' carrying K_eps - K_0', 'P(eps,infty) edge ' + cl.split('P edge ')[1]],
                 symmetry_factor='1 (the 1/2 of 1/2{I,I} cancels against the two assignments of v1, v2)',
                 eps_inv=f'{v.get("-1", "0")} u eps^-1 int (A|_V, dc)', log_eps='0',
                 finite=f'{finc} u int (d*A|_V, d*dc)', eps_inv_sq='0')
        graphs.append(g)
        chains[f'T2 {ids[cl]}'] = dict(
            amplitude='D-graph = -(hbar/2)[d_P d_(K_eps - K_0)(I I)]_cross = -d/d delta W_1[P + delta (K_eps - K_0)] (delta odd) '
                      '= (1/4) sum_{words, position of the special edge} (-1)^(odd parameters before it) STr[...], kernel rows carry the parity operator',
            numerator='STr[F Pi d^2I(v1) N(q - P) Pi d^2I(v2)] with F = -(fibre identity) on the field-antifield pair of the special edge, '
                      'N the P(eps,infty) numerator; colour tr_g(ad ad) = tr_g(ad_X ad_Y)',
            schwinger='special edge: (e^{-eps k^2} - 1) = pinned edge (t = eps) minus contact term (K_0); P edge: int_eps^infty dt e^{-t k^2}',
            eps_expansion='pinned part: (16pi^2)^-1 sum_k 2^-k eps^{-1-k} int_1^infty du (1+u)^{-2-k} C_k(1/(1+u), u/(1+u)) e^{-eps p^2 u/(1+u)}; '
                          'contact part: -(16pi^2)^-1 sum_k 2^-k C_k eps^{-1-k}/(1+k); only integer powers of eps occur (no log eps)',
            local_part=f'eps^-1: {v.get("-1", "0")} u int (A|_V,dc); eps^0: {finc} u int (d*A|_V, d*dc)')
    graphs.append(dict(id='QME2-CT', description='counterterm vertex -Q I^CT[eps] (quadratic I^CT agreed in round 1)', vertices=['-Q I^CT[eps]'],
                       internal_lines=[], symmetry_factor='1',
                       eps_inv='-(9/4) u eps^-1 Q int (A|_V,A|_V) = +(9/2) u eps^-1 int (A|_V, dc)',
                       log_eps='0  [Q int (dA,dA) = 0; Q int (dc,A*) restricted to (c, A|_V) is -1/2 int (dc, d*dA|_V) up to sign = 0 since d*d* = 0; '
                               'coefficient of int (d*A|_V,d*A|_V) in the log part is 0 (round 1)]', finite='0'))
    tot_fin = sp.nsimplify(raw['total_over_i_pa']['0'].replace('(p0**2 + p1**2 + p2**2 + p3**2)', '1'))
    return dict(
        sector='external legs c and A|_V (ghost number 1, quadratic)',
        defect_vertex='D_eps = 1/2({I,I}_eps - {I,I}_0) + hbar Delta_eps I - Q I^CT[eps] - {I,I^CT[eps]}_eps + O(hbar^2); in this sector hbar Delta_eps I = 0 '
                      '(every self-contraction of a vertex has colour factor sum_a ([X,t_a],t^a) = 0), {I,I^CT}_eps has at least three legs; '
                      'one-loop graphs: 1/2({I,I}_eps - {I,I}_0) closed by one P(eps,infty) edge (bubble with one special edge), and -Q I^CT[eps]',
        identity_used='Q e^{hbar d_P} e^{J/hbar} = e^{hbar d_P}(Q + hbar Delta_eps) e^{J/hbar}, J = I - I^CT[eps], so Q W(P(eps,infty),J) = connected graphs with one D_eps vertex',
        other_quadratic_sectors='ghost number 1 bilinears: (c, A) and (c, lambda) only (coordinate degrees c:+1, A,lambda:0, A*,lambda*:-1, c*:-2). '
                                '(c, A|_Wperp) vanishes: no Spin(6)-invariant pairing of a scalar with a vector of W_R^perp (regulator R-covariant); '
                                '(c, lambda) vanishes: no Spin(4)-invariant bilinear of a scalar and a spinor.',
        graphs=graphs,
        eps_inv_total='special-edge graphs: -(9/2) u eps^-1 int (A|_V, dc); counterterm: +(9/2) u eps^-1 int (A|_V, dc); total 0',
        eps_inv_sq_total='0 (each graph)',
        log_total='0 (each graph: special-edge graphs contain only integer powers of eps; counterterm log part has no (c, A|_V) component)',
        finite_total=f'{tot_fin} u int (d*A|_V, d*dc)   [plane-wave coefficient of r_c: {tot_fin} i (p.a) p^2 (X_0,X_1) u, c at +p, A|_V polarisation a at -p]',
        obstruction=f'O = {tot_fin} u int (d*A|_V, d*dc) = (27/16) u Q int (d*A|_V, d*A|_V)   [Q int (d*A|_V,d*A|_V) = -2 int (d*A|_V, d*dc), from Q c = -dc in the A|_V slot]',
        q_exact='yes: O = Q( (27/16) u int (d*A|_V, d*A|_V) ); in this sector every dimension-4 local functional is a multiple of int (d*A|_V, d*dc) '
                '(int (dc, d*dA|_V) = 0), so the class is fixed by one number',
        finite_counterterm='F = (27/16) u int (d*A|_V, d*A|_V) added to I^CT[eps] (renormalised interaction I - I^CT[eps] - F); then D changes by -Q F = -O and the '
                           'quadratic obstruction vanishes. In Lemma 12.0.4 form: lift S + hbar S\' with S\' = -(27/16)(tr_g(ad_X ad_Y)/((X,Y)16pi^2)) int (d*A|_V,d*A|_V).',
        cross_check=dict(method='independent: Q W_1|_(c,A|_V) = -V^AA(i p, a), the longitudinal part of the regulated A|_V self-energy (round-1 numerators, no special edges), '
                                'exact eps -> 0 expansion via incomplete Gamma functions in the Feynman parameter x',
                         result=chk['QW1 = -V_AA(ip, a): prediction for the D-graph sum'],
                         per_graph=dict((k, v) for k, v in chk.items() if k not in ('V_AA_longitudinal_total', 'QW1 = -V_AA(ip, a): prediction for the D-graph sum')),
                         agreement='exact: eps^-1 -9/2 and finite -27/8 from both routes; the gauge-loop and ghost-loop log(p^2 eps) and log 2 pieces cancel in the sum'),
        interpretation='the renormalised Feynman-gauge A|_V self-energy at scale infinity is transverse up to the local term (27/16) u int (d*A|_V,d*A|_V); '
                       'the QME obstruction in the (c,A|_V) sector is its Q-variation',
        cubic_sector='open (not computed): requires triangles and quartic-vertex bubbles with one special edge, the vertex {I, I^CT[eps]}_eps, and the '
                     'one-particle-reducible graphs -Q I^CT[eps] --P-- I. Structural statement only: the one-loop obstruction is closed under Q + {I,-}, so its cubic '
                     'part equals (27/16) u {I, int (d*A|_V,d*A|_V)} (cubic part) plus a (Q + {I,-})-closed cubic remainder.',
        verdict='eps^-1 cancels against -Q of (9/4) u eps^-1 int (A|_V,A|_V); no log eps; finite one-loop obstruction in the (c,A|_V) sector is nonzero and Q-exact, '
                'removed by the finite counterterm (27/16) u int (d*A|_V,d*A|_V).'), chains

def t3():
    sw3 = {'traceless symmetric scalar bilinear (iota_w1 A, iota_w2 A)': 0, 'trace scalar bilinear sum_w (iota_w A, iota_w A)': -24,
           'letter iota_w A': -4, 'letter A|_V (transverse part)': -4, 'letter A|_V (longitudinal part)': 0, 'letter lambda': -8,
           'letter c (gauge ghost, computed in round 2)': -2}
    rows = {k: dict(SW3_Delta_equals_minus_E=v, SW17_Delta_equals_2E=-2 * v) for k, v in sw3.items()}
    return dict(
        normalisations='SW (3): Delta = -E; SW (17): Delta = 2E = -2 x (SW (3) Delta); tabulated: the coefficient of u O in the anomalous part of '
                       '[(d_CE + Q) O[infty]](Delta), same sign convention as round 1 (DEFINITIONS_v2), units u',
        table=rows,
        konishi_check='trace: -24 u (SW (3)) = 3 g^2 N/(4 pi^2) for su(N) with u = -g^2 N/(32 pi^2); in SW (17): +48 u = -3 g^2 N/(2 pi^2) (the eigenvalue of Delta = 2E is -2 times)',
        ghost_letter_computation=dict(
            O='int omega (X, c), X a fixed element of g^vee (Q-closed: Q has no component into the c slot)',
            route='round-1 route (i): anomalous part = -2 hbar sum_k b_k d^cross_{P(0,infty)}(B_k, O) - Q O_1, O_1 = -2 hbar d^cross_{P(0,infty)}(R, O), '
                  'sum_k b_k B_k = Q R, R = 2 int (A,A*) - 4 int (lambda*,lambda) - int (c*,c)',
            step1='P(0,infty) joins c only to A*|_V: component A* (x) c = -(-,-)^-1 d_x (4 pi^2 |x-y|^2)^-1 (momentum space -i k/k^2); only B = int (dc, A*) '
                  '(log coefficient b = -1) contributes: d^cross_P(int (dc,A*), O) = -int (dc, X d (d*d)^-1 omega) = -O, so -2 hbar b d^cross = -2 u O',
            step2='d^cross_P(R, O) = 2 int (A, -X d (d*d)^-1 omega) => O_1 = 4 u int (A|_V, X d (d*d)^-1 omega); Q O_1 = 4 u int (-dc, X d (d*d)^-1 omega) = -4 u O; -Q O_1 = +4 u O',
            sum='-2 + 4 = +2 in the orientation of the identity anomalous part = 2 hbar {R, O}_0 (from [d_P, Q] = d_{QP}, Q P(0,infty) = K_0); '
                'the same computation gives +4 for iota_w A and +8 for lambda, which round 1 tabulates as -4 and -8; in that convention the ghost letter is -2 '
                '(SW (3)) and +4 (SW (17))',
            general_rule='Q-closed letter of field phi: coefficient = -2 x (weight of phi in {Psi,-}), weights A: 2, lambda: 4, c: 1'))

if __name__ == '__main__':
    T1, ch1 = t1()
    T2, ch2 = t2()
    res = dict(
        engine='B (round 2): exact Gaussian-integer arithmetic; method B (eps d/deps via pinned edges, Taylor expansion in the external momenta, '
               'exact radial integral and exact angular average on the 24 short roots of F4)',
        units=U, legend=LEGEND,
        T1=dict(external_sets=T1, coupling_direction=coupling_direction(),
                summary='every external set: computed total log eps coefficient = (A)-prediction x tree vertex; no eps^-1, no eps^-1/2 pieces; '
                        'non-Lie colour parts (symmetrised adjoint traces) and parity-odd parts cancel in each total',
                log_part_of_ICT_cubic_quartic='u log eps [6 (A|_V)^3 + 6 A|_V(A|_Wperp)^2 parts of int(-1/2 (dA, 1/2[A^A])) + 10 int 1/2(lambda, (dslash_A - dslash)lambda) '
                                              '+ int ([A,c],A*) + int ([c,lambda],lambda*) + int 1/2([c,c],c*) + 8 int(-1/4(1/2[A^A],1/2[A^A]))] '
                                              '= u log eps {I, 2 int (A,A*) - 4 int (lambda*,lambda) - int (c*,c)}'),
        T2=T2, T3=t3(), formula_chains={**ch1, **ch2},
        checks=['24 short roots of F4 reproduce all S^3 monomial averages up to degree 5 exactly',
                'round-1 omega=1 two-point log coefficients reproduced by the round-2 engine: lambda-lambda 4, A*-c -1',
                'colour reduction tr_g(M_n...M_1) -> {([X,Y],Z), ([.,.],[.,.]), symmetrised trace} checked against explicit su(3) traces '
                '(non-orthonormal basis, rescaled pairing): 20/20, Jacobi ([X0,X2],[X1,X3]) = ([X0,X1],[X2,X3]) + ([X0,X3],[X1,X2])',
                '0-dim Grassmann toy: word expansion with stripped odd parameters exact for triangles and quartic-vertex bubbles with odd legs; '
                'transposed Hessian order fails at n = 3 (passes at n <= 2); odd special edge requires the parity operator on the kernel row index',
                'T1 fits: overdetermined exact linear systems on 3-11 random integer configurations per set, consistent (zero residual)',
                'T1 per-graph sanity: (A|_V)^3 gauge+ghost+swordfish -4/3, scalars 2, fermions 16/3 (x tree) reproduce the textbook Feynman-gauge split of the '
                'three-gluon vertex renormalisation in the ratio -2',
                'T2: special-edge D-graphs equal -V^AA(ip, a) (independent longitudinal self-energy) exactly in eps^-1 and eps^0'],
        files=dict(engine='engine_b/r2_engine.py', T1='engine_b/r2_T1.py (+ r2_T1_*.json)', T2='engine_b/r2_T2.py, r2_T2_check.py (+ r2_T2_raw.json, r2_T2_check.json)',
                   toys='engine_b/toy_r2.py, toy_r2b.py, toy_r2c.py', colour='engine_b/r2_colour_check.py', regression='engine_b/r2_regress.py'))
    json.dump(res, open('results_r2.json', 'w'), indent=1)
    print(json.dumps(res['T1']['coupling_direction'], indent=1))
    for k, v in T1.items():
        print(k, '|', v['total'], '| predicted', v['predicted_from_A'].split('  ')[0], '|', v['verdict'])
    print(T2['finite_total']); print(T2['obstruction'])
