"""
exhaustive_alg.py -- exact basis-level verification of the psl(4|4) realization of algebra.py / structure.py.

For the 62-element basis e_1..e_62 of structure.basis() (30 even | 32 odd, grouped by the [SW (8),(17)] names):
  1. structure constants: [e_i, e_j] = sum_m C_ij^m e_m, computed from the realized bracket for ALL 62^2 ordered
     pairs, with exact reconstruction (closure of the span) asserted for every pair;
  2. graded antisymmetry  C_ij = -(-1)^{|i||j|} C_ji   for all ordered pairs;
  3. super-Jacobi on ALL 62^3 ordered basis triples, in the equivalent matrix form
         ad_{e_i} ad_{e_j} - (-1)^{|i||j|} ad_{e_j} ad_{e_i} = ad_{[e_i,e_j]}     for all (i,j)
     (the (n,k) entry of this identity is the Jacobiator on (e_i, e_j, e_k));
  4. Killing form  str(ad_{e_i} ad_{e_j}) = 0  for all 62^2 ordered pairs;
  5. [SW (17)] supermatrix block rule for all 62^2 ordered pairs, in a strict form: the diagonal block (3,3)
     admits only the R-symmetry subspace (representatives normalised with traceless blocks);
  6. ad_Delta eigenvalues on every basis element;
  7. parity bookkeeping: the bracket is even (C_ij^m = 0 unless |m| = |i| + |j| mod 2).
Arithmetic: sympy DomainMatrix over QQ_I (exact Gaussian rationals).
"""
import os, sys, json, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy import QQ_I
import algebra as Al
import structure as St

Al.set_cR(sp.Integer(1))


def strict_allowed(n1, n2):
    """block rule with representatives normalised so that every diagonal block of the (4|4) supermatrix
    is traceless except for the Delta direction diag(1_2,-1_2|0): target (1,1) -> {M, Delta},
    (2,2) -> {Mbar, Delta}, (3,3) -> {R}."""
    DIAG = {(1, 1): ['M', 'Delta'], (2, 2): ['Mbar', 'Delta'], (3, 3): ['R']}
    out = set()
    for (a, b) in St.BLOCK[n1]:
        for (c, d) in St.BLOCK[n2]:
            tgt = []
            if b == c:
                tgt.append((a, d))
            if d == a:
                tgt.append((c, b))
            for t in tgt:
                if t in DIAG:
                    out.update(DIAG[t])
                else:
                    out.add(St.OFFD[t])
    return out


def main(cR=sp.Integer(1), tag=''):
    Al.set_cR(cR)
    T0 = time.time()
    out = {'coefficient_of_so(W_perp)_part_in_odd_odd_bracket': str(cR)}
    B = St.basis()
    C = St.Coords(B)
    N = C.ne + len(C.odd)
    par = [C.parity(k) for k in range(N)]
    names = C.names
    out['basis'] = {'dimension_even': C.ne, 'dimension_odd': len(C.odd),
                    'subspaces': {k: len(v) for k, v in B.items()}}
    # ---- 1. structure constants for all ordered pairs (coords() asserts exact reconstruction) ----
    t = time.time()
    CC = [[None] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            br = Al.bracket(C.element(i), C.element(j))
            CC[i][j] = [sp.nsimplify(sp.expand(z)) for z in C.coords(br)]
    out['structure_constants'] = {'ordered_pairs': N * N, 'closure_exact_reconstruction': True,
                                  'seconds': round(time.time() - t, 1)}
    # ---- 7. parity ----
    bad_par = [(i, j) for i in range(N) for j in range(N)
               for m in range(N) if CC[i][j][m] != 0 and par[m] != (par[i] + par[j]) % 2]
    out['bracket_is_even'] = {'violations': len(bad_par)}
    # ---- 2. graded antisymmetry ----
    bad_as = []
    for i in range(N):
        for j in range(N):
            s = -(-1) ** (par[i] * par[j])
            if any(sp.expand(CC[i][j][m] - s * CC[j][i][m]) != 0 for m in range(N)):
                bad_as.append((i, j))
    out['graded_antisymmetry'] = {'ordered_pairs_checked': N * N, 'violations': len(bad_as)}
    # ---- ad matrices: (ad_i)_{m,k} = C_ik^m ----
    K = QQ_I

    def conv(z):
        return K.from_sympy(sp.nsimplify(z))
    AD = []
    for i in range(N):
        rows = [[conv(CC[i][k][m]) for k in range(N)] for m in range(N)]
        AD.append(DomainMatrix(rows, (N, N), K).to_sparse())
    # ---- 3. super-Jacobi on all ordered triples ----
    t = time.time()
    bad_j = []
    nz_comm = 0
    for i in range(N):
        for j in range(N):
            s = (-1) ** (par[i] * par[j])
            lhs = AD[i] * AD[j] - (AD[j] * AD[i] if s == 1 else -(AD[j] * AD[i]))
            rhs = DomainMatrix.zeros((N, N), K).to_sparse()
            for m in range(N):
                cm = CC[i][j][m]
                if cm != 0:
                    rhs = rhs + AD[m] * conv(cm)
            D = lhs - rhs
            if not D.is_zero_matrix:
                bad_j.append((i, j))
            if not rhs.is_zero_matrix:
                nz_comm += 1
    out['super_Jacobi'] = {'ordered_triples_checked': N ** 3, 'violations': len(bad_j),
                           'pairs_with_nonzero_bracket': nz_comm, 'seconds': round(time.time() - t, 1)}
    # ---- 4. Killing form ----
    t = time.time()
    sgn = [K.one if p == 0 else -K.one for p in par]
    bad_k = []
    for i in range(N):
        for j in range(N):
            P = AD[i] * AD[j]
            tr = K.zero
            Pd = P.rep.to_sdm() if hasattr(P.rep, 'to_sdm') else P.rep
            for k, row in Pd.items():
                if k in row:
                    tr += sgn[k] * row[k]
            if tr != K.zero:
                bad_k.append((i, j, str(K.to_sympy(tr))))
    out['Killing_form_str_ad_ad'] = {'ordered_pairs_checked': N * N, 'nonzero_entries': len(bad_k),
                                     'seconds': round(time.time() - t, 1)}
    # also: the even part restricted to itself carries nondegenerate Killing forms of its two sl(4) factors
    ev = list(range(C.ne))
    conf = [k for k in ev if names[k] in ('P', 'K', 'M', 'Mbar', 'Delta')]
    rsym = [k for k in ev if names[k] == 'R']

    def killing_sub(idx):
        """Killing form of the subalgebra spanned by idx (assumed closed), as a matrix."""
        sub = {k: n for n, k in enumerate(idx)}
        M = sp.zeros(len(idx), len(idx))
        for a, i in enumerate(idx):
            for b, j in enumerate(idx):
                tr = 0
                for k in idx:
                    for l in idx:
                        tr += CC[i][l][k] * CC[j][k][l]
                M[a, b] = tr
        return M
    clos = {}
    for nm, idx in (('conformal', conf), ('R', rsym)):
        S = set(idx)
        clos[nm] = all(CC[i][j][m] == 0 for i in idx for j in idx for m in range(N) if m not in S)
    KC, KR = killing_sub(conf), killing_sub(rsym)
    out['even_part'] = {'conformal_span_closed': clos['conformal'], 'R_span_closed': clos['R'],
                        'conformal_R_commute': all(CC[i][j][m] == 0 for i in conf for j in rsym for m in range(N)),
                        'dim_conformal': len(conf), 'dim_R': len(rsym),
                        'rank_Killing_conformal': KC.rank(), 'rank_Killing_R': KR.rank()}
    # ---- 5. block rule ----
    bad_b, bad_b_loose = [], []
    for i in range(N):
        for j in range(N):
            sup = set(names[m] for m in range(N) if CC[i][j][m] != 0)
            al = strict_allowed(names[i], names[j]) | strict_allowed(names[j], names[i])
            if not sup <= al:
                bad_b.append((i, j, names[i], names[j], sorted(sup)))
    out['SW17_block_rule_strict'] = {'ordered_pairs_checked': N * N, 'violations': len(bad_b),
                                     'examples': bad_b[:5]}
    # ---- 6. ad_Delta eigenvalues ----
    dI = [k for k in range(N) if names[k] == 'Delta'][0]
    w = {}
    eig_ok = True
    for k in range(N):
        col = [CC[dI][k][m] for m in range(N)]
        others = [m for m in range(N) if m != k and col[m] != 0]
        if others:
            eig_ok = False
        w.setdefault(names[k], set()).add(str(col[k]))
    out['ad_Delta_eigenvalues_on_basis'] = {'every_basis_element_is_an_eigenvector': eig_ok,
                                            'eigenvalues': {k: sorted(v) for k, v in w.items()}}
    # ---- derived: [g,g] = g, centre ----
    span = DomainMatrix([[conv(CC[i][j][m]) for m in range(N)] for i in range(N) for j in range(N)], (N * N, N), K)
    out['derived_algebra_dimension'] = int(span.rank())
    cen = DomainMatrix([[conv(CC[i][j][m]) for i in range(N)] for j in range(N) for m in range(N)], (N * N, N), K)
    out['centre_dimension'] = N - int(cen.rank())
    out['all_checks_pass'] = (not bad_par and not bad_as and not bad_j and not bad_k and not bad_b and eig_ok)
    out['_seconds_total'] = round(time.time() - T0, 1)
    # store structure constants (nonzero entries, as strings) for downstream reuse
    sc = {'names': names, 'parity': par,
          'C': [[i, j, m, str(CC[i][j][m])] for i in range(N) for j in range(N) for m in range(N) if CC[i][j][m] != 0]}
    if not tag:
        json.dump(sc, open(os.path.join(HERE, 'out', 'psl44_structure_constants.json'), 'w'))
    json.dump(out, open(os.path.join(HERE, 'out', 'exhaustive_alg%s.json' % tag), 'w'), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--controls':
        # negative controls: a wrong so(W_perp) coefficient in [odd, odd] must violate super-Jacobi
        main(sp.Rational(1, 2), '_control_cR_half')
        main(sp.Integer(0), '_control_cR_zero')
    else:
        main()
