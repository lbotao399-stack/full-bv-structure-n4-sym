"""Engine A, part Q, task TQ0 (gate G-TC0): the quadratic Q-closed directions along Sigma (Prop prop:qs-freexi (i)).
The five covariant structures, linear in xi in S_+ (x) W, degree -1, Weyl weight 9/2: T1..T4 of Table tab:qs-terms and
T5' = int ((d^* A|_V) xi, lam^*).  A quadratic functional int (X e, e^*) is Q-closed iff [X, Q] = 0 (graded commutator of
the two odd operators); X is read off the block (qa_e0.block_operator).  Exact linear algebra over Q(i) at several
momenta; the kernel is compared with span{S^V_2 = T1 + T3, S^perp_2 = T2 + T4}.  Exact directions: a Q-exact functional
would be Q of a degree -2 quadratic functional of weight 9/2; the degree -2 quadratic monomials (c^*, f), (A^*, A^*),
(A^*, lam^*), (lam^*, lam^*) have weights >= 4 + 1/2 only for (c^*, f) with f of weight 1/2, which does not exist."""
import numpy as np
import sympy as sp
import qa_s1 as S1, qa_e0 as E0


def t5prime(xi):
    xi = np.asarray(xi, complex)
    def fn(pols, moms):
        A, ls = pols
        k = moms[0]
        dstar = np.zeros(10, complex); dstar[:4] = -1j * np.asarray(k[:4])         # d^* A_V = -d.A -> -i k.A
        return 2 * np.einsum('M,iMe,a,jae->ij', dstar, A, xi, ls, optimize=True)
    return {'T5p': (['A', 'ls'], fn)}


def run(xi, log=print):
    blocks = [S1.quad_block(S1.odd_terms(xi, (t,))) for t in ('T1', 'T2', 'T3', 'T4')] + [S1.quad_block(t5prime(xi))]
    ps = [np.array(p, dtype=complex) for p in ([1, -2, 3, 1], [0, 1, 2, -1], [2, 1, 0, 3])]
    rows = []
    for p in ps:
        Qm = E0.op_Q(p)
        cols = []
        for b in blocks:
            X = E0.block_operator(b, p)
            Cm = X @ Qm + Qm @ X
            cols.append(Cm.ravel())
        rows.append(np.column_stack(cols))
    A = np.vstack(rows)
    As = sp.Matrix(A.shape[0], A.shape[1], lambda i, j: sp.nsimplify(round(A[i, j].real, 9)) + sp.I * sp.nsimplify(round(A[i, j].imag, 9)))
    ker = As.nullspace()
    kv = [[str(sp.nsimplify(x)) for x in (v / [y for y in v if y != 0][0])] for v in ker]
    target = sp.Matrix([[1, 0, 1, 0, 0], [0, 1, 0, 1, 0]]).T
    same = (sp.Matrix.hstack(*ker).rank() == 2) and (sp.Matrix.hstack(*(list(ker) + [target[:, 0], target[:, 1]])).rank() == 2)
    out = {'structures': ['T1', 'T2', 'T3', 'T4', "T5'=((d^*A_V) xi, lam^*)"], 'kernel_dim': len(ker), 'kernel_basis': kv,
           'kernel = span{T1+T3, T2+T4}': bool(same), 'exact_directions': 0}
    log('  TQ0: kernel of [., Q] on the five structures: dim %d, basis %s, equals span{S^V_2, S^perp_2}: %s' % (len(ker), kv, same))
    return out
