"""Route A: the S^(2)(Q1,Q2) vertex block for engine A's cycle code (qa_s1 conventions: canonical vertex functions
scaled x2, colour-diagonal; for the two identical lam* slots the canonical function is the full symmetrised vertex
function), Definition c2-S2 in engine A's realisation."""
import numpy as np
import ca_common as A
import ca_tree as T


def s2_terms(Q1, Q2, parts=('ll', 'Ac')):
    G12 = T.gam(Q1, Q2)
    Gs = sum(G12[M] * A.SS[M].T for M in range(10))       # Gamma(w1,w2).G12 = w1^T Gs w2
    terms = {}
    if 'll' in parts:
        def fll(pols, moms):
            w1, w2 = pols
            a = np.einsum('iAe,AB,jBe->ij', w1, Gs, w2, optimize=True)            # (G12, Gamma(w1,w2)) (symmetric)
            q1w1 = np.einsum('A,iAe->ie', Q1, w1); q2w2 = np.einsum('A,jAe->je', Q2, w2)
            q1w2 = np.einsum('A,jAe->je', Q1, w2); q2w1 = np.einsum('A,iAe->ie', Q2, w1)
            b = np.einsum('ie,je->ij', q1w1, q2w2) + np.einsum('ie,je->ij', q2w1, q1w2)
            # vertex function = 2 * [ 1/4 * 2 a - 1/2 b ] (x2 scaling)
            return a - b
        terms['S2ll'] = (['ls', 'ls'], fll)
    if 'Ac' in parts:
        def fAc(pols, moms):
            Ap, cs = pols
            return -2 * np.einsum('M,iMe,je->ij', G12, Ap, cs, optimize=True)
        terms['S2Ac'] = (['A', 'cs'], fAc)
    return terms


def gen_s2(Q1, Q2, name='S2'):
    return A.SA.gen_from_terms(name, s2_terms(Q1, Q2), 0)
