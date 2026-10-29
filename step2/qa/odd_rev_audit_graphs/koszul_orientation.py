# Audit (odd_rev_audit_graphs), script 2: relative Koszul sign of the two orientations of a two-insertion chain,
# and the consequence for the leading large-distance term.
#
# Parities of the functional coordinates (degree + fermion number mod 2; ch01/ch03 conventions, DESIGN.md degrees):
#   c: odd, A (A_V, A_W): even, lambda: odd, A*: odd, lambda*: even, c*: even.
# A complete Wick contraction of an ordered product of legs is evaluated by moving the two legs of every pair next to
# each other in the order (first, second) written in the pair; the sign is the Koszul sign of that permutation.
# Chain types (the two insertions u1, u2 and the two ghost legs c_z, c_z' of the rest of the cycle):
#   odd-odd  : u1 = (lambdabar_1, A*_1), u2 = (lambda_2, A*_2): pairs (lambdabar_1,lambda_2), (A*_1, c_.), (A*_2, c_.)
#   even-even: u1 = (A_1, A*_1),        u2 = (A_2, A*_2):       pairs (A_1, A_2),             (A*_1, c_.), (A*_2, c_.)
# pattern 1: (A*_1,c_z)(A*_2,c_z');  pattern 2: (A*_1,c_z')(A*_2,c_z).  Every pair is written (A*, c), so the same
# propagator component appears in both patterns.
import json, os, itertools
import numpy as np

def koszul_sign(order, parity, pairs):
    """order: list of leg names in the product; pairs: list of (a,b). Returns the sign of bringing each pair
    adjacent as (a,b), pairs in the listed order, with Koszul signs for odd legs."""
    target = [x for p in pairs for x in p]
    seq = list(order); sign = 1
    # bubble the sequence into target order, counting swaps of two odd legs
    for i, t in enumerate(target):
        j = seq.index(t)
        while j > i:
            a, b = seq[j - 1], seq[j]
            if parity[a] and parity[b]:
                sign = -sign
            seq[j - 1], seq[j] = b, a
            j -= 1
    return sign

out = {}
P = {'lamb1': 1, 'lam2': 1, 'As1': 1, 'As2': 1, 'cz': 1, 'czp': 1, 'A1': 0, 'A2': 0}
# product order: S^(1)(u1) S^(1)(u2) R, each insertion written (field leg, antifield leg), R = c_z c_z'
oo = ['lamb1', 'As1', 'lam2', 'As2', 'cz', 'czp']
ee = ['A1', 'As1', 'A2', 'As2', 'cz', 'czp']
for nm, order, first in (('odd-odd', oo, ('lamb1', 'lam2')), ('even-even', ee, ('A1', 'A2'))):
    s1 = koszul_sign(order, P, [first, ('As1', 'cz'), ('As2', 'czp')])
    s2 = koszul_sign(order, P, [first, ('As1', 'czp'), ('As2', 'cz')])
    # robustness: other orderings of the factors in the product give the same relative sign
    rel = set()
    for perm in itertools.permutations([0, 1, 2]):
        blocks = [order[0:2], order[2:4], order[4:6]]
        o2 = [x for k in perm for x in blocks[k]]
        a = koszul_sign(o2, P, [first, ('As1', 'cz'), ('As2', 'czp')])
        b = koszul_sign(o2, P, [first, ('As1', 'czp'), ('As2', 'cz')])
        rel.add(a * b)
    out[nm] = {'sign_pattern_1': s1, 'sign_pattern_2': s2, 'relative_sign': s1 * s2,
               'relative_sign_over_all_block_orderings': sorted(rel)}

# Consequence: total = sum_{z,z'} r(z,z') [C(z,z') - C(z',z)] with C the chain kernel (u1 at the z end).
# Model the large-distance expansion C(z,z') = C0 + C1.z + C2.z' + O(|y|^-1 more); r an arbitrary rest kernel on a grid.
rng = np.random.default_rng(29)
Z = rng.normal(size=(7, 4))                 # positions of the ghost ends of the rest
r = rng.normal(size=(7, 7))                 # arbitrary (not symmetrised) rest kernel r(z,z')
C0, C1, C2 = rng.normal(), rng.normal(size=4), rng.normal(size=4)
def Cker(z, zp, order):
    return C0 + (C1 @ z + C2 @ zp if order >= 1 else 0.0)
tot0 = sum(r[i, j] * (Cker(Z[i], Z[j], 0) - Cker(Z[j], Z[i], 0)) for i in range(7) for j in range(7))
tot1 = sum(r[i, j] * (Cker(Z[i], Z[j], 1) - Cker(Z[j], Z[i], 1)) for i in range(7) for j in range(7))
mom = sum((r[i, j] - r[j, i]) * Z[i] for i in range(7) for j in range(7))
out['leading_term_after_orientation_sum'] = float(tot0)
out['first_moment_term'] = float(tot1)
out['first_moment_prediction_(C1-C2).antisym_moment'] = float((C1 - C2) @ mom)
out['statement'] = ('relative orientation sign -1 for both chain types (odd-odd through lambda lambdabar, even-even through '
                    'A_V A_V); the z-independent leading term C0 cancels identically, the first surviving term is '
                    '(C1 - C2) . sum (r(z,z\') - r(z\',z)) z, one power of |y| lower: the orientation sum lowers the IR degree '
                    'by one for every chain with two ghost ends, independently of the generators.')
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'koszul_orientation.json'), 'w'), indent=1)
