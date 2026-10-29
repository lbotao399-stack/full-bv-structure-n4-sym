"""QA: recompute the P7 targets of engine B's two-site families independently: -1/3 x (derivative map) in B's cyclic
bases (graded cyclicity for fermion letters) and -2/3 N on the diagonal; compare with B's sumC matrices."""
import json, os, re, itertools
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
B = json.load(open(os.path.join(HERE, '..', '..', '..', 'engine_b', 'results.json')))
out = {}
def canon(k, l, X, Y, ferm):
    """tr(d^k X d^l Y) -> (sign, key) in a canonical cyclic form"""
    if X == Y:
        if k <= l:
            return 1, (k, l)
        return (-1 if ferm else 1), (l, k)
    return 1, (k, l)
def dmap(basis_in, basis_out, X, Y, ferm):
    idx = {b: i for i, b in enumerate(basis_out)}
    Mx = sp.zeros(len(basis_out), len(basis_in))
    for j, (k, l) in enumerate(basis_in):
        for (kk, ll) in ((k + 1, l), (k, l + 1)):
            if X == Y and kk == ll and ferm:
                continue
            s, key = canon(kk, ll, X, Y, ferm)
            Mx[idx[key], j] += s
    return Mx
def parse_pairs(basis):
    res = []
    for b in basis:
        nums = [int(x) for x in re.findall(r"(\d+)", b)]
        if 'd^' in b:
            res.append((nums[0], nums[-1]) if len(nums) == 2 else (nums[0], nums[2]))
        else:
            res.append((nums[0], nums[1]))
    return res
fams = {'F1_Z_sl2': ('Z', 'Z', False, 4, 'P7b_sumC_Pb', 'P7c_sumC_D0'), 'F2_psi_sl2': ('psi', 'psi', True, 4, 'P7b_sumC_Pb', 'P7c_sumC_D0'),
        'F3_gluon': ('F', 'F', False, 6, 'sumC_Pb_in_basis', 'sumC_D0_matrix'), 'F3_mixed': ('F', 'psi', False, 5, 'sumC_Pb_in_basis', 'sumC_D0_columns')}
for fam, (X, Y, ferm, N0, kP, kD) in fams.items():
    d = B[fam]; res = {}
    Ss = sorted(int(s) for s in d if s.isdigit())
    for S in Ss:
        e = d[str(S)]
        if e.get('dim', 0) == 0 or kP not in e:
            continue
        bin_ = [tuple(map(int, re.findall(r"d\^(\d+)", b))) if 'd^' in b else tuple(int(x) for x in re.findall(r"'X', (\d+)", b)) for b in e['basis']]
        nxt = d.get(str(S + 1))
        if nxt and 'basis' in nxt:
            bout = [tuple(map(int, re.findall(r"d\^(\d+)", b))) if 'd^' in b else tuple(int(x) for x in re.findall(r"'X', (\d+)", b)) for b in nxt['basis']]
        else:
            # construct the S+1 cyclic basis
            if X == Y:
                bout = [(k, S + 1 - k) for k in range(0, (S + 1) // 2 + 1) if not (ferm and 2 * k == S + 1)]
            else:
                bout = [(k, S + 1 - k) for k in range(S + 2)]
        D = dmap(bin_, bout, X, Y, ferm)
        tgtP = -sp.Rational(1, 3) * D
        mP = sp.Matrix(sp.sympify(e[kP]))
        mD = sp.Matrix(sp.sympify(e[kD]))
        tgtD = -sp.Rational(2, 3) * (N0 + S) * sp.eye(mD.shape[0])
        if kP == 'sumC_Pb_in_basis':
            mP = mP.T          # B stores rows = inputs for these families
        res[S] = {'P7b_B_equals_QA_target': mP.shape == tgtP.shape and sp.simplify(mP - tgtP) == sp.zeros(*tgtP.shape),
                  'P7c_B_equals_QA_target': sp.simplify(mD - tgtD) == sp.zeros(*tgtD.shape),
                  'QA_target_P': str(tgtP.tolist()), 'B_P': e[kP]}
    out[fam] = res
    print(fam, {S: (v['P7b_B_equals_QA_target'], v['P7c_B_equals_QA_target']) for S, v in res.items()})
json.dump(out, open(os.path.join(HERE, '..', 'qa_targets_B.json'), 'w'), indent=1)
