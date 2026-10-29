"""tree(Ob_2, I_3; P(0,oo)) on plane waves, own edge-rule contraction (rb_loop vertices/kernels), with
Ob_2 = -27/8 unit int(d^*dc, d^*A|_V).  Checks locality (1/|k|^2 cancels) and the claimed values
27/8 (G1 - G2 + G3) - 27/4 (G4 + G5) in all four sectors (cross term of d_P(Ob_2 I_3), weight 1).
"""
import sys, os, json, itertools, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
import rb_func as RF
import rb_loop as RL
from rb_alg import PARITY

rng = np.random.default_rng(321)


def rint(s):
    return rng.integers(-3, 4, size=s).astype(float)


def T_E2(ks):
    """int (d^*d c, d^*A|_V), slots (c, A): |k_c|^2 * (-i) (k_A)_M delta_ab, M in V."""
    kc, kA = ks
    T = np.zeros((3, 10, 3), dtype=complex)
    for M in range(4):
        T[:, M, :] = (kc @ kc) * (-1j) * kA[M] * np.eye(3)
    return T.reshape(3, 30)


def tree_value(legs):
    """legs: list of (slot, vec, p) in leg order (3 legs)."""
    total = 0j
    n = len(legs)
    for j_ob in range(n):                       # which leg sits on Ob_2
        others = [j for j in range(n) if j != j_ob]
        for ob_ext in (0, 1):                   # Ob_2 slot (c=0, A=1) carrying the external leg
            ob_slots = ['c', 'A']
            if ob_slots[ob_ext] != legs[j_ob][0]:
                continue
            ob_int = 1 - ob_ext
            for vt, (slots, Tf) in RL.VERTS.items():
                if len(slots) != 3:
                    continue
                for perm in itertools.permutations(range(3)):
                    e1, e2, it = perm
                    if slots[e1] != legs[others[0]][0] or slots[e2] != legs[others[1]][0]:
                        continue
                    if (ob_slots[ob_int], slots[it]) not in RL.PSET:
                        continue
                    # avoid double counting identical slot-type patterns of the same vertex
                    # (full multilinear vertex functions: one representative per role pattern)
                    key = (vt, j_ob, ob_ext, slots[e1], slots[e2], slots[it])
                    if key in seen:
                        continue
                    seen.add(key)
                    pob = np.array(legs[j_ob][2], dtype=float)
                    k = -pob                        # enters Ob_2 through its internal slot
                    moms_ob = [None, None]; moms_ob[ob_ext] = pob; moms_ob[ob_int] = k
                    Tob = T_E2(moms_ob)
                    Tob = np.transpose(Tob, [ob_ext, ob_int])
                    Mob = np.tensordot(np.asarray(legs[j_ob][1]).reshape(-1), Tob, axes=(0, 0))   # [int]
                    moms_v = [None] * 3
                    moms_v[e1] = np.array(legs[others[0]][2], dtype=float)
                    moms_v[e2] = np.array(legs[others[1]][2], dtype=float)
                    moms_v[it] = -k
                    Tv = np.transpose(Tf(moms_v), [e1, e2, it])
                    Tv = np.tensordot(np.asarray(legs[others[0]][1]).reshape(-1), Tv, axes=(0, 0))
                    Tv = np.tensordot(np.asarray(legs[others[1]][1]).reshape(-1), Tv, axes=(0, 0))   # [int]
                    K = RL.kernel('P', ob_slots[ob_int], slots[it], k, 'V' if ob_slots[ob_int] == 'A' else None)
                    if ob_slots[ob_int] == 'A':
                        K = RL.kernel('P', 'A', 'A', k, 'V') + RL.kernel('P', 'A', 'A', k, 'W')
                    val = Mob @ K @ Tv / (k @ k)
                    # Koszul: sequence [Ob.s0, Ob.s1, V.s0, V.s1, V.s2] -> [pair(Ob.int, V.int), ext in leg order]
                    seq = [('o', 0), ('o', 1), ('v', 0), ('v', 1), ('v', 2)]
                    par = {('o', 0): PARITY['c'], ('o', 1): PARITY['A']}
                    for s in range(3):
                        par[('v', s)] = PARITY[slots[s]]
                    pos = {j_ob: ('o', ob_ext), others[0]: ('v', e1), others[1]: ('v', e2)}
                    target = [('o', ob_int), ('v', it)] + [pos[j] for j in range(n)]
                    total += RL.koszul(seq, target, par) * val
    return total


SECT = {'cAVAV': (['c', 'A', 'A'], 'V', ['G1', 'G2']), 'cAWAW': (['c', 'A', 'A'], 'W', ['G3']),
        'cll': (['c', 'lam', 'lam'], None, ['G5']), 'ccAs': (['c', 'c', 'Ast'], 'V', ['G4'])}
GF = {'G1': RF.F_G1, 'G2': RF.F_G2, 'G3': RF.F_G3, 'G4': RF.F_G4, 'G5': RF.F_G5}
out = {}
for sec, (sl, part, Gn) in SECT.items():
    rows, vals = [], []
    for _ in range(5):
        while True:
            p1, p2 = rint(4), rint(4); p0 = -p1 - p2
            if p0.any() and p1.any() and p2.any() and (p1 + p2).any() and (p0 + p1).any() and (p0 + p2).any():
                break
        legs = []
        for s, p in zip(sl, [p0, p1, p2]):
            v = rint(RF.SHAPE[s])
            if s in ('A', 'Ast'):
                if part == 'V':
                    v[4:] = 0
                else:
                    v[:4] = 0
            legs.append((s, v, p))
        seen = set()
        vals.append(-27.0 / 8.0 * tree_value(legs))
        rows.append([RF.value(GF[g], legs) for g in Gn])
    A = np.array(rows); b = np.array(vals)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    out[sec] = {'coef_on_' + ','.join(Gn): [str(Fr(x.real).limit_denominator(1000)) for x in sol],
                'max_imag': float(np.max(np.abs(sol.imag))), 'resid': float(np.max(np.abs(A @ sol - b))),
                'max_value': float(np.max(np.abs(b)))}
    print(sec, out[sec], flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rb_tree_check.json'), 'w'), indent=1)
