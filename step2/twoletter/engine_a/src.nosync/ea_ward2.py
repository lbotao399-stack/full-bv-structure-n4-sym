"""Engine A, gate G3 continued: Ward identity W1 (Prop prop:tl-W1) for the antifield letters, explicit su(2) route
(ea_explicit.value_vertex_operator): for N = int omega <zeta, X^*>, {I_T, N} is the vertex T with one leg fixed to zeta
(common constant on both sides; the odd parameter of a fermionic zeta-leg stands leftmost on both sides), and
   lc({I,N}) = u sum_T m_T {I_T,N}
on every multi-field output leg set (the one-field outputs carry the kinetic term m QN and are not tested here).
Multipliers (proof of Thm c4-log): AAA 6, AAAA 8, llA 10, AcAs 1, cllS 1, cccS 1."""
import numpy as np
from fractions import Fraction
import ea_core as E, ea_explicit as X
C = E.C
MULT = {'AAA': 6, 'AAAA': 8, 'llA': 10, 'AcAs': 1, 'cllS': 1, 'cccS': 1}
CASES = {
    'A*_perp': ('A', 'P', [['A:V', 'A:P'], ['l:lamP', 'l:lamP'], ['l:lamM', 'l:lamM'], ['c:', 'As:P'],
                          ['A:V', 'A:V', 'A:P'], ['A:P', 'A:P', 'A:P']]),
    'A*_V': ('A', 'V', [['A:V', 'A:V'], ['A:P', 'A:P'], ['l:lamP', 'l:lamM'], ['c:', 'As:V'], ['A:V', 'A:V', 'A:V'],
                        ['A:V', 'A:P', 'A:P']]),
    'lambda*': ('l', None, [['l:lamP', 'A:V'], ['l:lamM', 'A:V'], ['l:lamP', 'A:P'], ['l:lamM', 'A:P'], ['c:', 'ls:lamsP'], ['c:', 'ls:lamsM']]),
    'c*': ('c', None, [['A:V', 'As:V'], ['A:P', 'As:P'], ['l:lamP', 'ls:lamsP'], ['l:lamM', 'ls:lamsM'], ['c:', 'cs:']]),
}


def rpolfull(rng, t, which=None):
    x = np.zeros(C.SHAPE[t], complex)
    if t in ('A', 'As'):
        Ms = range(4) if which == 'V' else range(4, 10)
        for m in Ms:
            x[m] = rng.integers(-1, 2, 3)
    elif t in ('l', 'ls'):
        if which is None:
            v = rng.integers(-1, 2, 16)
        else:
            v = sum(int(c) * b for c, b in zip(rng.integers(-1, 2, 8), E.BASIS[which]))
        x = np.einsum('i,a->ia', v, rng.integers(-1, 2, 3)).astype(complex)
    else:
        x = rng.integers(-1, 2, 3).astype(complex)
    return x


def gate(seed=3, ntrials=2):
    rng = np.random.default_rng(seed)
    res = {}
    for name, (zt, zw, legsets) in CASES.items():
        rows = []
        for ls in legsets:
            for it in range(ntrials):
                zeta = rpolfull(rng, zt, zw)
                outs = []
                for spec in ls:
                    t, w = spec.split(':')
                    outs.append((t, rpolfull(rng, t, w or None), rng.integers(-2, 3, 4)))
                pom = -sum((o[2] for o in outs), np.zeros(4, int))
                zl = (zt, zeta, pom)
                lhs = E.Q0; rhs = 0j
                for vt, m in MULT.items():
                    if zt not in C.CANON[vt]:
                        continue
                    r = X.value_vertex_operator(vt, zl, outs)
                    if 'log' in r:
                        lhs = E.q_add(lhs, r['log'])
                    rhs += m * X.tree_vertex_operator(vt, zl, outs)
                rhs *= C.C_ADJ
                l = complex(float(lhs[0]), float(lhs[1]))
                rows.append({'legs': ls, 'lc': E.qstr(lhs), 'u_sum_m_T': str(rhs), 'pass': abs(l - rhs) < 1e-9 * max(1, abs(rhs))})
        res[name] = rows
        res[name + '_pass'] = all(r['pass'] for r in rows)
        res[name + '_nontrivial'] = sum(1 for r in rows if r['lc'] != '0')
        print(name, res[name + '_pass'], 'nontrivial', res[name + '_nontrivial'], 'of', len(rows), flush=True)
    return res


if __name__ == '__main__':
    gate()
