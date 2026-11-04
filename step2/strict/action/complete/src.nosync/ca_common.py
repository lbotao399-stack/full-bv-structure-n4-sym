"""Route A (engine A code path) for the action-side completion: imports engine A's verified modules read-only
(step2/strict/action/engine_a/src.nosync: sa_core cycles with GSimplex exact simplex integrals, sa_ob position
injection, sa_mode exact mode model (Q, bracket of Def c1-bracket), sa_fit); realisation of cliff.py:
rho(e_M)|_Sigma = sb[M] (Sigma -> Sigma^*), rho(e_M)|_Sigma^* = s[M], pairing (u, w) = u^T w."""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
COMPLETE = os.path.dirname(HERE)
ASRC = os.path.abspath(os.path.join(COMPLETE, '..', 'engine_a', 'src.nosync'))
if ASRC not in sys.path:
    sys.path.insert(0, ASRC)
import numpy as np
import sa_core as SA
import sa_ob as SO
import sa_fit as SF
import sa_mode as MD
import qa_s1 as S1
LOGF = os.path.join(COMPLETE, 'PROGRESS.txt')
SB = [np.asarray(m, complex) for m in SA.E.SGB]
SS = [np.asarray(m, complex) for m in SA.E.SG]
MP = np.eye(16, dtype=complex)
TYPE2SP = {'A': 'A', 'As': 'As', 'l': 'lam', 'ls': 'lams', 'c': 'c', 'cs': 'cs'}


def log(s, tag='routeA'):
    line = time.strftime('%Y-%m-%d %H:%M:%S') + ' ' + tag + ': ' + s
    print(line, flush=True)
    with open(LOGF, 'a') as f:
        f.write(line + '\n')


def generic_legs(L):
    """engine A legs -> generic dicts (sp, u, p) for the candidate densities."""
    return [dict(sp=TYPE2SP[SA.E.SLOT[l['t']]], u=np.asarray(l['pol'], complex),
                 p=[int(round(complex(x).real)) for x in l['p0']]) for l in L]


def mm_legs(L):
    return [(SA.E.SLOT[l['t']], np.asarray(l['pol'], complex)) for l in L]


def qmat_A(k):
    """symbol of Q (Prop c1-Q) in engine A's realisation on the 54-slot space of engine B's ordering
    (A 10, lam 16, c, As 10, lams 16, cs): c -> -i k c (A_V); A -> -1/2(k^2 A_V - k k.A_V), -1/2 k^2 A_perp (A*);
    A* -> -i (k, A*_V) (c*); lam -> i rho(k) lam (lam*), rho(k) = sum k_mu sb[mu]."""
    OFF = {'A': 0, 'lam': 10, 'c': 26, 'As': 27, 'lams': 37, 'cs': 53}
    M = np.zeros((54, 54), complex)
    k = np.asarray(k, complex)
    k2 = k @ k
    for m in range(4):
        M[OFF['A'] + m, OFF['c']] = -1j * k[m]
        M[OFF['cs'], OFF['As'] + m] = -1j * k[m]
        for n in range(4):
            M[OFF['As'] + m, OFF['A'] + n] = -0.5 * (k2 * (m == n) - k[m] * k[n])
    for m in range(4, 10):
        M[OFF['As'] + m, OFF['A'] + m] = -0.5 * k2
    rho = sum(k[m] * SB[m] for m in range(4))
    M[OFF['lams']:OFF['lams'] + 16, OFF['lam']:OFF['lam'] + 16] = 1j * rho
    return M
