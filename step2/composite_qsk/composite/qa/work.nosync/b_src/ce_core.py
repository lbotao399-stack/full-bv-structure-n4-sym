"""Engine B (route B), PART C: slot space, edge matrices and their identities.

Slot space INT6 = A (10), lam (16), c (1), As (10), lams (16), cs (1): 54 components (realisation of the Step-1
engine_b/realization.py; colour stripped).  Edge matrices are the components E^{XY}(k) of a kernel in the slot basis,
read as in manuscript Def c3-P: E^{XY} = d_E(f_X f_Y), row end at x with plane wave e^{ikx}, column end at y with
e^{-iky}; (Y (x) 1)E has matrix y(k) E(k), (1 (x) Y)E has matrix (-1)^{|Y||X|} E(k) y(-k)^T.

  N(k)    numerator of P(eps,oo) = N(k) e^{-eps k^2}/k^2              (Step-1 engine_b/loops.prop_numerator)
  k0      K_0 = -(fibre pairing)^{-1}                                  (manuscript Def c1-bracket)
  qgf(k)  symbol of Q^GF                                               (manuscript Def c3-qgf)
  qq(k)   symbol of Q                                                  (manuscript Prop c1-Q)
  MX(k)   numerator of the square kernel X = ((P*)^2 (x) 1)K_0 = MX(k) e^{-2 eps k^2}/k^4, MX = qgf N
Identities checked by selftest(): N = qgf k0;  qq N + Pi N qq(-k)^T = k^2 k0  (Prop c3-QP);  MX = qgf^2 k0;
MX graded antisymmetric; qq MX - Pi MX qq(-k)^T... (Q (x) 1 + 1 (x) Q) K_{P_b} = 0 (Lemma lem:cc-kernels (v));
the components of MX lie in c (x) A|_V and A|_V (x) c (Lemma lem:cc-kernels (iii))."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
# reused verified Step-1 modules (engine_b/vertices.py, loops.py, realization.py) vendored in src.nosync/vendor
# (iCloud evicts files outside *.nosync folders; copies identical to the originals of 2026-10-29)
EB1 = os.path.join(HERE, 'vendor')
if EB1 not in sys.path:
    sys.path.insert(0, EB1)
_cwd = os.getcwd()
os.chdir(EB1)
import sympy as sp
import vertices as V
import loops as L
os.chdir(_cwd)

I = sp.I
INT6 = ['A', 'lam', 'c', 'As', 'lams', 'cs']
DIM = {s: V.DIM[s] for s in INT6}
PAR = {s: V.PAR[s] for s in INT6}
OFF = {}
_o = 0
for _s in INT6:
    OFF[_s] = _o
    _o += DIM[_s]
NT = _o
SPEC_OF = []
for _s in INT6:
    SPEC_OF += [_s] * DIM[_s]
PARITY = [PAR[s] for s in SPEC_OF]
assert all(OFF[s] == L.OFF[s] for s in L.INT)
Mp = V.Mp
Mpinv = V.Mpinv
Rm = __import__('realization').R      # rho(e_M): Sigma -> Sigma^*
Rbm = __import__('realization').Rb    # rho(e_M): Sigma^* -> Sigma

# letter slots (names as in the census / twoletter engine): (species, component list)
SLOT = {'AV': ('A', list(range(4))), 'AW': ('A', list(range(4, 10))),
        'lamP': ('lam', list(range(8))), 'lamM': ('lam', list(range(8, 16))), 'c': ('c', [0]),
        'AsV': ('As', list(range(4))), 'AsW': ('As', list(range(4, 10))),
        'lamsM': ('lams', list(range(8))), 'lamsP': ('lams', list(range(8, 16))), 'cs': ('cs', [0])}
WEYL = {'c': 0, 'AV': 1, 'AW': 1, 'lamP': sp.Rational(3, 2), 'lamM': sp.Rational(3, 2), 'AsV': 3, 'AsW': 3,
        'lamsP': sp.Rational(5, 2), 'lamsM': sp.Rational(5, 2), 'cs': 4}
DEG = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 0, 'lamM': 0, 'AsV': -1, 'AsW': -1, 'lamsP': -1, 'lamsM': -1, 'cs': -2}


def zero():
    return sp.zeros(NT, NT)


def N_matrix(k):
    """54 x 54 numerator of P(eps,oo) at momentum k (symbols or numbers)."""
    N37 = L.prop_numerator(list(k))
    M = zero()
    for i in range(N37.rows):
        for j in range(N37.cols):
            if N37[i, j] != 0:
                M[i, j] = N37[i, j]
    return M


def k0_matrix():
    M = zero()
    for m in range(10):
        M[OFF['A'] + m, OFF['As'] + m] = -1
        M[OFF['As'] + m, OFF['A'] + m] = -1
    for i in range(16):
        for j in range(16):
            M[OFF['lam'] + i, OFF['lams'] + j] = -Mpinv[j, i]
            M[OFF['lams'] + i, OFF['lam'] + j] = -Mpinv[i, j]
    M[OFF['c'], OFF['cs']] = -1
    M[OFF['cs'], OFF['c']] = -1
    return M


def qgf_matrix(k):
    """symbol of Q^GF (Def c3-qgf): A* -> -2 A* (in A); A -> i iota_k A|_V (in c); c* -> i k c* (in A*|_V);
    lam* -> -i rho(k) lam* (in lam)."""
    M = zero()
    for m in range(10):
        M[OFF['A'] + m, OFF['As'] + m] = -2
    for m in range(4):
        M[OFF['c'], OFF['A'] + m] = I * k[m]
        M[OFF['As'] + m, OFF['cs']] = I * k[m]
    rho = sum((k[m] * Rbm[m] for m in range(4)), sp.zeros(16, 16))
    for i in range(16):
        for j in range(16):
            if rho[i, j] != 0:
                M[OFF['lam'] + i, OFF['lams'] + j] = -I * rho[i, j]
    return M


def q_matrix(k):
    """symbol of Q (Prop c1-Q): c -> -i k c (in A|_V); A -> -1/2(|k|^2 A_V - k (k.A_V)), -1/2 |k|^2 A_perp (in A*);
    A* -> -i (k, A*|_V) (in c*); lam -> i rho(k) lam (in lam*)."""
    M = zero()
    k2 = sum(k[m] ** 2 for m in range(4))
    for m in range(4):
        M[OFF['A'] + m, OFF['c']] = -I * k[m]
        M[OFF['cs'], OFF['As'] + m] = -I * k[m]
        for n in range(4):
            M[OFF['As'] + m, OFF['A'] + n] = sp.Rational(-1, 2) * (k2 * (1 if m == n else 0) - k[m] * k[n])
    for m in range(4, 10):
        M[OFF['As'] + m, OFF['A'] + m] = sp.Rational(-1, 2) * k2
    rho = sum((k[m] * Rm[m] for m in range(4)), sp.zeros(16, 16))
    for i in range(16):
        for j in range(16):
            if rho[i, j] != 0:
                M[OFF['lams'] + i, OFF['lam'] + j] = I * rho[i, j]
    return M


def MX_matrix(k):
    return (qgf_matrix(k) * N_matrix(k)).applyfunc(sp.expand)


PI = sp.diag(*[(-1) ** p for p in PARITY])


def selftest():
    k = sp.symbols('k0:4')
    mk = [-x for x in k]
    k2 = sum(x * x for x in k)
    out = {}
    N = N_matrix(k); K0 = k0_matrix(); G = qgf_matrix(k); Qm = q_matrix(k); Qmm = q_matrix(mk)
    out['N = qgf k0'] = (G * K0 - N).applyfunc(sp.expand).is_zero_matrix
    out['qq N + Pi N qq(-k)^T = |k|^2 k0 (Prop c3-QP)'] = (Qm * N + PI * N * Qmm.T - k2 * K0).applyfunc(sp.expand).is_zero_matrix
    out['Q^2 = 0'] = (Qm * Qm).applyfunc(sp.expand).is_zero_matrix
    out['[Q,Q^GF] = |k|^2'] = (Qm * G + G * Qm - k2 * sp.eye(NT)).applyfunc(sp.expand).is_zero_matrix
    out['(Q (x) 1 + 1 (x) Q) K0 = 0'] = (Qm * K0 + PI * K0 * Qmm.T).applyfunc(sp.expand).is_zero_matrix
    out['K0 graded symmetric'] = all(sp.expand(K0[i, j] - (-1) ** (PARITY[i] * PARITY[j]) * K0[j, i]) == 0
                                    for i in range(NT) for j in range(NT))
    MX = MX_matrix(k)
    out['MX = qgf^2 k0'] = (G * G * K0 - MX).applyfunc(sp.expand).is_zero_matrix
    MXm = MX_matrix(mk)
    out['MX graded antisymmetric: MX^{XY}(k) = -(-1)^{|X||Y|} MX^{YX}(-k)'] = all(
        sp.expand(MX[i, j] + (-1) ** (PARITY[i] * PARITY[j]) * MXm[j, i]) == 0 for i in range(NT) for j in range(NT))
    nz = sorted(set((SPEC_OF[i] + ('V' if SPEC_OF[i] in ('A', 'As') and i - OFF[SPEC_OF[i]] < 4 else ''),
                     SPEC_OF[j] + ('V' if SPEC_OF[j] in ('A', 'As') and j - OFF[SPEC_OF[j]] < 4 else ''))
                    for i in range(NT) for j in range(NT) if MX[i, j] != 0))
    out['MX components'] = nz
    out['MX components only c(x)A_V, A_V(x)c'] = set(nz) == {('c', 'AV'), ('AV', 'c')}
    out['MX^{c,A_m} = 2 i k_m, MX^{A_m,c} = 2 i k_m'] = all(sp.expand(MX[OFF['c'], OFF['A'] + m] - 2 * I * k[m]) == 0 and
                                                          sp.expand(MX[OFF['A'] + m, OFF['c']] - 2 * I * k[m]) == 0 for m in range(4))
    # K_{P_b} = 2 i (b.k) MX ; (Q (x) 1 + 1 (x) Q) K_{P_b} = 0
    b = sp.symbols('b0:4')
    KP = (2 * I * sum(b[m] * k[m] for m in range(4))) * MX
    KPm = (2 * I * sum(b[m] * mk[m] for m in range(4))) * MXm
    out['(Q (x) 1 + 1 (x) Q) K_{P_b} = 0 (Lemma lem:cc-kernels (v))'] = (Qm * KP + PI * KP * Qmm.T).applyfunc(sp.expand).is_zero_matrix
    out['K_{P_b} graded symmetric'] = all(sp.expand(KP[i, j] - (-1) ** (PARITY[i] * PARITY[j]) * KPm[j, i]) == 0
                                          for i in range(NT) for j in range(NT))
    return out


if __name__ == '__main__':
    import json
    r = selftest()
    for kk, v in r.items():
        print(kk, ':', v)
