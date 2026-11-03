"""QA: engine A's gate G-R1 (read-only import): the term C = {S^(1)(xi), L_1^(0)(O)} on the (c, lam) sector, the two
output orderings v1 = (c, A_V), v2 = (A_V, c) separately.  Engine A combined them only with a relative minus sign
(both 'Koszul readings' tested are differences), while B = L_1^(0)(-F1) uses the plain sum over the orderings of the two
odd outputs (c, lam) (qa_util.both_orders)."""
import os, sys, json
import numpy as np
from fractions import Fraction
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_s1 as S1, qa_core as QC, qa_util as U
SB = S1.SB
resA = json.load(open(os.path.join(QSK, 'engine_a/results.json')))
xi = np.array([complex(x) for x in resA['tasks']['TQ2']['xi']])
R = U.Rng(77)
g0 = QC.zero_gen()
out = []
for ip in range(2):
    w1, w2 = R.gv(6), R.gv(6)
    for lt in ('lamM', 'lamP'):
        qc, ql = R.mom(), R.mom()
        l = R.vec(lt)
        if lt != 'lamM':
            continue
        pv = np.array([-(xi @ SB[m] @ l) for m in range(4)])
        Xa = ('AW', w1, None); Xb = ('AW', w2, None)
        v1 = QC.pair_value(Xa, Xb, [('c', np.array([1 + 0j]), qc), ('AV', pv, ql)], g0)['fin']
        v2 = QC.pair_value(Xa, Xb, [('AV', pv, ql), ('c', np.array([1 + 0j]), qc)], g0)['fin']
        # also L_1^(0)(O) itself on (c, A_V) at generic A_V polarisation
        pg = np.array([1, 2, -1, 1], dtype=complex)
        u1 = QC.pair_value(Xa, Xb, [('c', np.array([1 + 0j]), qc), ('AV', pg, ql)], g0)['fin']
        u2 = QC.pair_value(Xa, Xb, [('AV', pg, ql), ('c', np.array([1 + 0j]), qc)], g0)['fin']
        rec = {'point': ip, 'v1 (c,AV) at A_V = -Gamma(xi,l)|_V': U.vstr(v1), 'v2 (AV,c)': U.vstr(v2),
               'v1+v2': U.vstr(U.vadd(v1, v2)), 'v1-v2': U.vstr(U.vadd(v1, U.vneg(v2))),
               'L1^(0)(O) on (c,AV) generic pol: (c,AV)': U.vstr(u1), '(AV,c)': U.vstr(u2)}
        print(rec, flush=True)
        out.append(rec)
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/r1_check.json'), 'w'), indent=1)
