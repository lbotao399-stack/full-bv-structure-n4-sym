"""QA: engine A's G-R1 term B = L_1^(0)(-F1) on the (c, lam) sector, the two output orderings separately."""
import os, sys, json
import numpy as np
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_s1 as S1, qa_core as QC, qa_util as U
SB = S1.SB
resA = json.load(open(os.path.join(QSK, 'engine_a/results.json')))
xi = np.array([complex(x) for x in resA['tasks']['TQ2']['xi']])
R = U.Rng(77); g0 = QC.zero_gen()
out = []
for ip in range(2):
    w1, w2 = R.gv(6), R.gv(6)
    for lt in ('lamM', 'lamP'):
        qc, ql = R.mom(), R.mom()
        l = R.vec(lt)
        if lt != 'lamM':
            continue
        legs = [('c', np.array([1 + 0j]), qc), (lt, l, ql)]
        b1 = U.vzero(); b2 = U.vzero()
        for wa, wb in ((w1, w2), (w2, w1)):
            zl = sum(wb[i] * (SB[4 + i] @ xi) for i in range(6))
            Xa = ('AW', wa, None); Xb = ('lamP', zl, None)
            b1 = U.vadd(b1, QC.pair_value(Xa, Xb, legs, g0)['fin'])
            b2 = U.vadd(b2, QC.pair_value(Xa, Xb, legs[::-1], g0)['fin'])
        rec = {'point': ip, 'b1 (c,lam) of L1^(0)(F1)': U.vstr(b1), 'b2 (lam,c)': U.vstr(b2),
               'b1+b2': U.vstr(U.vadd(b1, b2)), 'b1-b2': U.vstr(U.vadd(b1, U.vneg(b2)))}
        print(rec, flush=True); out.append(rec)
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/r1_check2.json'), 'w'), indent=1)
