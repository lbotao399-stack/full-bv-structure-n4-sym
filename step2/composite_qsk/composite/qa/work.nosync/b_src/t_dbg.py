import time, sys
import sympy as sp
from sympy.polys.domains import QQ_I
import ce_job as CJ
I = sp.I
w1 = {4: 1, 5: 2, 6: -1, 7: 0, 8: 3, 9: 1}
w2 = {4: 2, 5: -1, 6: 1, 7: 1, 8: 0, 9: -2}
Y0 = {4: 1, 5: 0, 6: 2, 7: 0, 8: -1, 9: 1}; Y1 = {4: 0, 5: 1, 6: 1, 7: 2, 8: 0, 9: 1}
job = CJ.CJob('AW', 'AW', ['AW', 'AW'], pol=[Y0, Y1], za=w1, zb=w2)
R = job.R
def hpart(res):
    out = {}
    for m, d in res.items():
        for b, P in d.items():
            red = R.zero
            for mu in range(4):
                red += P.diff(job.av[mu]).diff(job.hv[mu])
            red = red * R(QQ_I(0, -1))
            # keep only h-free, a-free part of the result
            red = R.from_dict({k: v for k, v in red.items() if sum(k[i] for i in job.h_idx + job.a_idx) == 0 and k[job.hm_idx] == 0})
            out[b] = out.get(b, R.zero) + red
    return out
def hmpart(res):
    out = {}
    for m, d in res.items():
        for b, P in d.items():
            red = P.diff(job.g['hm'])
            red = R.from_dict({k: v for k, v in red.items() if sum(k[i] for i in job.h_idx + job.a_idx) == 0})
            out[b] = out.get(b, R.zero) + red
    return out
for seq in job.sequences():
    k = len(seq)
    for c in range(k + 1):
        for s in range(k + 1):
            if s == c: continue
            tA = job.config(seq, c, s, c if c >= 1 else None, 'D2')
            tB = job.config(seq, c, s, c + 1 if c + 1 <= k else None, 'D2')
            tX = job.config(seq, c, s, None, 'X')
            for w in set(tA[2]) | set(tB[2]):
                A = job.integrate(tA[2].get(w, R.zero), tA[1], tA[0], (0,)) if w in tA[2] else {}
                B = job.integrate(tB[2].get(w, R.zero), tB[1], tB[0], (0,)) if w in tB[2] else {}
                Xc = job.integrate(tX[2].get(w, R.zero), tX[1], tX[0], (1,), cut_override='Ksingle2') if w in tX[2] else {}
                # X graph (two segments) at eps^0, constant line factor 1:
                X0 = job.integrate(tX[2].get(w, R.zero), tX[1], tX[0], (0,)) if w in tX[2] else {}
                hA, hB = hpart(A), hpart(B)
                diff = {b: (hA.get(b, R.zero) - hB.get(b, R.zero)).as_expr() for b in set(hA) | set(hB)}
                pred = {}
                for b, P in X0.get(0, {}).items():
                    pred[b] = pred.get(b, 0) + (-2 * P).as_expr()
                for b, P in Xc.get(1, {}).items():
                    pred[b] = pred.get(b, 0) + (8 * P).as_expr()
                print(seq, 'c', c, 's', s, 'w', w, '\n   A-B h-part:', diff, '\n   -2X + 8 g1 :', pred, flush=True)
