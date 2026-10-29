# exact rank-one elements of S_+ (x) W in the Step-1 engine basis: simultaneous eigenvectors of the commuting
# operators sigbar[M] sig[N] for (M,N) = (0,1),(2,3),(4,5),(6,7),(8,9) (weight vectors), lying in Omega = +1.
import os, sys, itertools, json
import numpy as np
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
sys.path.insert(0, os.path.join(ROOT, 'superconformal'))
import model as Mo
np.seterr(all='ignore')
sig = [np.array([[complex(*map(float, Mo.SIG[M][a][b])) for b in range(16)] for a in range(16)]) for M in range(10)]
sigb = [np.array([[complex(*map(float, Mo.SIGB[M][a][b])) for b in range(16)] for a in range(16)]) for M in range(10)]
mm = lambda *A: __import__('functools').reduce(lambda X, Y: X.dot(Y), A)
Om = mm(sigb[0], sig[1], sigb[2], sig[3])
H = [mm(sigb[M], sig[N]) for (M, N) in [(0,1),(2,3),(4,5),(6,7),(8,9)]]
# H^2 = -1 on Sigma (M != N anticommute, square 1): eigenvalues +-i; projector (1 - i s H)/2 onto eigenvalue i s
res = []
for signs in itertools.product([1, -1], repeat=5):
    Pj = np.eye(16, dtype=complex)
    for s, h in zip(signs, H):
        Pj = Pj.dot((np.eye(16) - 1j*s*h)/2)
    for a in range(16):
        v = Pj[:, a]
        if np.linalg.norm(v) > 1e-9:
            break
    else:
        continue
    if np.linalg.norm(Om.dot(v) - v) > 1e-9:
        continue
    g = max(abs(v.dot(sig[M]).dot(v)) for M in range(10))
    rk = np.linalg.matrix_rank(np.array([sig[M].dot(v) for M in range(10)]), tol=1e-9)
    vv = [(round(z.real*4)/4, round(z.imag*4)/4) for z in v]
    res.append({'weights': signs, 'max_|Gamma(xi,xi)|': g, 'rank_Gamma(xi,-)': int(rk),
                'vector_quarter_rounded': vv, 'rounding_err': float(np.linalg.norm(v - np.array([x+1j*y for x, y in vv])))})
for r in res:
    print(r['weights'], 'maxGamma=%.2e' % r['max_|Gamma(xi,xi)|'], 'rank=', r['rank_Gamma(xi,-)'], 'roundErr=%.1e' % r['rounding_err'])
json.dump([{k: (v if k != 'weights' else list(v)) for k, v in r.items()} for r in res],
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rank_one_vectors.json'), 'w'), indent=0)
