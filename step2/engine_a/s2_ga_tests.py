"""Consistency tests of the mode model s2_ga (run for S2_ALG=su2 and gl3):
 (1) Q on mode variables: restricted-pair formula = full bracket {S_free, y};
 (2) Q^2 = 0 on every mode variable;
 (3) Q Pnum = p^2 Khat in the mode model: d_{Pnum(p)}( Q(y_u y_v) ) = p^2 d_{Khat}(y_u y_v) for all pairs (u at p, v at -p)
     ([C07 Lemma 10.2.1] at T = infty: Q P(eps,infty) = K_eps), i.e. the propagator components of s2_core are the
     ones fixed by Q and the pairing;
 (4) the derivation form of Q on a random cubic monomial set equals the bracket {S_free, F}."""
import json, sys, itertools, random
from fractions import Fraction
import numpy as np
import s2_core as C
from r2_core import Q2
import s2_ga as GAm
from grass import GA, PAR

rng = np.random.default_rng(4)
p = np.array([1, -2, 0, 1])
M = GAm.Model([p, -p])
pm = [(0, 1), (1, 0)]
Sf = M.S_free(pm)
res = {'algebra': C.ALG}
bad1 = 0; bad2 = 0; nv = 0
allv = [M.var[(t, i, m)] for t in C.TYPES for i in range(C.SIZE[t]) for m in (0, 1)]
sample = random.Random(1).sample(allv, min(60, len(allv)))
for v in sample:
    a = M.Q_var(v, Sf, pm); b = M.Q_var_full(v, Sf, pm)
    if (a - b).t:
        bad1 += 1
    nv += 1
for v in allv:
    q = M.Q_var(v, Sf, pm)
    qq = M.Q_apply(q, Sf, pm)
    if qq.t:
        bad2 += 1
res['Q_restricted_vs_full_bracket_mismatches'] = '%d / %d' % (bad1, nv)
res['Q^2_nonzero_on_variables'] = '%d / %d' % (bad2, len(allv))
# (3)
Pk = M.kernel('P', 0, 1); Kk = M.kernel('K', 0, 1)
p2 = int(p @ p)
bad3 = 0; n3 = 0; nonzeroK = 0
for t1 in C.TYPES:
    for i1 in range(C.SIZE[t1]):
        u = M.var[(t1, i1, 0)]
        for t2 in C.TYPES:
            for i2 in range(C.SIZE[t2]):
                if (i1 + i2) % 3:              # subsample
                    continue
                w = M.var[(t2, i2, 1)]
                yy = GA.var(u) * GA.var(w)
                Qyy = M.Q_apply(yy, Sf, pm)
                lhs = M.dX(Qyy, Pk); rhs = M.dX(yy, Kk) * p2
                if (lhs - rhs).t:
                    bad3 += 1
                if rhs.t:
                    nonzeroK += 1
                n3 += 1
res['QP=p^2 K mismatches'] = '%d / %d pairs (%d with nonzero K)' % (bad3, n3, nonzeroK)
# (4)
q0 = np.array([1, 0, 2, -1]); q1 = np.array([0, 1, -1, 1]); q2_ = -(q0 + q1)
M2 = GAm.Model([q0, q1, q2_, -q0, -q1, -q2_])
pm2 = [(0, 3), (1, 4), (2, 5), (3, 0), (4, 1), (5, 2)]
Sf2 = M2.S_free(pm2)
F = M2.cubic('AcAs', [(0, 1, 2), (2, 0, 1)]) + M2.cubic('cllS', [(1, 2, 0)])
keys = list(F.t)[:40]
Fs = GA({k: F.t[k] for k in keys})
der = M2.Q_apply(Fs, Sf2, pm2)
br = M2.bracket0(Sf2, Fs, 0, pm2)
res['Q_derivation_vs_bracket_on_cubic_terms_mismatch'] = len((der - br).t)
res['Q_derivation_terms'] = len(der.t)
print(json.dumps(res, indent=1))
json.dump(res, open('logs/s2_ga_tests_%s.json' % C.ALG, 'w'), indent=1)
