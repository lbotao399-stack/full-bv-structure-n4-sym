"""Brute-force Grassmann check (engine_a/check_vertices.build, grass.GA) of the S^(1)(u e^{iPy}) Hessian of kt_engine:
density of ch02 eq:c2-S1even written directly on the field arrays, gamma placed last; compares every component pair."""
import sys, json
import numpy as np
from kt_engine import *
sys.path.insert(0, ROOT + '/engine_a')
import check_vertices as cv
from grass import GA, newvar
rng = np.random.default_rng(3)
u = rng.integers(-2, 3, 4).astype(float); P = rng.integers(-2, 3, 4).astype(float)
kj = rng.integers(-2, 3, 4).astype(float); kl = rng.integers(-2, 3, 4).astype(float)
H = s1_hessian(u, kj, kl, P)
sP = sum(P[m] * C.SARR[m] for m in range(4)); sbu = sum(u[m] * C.SBARR[m] for m in range(4)); Ml = 1j / 4 * (sP @ sbu)
uP = u @ P
def density(fld, der):
    tot = GA()
    for a in range(D):
        for M in range(10):
            LvA = GA()
            for n in range(4):
                LvA = LvA + der['A'][n, M, a] * complex(u[n])
            if M < 4:
                for n in range(4):
                    LvA = LvA + fld['A'][n, a] * complex(1j * P[M] * u[n])
            else:
                LvA = LvA + fld['A'][M, a] * complex(0.25j * uP)
            tot = tot + LvA * fld['As'][M, a]
        for be in range(16):
            vl = GA()
            for n in range(4):
                vl = vl + der['l'][n, be, a] * complex(u[n])
            vl = vl + fld['l'][be, a] * complex(1j / 8 * uP)
            for al in range(16):
                if Ml[be, al] != 0:
                    vl = vl + fld['l'][al, a] * complex(Ml[be, al])
            tot = tot - vl * fld['ls'][be, a]
        vc = GA()
        for n in range(4):
            vc = vc + der['c'][n, a] * complex(u[n])
        tot = tot - vc * fld['cs'][a]
    return tot
TY = ['c', 'A', 'As', 'l', 'ls', 'cs']
maxdev = 0; nchk = 0; nnz = 0
for tj in TY:
    for tl in TY:
        for trial in range(3):
            ij = rng.integers(0, C.SIZE[tj]); il = rng.integers(0, C.SIZE[tl])
            polj = np.zeros(C.SIZE[tj], complex); polj[ij] = 1; polj = polj.reshape(C.SHAPE[tj])
            poll = np.zeros(C.SIZE[tl], complex); poll[il] = 1; poll = poll.reshape(C.SHAPE[tl])
            fld, der, vs = cv.build([(tj, polj, kj), (tl, poll, kl)])
            g = newvar(99999 + nchk, 1)
            val = (density(fld, der) * GA.var(g)).coeff(vs + [g])
            mine = (-1) ** (PARITY[tj] * PARITY[tl]) * H[C.OFF[tj] + ij, C.OFF[tl] + il]
            maxdev = max(maxdev, abs(complex(val) - mine)); nchk += 1; nnz += abs(mine) > 0
# forced nonzero pairs
for (tj, ij, tl, il) in [('A', 1 * D + 0, 'As', 2 * D + 0), ('A', 5 * D + 1, 'As', 5 * D + 1), ('l', 3 * D + 2, 'ls', 7 * D + 2),
                         ('ls', 7 * D + 2, 'l', 3 * D + 2), ('c', 1, 'cs', 1), ('As', 0, 'A', 3 * D)]:
    polj = np.zeros(C.SIZE[tj], complex); polj[ij] = 1; polj = polj.reshape(C.SHAPE[tj])
    poll = np.zeros(C.SIZE[tl], complex); poll[il] = 1; poll = poll.reshape(C.SHAPE[tl])
    fld, der, vs = cv.build([(tj, polj, kj), (tl, poll, kl)])
    g = newvar(199999 + nchk, 1)
    val = (density(fld, der) * GA.var(g)).coeff(vs + [g])
    mine = (-1) ** (PARITY[tj] * PARITY[tl]) * H[C.OFF[tj] + ij, C.OFF[tl] + il]
    maxdev = max(maxdev, abs(complex(val) - mine)); nchk += 1; nnz += abs(mine) > 0
out = {'checked_pairs': nchk, 'nonzero_entries_checked': int(nnz), 'max_abs_deviation': float(maxdev)}
print(out); json.dump(out, open('check_s1_vertex.json', 'w'), indent=1)
