"""Item (b): one-loop graphs with one insertion of O = int omega (iota_{w1}A, iota_{w2}A), omega = e^{iq.x} dvol
(background-field route [CG2 Lemma 10.1.1.3]); delta-linear part of W_1 = 1/2 sum_n str[(P S d2I)^n P S d2O].
Units of reported numbers: hbar c_adj/(16 pi^2) (vacuum term: hbar/(16 pi^2))."""
import json, sys
import numpy as np
from oneloop import Engine
from sector_a import rat, E4
gname = sys.argv[1] if len(sys.argv) > 1 else 'su2'
eng = Engine(gname); d = eng.d; o = eng.off
rng = np.random.default_rng(23)
S = np.diag(eng.S)
LINES = {'AV': [o['A'] + M * d + a for M in range(4) for a in range(d)],
         'Aperp': [o['A'] + M * d + a for M in range(4, 10) for a in range(d)]}
def Nm(k, L):
    N = eng.Nprop(k); m = np.zeros(eng.n, bool); m[LINES[L]] = True
    return N * np.outer(m, m)

def HO(w1, w2):
    """d_j d_l O (momentum-independent apart from conservation): (iota_w1 x, iota_w2 y) + (iota_w1 y, iota_w2 x)"""
    H = np.zeros((eng.n, eng.n), complex)
    W1 = np.zeros(10); W1[4:] = w1; W2 = np.zeros(10); W2[4:] = w2
    blk = np.kron(np.outer(W1, W2) + np.outer(W2, W1), np.eye(d))
    H[o['A']:o['A'] + 10 * d, o['A']:o['A'] + 10 * d] = blk
    return H

def tri_lap(b1, b2, p1, p2, w1, w2, L):
    """1/2[str(X1 X2 XO) + str(X2 X1 XO)] : Laplacian (constant) of the numerator, lines (perp, L, perp)"""
    hO = HO(w1, w2)
    tot = 0
    for (ba, pa), (bb, pb) in (((b1, p1), (b2, p2)), ((b2, p2), (b1, p1))):
        def F(k):
            ha = eng.hessian([('A', ba, pa)], -k, k - pa)
            hb = eng.hessian([('A', bb, pb)], -(k - pa), k - pa - pb)
            return eng.str_(Nm(k, 'Aperp') @ S @ ha @ Nm(k - pa, L) @ S @ hb @ Nm(k - pa - pb, 'Aperp') @ S @ hO)
        k0 = np.array([0.2, -0.1, 0.4, 0.3])
        f0 = F(k0); lap = sum(F(k0 + e) + F(k0 - e) - 2 * f0 for e in E4)
        cub = max(abs(F(k0 + 2 * e) - 3 * F(k0 + e) + 3 * f0 - F(k0 - e)) for e in E4)
        assert cub < 1e-8, cub
        tot += 0.5 * lap
    return tot

def bub_num(b1, b2, p1, p2, w1, w2):
    """1/2 str(X4 XO): constant numerator"""
    k = np.array([0.3, 0.2, -0.5, 0.1])
    h4 = eng.hessian([('A', b1, p1), ('A', b2, p2)], -k, k - p1 - p2)
    return 0.5 * eng.str_(Nm(k, 'Aperp') @ S @ h4 @ Nm(k - p1 - p2, 'Aperp') @ S @ HO(w1, w2))

def rpol(which):
    x = rng.integers(-2, 3, size=(10, d)).astype(complex)
    if which == 'V': x[4:] = 0
    else: x[:4] = 0
    return x

def structs(b1, b2, w1, w2):
    W = lambda w: np.concatenate([np.zeros(4), w])
    iw = lambda a, w: W(w) @ a     # iota_w a in g
    TO = iw(b1, w1) @ iw(b2, w2) + iw(b2, w1) @ iw(b1, w2)
    Ttr = (w1 @ w2) * 2 * np.einsum('ia,ia->', b1[4:], b2[4:])
    TV = (w1 @ w2) * 2 * np.einsum('ia,ia->', b1[:4], b2[:4])
    return np.array([TO, Ttr, TV])

def run():
    res = {'graphs': [], 'fits': {}}
    for ext in ('perp', 'V'):
        rows = {}; M = []
        for it in range(8):
            p1 = rng.integers(-2, 3, size=4).astype(float); p2 = rng.integers(-2, 3, size=4).astype(float)
            if not p1.any(): p1[1] = 1
            if not p2.any(): p2[2] = 1
            w1 = rng.integers(-2, 3, size=6).astype(float); w2 = rng.integers(-2, 3, size=6).astype(float)
            if it % 3 == 0: w2 = w1.copy()
            b1, b2 = rpol(ext), rpol(ext)
            M.append(structs(b1, b2, w1, w2))
            gl = {}
            gl[f'A{ext}-A{ext}:O+AAAA:Aperp,Aperp'] = -bub_num(b1, b2, p1, p2, w1, w2)          # log: -int_0^1 g
            if ext == 'perp':
                gl['Aperp-Aperp:O+AAA+AAA:Aperp,AV,Aperp'] = -tri_lap(b1, b2, p1, p2, w1, w2, 'AV') / 8
                gl['Aperp-Aperp:O+AAA+AAA:Aperp,Aperp,Aperp'] = -tri_lap(b1, b2, p1, p2, w1, w2, 'Aperp') / 8
            else:
                gl['AV-AV:O+AAA+AAA:Aperp,Aperp,Aperp'] = -tri_lap(b1, b2, p1, p2, w1, w2, 'Aperp') / 8
                gl['AV-AV:O+AAA+AAA:Aperp,AV,Aperp'] = -tri_lap(b1, b2, p1, p2, w1, w2, 'AV') / 8
            for k, v in gl.items():
                rows.setdefault(k, []).append(v)
        M = np.array(M)
        tot = 0
        for k, vals in rows.items():
            x, *_ = np.linalg.lstsq(M, np.array(vals), rcond=None)
            resid = float(np.linalg.norm(M @ x - np.array(vals)))
            ent = {'id': k, 'log_eps': {nm: str(rat(complex(xx).real / eng.cadj)) for nm, xx in zip(['O', 'trace', 'AV_trace'], x)},
                   'resid': resid, 'imag_max': float(np.max(np.abs(np.imag(x))))}
            res['graphs'].append(ent); print(json.dumps(ent), flush=True)
            tot = tot + x
        res['fits'][ext] = {nm: str(rat(complex(xx).real / eng.cadj)) for nm, xx in zip(['O', 'trace', 'AV_trace'], tot)}
        print('TOTAL', ext, res['fits'][ext], flush=True)
    # vacuum: one edge, q = 0 component of omega: 1/2 str(N S HO) * (1/eps)
    w1 = np.array([1., 2, 0, -1, 0, 1]); w2 = np.array([0., 1, 1, 3, -1, 0])
    vac = 0.5 * eng.str_(eng.Nprop(np.zeros(4)) @ S @ HO(w1, w2))
    res['vacuum_eps_inv_over_(w1.w2)dimg'] = float((vac / ((w1 @ w2) * d)).real)
    print('vacuum eps^-1 coefficient / ((w1,w2) dim g) =', res['vacuum_eps_inv_over_(w1.w2)dimg'])
    return res

if __name__ == '__main__':
    r = run()
    json.dump(r, open(f'logs/sector_b_{gname}.json', 'w'), indent=1)
