"""K_v route at finite transfer endpoint T: T d/dT of the one-loop 1PI graph sums in the quadratic sector.

mode 'defect' : arity-one defect along K_v, special edge K_eps - K_0 on the cycle, other edges P(eps,T);
                sectors (A|_V,A|_V), (A|_Wperp,A|_Wperp), (lambda,lambda), (c,A^*)  [degree 0].
mode 'flow'   : arity-one part of W(P(eps,T), I + S^(1)(K_v)), no special edge;
                sectors (A|_V,A^*|_V), (A|_Wperp,A^*|_Wperp), (lambda,lambda^*), (c,c^*)  [degree -1].
K_v: v(y) = |y|^2 b - 2 (b.y) y, b = e_0.  With g_n(P) = graph sum with S^(1)(e_n e^{iPy}), legs b1 at p1, b2 at
-p1-P, the K_v functional on test functions is  int phi1^(p1) [ phi2^(-p1) kap0 + d_a phi2^(-p1) kap1_a
+ d_a d_b phi2^(-p1) kap2_ab ],  C^n_ab = delta_ab b_n - b_a delta_nb - b_b delta_na,
   kap0 = -C^n_ab d_a d_b g_n(0),  kap1_a = 2 C^n_ab d_b g_n(0),  kap2_ab = -C^n_ab g_n(0).
T d/dT: each edge P(eps,T) in turn pinned at t = T with weight T.  Units hbar/(16 pi^2)."""
import sys, json, time, itertools
import numpy as np
from kt_engine import *

SECTORS = {
    'defect': {'AVAV': ('A', 'V', 'A', 'V'), 'AWAW': ('A', 'W', 'A', 'W'), 'll': ('l', None, 'l', None), 'cAs': ('c', None, 'As', 'V')},
    'flow': {'AVAsV': ('A', 'V', 'As', 'V'), 'AWAsW': ('A', 'W', 'As', 'W'), 'lls': ('l', None, 'ls', None), 'ccs': ('c', None, 'cs', None)},
}
P1 = np.array([1., 2., -1., 1.])


def pol(t, which, seed):
    rng = np.random.default_rng(seed)
    x = np.zeros(C.SHAPE[t], complex)
    if t in ('A', 'As'):
        Ms = range(4) if which == 'V' else range(4, 10)
        for M in Ms:
            x[M, 0] = rng.integers(-2, 3) or 1
    elif t in ('l', 'ls'):
        x[:, 0] = rng.integers(-2, 3, 16)
    else:
        x[0] = 1.0
    return x


class Sector:
    def __init__(self, mode, name, p1=P1, seed=0):
        self.mode, self.name, self.p1 = mode, name, np.asarray(p1, float)
        t1, w1, t2, w2 = SECTORS[mode][name]
        self.pol1, self.pol2 = pol(t1, w1, 10 + seed), pol(t2, w2, 20 + seed)
        self.t1, self.t2 = t1, t2
        self.B1 = Insertion('bg', t=t1, pol=self.pol1); self.B2 = Insertion('bg', t=t2, pol=self.pol2)
        self.B12 = Insertion('bg2', legs=[(t1, self.pol1), (t2, self.pol2)]) if (t1 == 'A' and t2 == 'A') else None
        self.lab1 = 0 if PARITY[t1] else None; self.lab2 = 1 if PARITY[t2] else None
        self.KV = {}
        self.structs = self._enumerate()

    def ins(self, u, d1, d2):
        key = tuple(np.round(u, 12))
        if key not in self.KV:
            self.KV[key] = Insertion('S1', u=u)
        z = np.zeros(4); p1 = self.p1
        d1 = np.asarray(d1, float); d2 = np.asarray(d2, float)
        out = {'b1': (self.B1, hd_vec(p1), self.lab1), 'b2': (self.B2, (-p1, -d1, -d2), self.lab2),
               'K': (self.KV[key], (z, d1, d2), 2)}
        if self.B12 is not None:
            out['b12'] = (self.B12, (z, -d1, -d2), None)
        return out

    def _enumerate(self):
        dg1 = np.array([.3, .7, -.2, .5]); dg2 = np.array([-.4, .1, .6, .3]); ug = np.array([.5, -.3, .8, .2])
        ins = self.ins(ug, dg1, dg2)
        first = SLINES if self.mode == 'defect' else PLINES
        out = []
        for order in itertools.permutations(['b1', 'b2', 'K']):
            for L in itertools.product(first, PLINES, PLINES):
                pref = 0.5 if self.mode == 'defect' else 1 / 6
                if Cycle([ins[o] for o in order], list(L), pref).nonzero():
                    out.append((order, L, pref))
        if self.B12 is not None:
            for order in itertools.permutations(['b12', 'K']):
                for L in itertools.product(first, PLINES):
                    pref = 0.5 if self.mode == 'defect' else 1 / 4
                    if Cycle([ins[o] for o in order], list(L), pref).nonzero():
                        out.append((order, L, pref))
        return out

    def graphpolys(self, u, d1, d2):
        ins = self.ins(u, d1, d2)
        return [GraphPoly(Cycle([ins[o] for o in order], list(L), pref), d1, d2) for (order, L, pref) in self.structs]


def tdt_points(gp, eps, T, per=10, nsig=8):
    """Schwinger points for T d/dT: every P edge pinned at T in turn"""
    n = gp.cyc.n; special = gp.cyc.special
    pedges = list(range(1, n)) if special else list(range(n))
    tg, tw = log_grid(eps, T, per)
    sg, sw = lin_grid(0.0, eps, nsig)
    pts, wts = [], []
    for e in pedges:
        others = [f for f in pedges if f != e]
        grids = []
        if special:
            grids.append((0, sg, sw))
        for f in others:
            grids.append((f, tg, tw))
        if not grids:
            t = np.zeros((1, n)); t[0, e] = T
            pts.append(t); wts.append(np.array([T]))
            continue
        idx = [g[0] for g in grids]
        mesh = np.meshgrid(*[g[1] for g in grids], indexing='ij')
        wm = np.meshgrid(*[g[2] for g in grids], indexing='ij')
        npt = mesh[0].size
        t = np.zeros((npt, n)); t[:, e] = T
        w = np.full(npt, T)
        for j, f in enumerate(idx):
            t[:, f] = mesh[j].ravel(); w = w * wm[j].ravel()
        pts.append(t); wts.append(w)
    return np.concatenate(pts), np.concatenate(wts)


RUNS_B0 = [('R', a) for a in range(4)] + [('Q', n) for n in range(1, 4)]   # (R,0) = (Q,0)


def run_dirs(n_or_a, kind):
    if kind == 'R':
        return E4[0], E4[n_or_a], E4[n_or_a]          # u = b = e0, d1 = d2 = e_a
    return E4[n_or_a], E4[0], E4[n_or_a]              # u = e_n, d1 = b = e0, d2 = e_n


def kernels(J):
    """J[(kind, idx)] = hyper-dual [g, d1 g, d2 g, d1 d2 g] (sum or one graph); b = e0"""
    R = {a: J[('R', a)] for a in range(4)}
    Q = {0: J[('R', 0)]}
    for n in range(1, 4):
        Q[n] = J[('Q', n)]
    g = {n: Q[n][0] for n in range(4)}
    d0g = {n: Q[n][1] for n in range(4)}            # d_0 g_n
    dng = {n: Q[n][2] for n in range(4)}            # d_n g_n
    d0ng = {n: Q[n][3] for n in range(4)}           # d_0 d_n g_n
    dag0 = {a: R[a][1] for a in range(4)}           # d_a g_0
    daag0 = {a: R[a][3] for a in range(4)}          # d_a^2 g_0
    kap0 = -(sum(daag0.values()) - 2 * sum(d0ng.values()))
    div = sum(dng.values())
    kap1 = [2 * (dag0[a] - (a == 0) * div - d0g[a]) for a in range(4)]
    kap2 = [[-((a == b) * g[0] - (a == 0) * g[b] - (b == 0) * g[a]) for b in range(4)] for a in range(4)]
    return complex(kap0), [complex(x) for x in kap1], [[complex(x) for x in r] for r in kap2]


def compute(mode, name, eps_list, T_list, per=10, nsig=8, seed=0, p1=P1, log=print):
    t0 = time.time()
    sec = Sector(mode, name, p1=p1, seed=seed)
    log('%s %s: %d nonzero cycle structures (%.1fs)' % (mode, name, len(sec.structs), time.time() - t0))
    res = {}       # (eps,T) -> struct index -> run -> HD
    for run in RUNS_B0:
        u, d1, d2 = run_dirs(run[1], run[0])
        gps = sec.graphpolys(u, d1, d2)
        for (eps, T) in itertools.product(eps_list, T_list):
            for si, gp in enumerate(gps):
                t, w = tdt_points(gp, eps, T, per, nsig)
                res.setdefault((eps, T), {}).setdefault(si, {})[run] = gp.integrate(t, w)
        log('  run %s done (%.1fs)' % (str(run), time.time() - t0))
    out = {'mode': mode, 'sector': name, 'p1': list(p1), 'pol1': np.asarray(sec.pol1).real.tolist(),
           'pol2': np.asarray(sec.pol2).real.tolist(), 'b': [1, 0, 0, 0],
           'structures': [[list(o), list(L), pr] for (o, L, pr) in sec.structs], 'values': {}}
    for (eps, T), per_s in res.items():
        tot = {run: sum(per_s[si][run] for si in per_s) for run in RUNS_B0}
        k0, k1, k2 = kernels(tot)
        pg = {}
        for si in per_s:
            a0, a1_, a2_ = kernels(per_s[si])
            pg[si] = {'kap0': [a0.real, a0.imag], 'max|kap1|': max(abs(x) for x in a1_), 'max|kap2|': max(abs(x) for r in a2_ for x in r)}
        out['values']['eps=%g,T=%g' % (eps, T)] = {
            'eps': eps, 'T': T, 'TdT_kap0': [k0.real, k0.imag], 'TdT_kap1': [[x.real, x.imag] for x in k1],
            'TdT_kap2': [[[x.real, x.imag] for x in r] for r in k2], 'per_structure_TdT': pg}
    out['seconds'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    mode, name = sys.argv[1], sys.argv[2]
    tag = sys.argv[3] if len(sys.argv) > 3 else 'main'
    if mode == 'defect':
        eps_list = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
        T_list = [1.0, 10.0, 1e2, 1e3, 1e4, 1e5, 1e6]
    else:
        eps_list = [1e-2, 1e-4, 1e-6]
        T_list = [1.0, 1e2, 1e4, 1e6]
    per, nsig = (10, 8) if tag == 'main' else (15, 12)
    if tag != 'main':
        eps_list = eps_list[::2]; T_list = T_list[::2]
    out = compute(mode, name, eps_list, T_list, per, nsig, log=lambda s: print(s, flush=True))
    json.dump(out, open('kv_%s_%s_%s.json' % (mode, name, tag), 'w'), indent=1)
    print('saved', out['seconds'], 's')
