"""Arity two along (eta1, eta2) in Sigma^* x Sigma^*, quadratic sector of degree -1, at finite transfer endpoint T.

Generating function along the even ghosts t1, t2 of eta1, eta2 (ch02 def:c2-cochains; twisted master equation
coefficient -1/2 of S^(2)(c,c), engine plan check E0a):  J = I + t1 S^(1)(eta1) + t2 S^(1)(eta2) - t1 t2 S^(2)(eta1,eta2);
the arity-two component is the coefficient of t1 t2.  Killing spinors eps_a(y) = rho(y) eta_a = sum_c y_c s[c] eta_a
(ch02 def:c2-odd), y_c -> -i d_{P_c} e^{iPy}:
   S^(1)(eta) = sum_c (-i d_{P_c}) S^(1)_lin[s[c] eta e^{iPy}] |_{P=0} + S^(1)_const[eta],
   S^(1)_lin[z e^{iPy}] = int -(Gamma(eps,lambda),A^*) + 1/2 (rho(dA) eps, lambda^*) + cubic,  S^(1)_const = -(rho(A|_W) eta, lambda^*).
Part 'linlin': both vertices position-linear (quadratic parts) and S^(2) (bilinear in the two Killing spinors):
   kap0 = -sum_{a,b} d_{P1 a} d_{P2 b} G[s[a] eta1, s[b] eta2](0)   (the factor (-i)^2 included),
   lower jets kap1_a = sum_b d_{P2 b} G_ab(0), kap1'_b = sum_a d_{P1 a} G_ab(0), kap2_ab = -G_ab(0)
(B^(2)(phi1,phi2) = int phi1^(p1)[phi2^(-p1) kap0 + d_a phi2^ kap1_a + d_b phi2^ kap1'_b + d_a d_b phi2^ kap2_ab]).
mode 'defect': special edge K_eps - K_0 first on the cycle; mode 'flow': no special edge.  T d/dT as in kv_driver."""
import sys, json, time, itertools
import numpy as np
from kt_engine import *
import kv_driver as KD

rng0 = np.random.default_rng(777)
ETA1 = rng0.integers(-2, 3, 16) + 1j * rng0.integers(-2, 3, 16)
ETA2 = rng0.integers(-2, 3, 16) + 1j * rng0.integers(-2, 3, 16)
V1 = np.array([ETA1 @ SARR[M] @ ETA2 for M in range(4)])        # Gamma(eta1, eta2)|_V  (pairing on Sigma^*)
# degree -1 sectors (defect) and degree -2 sectors (flow W^(2), two S^(1)(eta) vertices of degree -1)
ETA_SECTORS = {'AWAsW': ('A', 'W', 'As', 'W'), 'AVAsV': ('A', 'V', 'As', 'V'), 'lls': ('l', None, 'ls', None), 'ccs': ('c', None, 'cs', None),
               'AVcs': ('A', 'V', 'cs', None), 'AWcs': ('A', 'W', 'cs', None)}


class EtaSector:
    def __init__(self, mode, name, p1=KD.P1, seed=0, with_S2=True):
        self.mode, self.name, self.p1 = mode, name, np.asarray(p1, float)
        t1, w1, t2, w2 = ETA_SECTORS[name]
        self.pol1, self.pol2 = KD.pol(t1, w1, 10 + seed), KD.pol(t2, w2, 20 + seed)
        self.B1 = Insertion('bg', t=t1, pol=self.pol1); self.B2 = Insertion('bg', t=t2, pol=self.pol2)
        self.lab1 = 0 if PARITY[t1] else None; self.lab2 = 1 if PARITY[t2] else None
        self.with_S2 = with_S2
        self.cacheE = {}
        self.structs = self._enumerate()

    def eins(self, a, b):
        """insertions for the run (a, b): eta1 linear part with z1 = s[a] eta1 at P1 = e1 e_a, eta2 with z2 = s[b] eta2"""
        key = (a, b)
        if key not in self.cacheE:
            z1 = SARR[a] @ ETA1; z2 = SARR[b] @ ETA2
            E1 = InsertionGen(lambda kj, kl, P, z=z1: s1eta_hessian(z, kj, kl, P))
            E2 = InsertionGen(lambda kj, kl, P, z=z2: s1eta_hessian(z, kj, kl, P))
            S2 = InsertionGen(lambda kj, kl, P, z1=z1, z2=z2: -s2_hessian(z1, z2, kj, kl, P))
            self.cacheE[key] = (E1, E2, S2)
        E1, E2, S2 = self.cacheE[key]
        d1, d2 = E4[a], E4[b]; z = np.zeros(4); p1 = self.p1
        return {'b1': (self.B1, hd_vec(p1), self.lab1), 'b2': (self.B2, (-p1, -d1, -d2), self.lab2),
                'e1': (E1, (z, d1, None), None), 'e2': (E2, (z, None, d2), None), 'S2': (S2, (z, d1, d2), None)}

    def _enumerate(self):
        first = SLINES if self.mode == 'defect' else PLINES
        out = []
        insg = self.eins(1, 2)
        for order in itertools.permutations(['b1', 'b2', 'e1', 'e2']):
            for L in itertools.product(first, PLINES, PLINES, PLINES):
                pref = 0.5 if self.mode == 'defect' else 1 / 8
                c = Cycle([insg[o] for o in order], list(L), pref)
                if c.nonzero():
                    out.append((order, L, pref))
        if self.with_S2:
            for order in itertools.permutations(['b1', 'b2', 'S2']):
                for L in itertools.product(first, PLINES, PLINES):
                    pref = 0.5 if self.mode == 'defect' else 1 / 6
                    c = Cycle([insg[o] for o in order], list(L), pref)
                    if c.nonzero():
                        out.append((order, L, pref))
        return out


def classify(order, L):
    """'chain' if e1, e2 are adjacent, joined by a lambda line, both other lines c-A^* (floating eta-eta chain)"""
    if 'S2' in order:
        return 'S2'
    n = len(order); i1, i2 = order.index('e1'), order.index('e2')
    for (ia, ib) in ((i1, i2), (i2, i1)):
        if (ia + 1) % n == ib:
            mid = L[ib]; before = L[ia]; after = L[(ib + 1) % n]
            if mid == 'l' and before == 'gh' and after == 'gh':
                return 'chain'
    return 'other'


def compute(mode, name, eps_list, T_list, per=10, nsig=8, log=print, with_S2=True):
    t0 = time.time()
    sec = EtaSector(mode, name, with_S2=with_S2)
    cls = [classify(o, L) for (o, L, _) in sec.structs]
    log('%s %s: %d cycle structures (chain %d, S2 %d, other %d), |v1| = %.3f (%.1fs)' % (
        mode, name, len(sec.structs), cls.count('chain'), cls.count('S2'), cls.count('other'), np.abs(V1).max(), time.time() - t0))
    res = {}
    for a in range(4):
        for b in range(4):
            ins = sec.eins(a, b)
            for si, (order, L, pref) in enumerate(sec.structs):
                gp = GraphPoly2(Cycle([ins[o] for o in order], list(L), pref), E4[a], E4[b])
                for (eps, T) in itertools.product(eps_list, T_list):
                    t, w = KD.tdt_points(gp, eps, T, per, nsig)
                    res.setdefault((eps, T), {}).setdefault(si, {})[(a, b)] = gp.integrate(t, w)
            log('  run (%d,%d) done (%.1fs)' % (a, b, time.time() - t0))
    out = {'mode': mode, 'sector': name, 'p1': list(sec.p1), 'eta1': [[x.real, x.imag] for x in ETA1],
           'eta2': [[x.real, x.imag] for x in ETA2], 'Gamma(eta1,eta2)|_V': [[x.real, x.imag] for x in V1],
           'structures': [[list(o), list(L), pr, c] for (o, L, pr), c in zip(sec.structs, cls)], 'values': {}}
    def kern(J):
        k0 = -sum(J[(a, b)][3] for a in range(4) for b in range(4))
        k1 = [sum(J[(a, b)][2] for b in range(4)) for a in range(4)]
        k1p = [sum(J[(a, b)][1] for a in range(4)) for b in range(4)]
        k2 = [[-J[(a, b)][0] for b in range(4)] for a in range(4)]
        return complex(k0), k1, k1p, k2
    for (eps, T), per_s in res.items():
        tot = {ab: sum(per_s[si][ab] for si in per_s) for ab in per_s[0]}
        k0, k1, k1p, k2 = kern(tot)
        bycls = {}
        for si in per_s:
            c = cls[si]
            bycls[c] = bycls.get(c, 0) + kern(per_s[si])[0]
        pg = {str(si): [kern(per_s[si])[0].real, kern(per_s[si])[0].imag] for si in per_s}
        out['values']['eps=%g,T=%g' % (eps, T)] = {
            'eps': eps, 'T': T, 'TdT_kap0': [k0.real, k0.imag], 'TdT_kap0_by_class': {c: [v.real, v.imag] for c, v in bycls.items()},
            'max|TdT kap1|': float(max(abs(complex(x)) for x in k1 + k1p)), 'max|TdT kap2|': float(max(abs(complex(x)) for r in k2 for x in r)),
            'per_structure_TdT_kap0': pg}
    out['seconds'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    mode, name = sys.argv[1], sys.argv[2]
    eps_list = [1e-2, 1e-4, 1e-6]
    T_list = [1.0, 1e2, 1e4, 1e6]
    out = compute(mode, name, eps_list, T_list, per=8, nsig=6, log=lambda s: print(s, flush=True))
    json.dump(out, open('eta_%s_%s.json' % (mode, name), 'w'), indent=1)
    print('saved', out['seconds'], 's')
