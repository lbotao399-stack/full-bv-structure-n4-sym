"""Item (a) graph by graph: masks on the propagator select the internal lines.  For every graph:
integrand after the Gaussian loop integral = (16pi^2)^{-1} [ T^{-2} g(u) + T^{-3} c ] e^{-T u(1-u) p^2}, u = t_2/T,
g and c fitted (exactly, rational) on index-free structures; singular parts: eps^{-1}: c/4 (+ seagull h);
log eps: -int_0^1 g + p^2 c/6.  Units of the reported numbers: hbar c_adj/(16 pi^2)."""
import json, sys
import numpy as np
from fractions import Fraction
from oneloop import Engine, PARITY
import check_vertices as cv
from sector_a import Bcross, rat, E4
gname = sys.argv[1] if len(sys.argv) > 1 else 'su2'
eng = Engine(gname); d = eng.d
cv.d = d; cv.SH = {'c': (d,), 'A': (10, d), 'l': (16, d), 'As': (10, d), 'ls': (16, d), 'cs': (d,)}
rng = np.random.default_rng(5)
o = eng.off
LINES = {'AV': [o['A'] + M * d + a for M in range(4) for a in range(d)],
         'Aperp': [o['A'] + M * d + a for M in range(4, 10) for a in range(d)],
         'gh': list(range(o['c'], o['c'] + d)) + [o['As'] + M * d + a for M in range(4) for a in range(d)],
         'l': list(range(o['l'], o['l'] + 16 * d))}
def Nmask(k, lines):
    N = eng.Nprop(k); m = np.zeros(eng.n, bool)
    for L in lines:
        m[LINES[L]] = True
    return N * np.outer(m, m)
S = np.diag(eng.S)

def bubble_data(b1, b2, p, lines1, lines2):
    """returns (g-coefficients [1,u,u^2], c) of 1/4[str(X1X2)+str(X2X1)] with propagator masks"""
    t1, pol1, _ = b1; t2, pol2, _ = b2
    sgn21 = (-1) ** (PARITY[t1] * PARITY[t2])
    def F12(k):
        h1 = eng.hessian([(t1, pol1, p)], -k, k - p); h2 = eng.hessian([(t2, pol2, -p)], p - k, k)
        return eng.str_(Nmask(k, lines1) @ S @ h1 @ Nmask(k - p, lines2) @ S @ h2)
    def F21(k):
        h2 = eng.hessian([(t2, pol2, -p)], -k, k + p); h1 = eng.hessian([(t1, pol1, p)], -k - p, k)
        return sgn21 * eng.str_(Nmask(k, lines2) @ S @ h2 @ Nmask(k + p, lines1) @ S @ h1)
    gc = np.zeros(3, complex); c = 0
    for F, sg in ((F12, 1), (F21, -1)):
        v0, vh, v1 = F(0 * p), F(sg * .5 * p), F(sg * p)
        # quadratic through (0,v0),(1/2,vh),(1,v1)
        a0 = v0; a2 = 2 * (v1 - 2 * vh + v0); a1 = v1 - v0 - a2
        gc += 0.25 * np.array([a0, a1, a2])
        lap = sum(F(e) + F(-e) - 2 * v0 for e in E4)
        c += 0.25 * lap / 4
    return gc, c

def seagull_data(b1, b2, p, lines):
    k = np.array([0.3, -0.7, 1.1, 0.2])
    h4 = eng.hessian([(b1[0], b1[1], p), (b2[0], b2[1], -p)], -k, k)
    return 0.5 * eng.str_(Nmask(k, lines) @ S @ h4)

def rpol(t, which=None):
    x = rng.integers(-2, 3, size=eng.shape[t]).astype(complex)
    if which == 'V': x[4:] = 0
    if which == 'perp': x[:4] = 0
    return x

GRAPHS = {
  'AV-AV': (('A', 'V'), ('A', 'V'), ['B1', 'B2', 'B3'], [
      ('AV-AV:AAA+AAA:AV,AV', 'b', ['AV'], ['AV']), ('AV-AV:AAA+AAA:Aperp,Aperp', 'b', ['Aperp'], ['Aperp']),
      ('AV-AV:AcAs+AcAs:c,As', 'b', ['gh'], ['gh']), ('AV-AV:llA+llA:l,l', 'b', ['l'], ['l']),
      ('AV-AV:AAAA:AV', 's', ['AV'], None), ('AV-AV:AAAA:Aperp', 's', ['Aperp'], None)]),
  'Aperp-Aperp': (('A', 'perp'), ('A', 'perp'), ['B4', 'B5'], [
      ('Aperp-Aperp:AAA+AAA:AV,Aperp', 'b', ['AV', 'Aperp'], ['AV', 'Aperp']),
      ('Aperp-Aperp:AcAs+AcAs:c,As', 'b', ['gh'], ['gh']), ('Aperp-Aperp:llA+llA:l,l', 'b', ['l'], ['l']),
      ('Aperp-Aperp:AAAA:AV', 's', ['AV'], None), ('Aperp-Aperp:AAAA:Aperp', 's', ['Aperp'], None)]),
  'l-l': (('l', None), ('l', None), ['B6'], [
      ('l-l:llA+llA:l,AV', 'b', ['l', 'AV'], ['l', 'AV']), ('l-l:llA+llA:l,Aperp', 'b', ['l', 'Aperp'], ['l', 'Aperp'])]),
  'As-c': (('As', None), ('c', None), ['B7'], [
      ('As-c:AcAs+AcAs:AV,c-As', 'b', ['AV', 'gh'], ['AV', 'gh']), ('As-c:AcAs+AcAs:Aperp,c-As', 'b', ['Aperp', 'gh'], ['Aperp', 'gh'])]),
}

def run():
    out = []
    totals = {}
    for sname, (t1, t2, Bn, glist) in GRAPHS.items():
        samples = []
        for it in range(6):
            p = rng.integers(-3, 4, size=4).astype(float)
            if not p.any(): p[0] = 1
            b1 = (t1[0], rpol(*t1), p); b2 = (t2[0], rpol(*t2), -p)
            Bc = np.array([Bcross(nm, b1, b2, p) for nm in Bn])
            row = {}
            for (gid, kind, l1, l2) in glist:
                if kind == 'b':
                    gc, c = bubble_data(b1, b2, p, l1, l2)
                    s2 = p @ p
                    row[gid] = (c / 4, -(gc[0] + gc[1] / 2 + gc[2] / 3) + s2 * c / 6, gc, c)
                else:
                    h = seagull_data(b1, b2, p, l1)
                    row[gid] = (h, 0, None, None)
            samples.append((Bc, row))
        M = np.array([smp[0] for smp in samples])
        for (gid, kind, l1, l2) in glist:
            ent = {'id': gid, 'sector': sname, 'kind': 'bubble' if kind == 'b' else 'seagull'}
            for idx, key in ((0, 'eps_inv'), (1, 'log_eps')):
                rhs = np.array([smp[1][gid][idx] for smp in samples])
                x, *_ = np.linalg.lstsq(M, rhs, rcond=None)
                ent[key] = {nm: str(rat(complex(xx).real / eng.cadj)) for nm, xx in zip(Bn, x)}
                ent[key + '_imag_max'] = float(np.max(np.abs(np.imag(x))))
                ent[key + '_resid'] = float(np.linalg.norm(M @ x - rhs))
                for nm, xx in zip(Bn, x):
                    totals.setdefault((sname, key, nm), 0)
                    totals[(sname, key, nm)] += complex(xx).real / eng.cadj
            out.append(ent)
            print(json.dumps(ent), flush=True)
    tot = {f'{a}|{b}|{c}': str(rat(v)) for (a, b, c), v in totals.items()}
    print('TOTALS', json.dumps(tot))
    return out, tot

if __name__ == '__main__':
    out, tot = run()
    json.dump({'graphs': out, 'totals': tot, 'c_adj': eng.cadj}, open(f'logs/graphs_a_{gname}.json', 'w'), indent=1)
