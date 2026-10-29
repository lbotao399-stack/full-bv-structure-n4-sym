"""Route (ii): UV-boundary (special-edge) graphs.  d/dlambda P(eps/lambda^2,infty) = 2 eps (Q^GF x 1) K_eps, so the
dilatation acting on an internal line is the special edge 2 eps (Q^GF x 1)K_eps at t_e = eps.  For each two-edge graph
we fit the Gaussian-integrated integrand T^{-2} g(u) + T^{-3} c on index-free structures and evaluate each face
(t_1 = eps or t_2 = eps) exactly: face-1 log part = -int_{1/2}^1 g + p^2 c/12, face-2 = -int_0^{1/2} g + p^2 c/12,
eps^{-1} part per face c/8.  Units hbar c_adj/(16 pi^2)."""
import json, sys
import numpy as np
import sympy as sp
import graphs_a as G
from sector_a import Bcross, rat
from cliff import sb
eng = G.eng; d = eng.d
u = sp.symbols('u')

def structs(sector, b1, b2, p):
    a1, a2 = b1[1], b2[1]
    if sector == 'AV-AV':
        return [p @ p * np.sum(a1[:4] * a2[:4]), np.einsum('m,ma->a', p, a1[:4]) @ np.einsum('m,ma->a', p, a2[:4])], [np.sum(a1[:4] * a2[:4])]
    if sector == 'Aperp-Aperp':
        return [p @ p * np.sum(a1[4:] * a2[4:])], [np.sum(a1[4:] * a2[4:])]
    if sector == 'l-l':
        rp = sum(p[m] * sb[m] for m in range(4))
        return [np.einsum('ia,ij,ja->', a1, rp, a2)], []
    if sector == 'As-c':
        return [np.einsum('m,ma,a->', p, a1[:4], a2)], []

def run():
    out = {}
    for sname, (t1, t2, Bn, glist) in G.GRAPHS.items():
        rows = []
        for it in range(6):
            p = G.rng.integers(-3, 4, size=4).astype(float)
            if not p.any(): p[0] = 1
            b1 = (t1[0], G.rpol(*t1), p); b2 = (t2[0], G.rpol(*t2), -p)
            gs, cs = structs(sname, b1, b2, p)
            Bc = [Bcross(nm, b1, b2, p) for nm in Bn]
            data = {}
            for (gid, kind, l1, l2) in glist:
                if kind == 'b':
                    data[gid] = G.bubble_data(b1, b2, p, l1, l2)
            rows.append((gs, cs, Bc, data))
        # conversion: B crosses in terms of g-structures (log parts are quadratic in p) and c-structures
        Mg = np.array([r[0] for r in rows]); Mc = np.array([r[1] for r in rows]) if rows[0][1] else None
        MB = np.array([r[2] for r in rows])
        for (gid, kind, l1, l2) in glist:
            if kind != 'b':
                continue
            ent = {}
            # fit g coefficients (1,u,u^2) on g-structures, c on c-structures
            gpoly = []
            for j in range(3):
                x, *_ = np.linalg.lstsq(Mg, np.array([r[3][gid][0][j] for r in rows]), rcond=None)
                gpoly.append([rat(complex(xx).real / eng.cadj) if abs(complex(xx).imag) < 1e-9 else complex(xx) / eng.cadj for xx in x])
            if Mc is not None:
                xc, *_ = np.linalg.lstsq(Mc, np.array([r[3][gid][1] for r in rows]), rcond=None)
                cval = [rat(complex(xx).real / eng.cadj) for xx in xc]
            else:
                cval = []
            ent['g(u)'] = {f'struct{k}': str(sum(sp.nsimplify(gpoly[j][k]) * u ** j for j in range(3))) for k in range(len(gpoly[0]))}
            ent['c'] = {f'cstruct{k}': str(v) for k, v in enumerate(cval)}
            # faces
            for face, (lo, hi) in (('t1=eps', (sp.Rational(1, 2), 1)), ('t2=eps', (0, sp.Rational(1, 2)))):
                vals = {}
                for k in range(len(gpoly[0])):
                    gk = sum(sp.nsimplify(gpoly[j][k]) * u ** j for j in range(3))
                    vals[f'struct{k}'] = str(-sp.integrate(gk, (u, lo, hi)))
                vals['p^2*cstruct'] = [str(sp.nsimplify(cv) / 12) for cv in cval]
                vals['eps_inv_cstruct'] = [str(sp.nsimplify(cv) / 8) for cv in cval]
                ent[face] = vals
            out[gid] = ent
            print(gid, json.dumps(ent), flush=True)
        # structure dictionary: B crosses in the structure basis (for conversion in the report)
        allS = np.array([list(r[0]) + list(r[1]) for r in rows])
        conv = {}
        for j, nm in enumerate(Bn):
            x, *_ = np.linalg.lstsq(allS, MB[:, j], rcond=None)
            conv[nm] = [str(rat(complex(xx).real)) + ('+i*' + str(rat(complex(xx).imag)) if abs(complex(xx).imag) > 1e-9 else '') for xx in x]
        out['conversion_' + sname] = conv
        print('conversion', sname, conv, flush=True)
    return out

if __name__ == '__main__':
    r = run()
    json.dump(r, open('logs/faces_su2.json', 'w'), indent=1)
