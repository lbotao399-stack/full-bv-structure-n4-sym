"""QA: su(1|1) sector from the independent harmonic action (qa_harmonic.py): two-site H_12 on {Z = c3 c4, psi = a1 c4},
then graded cyclic spectra of sum_l H_{l,l+1} on single traces tr(Z^{L-M} psi^M), L <= 6."""
import itertools, json, os
import sympy as sp
import qa_harmonic as HA

Z, PSI = HA.Z, HA.PSI[4]
L1 = {'Z': Z, 'psi': PSI}
PAR = {'Z': 0, 'psi': 1}
# two-site density in the letter basis
names = ['Z', 'psi']
dens = {}
for x, y in itertools.product(names, names):
    v = HA.two(L1[x], L1[y])
    h = HA.H12(v)
    # decompose h on the basis two(L1[u], L1[w])
    row = {}
    for u, w in itertools.product(names, names):
        b = HA.two(L1[u], L1[w])
        if not b:
            continue
        (k, cb), = b.items()
        if k in h:
            row[(u, w)] = sp.nsimplify(h[k] / cb)
    # check that h is spanned
    rec = {}
    for (u, w), c in row.items():
        for kk, cc in HA.two(L1[u], L1[w]).items():
            rec[kk] = rec.get(kk, 0) + c * cc
    assert all(sp.simplify(rec.get(k, 0) - h.get(k, 0)) == 0 for k in set(rec) | set(h)), 'not closed'
    (k0, c0), = v.items()
    dens[(x, y)] = {uw: sp.nsimplify(c / 1) for uw, c in row.items()}
# expected 1 - Pi (graded): (x,y) -> (x,y) - (-1)^{|x||y|} (y,x)
expect = {}
for x, y in itertools.product(names, names):
    e = {(x, y): sp.Integer(1)}
    sg = (-1) ** (PAR[x] * PAR[y])
    e[(y, x)] = e.get((y, x), 0) - sg
    expect[(x, y)] = {k: v for k, v in e.items() if v != 0}
ok_density = all({k: v for k, v in dens[p].items() if v != 0} == expect[p] for p in dens)


def words(L, M):
    out = []
    for pos in itertools.combinations(range(L), M):
        w = tuple('psi' if i in pos else 'Z' for i in range(L))
        out.append(w)
    return out


def cyc_sign(w):
    """sign of moving the last letter to the front in a graded trace"""
    last = PAR[w[-1]]
    rest = sum(PAR[x] for x in w[:-1])
    return (-1) ** (last * rest)


def apply_H(vec, L):
    out = {}
    for w, c in vec.items():
        for l in range(L):
            m = (l + 1) % L
            if m == 0:
                # pair (w[L-1], w[0]): rotate so that the pair is adjacent, act, rotate back (graded)
                w2 = (w[-1],) + w[:-1]; s2 = cyc_sign(w)
                for (u, x), cf in dens[(w2[0], w2[1])].items():
                    nw = (u, x) + w2[2:]
                    # rotate back: move first letter to the end
                    back = nw[1:] + (nw[0],)
                    sb = (-1) ** (PAR[nw[0]] * sum(PAR[t] for t in nw[1:]))
                    out[back] = out.get(back, 0) + c * cf * s2 * sb
            else:
                for (u, x), cf in dens[(w[l], w[m])].items():
                    nw = w[:l] + (u, x) + w[m + 1:]
                    out[nw] = out.get(nw, 0) + c * cf
    return {k: sp.nsimplify(v) for k, v in out.items() if v != 0}


def cyclic_basis(L, M):
    seen = set(); basis = []
    for w in words(L, M):
        if w in seen:
            continue
        # graded cyclic orbit with signs
        vec = {}; cur = w; s = 1
        for _ in range(L):
            vec[cur] = vec.get(cur, 0) + s
            s = s * cyc_sign(cur)
            cur = (cur[-1],) + cur[:-1]
            seen.add(cur)
        vec = {k: v for k, v in vec.items() if v != 0}
        if vec:
            basis.append(vec)
    return basis


R = {'two_site_density_equals_1_minus_graded_Pi': ok_density, 'density': {str(k): {str(kk): str(vv) for kk, vv in v.items()} for k, v in dens.items()}, 'spectra': {}}
for L in range(2, 7):
    for M in range(0, L + 1):
        B = cyclic_basis(L, M)
        if not B:
            R['spectra']['L%d_M%d' % (L, M)] = []
            continue
        keys = sorted(set(k for b in B for k in b))
        Bm = sp.Matrix([[b.get(k, 0) for b in B] for k in keys])
        imgs = [apply_H(b, L) for b in B]
        Im = sp.Matrix([[v.get(k, 0) for v in imgs] for k in keys])
        X = (Bm.T * Bm).inv() * (Bm.T * Im)
        assert sp.simplify(Bm * X - Im) == sp.zeros(*Im.shape)
        ev = []
        for e, m in X.eigenvals().items():
            ev += [sp.nsimplify(sp.simplify(e))] * m
        R['spectra']['L%d_M%d' % (L, M)] = sorted([str(x) for x in ev], key=lambda s: float(sp.sympify(s)))
json.dump(R, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'qa_su11.json'), 'w'), indent=1)
print('density 1-Pi:', ok_density)
for k, v in R['spectra'].items():
    print(k, v)
