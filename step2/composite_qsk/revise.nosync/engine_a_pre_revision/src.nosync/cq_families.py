"""Engine A, part C: families F0, F2 (su(1|1) cyclic spectra) and the on-shell universality of the planar densities.

universality(): for the field pairs (X_a, X_b) with field letters A|_perp, lambda_+, lambda_- (gauge-invariant letters),
random letter vectors zeta_a, zeta_b and random on-shell legs (null Gaussian momenta; A|_perp any polarisation,
lambda legs in the Weyl kernel of the right chirality), the planar density of an adjacent pair (factor 1/2) is compared with
   d0 = -4/3 (X_a X_b)  [same order],  0 for the exchanged word;   dP_beta = -1/6 i beta.(q_1 + q_2) (X_a X_b)
(Conjecture P7 (b), (c) at the level of the two-site density; the letter-on-leg values zeta.pol form the normalisation).
su11_spectra(): graded cyclic spectra of D_1 = 4 lc^natural + Xi_Delta from the computed two-site densities (log kernel
from cq_gval.logkernel, composite from cq_su11), E = D_1 / (-8u), against ha_su11.json (harmonic action [B03c])."""
import itertools, json, os, time
from fractions import Fraction
import numpy as np
import sympy as sp
import cq_core as Q, cq_gval as GV, cq_eval as V, cq_su11 as S11, cq_fsl2 as F2, cq_ops as CO
import ea_core as E

ONE = {(0, 0, 0, 0): 1}


def _onshell(t, q, rng):
    if t == 'AW':
        return np.array([complex(int(a), int(b)) for a, b in zip(rng.integers(-2, 3, 6), rng.integers(-2, 3, 6))])
    B = np.column_stack(E.BASIS['lamP' if t == 'lamP' else 'lamM']).astype(complex)
    W = sum(complex(q[mu]) * np.asarray(E.SGB[mu]) for mu in range(4))
    M = W @ B
    Ms = sp.Matrix(M.shape[0], M.shape[1], lambda i, j: sp.Integer(int(round(M[i, j].real))) + sp.I * int(round(M[i, j].imag)))
    ns = Ms.nullspace()
    u = np.zeros(16, complex)
    for v in ns:
        den = sp.ilcm(*[sp.fraction(sp.nsimplify(sp.re(x)))[1] for x in v] + [sp.fraction(sp.nsimplify(sp.im(x)))[1] for x in v])
        x = np.array([complex(sp.N(c * den)) for c in v])
        u = u + int(rng.integers(1, 4)) * (B @ x)
    u = np.array([complex(round(z.real), round(z.imag)) for z in u])
    assert np.allclose(W @ u, 0)
    return u


def _rv(n, rng):
    return np.array([complex(int(a), int(b)) for a, b in zip(rng.integers(-2, 3, n), rng.integers(-2, 3, n))])


ROWN = {'AW': 6, 'lamP': 16, 'lamM': 16}


def _proj(t, z):
    """letter vector of a definite chirality: the plain dot product pairs S_+ (x) W with S_+ (x) W only (G5S symmetric),
    so a lambda_+ letter is (1 + G5S) zeta, a lambda_- letter (1 - G5S) zeta"""
    if t == 'lamP':
        return (np.eye(16, dtype=int) + np.asarray(E.G5S)) @ z
    if t == 'lamM':
        return (np.eye(16, dtype=int) - np.asarray(E.G5S)) @ z
    return z
FIELD_PAIRS = [('AW', 'AW'), ('AW', 'lamP'), ('AW', 'lamM'), ('lamP', 'AW'), ('lamP', 'lamP'), ('lamM', 'lamM'),
               ('lamP', 'lamM'), ('lamM', 'lamP')]


def universality(pairs=FIELD_PAIRS, npts=2, seed=5, beta=(1, 2, -1, 3), log=print):
    rng = np.random.default_rng(seed)
    nulls = [q for q in S11.NULLS]
    out = {}
    for ta, tb in pairs:
        rows = []
        ok = True
        for it in range(npts):
            za, zb = _proj(ta, _rv(ROWN[ta], rng)), _proj(tb, _rv(ROWN[tb], rng))
            i, j = rng.choice(len(nulls), 2, replace=False)
            q1, q2 = nulls[i], nulls[j]
            for word in ((ta, tb), (tb, ta)):
                if word == (tb, ta) and ta == tb:
                    continue
                u1, u2 = _onshell(word[0], q1, rng), _onshell(word[1], q2, rng)
                legs = [(word[0], u1, q1), (word[1], u2, q2)]
                Xa, Xb = (ta, za, ONE), (tb, zb, ONE)
                same = word == (ta, tb)
                norm = complex(np.dot(za, u1) * np.dot(zb, u2)) if same else complex(np.dot(zb, u1) * np.dot(za, u2))
                d0 = GV.combine(GV.kernel(Xa, Xb, legs, 'D0'), 1)
                dP = GV.combine(GV.kernel(Xa, Xb, legs, ('Pv', tuple(beta))), 1)
                d0s = F2.symv(d0) / 2; dPs = F2.symv(dP) / 2
                t0 = (sp.Rational(-4, 3) * F2.symc(norm)) if same else 0
                tP = (sp.Rational(-1, 6) * F2.symc(1j * np.dot(beta, q1 + q2)) * F2.symc(norm)) if same else 0
                r0 = sp.simplify(sp.expand(d0s - t0)); rP = sp.simplify(sp.expand(dPs - tP))
                ok &= (r0 == 0 and rP == 0)
                rows.append({'word': list(word), 'd0': str(sp.expand(d0s)), 'target_d0': str(t0), 'dP': str(sp.expand(dPs)),
                             'target_dP': str(tP), 'pass': bool(r0 == 0 and rP == 0)})
        out['%s|%s' % (ta, tb)] = {'rows': rows, 'pass': bool(ok)}
        log('  universality %s|%s: %s' % (ta, tb, ok))
    return out


# ------------------------------------------------------------------------------------------------ su(1|1) spectra
def su11_two_site(log=print):
    """two-site densities in units u on the words over {Z, psi}: log (planar factor 1/2), composite D0"""
    Zl, Pl = S11.LET['Z'], S11.LET['psi']
    good = [q for q in S11.NULLS if abs(S11.lol('psi', q)) > 1e-9]
    q1, q2 = good[0], good[1]
    lets = {'Z': Zl, 'psi': Pl}
    dens = {'log': {}, 'D0': {}}
    for pa, pb in itertools.product(['Z', 'psi'], repeat=2):
        words = [(pa, pb)] if pa == pb else [(pa, pb), (pb, pa)]
        for w in words:
            legs = [S11.leg(w[0], q1), S11.leg(w[1], q2)]
            norm = S11.lol(w[0], q1) * S11.lol(w[1], q2)
            lv = GV.logkernel(lets[pa], lets[pb], legs)
            cv = GV.combine(GV.kernel(lets[pa], lets[pb], legs, 'D0'), 1)
            dens['log'][(pa, pb, w)] = sp.nsimplify(F2.symv(lv) / 2 / F2.symc(norm))
            dens['D0'][(pa, pb, w)] = sp.nsimplify(F2.symv(cv) / 2 / F2.symc(norm))
    return dens


def graded_cyclic_spectrum(h2, L, M):
    """h2[(x,y)] = {(x',y'): coeff}: two-site operator; H = sum over the L cyclic links on graded cyclic states of
    tr(X_1..X_L) with M letters psi.  Returns the eigenvalues (sympy) on the graded-cyclic subspace."""
    par = {'Z': 0, 'psi': 1}
    words = [w for w in itertools.product(['Z', 'psi'], repeat=L) if w.count('psi') == M]
    idx = {w: i for i, w in enumerate(words)}
    n = len(words)
    def shift(w):            # tr(X_1 X_2 .. X_L) = (-1)^{|X_1|(|X_2|+..)} tr(X_2 .. X_L X_1)
        s = (-1) ** (par[w[0]] * sum(par[x] for x in w[1:]))
        return s, w[1:] + w[:1]
    T = sp.zeros(n, n)
    for w in words:
        s, w2 = shift(w)
        T[idx[w2], idx[w]] += s
    H0 = sp.zeros(n, n)                       # link (1,2) acting in place
    for w in words:
        for (o1, o2), c in h2.get((w[0], w[1]), {}).items():
            w2 = (o1, o2) + w[2:]
            if w2 in idx:
                H0[idx[w2], idx[w]] += c
    Hs = sp.zeros(n, n)
    Tk = sp.eye(n); Tinv = T.inv()
    for l in range(L):
        Hs += Tk.inv() * H0 * Tk
        Tk = T * Tk
    # graded-cyclic subspace: T v = v
    ns = (T - sp.eye(n)).nullspace()
    if not ns:
        return []
    B = sp.Matrix.hstack(*ns)
    Hr = (B.H * B).inv() * B.H * Hs * B
    ev = Hr.eigenvals()
    return sorted([sp.nsimplify(sp.simplify(e)) for e, m in ev.items() for _ in range(m)], key=lambda z: (float(sp.re(z)), float(sp.im(z))))


def su11_spectra(Lmax=6, theta=sp.Rational(-1, 3), log=print):
    dens = su11_two_site(log)
    wF = {'Z': sp.Integer(1), 'psi': sp.Rational(3, 2)}
    cphi = {'Z': sp.Rational(2, 3), 'psi': sp.Rational(1, 3)}
    nF = {'Z': 2, 'psi': 4}
    # D_1 two-site density: 4 (log - u n_F/2 per letter share) + composite + one-body share (Lemma lem:cc-onebody)
    h2 = {}; hlc = {}; hxi = {}
    for (pa, pb, w), lv in dens['log'].items():
        cv = dens['D0'][(pa, pb, w)]
        d = h2.setdefault((pa, pb), {}); dl = hlc.setdefault((pa, pb), {}); dx = hxi.setdefault((pa, pb), {})
        lcv = lv - (sp.Rational(nF[pa] + nF[pb], 2) if w == (pa, pb) else 0)
        ob = ((cphi[pa] - 2 * theta * wF[pa]) / 2 + (cphi[pb] - 2 * theta * wF[pb]) / 2) if w == (pa, pb) else 0
        d[w] = 4 * lcv + cv + ob; dl[w] = lcv; dx[w] = cv + ob
    target = json.load(open(os.path.join(Q.ROOT, 'step2', 'composite_qsk', 'composite', 'theory', 'src.nosync', 'ha_su11.json')))
    out = {'two_site_lc_natural': {str(k): {str(w): str(v) for w, v in d.items()} for k, d in hlc.items()},
           'two_site_Xi_Delta_translation_invariant': {str(k): {str(w): str(v) for w, v in d.items()} for k, d in hxi.items()},
           'spectra': {}}
    allok = True
    for L in range(2, Lmax + 1):
        for M in range(0, L + 1):
            ev = graded_cyclic_spectrum(h2, L, M)
            E_ = [sp.nsimplify(sp.simplify(-e / 8)) for e in ev]
            tgt = []
            for val, mult in target.get('L%d_M%d' % (L, M), []):
                tgt += [sp.sympify(val)] * mult
            ok = sorted([float(sp.N(x)) for x in E_]) == sorted([float(sp.N(x)) for x in tgt]) or \
                np.allclose(sorted([float(sp.N(x)) for x in E_]), sorted([float(sp.N(x)) for x in tgt]))
            allok &= bool(ok)
            out['spectra']['L%d_M%d' % (L, M)] = {'E_H_from_D1': [str(x) for x in E_], 'target': [str(x) for x in tgt], 'pass': bool(ok)}
    out['pass'] = bool(allok)
    log('  su(1|1) cyclic spectra L<=%d: %s' % (Lmax, allok))
    return out


def so6_density(seed=3, log=print):
    """planar (A_perp, A_perp) two-site composite density as an SO(6) tensor (off shell, real momenta):
    C(e_i, e_j) -> lambda (Y_i Y_j) + mu delta_ij (Y.Y) + nu (Y_j Y_i): fitted exactly, D0 and P_beta"""
    p1 = np.array([1, 0, 2, -1]); p2 = np.array([0, 1, -1, 1]); beta = (1, 2, -1, 3)
    def u6(i):
        v = np.zeros(6, complex); v[i] = 1; return v
    out = {}
    for kd in ('D0', ('Pv', beta)):
        vals = {}
        for (a, b, c, d) in [(0, 1, 0, 1), (0, 1, 1, 0), (0, 0, 1, 1), (0, 0, 0, 0)]:
            legs = [('AW', u6(c), p1), ('AW', u6(d), p2)]
            v = GV.combine(GV.kernel(('AW', u6(a), ONE), ('AW', u6(b), ONE), legs, kd), 1)
            vals[(a, b, c, d)] = F2.symv(v) / 2
        out[str(kd)] = {'lambda (same order)': str(vals[(0, 1, 0, 1)]), 'nu (exchanged)': str(vals[(0, 1, 1, 0)]),
                        'mu (trace)': str(vals[(0, 0, 1, 1)]), 'diagonal = lambda+mu+nu': str(vals[(0, 0, 0, 0)])}
    bP = sp.Rational(-1, 6) * F2.symc(1j * np.dot(beta, p1 + p2))
    ok = (sp.simplify(sp.sympify(out['D0']['lambda (same order)']) + sp.Rational(4, 3)) == 0 and
          sp.sympify(out['D0']['nu (exchanged)']) == 0 and sp.sympify(out['D0']['mu (trace)']) == 0 and
          sp.simplify(sp.sympify(out[str(('Pv', beta))]['lambda (same order)']) - bP) == 0 and
          sp.sympify(out[str(('Pv', beta))]['nu (exchanged)']) == 0 and sp.sympify(out[str(('Pv', beta))]['mu (trace)']) == 0)
    out['pass_universal_-4/3_and_-1/6'] = bool(ok)
    log('  SO(6) scalar density: %s' % ok)
    return out
