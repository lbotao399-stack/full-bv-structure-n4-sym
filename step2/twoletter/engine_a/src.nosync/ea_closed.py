"""Engine A: closed forms in the derivative numbers, verified against explicit evaluations with derivative letters.
(1) k = 1 with constant chain numerator N_V (quartic bubbles): eq (eq:tl-bubbleclosed),
    16 pi^2 [log] = -N_V sum_{a,b,j} (i(s,R))^a (i(u,R))^b |u-s|^{2j} (|R|^2/4)^j (a+j)!(b+j)!/(a!b!(j!)^2(a+b+2j+1)!),
    for the point splitting X_a(x+s) X_b(x+u), R the total output momentum (letter momenta p_a = -l, p_b = l+R).
(2) general chains (Prop prop:tl-genfun), in engine variables (edges k+q_i, p_a = k+R, p_b = -k, sigma = u-s):
    sum s^al u^be/(al! be!) [log]K = -(1/16pi^2) e^{i(s,R)} int_Delta N_V(i d_sigma)[ e^{-i(sigma,c)} (-1)^E
         sum_m (|sigma|^2/4)^{m+E-2} Phi^m/(m!(m+E-2)!) ] d alpha,   c = -sum alpha_i q_i,
    an independent symbolic implementation (sympy) of the Bessel resummation, compared order by order in (s,u)."""
import math, itertools
from fractions import Fraction
import numpy as np
import sympy as sp
import ea_core as E, ea_kern as K


def letter_genpoly(vec, n):
    """(i p.vec)^n as {exp: Gaussian int} (times n! removed by the caller)"""
    out = {}
    for c in itertools.product(range(4), repeat=n):
        e = [0] * 4
        coef = 1j ** n
        for mu in c:
            e[mu] += 1; coef *= int(vec[mu])
        if coef != 0:
            out[tuple(e)] = out.get(tuple(e), 0) + coef
    return {k: complex(round(v.real), round(v.imag)) for k, v in out.items() if v != 0}


def engine_degree(blocks, ta, tb, za, zb, s, u, n):
    """sum_{A+B=n} K^log(letters with (i p_a.s)^A/A!, (i p_b.u)^B/B!) contracted with za, zb"""
    tot = E.Q0
    for A in range(n + 1):
        B = n - A
        fa = letter_genpoly(s, A); fb = letter_genpoly(u, B)
        if not fa or not fb:
            continue
        r = K.kernel(blocks, ta, tb, kinds=('log',), fa=fa, fb=fb)
        if not r:
            continue
        v = K.contract(r['log'], za, zb)
        f = Fraction(1, math.factorial(A) * math.factorial(B))
        tot = E.q_add(tot, (v[0] * f, v[1] * f))
    return tot


def bubble_closed_degree(Nv, R, s, u, n):
    """eq:tl-bubbleclosed, degree-n part at numeric s, u"""
    sR = int(np.dot(s, R)); uR = int(np.dot(u, R)); us2 = int(np.dot(u - s, u - s)); R2 = int(np.dot(R, R))
    tot = 0
    for j in range(n // 2 + 1):
        for a in range(n - 2 * j + 1):
            b = n - 2 * j - a
            c = Fraction(math.factorial(a + j) * math.factorial(b + j),
                        math.factorial(a) * math.factorial(b) * math.factorial(j) ** 2 * math.factorial(a + b + 2 * j + 1))
            tot += (1j * sR) ** a * (1j * uR) ** b * us2 ** j * Fraction(R2, 4) ** j * c
    return -Nv * tot


def check_bubbles(nmax=4, seed=7):
    rng = np.random.default_rng(seed)
    out = []
    for (ta, tb, outt) in (('AW', 'AW', 'AW'), ('AV', 'AV', 'AV'), ('AV', 'AW', None)):
        for n in range(nmax + 1):
            za = rng.integers(-2, 3, 6 if ta == 'AW' else 4).astype(complex)
            zb = rng.integers(-2, 3, 6 if tb == 'AW' else 4).astype(complex)
            def pol(t):
                v = np.zeros(10, complex)
                if t == 'AW': v[4:] = rng.integers(-2, 3, 6)
                else: v[:4] = rng.integers(-2, 3, 4)
                return v
            if outt is None:
                p1, p2 = pol('AV'), pol('AW')
            else:
                p1, p2 = pol(outt), pol(outt)
            q1, q2 = rng.integers(-2, 3, 4), rng.integers(-2, 3, 4)
            s, u = rng.integers(-2, 3, 4), rng.integers(-2, 3, 4)
            blocks = [('qua', (p1, q1), (p2, q2))]
            base = K.kernel(blocks, ta, tb, kinds=('log',))
            Nv = -complex(*map(float, K.contract(base['log'], za, zb)))      # 16pi^2 [log] at s=u=0 is -N_V
            got = engine_degree(blocks, ta, tb, za, zb, s, u, n)
            pred = bubble_closed_degree(Nv, q1 + q2, s, u, n)
            g = complex(float(got[0]), float(got[1]))
            ok = abs(g - complex(pred)) < 1e-9 * max(1, abs(g))
            out.append({'pair': ta + ',' + tb, 'n': n, 'engine': E.qstr(got), 'closed_form': str(pred), 'pass': bool(ok)})
    return out


def genfun_degree(Npoly_k, qs, R, s, u, n):
    """independent sympy evaluation of the generating function (degree n part) at numeric s, u"""
    sg = sp.symbols('sg0:4')
    Eg = len(qs)
    al = sp.symbols('a0:%d' % Eg)
    c = [-sum(al[i] * int(qs[i][mu]) for i in range(Eg)) for mu in range(4)]
    Phi = sp.expand(sum(al[i] * int(qs[i] @ qs[i]) for i in range(Eg)) - sum(x * x for x in c))
    degN = max(sum(e) for e in Npoly_k)
    order = n + degN
    s2 = sum(x * x for x in sg)
    sc = sum(sg[mu] * c[mu] for mu in range(4))
    expo = sum((-sp.I * sc) ** k / sp.factorial(k) for k in range(order + 1))
    bes = 0
    m = 0
    while 2 * (m + Eg - 2) <= order:
        if m + Eg - 2 >= 0:
            bes += (s2 / 4) ** (m + Eg - 2) * Phi ** m / (sp.factorial(m) * sp.factorial(m + Eg - 2))
        m += 1
    G = sp.expand((-1) ** Eg * expo * bes)
    # truncate to total sigma degree <= order
    Pg = sp.Poly(G, *sg)
    G = sum(cf * sp.prod([sg[i] ** e[i] for i in range(4)]) for e, cf in Pg.terms() if sum(e) <= order)
    # apply N_V(i d_sigma)
    H = 0
    for e, cf in Npoly_k.items():
        t = G
        for mu in range(4):
            for _ in range(e[mu]):
                t = sp.I * sp.diff(t, sg[mu])
        H += (sp.Rational(int(round(cf.real))) + sp.I * int(round(cf.imag))) * t
    H = sp.expand(H)
    # simplex integral
    P = sp.Poly(H, *al)
    Hint = 0
    for e, cf in P.terms():
        Hint += cf * sp.Rational(math.prod(math.factorial(x) for x in e), math.factorial(sum(e) + Eg - 1))
    # multiply e^{i s.R}, substitute sigma = u - s, extract degree n in (s,u) via scaling lam
    lam = sp.Symbol('lam')
    sub = {sg[mu]: lam * int(u[mu] - s[mu]) for mu in range(4)}
    Hs = sp.expand(sp.sympify(Hint).subs(sub))
    sR = int(np.dot(s, R))
    ex = sum((sp.I * lam * sR) ** k / sp.factorial(k) for k in range(n + 1))
    tot = sp.expand(-Hs * ex)
    return complex(sp.N(tot.coeff(lam, n), 30))


def check_genfun(nmax=3, seed=9):
    rng = np.random.default_rng(seed)
    out = []
    cases = []
    eW = np.zeros(10, complex); eW[4] = 1; eW[6] = -1; eW[8] = 2
    eV = np.zeros(10, complex); eV[0] = 1; eV[1] = 2; eV[3] = -1
    cases.append(('k=1 cubic AAA', 'AV', 'AV', [('cub', 'AV', eV, np.array([1, 0, -1, 1]))], np.array([1, 2, 0, -1], complex), np.array([0, 1, 1, 2], complex)))
    cases.append(('k=2 gluon exchange', 'AW', 'AW', [('cub', 'AW', eW, np.array([1, 0, 1, 0])), ('cub', 'AW', eW, np.array([0, 1, -1, 1]))],
                  np.array([1, 0, -1, 2, 0, 1], complex), np.array([0, 2, 1, -1, 1, 0], complex)))
    cases.append(('k=1 lamP lamM -> AV', 'lamP', 'lamM', [('cub', 'AV', eV, np.array([1, -1, 0, 2]))],
                  E.BASIS['lamP'][0].astype(complex) + E.BASIS['lamP'][3], E.BASIS['lamM'][1].astype(complex)))
    for (name, ta, tb, blocks, za, zb) in cases:
        Kp, qs, scale = K.chain_poly(blocks, ta, tb)
        sgnp = -1 if (E.PARITY_T[ta] + E.PARITY_T[tb]) % 2 else 1
        Nsc = {}
        for ex, M in Kp.items():
            v = complex(za @ M @ zb)
            if v:
                Nsc[ex] = v
        R = sum((b[3] for b in blocks), np.zeros(4, dtype=int))
        for n in range(nmax + 1):
            s, u = rng.integers(-2, 3, 4), rng.integers(-2, 3, 4)
            got = engine_degree(blocks, ta, tb, za, zb, s, u, n)
            pred = genfun_degree(Nsc, qs, R, s, u, n) / scale * sgnp
            g = complex(float(got[0]), float(got[1]))
            ok = abs(g - pred) < 1e-9 * max(1, abs(g))
            out.append({'chain': name, 'n': n, 'engine': E.qstr(got), 'genfun': '%.12g%+.12gi' % (pred.real, pred.imag), 'pass': bool(ok)})
    return out


if __name__ == '__main__':
    for r in check_bubbles(3):
        print(r)
    for r in check_genfun(2):
        print(r)
