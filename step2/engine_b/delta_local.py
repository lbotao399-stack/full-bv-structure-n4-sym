"""Local functionals with explicit position dependence (at most linear in x), in the canonical form used throughout Part D4.

A two-leg functional on (leg 1: u1 e^{i p1 x}, leg 2: u2 e^{i p2 x}) has the value
    alpha(p1) delta(p1+p2) + beta^mu(p1) d/dp2^mu delta(p1+p2)   (delta = (2 pi)^4 delta, leg 2 = eliminated leg)
<-> int dvol [ (leg2, alpha(-i d) leg1) + x^mu (leg2, i beta_mu(-i d) leg1) ]   (x on leg 2).
Injection representation (int x^mu f = -i d/dkappa_mu int e^{i kappa.x} f |_0): the x-part is the translation-invariant
functional int e^{i kappa.x} (leg2, i beta_mu(-i d) leg1) carried by an extra even leg of species 'k<mu>' (momentum kappa).
Brackets and Q are the step-2 tree contractions of fa.py ({f,g} = -B_{K_0}(f,g), Q = {S_2,-}), which act on these
injection legs as on any other leg."""
import sympy as sp
import delta_engine as D
import fa

V = D.V
I = sp.I
for mu in range(4):
    V.DIM[f'k{mu}'] = 1; V.PAR[f'k{mu}'] = 0
    fa.PAR[f'k{mu}'] = 0; fa.DIM[f'k{mu}'] = 1
fa.PAR['xi'] = 1; fa.DIM['xi'] = 1; fa.PAR['om'] = 0; fa.DIM['om'] = 1

_cnt = [0]
def canon(s1, s2, alpha, beta, name=None, parity=None):
    """alpha(u1, u2, p1) -> number; beta(u1, u2, p1, mu) -> number (beta^mu). Returns an fa.Functional."""
    _cnt[0] += 1
    nm = name or f'canon{_cnt[0]}'
    terms = []
    if alpha is not None:
        V.TERMS[nm + '_a'] = (sp.Integer(1), (s1, s2), (lambda f: (lambda s: f(s[0]['u'], s[1]['u'], s[0]['p'])))(alpha), ('k', (0, 1)))
        terms.append(nm + '_a')
    if beta is not None:
        for mu in range(4):
            V.TERMS[nm + f'_x{mu}'] = (sp.Integer(1), (s1, s2, f'k{mu}'),
                                       (lambda f, m: (lambda s: I * f(s[0]['u'], s[1]['u'], s[0]['p'], m) * s[2]['u'][0]))(beta, mu),
                                       ('k', (0, 1)))
            terms.append(nm + f'_x{mu}')
    par = parity if parity is not None else (V.PAR[s1] + V.PAR[s2]) % 2
    return fa.local([(1, t) for t in terms], nm, par)

def to_canon(F, s1, s2, u1, u2, P, KAP):
    """evaluate a two-leg functional (possibly with injection legs) in canonical form on legs (s1,u1,P), (s2,u2,-P-kappa).
    Returns (alpha, [beta^mu]) as sympy expressions in P."""
    leg1 = dict(sp=s1, u=u1, p=list(P))
    leg2_0 = dict(sp=s2, u=u2, p=[-x for x in P])
    a0 = F([leg1, leg2_0]).get('k', 0)
    alpha = sp.expand(a0); beta = []
    for mu in range(4):
        leg2 = dict(sp=s2, u=u2, p=[-P[m] - KAP[m] for m in range(4)])
        legk = dict(sp=f'k{mu}', u=[sp.Integer(1)], p=list(KAP))
        v = F([leg1, leg2, legk]).get('k', 0)
        v = sp.expand(v)
        alpha += -I * sp.diff(v, KAP[mu]).xreplace({k: 0 for k in KAP})
        beta.append(sp.expand(-I * v.xreplace({k: 0 for k in KAP})))
    return sp.expand(alpha), beta

# ------------------------------------------------------------------ S_M and the class-(i) log coefficient X
def dot(a, b, comps):
    return sum(a[i] * b[i] for i in comps)
AVc = list(range(4)); APc = list(range(4, 10)); ALL10 = list(range(10))

def lam_pair(u1, u2):
    return (sp.Matrix(u1).T * V.Mp * sp.Matrix(u2))[0]

def S_M_functional():
    """S_M = int(v.A, A*) - int(v.lambda, lambda*) - int(v.c, c*) in canonical form (x on the antifield leg)."""
    fA = canon('A', 'As', lambda u, w, p: 2 * dot(u, w, ALL10), lambda u, w, p, m: 2 * p[m] * dot(u, w, ALL10), 'SM_A')
    fl = canon('lam', 'lams', lambda u, w, p: -3 * lam_pair(u, w), lambda u, w, p, m: -2 * p[m] * lam_pair(u, w), 'SM_l')
    fc = canon('c', 'cs', None, lambda u, w, p, m: -2 * p[m] * u[0] * w[0], 'SM_c')
    return fA + fl + fc

def X_functional():
    """log eps coefficient of class (i) (Part 1, exact): sector A|_V: alpha = (u_f,u_a)_V, beta^mu = u_f^mu (p,u_a)_V;
    sector c: beta^mu = -p^mu u_c u_c*; sector lambda: beta^mu = -(rho(e_mu) rho(p) u_f, u_a); sector A|_Wperp: 0."""
    from realization import R, Rb
    fV = canon('A', 'As', lambda u, w, p: dot(u, w, AVc), lambda u, w, p, m: (u[m] if m < 4 else 0) * dot(p, w, AVc), 'X_V')
    fc = canon('c', 'cs', None, lambda u, w, p, m: -p[m] * u[0] * w[0], 'X_c')
    def bl(u, w, p, m):
        Rp = sum((p[n] * R[n] for n in range(4)), sp.zeros(16, 16))
        return -lam_pair(list(Rb[m] * Rp * sp.Matrix(u)), w)
    fl = canon('lam', 'lams', None, bl, 'X_l')
    return fV + fc + fl
