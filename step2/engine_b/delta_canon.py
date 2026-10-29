"""Canonical-form evaluation of graph sums with one S_M letter (and optionally one O letter) and one special edge.

S_M = S_w + S_E, S_E = sum_mu (-i d/dkappa_mu) T_mu[kappa].  One pass with a symbolic direction a (T_a = sum a_mu T_mu) and the
injection kappa = eta * a (eta^2 = 0):  G(a, eta a) = sum_mu a_mu G_mu(0) + eta sum_{mu nu} a_mu a_nu d_nu G_mu(0), so
   beta^mu  = -i [a_mu] G|_{eta=0},     alpha_E = -i sum_mu (1/2) d^2/da_mu^2 [eta] G,     alpha = G_w + alpha_E.
Canonical form: value = alpha delta(sum p) + beta^mu d/dp_2^mu delta(sum p), leg 2 = eliminated leg (momentum -p1 - k - kappa)."""
import sympy as sp
import delta_engine as D

I = sp.I
ETA = sp.Symbol("eta_inj")
q = D.q
AD = D.AD

def trunc_eta(expr):
    expr = sp.expand(expr)
    return sp.expand(expr.coeff(ETA, 0) + ETA * expr.coeff(ETA, 1))

def build_legs(leg1, s2, u2, p1, kvec, om, kap):
    legs = [dict(sp=leg1['sp'], u=leg1['u'], p=list(p1)),
            dict(sp=s2, u=u2, p=[-p1[m] - kvec[m] - kap[m] for m in range(4)]),
            dict(sp='xi', u=[sp.Integer(1)], p=list(kap))]
    if om:
        legs.append(dict(sp='om', u=[sp.Integer(1)], p=list(kvec)))
    return legs

def partitions(legs, sm_terms, om):
    s1, s2 = legs[0]['sp'], legs[1]['sp']
    base = [('SM', (2,), sm_terms)] + ([('O', (3,), ['O'])] if om else [])
    parts = [[('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC)] + base]
    if s1 == 'A' and s2 == 'A':
        parts.append([('quartic', (0, 1), D.V.QUARTIC)] + base)
    return parts

def ops_T():
    ops = {}
    for mu in range(4):
        ops[f'b{mu}'] = (lambda m: (lambda e: sp.expand(-I * sp.diff(e, AD[m]).xreplace({ETA: 0, **{a: 0 for a in AD}}))))(mu)
    def aE(e):
        c1 = sp.expand(e).coeff(ETA, 1)
        return sp.expand(-I * sum(sp.Rational(1, 2) * sp.diff(c1, AD[m], 2) for m in range(4)))
    ops['aE'] = aE
    return ops

LAMB = sp.Symbol('lam_ext')

def trunc_lam(expr, nmax):
    expr = sp.expand(expr)
    return sp.expand(sum(LAMB ** j * expr.coeff(LAMB, j) for j in range(nmax + 1)))

def special_scaled(Num, curs, s, ops, omega):
    """one special edge; all external momenta carry the factor LAMB.  The eps^n coefficient (n = -1, 0) has momentum degree
    omega + 2n, i.e. it is the LAMB^(omega+2n) Taylor coefficient: only that coefficient is computed."""
    E = len(curs)
    ks = [q[m] - curs[s][m] for m in range(4)]
    Nt = trunc_lam(-Num * sum(x * x for x in ks), omega)
    shift = [sum(D.BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
    Cm = D.L.gaussian_moments(Nt, shift)
    Cm = {k: trunc_lam(v, omega) for k, v in Cm.items()}
    rho = sp.expand(sum(D.BET[e] * sum(x * x for x in curs[e]) for e in range(E)) - sum(x * x for x in shift))
    order = [s] + [e for e in range(E) if e != s]
    res = {nm: {} for nm in ops}
    for n in (-1, 0):
        deg = omega + 2 * n
        base = 0
        if deg >= 0:
            for k, Ck in Cm.items():
                j = k + 2 - E + n
                if j < 0:
                    continue
                if 2 * j > deg:
                    continue
                base += sp.expand(Ck * rho ** j).coeff(LAMB, deg) * sp.Rational((-1) ** j, 2 ** k * sp.factorial(j))
        base = sp.expand(base)
        for nm, op in ops.items():
            expr = sp.expand(op(base)) if base != 0 else 0
            acc = {'1': 0, 'L': 0}
            if expr != 0:
                poly = sp.Poly(expr, *D.BET[:E])
                for mon, cf in poly.terms():
                    al = tuple(mon[e] for e in order)
                    J = _mcache(al, E - n + sum(al))
                    acc['1'] += cf * J[0]; acc['L'] += cf * J[1]
            res[nm][n] = {kk: sp.expand(v) for kk, v in acc.items()}
    return res

_MC = {}
def _mcache(al, mm):
    key = (al, mm)
    if key not in _MC:
        J = D.master_special(al[0], list(al[1:]), mm)[0]
        _MC[key] = (sp.Rational(J['1'].numerator, J['1'].denominator), sp.Rational(J['L'].numerator, J['L'].denominator))
    return _MC[key]

def evaluate_scaled(leg1, s2, u2, p1, omega_w, kvec=(0, 0, 0, 0), om=False, filt=None, collect=None, passes=('w', 'T')):
    """as evaluate (mode K), with all external momenta scaled by LAMB and only the needed Taylor coefficients kept.
    omega_w: superficial degree of the graphs with the weight letter (the translation letter adds one)."""
    kvec = [sp.Integer(x) for x in kvec]
    p1s = [LAMB * x for x in p1]; ks_ = [LAMB * x for x in kvec]
    colour_legs = [0, 1]
    res = {n: dict(alpha={}, beta=[{} for _ in range(4)]) for n in (-1, 0)}
    def add(d, key, v):
        d[key] = sp.expand(d.get(key, 0) + v)
    if 'w' in passes:
        legs = build_legs(leg1, s2, u2, p1s, ks_, om, [0] * 4)
        for part in partitions(legs, D.SW_TERMS, om):
            for r in D.enumerate_terms(legs, part, 'K', colour_legs, filt):
                v = special_scaled(trunc_lam(r['Num'], omega_w), r['curs'], r['s'], {'w': lambda e: e}, omega_w)
                Lk = f'L{len(r["curs"])}'
                for n in (-1, 0):
                    add(res[n]['alpha'], '1', r['fac'] * v['w'][n]['1']); add(res[n]['alpha'], Lk, r['fac'] * v['w'][n]['L'])
                if collect is not None:
                    collect.append(('w', r['word'], r['edges'], r['s'], {n: {k_: str(r['fac'] * x) for k_, x in v['w'][n].items()} for n in (-1, 0)}))
    if 'T' in passes:
        kap = [LAMB * ETA * AD[m] for m in range(4)]
        legs = build_legs(leg1, s2, u2, p1s, ks_, om, kap)
        ops = ops_T()
        for part in partitions(legs, D.TA_TERMS, om):
            for r in D.enumerate_terms(legs, part, 'K', colour_legs, filt):
                Num = trunc_lam(trunc_eta(r['Num']), omega_w + 1)
                curs = [[trunc_eta(x) for x in c] for c in r['curs']]
                v = special_scaled(Num, curs, r['s'], ops, omega_w + 1)
                v = {nm: {n: {kk: trunc_eta(x) for kk, x in v[nm][n].items()} for n in v[nm]} for nm in v}
                Lk = f'L{len(r["curs"])}'
                for n in (-1, 0):
                    add(res[n]['alpha'], '1', r['fac'] * v['aE'][n]['1']); add(res[n]['alpha'], Lk, r['fac'] * v['aE'][n]['L'])
                    for mu in range(4):
                        add(res[n]['beta'][mu], '1', r['fac'] * v[f'b{mu}'][n]['1'])
                        add(res[n]['beta'][mu], Lk, r['fac'] * v[f'b{mu}'][n]['L'])
                if collect is not None:
                    collect.append(('T', r['word'], r['edges'], r['s'], {nm: {n: {k_: str(r['fac'] * x) for k_, x in v[nm][n].items()} for n in (-1, 0)} for nm in v}))
    return res

def evaluate(leg1, s2, u2, p1, kvec=(0, 0, 0, 0), om=False, mode='K', filt=None, collect=None):
    """returns {n: {'alpha': {'1','L3','L4','L2'..}, 'beta': [...]}} for n = -1, 0 (mode 'K')."""
    kvec = [sp.Integer(x) for x in kvec]
    colour_legs = [0, 1]
    res = {n: dict(alpha={}, beta=[{} for _ in range(4)]) for n in (-1, 0)}
    def add(d, key, v):
        d[key] = sp.expand(d.get(key, 0) + v)
    # weight letter
    legs = build_legs(leg1, s2, u2, p1, kvec, om, [0] * 4)
    for part in partitions(legs, D.SW_TERMS, om):
        for r in D.enumerate_terms(legs, part, mode, colour_legs, filt):
            v = D.special_coeffs_multi(r['Num'], r['curs'], r['s'], {'w': lambda e: e})
            Lk = f'L{len(r["curs"])}'
            for n in (-1, 0):
                add(res[n]['alpha'], '1', r['fac'] * v['w'][n]['1']); add(res[n]['alpha'], Lk, r['fac'] * v['w'][n]['L'])
            if collect is not None:
                collect.append(('w', r['word'], r['edges'], r['s'], {n: {k_: str(r['fac'] * x) for k_, x in v['w'][n].items()} for n in (-1, 0)}))
    # translation letter, kappa = eta * a
    kap = [ETA * AD[m] for m in range(4)]
    legs = build_legs(leg1, s2, u2, p1, kvec, om, kap)
    ops = ops_T()
    for part in partitions(legs, D.TA_TERMS, om):
        for r in D.enumerate_terms(legs, part, mode, colour_legs, filt):
            Num = trunc_eta(r['Num'])
            curs = [[trunc_eta(x) for x in c] for c in r['curs']]
            v = D.special_coeffs_multi(Num, curs, r['s'], ops)
            Lk = f'L{len(r["curs"])}'
            for n in (-1, 0):
                add(res[n]['alpha'], '1', r['fac'] * v['aE'][n]['1']); add(res[n]['alpha'], Lk, r['fac'] * v['aE'][n]['L'])
                for mu in range(4):
                    add(res[n]['beta'][mu], '1', r['fac'] * v[f'b{mu}'][n]['1'])
                    add(res[n]['beta'][mu], Lk, r['fac'] * v[f'b{mu}'][n]['L'])
            if collect is not None:
                collect.append(('T', r['word'], r['edges'], r['s'], {nm: {n: {k_: str(r['fac'] * x) for k_, x in v[nm][n].items()} for n in (-1, 0)} for nm in v}))
    return res

def evaluate_scaled4(leg1, s2, u2, p1, omega_w, kvec=(0, 0, 0, 0), om=False, filt=None, collect=None, passes=('w', 'T')):
    """as evaluate_scaled, but the translation letter is evaluated in four passes with numeric directions a = e_mu and
    injection kappa = eta e_mu (only eta symbolic):  beta^mu = -i G_mu(0),  alpha_E = -i sum_mu [eta] G_mu(eta e_mu)."""
    kvec = [sp.Integer(x) for x in kvec]
    p1s = [LAMB * x for x in p1]; ks_ = [LAMB * x for x in kvec]
    colour_legs = [0, 1]
    res = {n: dict(alpha={}, beta=[{} for _ in range(4)]) for n in (-1, 0)}
    def add(d, key, v):
        d[key] = sp.expand(d.get(key, 0) + v)
    if 'w' in passes:
        legs = build_legs(leg1, s2, u2, p1s, ks_, om, [0] * 4)
        for part in partitions(legs, D.SW_TERMS, om):
            for r in D.enumerate_terms(legs, part, 'K', colour_legs, filt):
                v = special_scaled(trunc_lam(r['Num'], omega_w), r['curs'], r['s'], {'w': lambda e: e}, omega_w)
                Lk = f'L{len(r["curs"])}'
                for n in (-1, 0):
                    add(res[n]['alpha'], '1', r['fac'] * v['w'][n]['1']); add(res[n]['alpha'], Lk, r['fac'] * v['w'][n]['L'])
                if collect is not None:
                    collect.append(('w', r['word'], r['edges'], r['s'], {n: {k_: str(r['fac'] * x) for k_, x in v['w'][n].items()} for n in (-1, 0)}))
    if 'T' in passes:
        ops = {'b': lambda e: sp.expand(-I * e.xreplace({ETA: 0})),
               'aE': lambda e: sp.expand(-I * sp.expand(e).coeff(ETA, 1))}
        for mu in range(4):
            kap = [LAMB * ETA * (1 if m == mu else 0) for m in range(4)]
            legs = build_legs(leg1, s2, u2, p1s, ks_, om, kap)
            for part in partitions(legs, D.T_TERMS(mu), om):
                for r in D.enumerate_terms(legs, part, 'K', colour_legs, filt):
                    Num = trunc_lam(trunc_eta(r['Num']), omega_w + 1)
                    curs = [[trunc_eta(x) for x in c] for c in r['curs']]
                    v = special_scaled(Num, curs, r['s'], ops, omega_w + 1)
                    Lk = f'L{len(r["curs"])}'
                    for n in (-1, 0):
                        add(res[n]['alpha'], '1', r['fac'] * v['aE'][n]['1']); add(res[n]['alpha'], Lk, r['fac'] * v['aE'][n]['L'])
                        add(res[n]['beta'][mu], '1', r['fac'] * v['b'][n]['1'])
                        add(res[n]['beta'][mu], Lk, r['fac'] * v['b'][n]['L'])
                    if collect is not None:
                        collect.append((f'T{mu}', r['word'], r['edges'], r['s'], {nm: {n: {k_: str(r['fac'] * x) for k_, x in v[nm][n].items()} for n in (-1, 0)} for nm in v}))
    return res
