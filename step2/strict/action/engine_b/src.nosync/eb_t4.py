"""Task T4: Ob^(1)(xi) on two legs (former TQ1), both chiralities, basis g1..g4 of Prop at-xi2; poles (G-pole).
Normalisation: odd generator, CE coordinate even -> theory component = raw left coefficient (eb_t1 notes)."""
import json, os, random, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_fit as F
import ce_core as C

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')
I = sp.I
Mp = C.Mp
MR = [Mp * C.Rm[M] for M in range(10)]
LAMP = list(range(8)); LAMM = list(range(8, 16))            # lam components: S_+ (x) W, S_- (x) W^*
LAMS_M = list(range(8)); LAMS_P = list(range(8, 16))        # lam* components: S_- (x) W, S_+ (x) W^*


def rv(n, lo=-3, hi=3):
    return [Fr(random.randint(lo, hi)) for _ in range(n)]


def symvec(v):
    return [sp.nsimplify(complex(x).real) + I * sp.nsimplify(complex(x).imag) for x in v]


def gamma_xi(xi_s, lam_u, M):
    return (sp.Matrix(xi_s).T * MR[M] * sp.Matrix(lam_u))[0]


def basis_xi(xi_s, tag):
    """g1..g4 and T0..T4 (Prop at-xi2) as plane-wave functionals at the numeric xi (sympy list)."""
    def g1(s):      # slots (lam, A): int (Gamma_V, d*d A_V); d*dA -> p^2 a - p (p.a)
        l, a = s
        p = a['p']; p2 = sum(x * x for x in p); pa = sum(p[m] * a['u'][m] for m in range(4))
        return sum(gamma_xi(xi_s, l['u'], n) * (p2 * a['u'][n] - p[n] * pa) for n in range(4))

    def g2(s):      # int (Gamma_V, d d* A_V); dd*A -> p (p.a)
        l, a = s
        p = a['p']; pa = sum(p[m] * a['u'][m] for m in range(4))
        return sum(gamma_xi(xi_s, l['u'], n) * p[n] * pa for n in range(4))

    def g3(s):      # int (Gamma_perp, d*d A_perp) -> p^2 w
        l, a = s
        p = a['p']; p2 = sum(x * x for x in p)
        return sum(gamma_xi(xi_s, l['u'], M) * p2 * a['u'][M] for M in range(4, 10))

    def g4(s):      # int (d*d c, (xi, lam*)) -> p_c^2 xi^T Mp lam*
        c, ls = s
        p = c['p']; p2 = sum(x * x for x in p)
        return p2 * c['u'][0] * (sp.Matrix(xi_s).T * Mp * sp.Matrix(ls['u']))[0]

    def gAG(rng):   # int (A, Gamma(xi,lam)) restricted (pole structure)
        def f(s):
            l, a = s
            return sum(gamma_xi(xi_s, l['u'], M) * a['u'][M] for M in rng)
        return f
    out = {
        'g1': [F.reg(f'g1_{tag}', 1, ('lam', 'A'), g1)],
        'g2': [F.reg(f'g2_{tag}', 1, ('lam', 'A'), g2)],
        'g3': [F.reg(f'g3_{tag}', 1, ('lam', 'A'), g3)],
        'g4': [F.reg(f'g4_{tag}', 1, ('c', 'lams'), g4)],
        'int(Gamma_V,A_V)': [F.reg(f'gAGV_{tag}', 1, ('lam', 'A'), gAG(range(4)))],
        'int(Gamma_perp,A_perp)': [F.reg(f'gAGp_{tag}', 1, ('lam', 'A'), gAG(range(4, 10)))],
    }
    return out


SECTORS = {
    # sector: (leg makers, basis labels, pole labels)
    '+': {'c_lams+': ('c', 'lams', LAMS_P, ['g4'], []),
          'AV_lam-': ('A', 'lam', LAMM, ['g1', 'g2'], ['int(Gamma_V,A_V)']),
          'Ap_lam+': ('A', 'lam', LAMP, ['g3'], ['int(Gamma_perp,A_perp)'])},
    '-': {'c_lams-': ('c', 'lams', LAMS_M, ['g4'], []),
          'AV_lam+': ('A', 'lam', LAMP, ['g1', 'g2'], ['int(Gamma_V,A_V)']),
          'Ap_lam-': ('A', 'lam', LAMM, ['g3'], ['int(Gamma_perp,A_perp)'])},
}


def make_legs(spec):
    sp0, sp1, comps, _, _ = spec
    p = rv(4)
    if sp0 == 'c':
        u0 = np.array([1], complex)
    else:
        u0 = np.zeros(10, complex)
        rng = range(4) if spec[3][0] in ('g1', 'g2') else range(4, 10)
        for M in rng:
            u0[M] = random.randint(-2, 2)
    u1 = np.zeros(16, complex)
    for c in comps:
        u1[c] = random.randint(-2, 2)
    # leg order: odd-parity conventions follow the leg list; put the spinor/antifield second
    return [dict(sp=sp0, u=u0, p=p), dict(sp=sp1, u=u1, p=[-x for x in p])]


def run(chir='+', seed=3, extra=2, log=print):
    random.seed(seed)
    tq, t3 = G.tmpl_xi()
    res = {}
    for sec, spec in SECTORS[chir].items():
        labels = spec[3]
        nconf = len(labels) + extra
        vals0, vals1, cands, candp = [], [], [], []
        st = {}
        t0 = time.time()
        for ic in range(nconf):
            xi = np.zeros(16, complex)
            for c in (range(8) if chir == '+' else range(8, 16)):
                xi[c] = random.randint(-2, 2)
            gen = Y.Gen('xi', 0, tq, gvec=xi, t3=t3)
            legs = make_legs(spec)
            r = Y.evaluate(legs, [gen], mlist=(0, 1), stats=st)
            vals0.append(r.get('k', {}).get(0, E.CLog()))
            vals1.append(r.get('k', {}).get(1, E.CLog()))
            Bs = basis_xi(symvec(xi), F.newtag())
            cands.append([F.value(Bs[l], legs) for l in labels])
            if spec[4]:
                candp.append([F.value(Bs[l], legs) for l in spec[4]])
        fit0, ok0 = F.fit(vals0, cands, labels)
        if spec[4]:
            fit1, ok1 = F.fit(vals1, candp, spec[4])
        else:
            fit1, ok1 = {}, all(v.is_zero() for v in vals1)
        res[sec] = dict(eps0=jsonable(fit0), eps0_consistent=ok0, pole=jsonable(fit1), pole_consistent=ok1,
                        nconf=nconf, chains=st.get('chains', 0), seconds=round(time.time() - t0, 1))
        log(f'T4 xi{chir} {sec}: {res[sec]}')
    return res


def jsonable(fitres):
    out = {}
    for k, d in (fitres or {}).items():
        out[k] = None if d is None else {lab: str(x) for lab, x in d.items()}
    return out


if __name__ == '__main__':
    out = {chir: run(chir) for chir in ('+', '-')}
    with open(os.path.join(OUT, 'T4.json'), 'w') as f:
        json.dump(out, f, indent=1)
