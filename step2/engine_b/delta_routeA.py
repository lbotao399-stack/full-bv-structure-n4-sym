"""Numerical class (i) (S_M on the cycle, all edges P(eps,infty)) at finite eps, canonical form, direct Schwinger quadrature:
  - sector (A|_Wperp, A*|_Wperp): finite limits alpha_i -> a0 (u_f,u_a), beta_i -> b0 p (u_f,u_a) (no log, X_perp = 0);
    route A of the arity-one obstruction: Q(class i) + D_M Sigma_ren = -(p^2 (a0 + 2 b0)) (u,u') - 8 p^2 (u,u')
    must equal the UV-boundary value -(22/3) p^2 (u,u');
  - sector (c, c*): eps d/deps of beta_i -> the exact log coefficient -p (u_c u_c*);
  - scaling (Ward) check: alpha_i(p, eps) = alpha_i(s p, eps/s^2), beta_i(s p, eps/s^2) = s beta_i(p, eps)."""
import json, sys, time, math
import numpy as np
import sympy as sp
import delta_engine as D
import delta_canon as C
import delta_numeric_i as NI

I = sp.I
def recs_for(sec_legs, terms, kap):
    leg1, s2, u2, p1 = sec_legs
    legs = [dict(sp=leg1['sp'], u=leg1['u'], p=list(p1)), dict(sp=s2, u=u2, p=[-p1[m] - kap[m] for m in range(4)]),
            dict(sp='xi', u=[sp.Integer(1)], p=list(kap))]
    part = [('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC), ('SM', (2,), terms)]
    return D.enumerate_terms(legs, part, 'P', [0, 1])

def canonical_numeric(leg1, s2, u2, p1, eps, n):
    """alpha, beta (numbers, units (16 pi^2)^-1) at finite eps."""
    alpha = 0j; beta = [0j] * 4
    recs = recs_for((leg1, s2, u2, p1), D.SW_TERMS, [0] * 4)
    v = NI.evaluate_graphs(recs, {'w': lambda e: e}, eps, n=n)
    alpha += v['w']
    for mu in range(4):
        kap = [C.ETA * (1 if m == mu else 0) for m in range(4)]
        recs = recs_for((leg1, s2, u2, p1), D.T_TERMS(mu), kap)
        ops = {'b': lambda e: sp.expand(-I * sp.expand(e).xreplace({C.ETA: 0})),
               'aE': lambda e: sp.expand(-I * sp.expand(e).coeff(C.ETA, 1))}
        v = NI.evaluate_graphs(recs, ops, eps, n=n)
        alpha += v['aE']; beta[mu] = v['b']
    return alpha, beta

if __name__ == '__main__':
    out = {}
    sec = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 48
    epss = [float(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else [1e-3, 1e-4, 1e-5]
    if sec == 'Ap':
        u1 = [0] * 4 + [-1, -3, -3, 2, -1, 1]; u2 = [0] * 4 + [1, 2, -1, -2, 3, -2]
        leg1 = dict(sp='A', u=[sp.Integer(x) for x in u1]); s2 = 'As'; u2s = [sp.Integer(x) for x in u2]
        uu = sum(a * b for a, b in zip(u1, u2))
    else:
        leg1 = dict(sp='c', u=[sp.Integer(2)]); s2 = 'cs'; u2s = [sp.Integer(-3)]; uu = -6
    p1 = [sp.Integer(x) for x in (-1, 1, -3, -2)]
    pn = np.array([float(x) for x in p1]); p2 = float(pn @ pn)
    for eps in epss:
        t0 = time.time()
        a, b = canonical_numeric(leg1, s2, u2s, p1, eps, n)
        bp = [complex(x) / (uu * pn[m]) for m, x in enumerate(b)]
        out[str(eps)] = dict(alpha_over_uu=[complex(a).real / uu, complex(a).imag / uu], beta_over_p_uu=[[x.real, x.imag] for x in bp],
                             time=round(time.time() - t0, 1))
        print(sec, eps, out[str(eps)], flush=True)
        json.dump(out, open(f'delta_routeA_{sec}_n{n}.json', 'w'), indent=1)
    # scaling (Ward) check: alpha(2p, eps/4) = alpha(p, eps), beta(2p, eps/4) = 2 beta(p, eps)
    eps = epss[0]
    a1, b1 = canonical_numeric(leg1, s2, u2s, p1, eps, n)
    a2, b2 = canonical_numeric(leg1, s2, u2s, [2 * x for x in p1], eps / 4, n)
    out['scaling_check'] = dict(eps=eps, alpha=[str(complex(a1)), str(complex(a2))],
                                beta=[[str(complex(x)) for x in b1], [str(complex(x) / 2) for x in b2]])
    print('scaling', out['scaling_check'], flush=True)
    json.dump(out, open(f'delta_routeA_{sec}_n{n}.json', 'w'), indent=1)
