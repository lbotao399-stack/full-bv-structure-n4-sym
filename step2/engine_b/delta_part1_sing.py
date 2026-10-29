"""Part D4 (1): the log eps coefficients of the 1PI one-loop graphs with S_M on the cycle (class (i)), i.e. the arity-one
counterterm I^CT[eps](Delta) on the quadratic sector, in the canonical form

   value on (field u_f e^{ipx}, antifield u_a e^{i p_a x}, xi) = alpha(p) delta(p+p_a) + beta^mu(p) d/dp_a^mu delta(p+p_a)
   <->  int dvol [ (antifield, alpha(-i d) field) + x^mu (antifield, i beta_mu(-i d) field) ]   (x on the antifield leg),

alpha = G_w - i sum_mu d/dkappa_mu G_mu |_{kappa=0},  beta^mu = -i G_mu |_{kappa=0}, G_w / G_mu = graphs with the S_w / T_mu[kappa]
letter, p_a = -p - kappa.  Units: unit = hbar tr_g(ad ad)/(16 pi^2 (X,Y)), colour (Z_f, Z_a)."""
import json, random, time, sys
import sympy as sp
import delta_engine as D

I = sp.I
P = list(sp.symbols('P0:4'))
KAP = list(sp.symbols('K0:4'))
random.seed(11)

SECT = {
    'Aperp': ('A', 'As', list(range(4, 10))),
    'AV': ('A', 'As', list(range(4))),
    'lam': ('lam', 'lams', list(range(16))),
    'c': ('c', 'cs', [0]),
}

def legs_for(sec, uf, ua, kap):
    sf, sa, comps = SECT[sec]
    pa = [-P[m] - kap[m] for m in range(4)]
    return [dict(sp=sf, u=uf, p=list(P)), dict(sp=sa, u=ua, p=pa), dict(sp='xi', u=[sp.Integer(1)], p=list(kap))]

def rand_pol(sp_, comps):
    u = [sp.Integer(0)] * D.V.DIM[sp_]
    for c in comps:
        u[c] = sp.Integer(random.choice([-3, -2, -1, 1, 2, 3]))
    return u

def class_i_log(sec, uf, ua, verbose=False):
    t0 = time.time()
    res = {}
    # weight part, kappa = 0
    legs = legs_for(sec, uf, ua, [0] * 4)
    part = [('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC), ('SM', (2,), D.SW_TERMS)]
    Gw = 0; ncl = {}
    for r in D.enumerate_terms(legs, part, 'P', [0, 1]):
        v = r['fac'] * D.log_coeff(r['Num'], r['curs'], 0, P)
        Gw += v
        key = str([w[0] for w in r['word']]) + str(r['edges'])
        ncl[key] = ncl.get(key, 0) + v
    res['Gw'] = sp.expand(Gw)
    # translation part
    legs = legs_for(sec, uf, ua, KAP)
    alphaE = 0; beta = [0] * 4
    for mu in range(4):
        part = [('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC), ('SM', (2,), D.T_TERMS(mu))]
        Gmu = 0
        for r in D.enumerate_terms(legs, part, 'P', [0, 1]):
            Gmu += r['fac'] * D.log_coeff(r['Num'], r['curs'], 1, P + KAP)
        Gmu = sp.expand(Gmu)
        alphaE += -I * sp.diff(Gmu, KAP[mu])
        beta[mu] = sp.expand(-I * Gmu.xreplace({k: 0 for k in KAP}))
    res['alphaE'] = sp.expand(alphaE.xreplace({k: 0 for k in KAP}))
    res['alpha'] = sp.expand(res['Gw'] + res['alphaE'])
    res['beta'] = beta
    res['time'] = round(time.time() - t0, 1)
    return res

def tree_SM(sec, uf, ua):
    """the tree S_M itself in the same canonical form (check of the TERMS): alpha = Sw value, beta^mu = -i T_mu value."""
    legs0 = legs_for(sec, uf, ua, [0] * 4)
    al = sum(D.V.vertex(t, legs0).get('k', 0) for t in D.SW_TERMS)
    be = []
    for mu in range(4):
        be.append(sp.expand(-I * sum(D.V.vertex(t, legs0).get('k', 0) for t in D.T_TERMS(mu))))
    return sp.expand(al), be

if __name__ == '__main__':
    secs = sys.argv[1].split(',') if len(sys.argv) > 1 else list(SECT)
    nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    out = {}
    for sec in secs:
        sf, sa, comps = SECT[sec]
        samples = []
        for k in range(nsamp):
            uf = rand_pol(sf, comps); ua = rand_pol(sa, comps if sa != 'lams' else list(range(16)))
            r = class_i_log(sec, uf, ua)
            ta, tb = tree_SM(sec, uf, ua)
            samples.append(dict(uf=[str(x) for x in uf], ua=[str(x) for x in ua], alpha=str(r['alpha']), Gw=str(r['Gw']),
                                alphaE=str(r['alphaE']), beta=[str(x) for x in r['beta']], tree_alpha=str(ta),
                                tree_beta=[str(x) for x in tb], time=r['time']))
            print(sec, k, 'alpha', r['alpha'], 'beta', r['beta'], 'tree', ta, tb, r['time'], 's', flush=True)
        out[sec] = samples
        json.dump(out, open(f'delta_part1_sing_{"_".join(secs)}.json', 'w'), indent=1)
