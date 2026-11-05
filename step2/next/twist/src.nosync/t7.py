import env, twdriver as D, sector as S, sympy as sp, json
T = S.T
def ci_of(key): return [i for i, c in enumerate(S.COMPS) if c[1] == key][0]
def gam_key(i):
    t = T.taus[i]; wt = tuple(a + b for a, b in zip(t, T.z)); k = [j for j, x in enumerate(wt) if x != 0][0]
    return ('A', 2 + k, wt[k] // 2)
def deriv_coord(sec, ci, d):
    """class coordinate of (d.partial) X = 2i y_{-d} X for the component functional ci"""
    B0, pi, R = D.projection(*sec)
    md = (-d[0], -d[1])
    for i, b in enumerate(B0):
        if b[0] == ci and tuple(b[1]) == (md,):
            return sp.nsimplify(2 * sp.I * pi[i]), [S.name(z) for z in B0]
    raise KeyError
out = {}
for d, lab in ((D.D1, 'd1'), (D.D2, 'd2')):
    out['c_' + lab] = str(deriv_coord(D.sec_dc(d), ci_of(('c',)), d)[0])
    for i in range(3):
        out['gamma%d_%s' % (i, lab)] = str(deriv_coord(D.sec_dgamma(i, d), ci_of(gam_key(i)), d)[0])
        out['beta%d_%s' % (i, lab)] = str(deriv_coord(D.sec_dbeta(i, d), ci_of(('l', T.u + T.taus[i])), d)[0])
print(json.dumps(out, indent=1))
json.dump(out, open('deriv_coords.json', 'w'), indent=1)
# check the class bases are 1-dim and print the hz vector for one gamma and beta sector
for sec, nm in ((D.sec_dgamma(2, D.D2), 'dgamma2_d2'), (D.sec_dbeta(0, D.D1), 'dbeta0_d1')):
    B0, pi, R = D.projection(*sec)
    print(nm, 'basis', [S.name(b) for b in B0], 'pi', pi, 'dimZ', len(R['Z']), 'dimB', len(R['Bim']))
