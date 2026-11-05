import env, twdriver as D, sympy as sp, time
t0 = time.time()
r = D.assemble(D.L_beta(0), D.L_gamma(0), sp.Rational(5, 2), 0, D.sec_dc(D.D1), D.sec_dc(D.D2))
print('Konishi beta0 gamma0 -> [d1c][d2c]:  UX =', r['UX'], '  VX =', r['VX'])
for k, v in r['parts'].items(): print('  ', k, v)
print(time.time() - t0)
