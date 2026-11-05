import env, json, time, twdriver as D, sympy as sp
out = {}
specs = [('betagamma_11', D.L_beta(1), D.L_gamma(1), D.sec_dc(D.D1), D.sec_dc(D.D2)),
         ('betagamma_22', D.L_beta(2), D.L_gamma(2), D.sec_dc(D.D1), D.sec_dc(D.D2)),
         ('betagamma_01', D.L_beta(0), D.L_gamma(1), D.sec_dc(D.D1), D.sec_dc(D.D2)),
         ('betabeta_12_c_g0', D.L_beta(1), D.L_beta(2), D.sec_dc(D.D1), D.sec_dgamma(0, D.D2)),
         ('betabeta_12_c_g0_d2d1', D.L_beta(1), D.L_beta(2), D.sec_dc(D.D2), D.sec_dgamma(0, D.D1)),
         ('bgamma_1', D.L_b(), D.L_gamma(1), D.sec_dc(D.D1), D.sec_dgamma(1, D.D2)),
         ('gammagamma_01_c_c', D.L_gamma(0), D.L_gamma(1), D.sec_dc(D.D1), D.sec_dc(D.D2)),
         ('betac_0', D.L_beta(0), D.L_c(), D.sec_dc(D.D1), D.sec_dgamma(0, D.D2))]
for nm, La, Lb, sU, sV in specs:
    t0 = time.time()
    try:
        r = D.assemble3(La, Lb, sU, sV, log=lambda s: None)
        out[nm] = {'UX': str(r['UX']), 'VX': str(r['VX'])}
    except Exception as e:
        out[nm] = {'error': repr(e)[:300]}
    print(nm, out[nm], round(time.time() - t0, 1), flush=True)
json.dump(out, open('entries_flavour.json', 'w'), indent=1)
