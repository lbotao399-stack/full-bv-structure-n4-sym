"""Numerical I-a at finite eps (no Taylor expansion), sector (c, A|_Wperp, A|_Wperp), config 0 and 1: approach to the exact
eps^0 value (11/4 + (19/4) log(3/2)) G3 (eps^-1 part 0)."""
import json, time, math, sympy as sp, mpmath
import numeric_check as NC, ia_engine as IA, sectors3 as S
out = {}
alpha3 = 11 / 4 + 19 / 4 * math.log(1.5)
for cid in (0, 1):
    legs = S.config('cApAp', cid)
    terms = []
    IA.evaluate(legs, collect=terms)
    G3 = complex(S.G_value('G3', legs))
    raw = json.load(open('ia_raw_cApAp.json'))[str(cid)]['classes']
    res = {}
    for eps in (1e-4, 1e-5, 1e-6):
        t0 = time.time()
        tot, pc = NC.ia_numeric(terms, eps, 20, 200)
        res[str(eps)] = dict(I_a_over_G3=str(tot / G3), seconds=round(time.time() - t0, 1))
        print(cid, eps, tot / G3, 'exact limit', alpha3, flush=True)
    out[str(cid)] = dict(G3_value=str(G3), exact_limit_over_G3=f'11/4 + 19/4*log(3/2) = {alpha3:.12f}', runs=res)
json.dump(out, open('numeric_check_ia_cApAp.json', 'w'), indent=1)
