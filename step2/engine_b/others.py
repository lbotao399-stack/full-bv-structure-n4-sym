"""Engine B, Step 2, Part C: the tree-level contributions of the three-leg defect flow, sector by sector.
  II-a      = -Q I^CT_3                          = -log eps * unit * Q Y_3
  III_K0    = -{I_3, I^CT_2}                     = -unit * [ (9/4) eps^-1 {I_3, int(A,A)} + log eps {I_3, Y_2} ]
  I-b + III_Khat = -B_Khat(I_3, W_1^(2) - I^CT_2),  Khat = (e^{-eps k^2} - 1) x fibre of K_0:
      eps^0     : I-b  +(9/4) unit B_{k^2 F}(I_3, int(A,A)),   III_Khat  -(9/4) unit B_{k^2 F}(I_3, int(A,A))
      eps log eps: I-b  + unit B_{k^2 F}(I_3, Y_2),            III_Khat  - unit B_{k^2 F}(I_3, Y_2)
      sum = -eps B_{k^2 F}(I_3, renormalised self-energy) + O(eps^2 log eps)  -> 0
  I-d + II-b = B_{P(eps,infty)}(D_2(eps), I_3),  D_2(eps) = special bubbles - Q I^CT_2 = Ob_2 + O(eps)  ->  B_{P(0,infty)}(Ob_2, I_3)
with Y_3 = 6 AAA + 10 lambda lambda A + AcA* + c lambda lambda* + ccc* and Y_2 = 4 kinA + 8 kinlambda - kingh (Step 1,
Theorem c4-log; checked here as {Z, S}), Ob_2 = -(27/8) int(d^*dc, d^*A|_V) (Theorem c4-obstruction).
All values: coefficient of unit * ([c_0,c_1],c_2), projected on G1..G5 with residual."""
import json, os, sys, time
import sympy as sp
import fa
import sectors3 as S
import vertices as V

I = sp.I
HERE = os.path.dirname(os.path.abspath(__file__))
Y3 = fa.local([(6, 'AAA'), (10, 'lamlamA'), (1, 'AcAs'), (1, 'clamlams'), (1, 'cccs')], 'Y_3', 0)
Y2 = fa.local([(4, 'kinA'), (8, 'kinlam'), (-1, 'kingh')], 'Y_2', 0)
Zfun = fa.local([(2, 'ZA'), (-4, 'Zlam'), (-1, 'Zc')], 'Z', 1)
fa.add_term('dsdc_dsA', 1, ('c', 'A'), lambda s: V.dot4(s[0]['p'], s[0]['p']) * s[0]['u'][0] * (-I * V.dot4(s[1]['p'], s[1]['u'])), ('k', (0, 1)))
OB2 = fa.local([(sp.Rational(-27, 8), 'dsdc_dsA')], 'Ob_2', 1)
K0 = fa.kernel_K0()
Kk2 = fa.kernel_scaled(K0, lambda q2: q2)
P0 = fa.kernel_P(0)

def fval(F, legs):
    d = F(legs)
    assert set(d) <= {'f'}, d
    return sp.nsimplify(sp.expand(d.get('f', 0)))

ITEMS = {
    'QY3': lambda legs: fval(fa.Q(Y3), legs),
    '{I3,Y2}': lambda legs: fval(fa.bracket(fa.I3, Y2), legs),
    '{I3,int(A,A)}': lambda legs: fval(fa.bracket(fa.I3, fa.MASS), legs),
    'B_k2F(I3,int(A,A))': lambda legs: fval(fa.Bfun(Kk2, fa.I3, fa.MASS), legs),
    'B_k2F(I3,Y2)': lambda legs: fval(fa.Bfun(Kk2, fa.I3, Y2), legs),
    'tree B_P0(Ob2,I3)': lambda legs: fval(fa.Bfun(P0, OB2, fa.I3), legs),
    '{Z,I3}-Y3': lambda legs: fval(fa.bracket(Zfun, fa.I3), legs) - fval(Y3, legs),
}

NCONF = {"cAVAV": 10, "cApAp": 5, "clamlam": 6, "ccAsV": 4}
def run(sector, nconf=None):
    nconf = nconf or NCONF[sector]
    vals = {k: [] for k in ITEMS}
    rows_struct = []
    for c in range(nconf):
        legs = S.config(sector, c)
        st = S.structures(sector, legs, 0)
        gv = {g: S.G_value(g, legs) for g in S.SECTOR_G[sector]}
        for k, f in ITEMS.items():
            vals[k].append((f(legs), st, gv))
    out = {}
    for k, rows in vals.items():
        names = sorted(rows[0][1])
        sol, res = S.solve_exact([(r[0], r[1]) for r in rows], names)
        gsol, gres = S.solve_exact([(r[0], r[2]) for r in rows], S.SECTOR_G[sector])
        out[k] = dict(G={g: str(v) for g, v in gsol.items()}, residual_G=str(max(abs(sp.N(x)) for x in gres)),
                      structures={n: str(v) for n, v in sol.items() if v != 0}, residual_full=str(max(abs(sp.N(x)) for x in res)))
    return out

if __name__ == '__main__':
    res = {}
    t0 = time.time()
    for sec in S.SECTORS:
        res[sec] = run(sec)
        print(sec, json.dumps({k: (v['G'], v['residual_G']) for k, v in res[sec].items()}), round(time.time() - t0, 1), flush=True)
    res['seconds'] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(HERE, 'others.json'), 'w'), indent=1)
