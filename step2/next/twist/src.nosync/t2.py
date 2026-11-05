import env, sector as S, sympy as sp, time
t0=time.time()
# test: Q of iota_{d1} A
ci = [i for i,c in enumerate(S.COMPS) if c[1]==('A',1,1)][0]   # plane 1 (so4 second), s=+1 -> weight (0,1) = d1
print(S.COMPS[ci][1], S.COMPS[ci][3])
r = S.apply_M(S.comp_functional(S.COMPS[ci]))
d = S.decompose(r)
for k,v in d.items(): print(S.COMPS[k][1], S.COMPS[k][3], sp.factor(v))
print(time.time()-t0)
