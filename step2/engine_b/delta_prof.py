import cProfile, pstats, time, sympy as sp
import delta_canon as C
D = C.D
u1 = [0]*4 + [-1,-3,-3,2,-1,1]; u2 = [0]*4 + [1,2,-1,-2,3,-2]
u1 = [sp.Integer(x) for x in u1]; u2 = [sp.Integer(x) for x in u2]
p1 = [sp.Integer(x) for x in (-1,1,-3,-2)]
t0 = time.time()
legs = C.build_legs(dict(sp='A', u=u1), 'A', u2, p1, [0]*4, False, [0]*4)
recs = []
for part in C.partitions(legs, D.SW_TERMS, False):
    recs += D.enumerate_terms(legs, part, 'K', [0, 1])
t1 = time.time()
print('enumerate w:', len(recs), 'records', round(t1-t0,1), 's')
for r in recs[:40]:
    D.special_coeffs_multi(r['Num'], r['curs'], r['s'], {'w': lambda e: e})
print('special 40:', round(time.time()-t1,1))
