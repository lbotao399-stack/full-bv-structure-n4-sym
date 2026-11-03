import numpy as np
exec(open('t_gs.py').read().split("cases=")[0])
import cq_gcore as G
from fractions import Fraction
def dump(ta,tb,louts):
    Xa=(ta,letvec(ta),None); Xb=(tb,letvec(tb),None)
    legs=[(t,legpol(t),rng.integers(-2,3,4)) for t in louts]
    for perm in [legs, legs[::-1]]:
        for st, blks, P0s, P1s, scale, odd, spos, sig, lines in QC.configs(perm, gx):
            K, qs, P0 = G.chain(blks, P0s, P1s, lines, ta, tb)
            if not K: continue
            Sx=G.GSimplex(qs); r=Sx.evaluate(K,'fin',s=sig)[0]
            v=QC._from_rl([None if x is None else (x[0]*Fraction(1,scale),x[1]*Fraction(1,scale)) for x in r], len(qs[0]), Xa[1], Xb[1])
            sg=QC.config_sign(odd,spos,sig,len(blks),gx.par)
            print([b[0] for b in st], [p[0] for p in perm], 'spos',spos,'sig',sig,'sign',sg,'val',CV.vstr(v))
print('--- AW AW -> c lamsP'); dump('AW','AW',['c','lamsP'])
print('--- lamP AW -> lamM lamM'); dump('lamP','AW',['lamM','lamM'])
