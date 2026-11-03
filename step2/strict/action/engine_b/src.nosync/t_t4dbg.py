import random, numpy as np
from fractions import Fraction as Fr
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_fit as F, eb_t4 as T
random.seed(3)
tq,t3=G.tmpl_xi()
spec=T.SECTORS['+']['AV_lam-']
for ic in range(4):
    xi=np.zeros(16,complex)
    for c in range(8): xi[c]=random.randint(-2,2)
    legs=T.make_legs(spec)
    gen=Y.Gen('xi',0,tq,gvec=xi,t3=t3)
    st={}
    r=Y.evaluate(legs,[gen],mlist=(0,1),stats=st)
    Bs=T.basis_xi(T.symvec(xi),f'd{ic}')
    c=[F.value(Bs[l],legs) for l in ['g1','g2','int(Gamma_V,A_V)','int(Gamma_perp,A_perp)']]
    print(st, {k:{m:complex(v) for m,v in d.items()} for k,d in r.items()}, [complex(float(a),float(b)) for a,b in c])
    print('  legs', legs[0]['u'][:4], legs[1]['u'][8:], legs[0]['p'])
