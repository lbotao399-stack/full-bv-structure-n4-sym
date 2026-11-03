import time, numpy as np, itertools
from fractions import Fraction
t0=time.time()
import qa_s1 as S1, qa_core as QC
import cq_eval as CV
rng = np.random.default_rng(7)
G5S = S1.G5S; G5SS = S1.G5SS
def rv(n): return np.array([complex(int(a),int(b)) for a,b in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
def letvec(t):
    if t=='AW': return rv(6)
    if t=='AV': return rv(4)
    if t=='c': return np.array([1])
    if t=='lamP': return (np.eye(16)+G5S)@rv(16)
    if t=='lamM': return (np.eye(16)-G5S)@rv(16)
    if t=='lamsP': return (np.eye(16)+G5SS)@rv(16)
    if t=='lamsM': return (np.eye(16)-G5SS)@rv(16)
    if t=='AsV': return rv(4)
    if t=='AsW': return rv(6)
def legpol(t): return letvec(t)
B = S1.chir_basis('lamP')
xi = sum(int(c)*v for c,v in zip([1,2,-1,1,0,1,-2,1],B))
gx = QC.odd_gen(xi)
def gs(ta,tb,louts, npts=3):
    out=[]
    for _ in range(npts):
        Xa=(ta,letvec(ta),None); Xb=(tb,letvec(tb),None)
        legs=[(t,legpol(t),rng.integers(-2,3,4)) for t in louts]
        vab=QC.pair_value(Xa,Xb,legs,gx)['fin']
        vba=QC.pair_value(Xb,Xa,legs[::-1],gx)['fin']
        out.append((CV.vstr(vab),CV.vstr(vba)))
    return out
cases=[('AW','lamP',['AW']),('AW','AW',['AW','lamP']),('AW','AW',['AV','lamM']),('AW','AW',['c','lamsP']),
       ('lamP','lamP',['AW','lamM']),('lamP','AW',['lamM','lamM'])]
for ta,tb,lo in cases:
    r=gs(ta,tb,lo)
    print(ta,tb,lo)
    for a,b in r: print('    ab:',a,' | ba:',b)
    print('   t',round(time.time()-t0,1))
