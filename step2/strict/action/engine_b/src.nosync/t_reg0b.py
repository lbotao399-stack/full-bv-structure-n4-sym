import time, sys, random
from fractions import Fraction as Fr
import numpy as np
import eb_core as E, eb_cycle as Y
import ce_core as C
R=[np.array(C.Rm[M].tolist(),dtype=complex) for M in range(10)]; Rb=[np.array(C.Rbm[M].tolist(),dtype=complex) for M in range(10)]
Mp=np.array(C.Mp.tolist(),dtype=complex)
random.seed(5)
def rvec(n): return [Fr(random.randint(-3,3)) for _ in range(n)]
def dot(a,b): return sum(complex(a[m])*complex(b[m]) for m in range(4))
def mom3():
    p1=rvec(4); p2=rvec(4); p0=[-(p1[m]+p2[m]) for m in range(4)]; return p0,p1,p2
L=np.log(1.5)
def run(legs):
    t=time.time(); st={}
    r=Y.evaluate(legs,[],mlist=(0,1),cyclic=True,stats=st)
    return r, time.time()-t, st
# (c, A_V, A_V)
rows=[]
for trial in range(3):
    p0,p1,p2=mom3()
    a1=np.zeros(10,complex); a2=np.zeros(10,complex); a1[:4]=[random.randint(-2,2) for _ in range(4)]; a2[:4]=[random.randint(-2,2) for _ in range(4)]
    legs=[dict(sp='c',u=np.array([1],complex),p=p0),dict(sp='A',u=a1,p=p1),dict(sp='A',u=a2,p=p2)]
    r,tt,st=run(legs)
    v=r.get('f',{}); v0=complex(v.get(0,E.CLog())); v1=complex(v.get(1,E.CLog()))
    a1p1=dot(a1,p1);a1p2=dot(a1,p2);a2p1=dot(a2,p1);a2p2=dot(a2,p2);a12=dot(a1,a2)
    G1=(dot(p2,p2)-dot(p1,p1))*a12-a1p2*a2p2+a1p1*a2p1; G2=a1p2*a2p2-a1p1*a2p1
    pred=(11/4+19/4*L)*G1+27/8*G2
    print('cAVAV', 'eps0', v0, 'pred', pred, 'pole', v1, round(tt,2), st, flush=True)
for trial in range(2):
    p0,p1,p2=mom3()
    w1=np.zeros(10,complex); w2=np.zeros(10,complex); w1[4:]=[random.randint(-2,2) for _ in range(6)]; w2[4:]=[random.randint(-2,2) for _ in range(6)]
    legs=[dict(sp='c',u=np.array([1],complex),p=p0),dict(sp='A',u=w1,p=p1),dict(sp='A',u=w2,p=p2)]
    r,tt,st=run(legs)
    v=r.get('f',{}); v0=complex(v.get(0,E.CLog())); v1=complex(v.get(1,E.CLog()))
    G3=np.dot(w1,w2)*(dot(p2,p2)-dot(p1,p1))
    print('cApAp', 'eps0', v0, 'pred', (11/4+19/4*L)*G3, 'pole', v1, round(tt,2), st, flush=True)
for trial in range(2):
    p0,p1,p2=mom3()
    l1=np.array([random.randint(-2,2) for _ in range(16)],complex); l2=np.array([random.randint(-2,2) for _ in range(16)],complex)
    legs=[dict(sp='c',u=np.array([1],complex),p=p0),dict(sp='lam',u=l1,p=p1),dict(sp='lam',u=l2,p=p2)]
    r,tt,st=run(legs)
    v=r.get('f',{}); v0=complex(v.get(0,E.CLog())); v1=complex(v.get(1,E.CLog()))
    rho=lambda p: sum(complex(p[m])*R[m] for m in range(4))
    G5=-1j*(l1@Mp@rho(p0)@l2)
    odd=1j*(l1@Mp@rho([p2[m]-p1[m] for m in range(4)])@Rb[0]@R[1]@Rb[2]@R[3]@l2)
    print('clamlam', 'eps0', v0, 'pred', -9*L*G5, 'odd', odd, 'pole', v1, round(tt,2), st, flush=True)
for trial in range(2):
    p0,p1,p2=mom3()
    b=np.zeros(10,complex); b[:4]=[random.randint(-2,2) for _ in range(4)]
    legs=[dict(sp='c',u=np.array([1],complex),p=p0),dict(sp='c',u=np.array([1],complex),p=p1),dict(sp='As',u=b,p=p2)]
    r,tt,st=run(legs)
    v=r.get('f',{}); v0=complex(v.get(0,E.CLog())); v1=complex(v.get(1,E.CLog()))
    print('ccAs', 'eps0', v0, 'pred 0', 'G4', -1j*dot(p2,b), 'pole', v1, round(tt,2), st, flush=True)
