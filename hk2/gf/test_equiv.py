import time, random, sympy as sp
import bvjet, models as Mo
X = bvjet.X
random.seed(7)
def rspin(nz=16):
    v=[0]*16
    for k in random.sample(range(16), nz):
        v[k]=random.choice([-2,-1,1,2])
    return v
b = [1, -2, 1, 3]
b2 = [2, 1, -1, 0]
om = sp.Matrix(4,4,lambda i,j: 0)
om[0,1]=1; om[1,0]=-1; om[2,3]=2; om[3,2]=-2; om[0,3]=1; om[3,0]=-1
xx = sum(t*t for t in X)
gens_even = {
  'P': [sp.sympify(t) for t in b],
  'M': [sp.expand(sum(om[m,n]*X[n] for n in range(4))) for m in range(4)],
  'D': [2*X[m] for m in range(4)],
  'K': [sp.expand(xx*b[m]-2*sum(b[n]*X[n] for n in range(4))*X[m]) for m in range(4)],
  'K2': [sp.expand(xx*b2[m]-2*sum(b2[n]*X[n] for n in range(4))*X[m]) for m in range(4)],
}
r = sp.zeros(6,6); r[0,1]=1; r[1,0]=-1; r[2,4]=2; r[4,2]=-2
r2 = sp.zeros(6,6); r2[1,2]=1; r2[2,1]=-1; r2[0,5]=1; r2[5,0]=-1
xi1, xi2 = rspin(), rspin()
eta1, eta2 = rspin(), rspin()
Z=[0]*16
odd = {'xi1':(xi1,Z), 'xi2':(xi2,Z), 'eta1':(Z,eta1), 'eta2':(Z,eta2)}

def S1(F, key, **kw):
    if key in gens_even: return Mo.S1_even(F, v=gens_even[key], **kw)
    if key in ('R','R2'): return Mo.S1_even(F, r=(r if key=='R' else r2), **kw)
    xi, eta = odd[key]; return Mo.S1_odd(F, xi, eta)

def run(frame, gamma2=None, cbar_transforms=True, verbose=True):
    M = Mo.build_model(frame != 'min')
    F = Mo.Fields(M, frame)
    S0 = Mo.S0(F)
    res = {}
    kw = {'cbar_transforms': cbar_transforms}
    t0=time.time()
    for key in list(gens_even)+['R']+list(odd):
        s1 = S1(F,key,**kw) if (key in gens_even or key=='R') else S1(F,key)
        res['ar1_'+key] = M.is_zero_functional(M.bracket(S0, s1))
    # arity 2, even-even
    def S1v(v): return Mo.S1_even(F, v=v, **kw)
    for (k1,k2) in [('P','K'),('D','K'),('K','K2'),('M','K'),('P','D')]:
        v,w = gens_even[k1], gens_even[k2]
        lhs = M.plus(M.bracket(S1v(v),S1v(w)), S1v(Mo.bracket_even_even(v,w)))
        res['ar2_%s_%s'%(k1,k2)] = M.is_zero_functional(lhs)
    # R,R'
    lhs = M.plus(M.bracket(Mo.S1_even(F,r=r), Mo.S1_even(F,r=r2)), Mo.S1_even(F, r=(r2*r - r*r2)))
    res['ar2_R_R2'] = M.is_zero_functional(lhs)
    # arity 2, even-odd
    for (k1,k2) in [('P','eta1'),('K','xi1'),('D','eta1'),('M','xi1'),('K','eta1')]:
        v = gens_even[k1]; xi,eta = odd[k2]
        xn, en = Mo.bracket_even_odd(v, xi, eta)
        lhs = M.plus(M.bracket(S1v(v), Mo.S1_odd(F,xi,eta)), Mo.S1_odd(F,xn,en))
        res['ar2_%s_%s'%(k1,k2)] = M.is_zero_functional(lhs)
    xi,eta = odd['xi1']
    xn,en = Mo.bracket_even_odd(None, xi, eta, r=r)
    lhs = M.plus(M.bracket(Mo.S1_even(F,r=r), Mo.S1_odd(F,xi,eta)), Mo.S1_odd(F,xn,en))
    res['ar2_R_xi1'] = M.is_zero_functional(lhs)
    # arity 2, odd-odd
    for (k1,k2) in [('xi1','xi2'),('xi1','eta1'),('eta1','eta2'),('xi1','xi1')]:
        (a1,e1),(a2,e2) = odd[k1], odd[k2]
        v, rot = Mo.bracket_odd_odd(a1,e1,a2,e2)
        s2 = Mo.S2_odd(F,a1,e1,a2,e2,gamma2=gamma2)
        lhs = M.plus(M.bracket(S0,s2), S1v(v), Mo.S1_even(F,r=rot),
                     M.scale(M.bracket(Mo.S1_odd(F,a1,e1), Mo.S1_odd(F,a2,e2)), -1))
        res['ar2_%s_%s'%(k1,k2)] = M.is_zero_functional(lhs)
    if verbose: print(frame, 'gamma2=',gamma2, 'cbar_transforms=',cbar_transforms, ' %.1fs'%(time.time()-t0))
    return res

if __name__=='__main__':
    import sys
    r0 = run('min'); print({k:v for k,v in r0.items() if not v} or 'all PASS', len(r0))
