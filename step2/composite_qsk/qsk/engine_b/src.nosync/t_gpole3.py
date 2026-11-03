"""eq qs-pole at m = 1 for (A, lam_+) -> A with general colour polarisations: all terms.
RHS = -ct_1({S,X_a} X_b) - ct_1(X_a {S,X_b}_quad) - ct_1(X_a {S,X_b}_cubic) - ct1_{xi,1}(sO)
 {S, A letter} = Gamma(xi, lam) (lam letter z'), {S, lam letter} = (1/2) rho(F) xi (quadratic: c_MN d_M A_N; cubic:
 (1/4) sum [A_M,A_N] zb^T rb_M r_N xi), vacuum ct_1(A,A) = V (za.zb), V = 2 (computed below),
 vacuum ct1_xi(c, lam_+) = Z (xi.zb), Z = 1/2 (counterterm cycle, computed in t_gpole2.py);
 colour: word (0,) <-> (e_b,[Y,e_a]); cubic vacuum -> -ct_1 = +V/2 Acub; sO -> -Z (xi.zb)(za.Y)."""
import sympy as sp, sys, qb_core as B, qb_job as QJ, qb_fsp as F
XIV=F.XIV; xv=sp.Matrix(XIV)
Aslot=sys.argv[1]   # 'AW' or 'AV'
comps=list(B.SLOT[Aslot][1])
# V from the twoletter engine: vacuum kernel of (A,A)
job=F.T.Job(Aslot,Aslot,[]);
try:
    r=job.kernel(keep_sequences=False); print('tl vacuum', {w:{k:str(v) for k,v in d.items()} for w,d in r['total'].items()})
except Exception as e:
    print('tl vacuum not available:', repr(e)[:200])
V=2
dj=QJ.DJob(Aslot,'lamP',[Aslot],'xi',XIV,tag='g3'); dk=dj.kernel(mlist=(1,))
D=F.to_expr(dj, dk['total'][(0,)][1]['1'])
# {S, X_a}: lam letter z'_c = sum_M za_M (xi^T MR_M)_c ; Gamma|_V -> lamM, Gamma|_perp -> lamP
lamslot='lamP' if Aslot=='AW' else 'lamM'
lc=list(B.SLOT[lamslot][1])
zprime={c: sp.expand(sum(sp.Symbol(f'za_{M}')*(xv.T*B.MR[M])[c] for M in comps)) for c in lc}
job=F.T.Job(lamslot,'lamP',[Aslot]); r=job.kernel(keep_sequences=False)
e=F.to_expr(job, r['total'][(0,)][1], {f'za_{c}':f'zl_{c}' for c in lc})
k1=sp.expand(e.subs({sp.Symbol(f'zl_{c}'):zprime[c] for c in lc},simultaneous=True))
kb,_=F.ct1_deriv_pair(Aslot,'AV',Aslot,'b')
cb=F.cMN([sp.Symbol(f'zb_{c}') for c in range(8)],'lamP')
k2=sp.expand(sum(kb.get((0,),{}).get(mn,0)*c for mn,c in cb.items()))
zbv=sp.zeros(16,1)
for c in range(8): zbv[c]=sp.Symbol(f'zb_{c}')
Acub=0
for M in comps:
    for N in comps:
        if M==N: continue
        c3=(zbv.T*B.Rbm[M]*B.Rm[N]*xv)[0]; c3r=(zbv.T*B.Rbm[N]*B.Rm[M]*xv)[0]
        Acub+=sp.Rational(1,2)*(c3-c3r)*sp.Symbol(f'za_{M}')*sp.Symbol(f'o0_{N}')
Acub=sp.expand(Acub)
s1=sp.expand(sum(sp.Symbol(f'za_{M}')*sp.Symbol(f'o0_{M}') for M in comps)*sum(XIV[c]*sp.Symbol(f'zb_{c}') for c in range(8)))
Z=sp.Rational(1,2)
rhs=sp.expand(-k1-k2+sp.Rational(V,2)*Acub-Z*s1)
print(Aslot,': terms D',len(sp.Add.make_args(D)),'k1',len(sp.Add.make_args(k1)),'k2',len(sp.Add.make_args(k2)),'Acub',len(sp.Add.make_args(Acub)))
print('D - rhs == 0 :', sp.expand(D-rhs)==0)
print('D + rhs == 0 :', sp.expand(D+rhs)==0)
a,b,c,d=sp.symbols('a b c d')
E=sp.expand(D-(a*k1+b*k2+c*Acub+d*s1)); gens=sorted(E.free_symbols-{a,b,c,d},key=str)
print('D = a k1 + b k2 + c Acub + d s1:', sp.solve(sp.Poly(E,*gens).coeffs(),[a,b,c,d],dict=True))
