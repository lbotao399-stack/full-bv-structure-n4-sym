"""vacuum-pole x spectator terms of eq qs-pole for general colour polarisations, pair (A_perp, lam_+), leg set A_perp."""
import sympy as sp, qb_core as B, qb_job as QJ, qb_fsp as F
from sympy.polys.domains import QQ_I
XIV=F.XIV
q=sp.symbols('q0:4')
# Z: vacuum counterterm cycle ct1_xi(c, lam_+) (all edges P), coefficient of eps^-1 in units hbar/16pi^2
j=QJ.DJob('c','lamP',[],'xi',XIV,tag='z')
tot=0
for seq in j.sequences():
    types,curs,words=j.config(seq,-1)
    for w,P in words.items():
        e=-sp.expand(P.as_expr(*[sp.Symbol(n) for n in j.names]))   # remove the special-edge minus
        tot+=e
tot=sp.expand(tot)
# Num(q) = kappa |q|^2 ; [eps^-1] of int e^{-2 eps q^2} Num/q^4 = kappa/2
q2=sum(sp.Symbol(f'q{m}')**2 for m in range(4))
kappa=sp.cancel(tot/q2)
Z=sp.expand(kappa/2)
print('ct1_xi(c,lam+)|_0 [eps^-1] =', Z)
xi=XIV[:8]; zb=[sp.Symbol(f'zb_{c}') for c in range(8)]
print('  / (xi.zb) =', sp.simplify(Z/sum(a*b for a,b in zip(xi,zb))))
# cubic-vacuum structure: sum_{M,N in perp} c3^[MN] za_M Y_N, c3_MN = zb^T Rb_M R_N xi
xv=sp.Matrix(XIV); zbv=sp.zeros(16,1)
for c in range(8): zbv[c]=zb[c]
Acub=0
for M in range(4,10):
    for N in range(4,10):
        if M==N: continue
        c3=(zbv.T*B.Rbm[M]*B.Rm[N]*xv)[0]; c3r=(zbv.T*B.Rbm[N]*B.Rm[M]*xv)[0]
        Acub+=sp.Rational(1,2)*(c3-c3r)*sp.Symbol(f'za_{M}')*sp.Symbol(f'o0_{N}')
Acub=sp.expand(Acub)
# k1 and D from t_gpole
dj=QJ.DJob('AW','lamP',['AW'],'xi',XIV,tag='tg'); dk=dj.kernel(mlist=(1,))
D=F.to_expr(dj, dk['total'][(0,)][1]['1'])
zprime={c: sp.expand(sum(sp.Symbol(f'za_{4+K}')*(xv.T*B.MR[4+K])[c] for K in range(6))) for c in range(8)}
job=F.T.Job('lamP','lamP',['AW']); r=job.kernel(keep_sequences=False)
e=F.to_expr(job, r['total'][(0,)][1], {f'za_{c}':f'zl_{c}' for c in range(8)})
k1=sp.expand(e.subs({sp.Symbol(f'zl_{c}'):zprime[c] for c in range(8)},simultaneous=True))
za=[sp.Symbol(f'za_{c}') for c in range(4,10)]; Y=[sp.Symbol(f'o0_{c}') for c in range(4,10)]
s1=sp.expand(sum(a*b for a,b in zip(za,Y))*sum(a*b for a,b in zip(xi,zb)))
a,b,c=sp.symbols('a b c')
E=sp.expand(D-(a*k1+b*Acub+c*s1))
gens=sorted(E.free_symbols-{a,b,c},key=str)
print('D = a k1 + b Acub + c (za.Y)(xi.zb):', sp.solve(sp.Poly(E,*gens).coeffs(),[a,b,c],dict=True))
E2=sp.expand(k1-(a*Acub+c*s1)); gens=sorted(E2.free_symbols-{a,c},key=str)
print('k1 = a Acub + c (za.Y)(xi.zb):', sp.solve(sp.Poly(E2,*gens).coeffs(),[a,c],dict=True))
