import sympy as sp, qb_core as B, qb_job as QJ, qb_fsp as F
XIV=F.XIV
dj=QJ.DJob('AW','lamP',['AW'],'xi',XIV,tag='tg'); dk=dj.kernel(mlist=(1,),by_chain=True)
D=F.to_expr(dj, dk['total'][(0,)][1]['1'])
xv=sp.Matrix(XIV)
zprime={c: sp.expand(sum(sp.Symbol(f'za_{4+K}')*(xv.T*B.MR[4+K])[c] for K in range(6))) for c in range(8)}
job=F.T.Job('lamP','lamP',['AW']); r=job.kernel(keep_sequences=False)
e=F.to_expr(job, r['total'][(0,)][1], {f'za_{c}':f'zl_{c}' for c in range(8)})
k1=sp.expand(e.subs({sp.Symbol(f'zl_{c}'):zprime[c] for c in range(8)},simultaneous=True))
kb,_=F.ct1_deriv_pair('AW','AV','AW','b')
cb=F.cMN([sp.Symbol(f'zb_{c}') for c in range(8)],'lamP')
k2=sp.expand(sum(kb.get((0,),{}).get(mn,0)*c for mn,c in cb.items()))
za=[sp.Symbol(f'za_{c}') for c in range(4,10)]; Y=[sp.Symbol(f'o0_{c}') for c in range(4,10)]
zb=[sp.Symbol(f'zb_{c}') for c in range(8)]; xi=XIV[:8]
s1=sp.expand(sum(a*b for a,b in zip(za,Y))*sum(a*b for a,b in zip(xi,zb)))
def fit(E,name):
    c=sp.Symbol('c'); r=sp.expand(E-c*s1)
    sol=sp.solve(sp.Poly(r,*(za+Y+zb+xi)).coeffs(),c,dict=True)
    print(name,'= c*(za.Y)(xi.zb) with', sol if sol else 'not of that form')
print('nterms D',len(sp.Add.make_args(D)),'k1',len(sp.Add.make_args(k1)),'k2',len(sp.Add.make_args(k2)))
fit(sp.expand(D+k1+k2),'D - rhs')
fit(sp.expand(D+k1),'D + k1')
fit(sp.expand(D+k2),'D + k2')
fit(k1,'k1'); fit(k2,'k2')
# chain-by-chain of D: which configurations carry the (za.Y)(xi.zb) structure
for seq,sig,cres in dk['chains']:
    P=cres.get((0,),{}).get(1,{}).get('1')
    if P:
        Ee=F.to_expr(dj,P); print(seq,sig,len(sp.Add.make_args(Ee))); fit(Ee,'   chain')
