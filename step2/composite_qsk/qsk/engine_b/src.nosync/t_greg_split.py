import sympy as sp, qb_core as B, qb_job as QJ, qb_greg as G
from sympy.polys.domains import QQ_I
bsym=list(sp.symbols('b0:4'))
job=QJ.DJob('AW','AW',['AW','AW'],'Pb',bsym,tag='greg2')
res=job.kernel(mlist=(0,),by_chain=True)
R=job.R; g=job.g
comps=range(4,10)
za=[g[f'za_{c}'] for c in comps]; zb=[g[f'zb_{c}'] for c in comps]
Y0=[g[f'o0_{c}'] for c in comps]; Y1=[g[f'o1_{c}'] for c in comps]
dot=lambda x,y: sum((a*b for a,b in zip(x,y)),R.zero)
p0=[g[f'r0_{m}'] for m in range(4)]; p1=[g[f'r1_{m}'] for m in range(4)]
bb=[g[f'b{m}'] for m in range(4)]
S={'(a.Y0)(b.Y1)':dot(za,Y0)*dot(zb,Y1),'(a.Y1)(b.Y0)':dot(za,Y1)*dot(zb,Y0),'(a.b)(Y0.Y1)':dot(za,zb)*dot(Y0,Y1)}
basis={}
for kn,sv in S.items():
    for j,pj in ((0,p0),(1,p1)):
        basis[f'i(b.p{j}) {kn}']=sv*dot(bb,pj)*R(QQ_I(0,1))
acc={}
for seq,sig,cres in res['chains']:
    spos=[i for i,(kd,bl) in enumerate(seq) if kd.startswith('S')][0]
    key=('s<sig' if spos<sig else 's>sig', 'k_I=%d'%(len(seq)-1), 'adj' if sig in (spos,spos+1) else 'line')
    for w,dm in cres.items():
        P=dm.get(0,{}).get('1')
        if P:
            acc.setdefault((key,w),R.zero)
            acc[(key,w)]+=P
for (key,w),P in sorted(acc.items(), key=str):
    print(key,w,G.decompose(job,P,basis))
