import sympy as sp, qb_core as B, qb_job as QJ, sys, time
xi=[0]*16; xi[3]=1
sa,sb,outs = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
za={c:(1 if c==0 else 0) for c in B.SLOT[sa][1]} if len(sys.argv)<5 else None
t0=time.time()
j=QJ.DJob(sa,sb,outs,'xi',xi,tag='bg')
r=j.kernel(mlist=(0,),by_chain=True)
print('total words', {w:{m:{b:len(P.terms()) for b,P in d.items()} for m,d in dm.items()} for w,dm in r['total'].items()}, round(time.time()-t0,1),'s', round(QJ.peak_mb()),'MB')
for seq,sig,cres in r['chains']:
    print(seq,sig,{w:{m:{b:len(P.terms()) for b,P in d.items()} for m,d in dm.items()} for w,dm in cres.items()})
