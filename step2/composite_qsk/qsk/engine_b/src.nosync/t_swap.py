import sympy as sp, qb_core as B, qb_job as QJ, sys
xi=[0]*16; xi[3]=1; xi[4]=sp.Rational(1,3)
sa,sb,outs = sys.argv[1], sys.argv[2], sys.argv[3].split(',')
def ex(job,P,swap):
    syms=[]
    for n in job.names:
        if swap and n.startswith('za_'): n='zb_'+n[3:]
        elif swap and n.startswith('zb_'): n='za_'+n[3:]
        syms.append(sp.Symbol(n))
    return sp.expand(P.as_expr(*syms)) if P else 0
j1=QJ.DJob(sa,sb,outs,'xi',xi,tag='sw1'); r1=j1.kernel(mlist=(0,1))['total']
j2=QJ.DJob(sb,sa,outs,'xi',xi,tag='sw1'); r2=j2.kernel(mlist=(0,1))['total']
pa=B.V.PAR[B.SLOT[sa][0]]; pb=B.V.PAR[B.SLOT[sb][0]]
ok=True
for w,dm in r1.items():
    for m,d in dm.items():
        for b,P in d.items():
            wr=tuple(reversed(w))
            Q=r2.get(wr,{}).get(m,{}).get(b)
            e1=ex(j1,P,False); e2=ex(j2,Q,True) if Q else 0
            s=(-1)**(pa*pb+len(w))
            d1=sp.expand(e2-s*e1); d2=sp.expand(e2+s*e1)
            print(w,m,b,'swap-identity holds:',d1==0,' opposite:',d2==0, 'nterms',len(P.terms()))
            ok = ok and d1==0
print('G-swap', ok)
