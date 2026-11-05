import json, time
from test_equiv import run
out = {}
t0=time.time()
for (frame,g2,ct,label) in [('C',-1,True,'C_transforming_gamma2=-1'),
                            ('Psi',-1,True,'Psi_transforming_gamma2=-1'),
                            ('Psi',0,True,'NEG_Psi_transforming_gamma2=0'),
                            ('Psi',0,False,'Psi_inert_gamma2=0'),
                            ('Psi',-1,False,'NEG_Psi_inert_gamma2=-1')]:
    res = run(frame, gamma2=g2, cbar_transforms=ct)
    out[label] = res
    print(label, 'failures:', [k for k,v in res.items() if not v], flush=True)
    json.dump(out, open('equiv_results.json','w'), indent=1)
print('total %.0fs'%(time.time()-t0))
