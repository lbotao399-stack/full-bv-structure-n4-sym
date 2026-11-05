import env, time, resource, json
import numpy as np
import qa_twist as TW, qa_core as QC, qa_util as U
t0=time.time()
T=TW.Twist()
g0=QC.zero_gen()
La=TW.letters(T,'beta',0); Lb=TW.letters(T,'gamma',0)
print('taus', T.taus, 'hol', T.hol, 'exa', T.exa)
# outputs: lambar (dc) and ghost c with one derivative
outs_c=[('c', np.array([1+0j]))]
rows=TW.entry(T, g0, La, Lb, ('lamM','c'), [TW.outs(T,'dc'), outs_c], 1, None, log=lambda s: print(s, round(time.time()-t0,1), flush=True), label='betagamma_ar0')
rows2=TW.entry(T, g0, La, Lb, ('c','lamM'), [outs_c, TW.outs(T,'dc')], 1, None, log=lambda s: print(s, round(time.time()-t0,1), flush=True), label='betagamma_ar0b')
json.dump({'r1':rows,'r2':rows2}, open('out_t1.json','w'), indent=1, default=str)
print('rss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6)
