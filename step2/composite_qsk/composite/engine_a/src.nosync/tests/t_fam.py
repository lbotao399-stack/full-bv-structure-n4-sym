import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cq_core, cq_families as FM
t0=time.time()
r=FM.so6_density(); print(json.dumps(r)[:600], '%.1fs'%(time.time()-t0))
r=FM.su11_spectra(Lmax=6); print(json.dumps(r['two_site_lc_natural']), json.dumps(r['two_site_Xi_Delta_translation_invariant'])); print({k:v['pass'] for k,v in r['spectra'].items()}, '%.1fs'%(time.time()-t0))
r=FM.universality(npts=1); print('%.1fs'%(time.time()-t0))
