import time
import qa_twist as TW
t0=time.time()
T=TW.Twist()
print('xi', T.xi); print('gam_w', T.gam_w); print('dc_out', T.dc_out)
res,T=TW.run(which=('betagamma','bstar_c','betabeta'), log=lambda s: print(s, round(time.time()-t0,1), flush=True))
