import env, time, resource
t0=time.time()
import qa_twist as TW, qa_core as QC
T=TW.Twist()
gx=QC.odd_gen(T.xi)
res,T=TW.run(which=('betagamma',), log=lambda s: print(s, round(time.time()-t0,1), flush=True))
print('rss MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6)
