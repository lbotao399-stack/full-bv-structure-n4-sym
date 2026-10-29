import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json, cq_twist as TW
for N in [0, 1, 2]:
    r = TW.run(N)
    print(json.dumps(r, indent=0)[:3000])
