import sys, json, os
import pairfit_core as PC
route = __import__(sys.argv[1])
import fitlib
for pair in sys.argv[2:]:
    log = route.B.log if hasattr(route, 'B') else print
    r = PC.analyse(route, pair, log=lambda s: (route.A.log(s, 'pairfit') if route.name == 'routeA' else route.B.log(s, 'pairfit')))
    json.dump(r, open(os.path.join(route.OUT, f'{route.name}_pair_{pair}_fit.json'), 'w'), indent=1)
