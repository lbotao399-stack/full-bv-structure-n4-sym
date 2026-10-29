import sympy as sp
import toy_verify as tv
from grassmann import G
def orders(F):
    out = {}
    for k, v in F.t.items():
        p = sp.Poly(sp.expand(v), tv.g, tv.h, tv.d)
        for mon, cf in p.terms():
            out.setdefault(mon, G())
            out[mon] = out[mon] + G({k: cf})
    return out
Wo = orders(tv.W1)
for S_par in (True, False):
    for order in ('ij', 'ji'):
        for strace in (True, False):
            C = orders(tv.candidate(S_par, order, strace))
            bad = []
            for mon in sorted(set(Wo) | set(C)):
                if sum(mon) > 3 or mon == (2, 1, 1):
                    continue
                dd = Wo.get(mon, G()) - C.get(mon, G())
                if not dd.iszero():
                    bad.append(mon)
            print(S_par, order, strace, "mismatch orders (g,h,d):", bad)
