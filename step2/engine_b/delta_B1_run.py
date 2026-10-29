"""B1: the flowed arity-one defect along Delta at two legs (1PI special-edge graphs with S_M on the cycle), canonical form."""
import json, sys, time, random
import sympy as sp
import delta_canon as C
SECT = {'AVAV': ('A', 'A', list(range(4)), list(range(4))), 'ApAp': ('A', 'A', list(range(4, 10)), list(range(4, 10))),
        'lamlam': ('lam', 'lam', list(range(16)), list(range(16))), 'cAsV': ('c', 'As', [0], list(range(4)))}
def rp(sp_, comps, rnd):
    u = [sp.Integer(0)] * C.D.V.DIM[sp_]
    for c in comps: u[c] = sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3]))
    return u
if __name__ == '__main__':
    secs = sys.argv[1].split(','); nsamp = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 23
    rnd = random.Random(seed)
    fn = f'delta_B1_{"_".join(secs)}_s{seed}.json'
    out = {}
    for sec in secs:
        s1, s2, c1, c2 = SECT[sec]; out[sec] = []
        for k in range(nsamp):
            t0 = time.time()
            u1 = rp(s1, c1, rnd); u2 = rp(s2, c2, rnd); p1 = [sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3])) for _ in range(4)]
            res = C.evaluate(dict(sp=s1, u=u1), s2, u2, p1)
            rec = dict(u1=[str(x) for x in u1], u2=[str(x) for x in u2], p1=[str(x) for x in p1],
                       res={str(n): dict(alpha={kk: str(v) for kk, v in res[n]['alpha'].items()},
                                         beta=[{kk: str(v) for kk, v in b.items()} for b in res[n]['beta']]) for n in res},
                       time=round(time.time() - t0, 1))
            out[sec].append(rec); print(sec, k, json.dumps(rec), flush=True)
            json.dump(out, open(fn, 'w'), indent=1)
