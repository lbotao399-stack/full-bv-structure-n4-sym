"""K5: 1PI one-loop graphs with one O vertex, one S_M vertex and one special edge Khat = K_eps - K_0 on the cycle,
two external legs (sector given), canonical form with leg 2 eliminated (momentum -p1 - k - kappa), omega momentum k."""
import json, sys, time, random
import sympy as sp
import delta_canon as C
D = C.D
SECT = {'ApAp': ('A', 'A', list(range(4, 10)), list(range(4, 10))), 'AVAV': ('A', 'A', list(range(4)), list(range(4)))}
def rp(sp_, comps, rnd):
    u = [sp.Integer(0)] * D.V.DIM[sp_]
    for c in comps: u[c] = sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3]))
    return u
if __name__ == '__main__':
    secs = sys.argv[1].split(','); nsamp = int(sys.argv[2]); seed = int(sys.argv[3]) if len(sys.argv) > 3 else 31
    rnd = random.Random(seed)
    fn = f'delta_K5_{"_".join(secs)}_s{seed}.json'
    out = {}
    for sec in secs:
        s1, s2, c1, c2 = SECT[sec]; out[sec] = []
        for k in range(nsamp):
            t0 = time.time()
            w1 = rp('A', range(4, 10), rnd); w2 = rp('A', range(4, 10), rnd)
            D.set_O(w1, w2)
            u1 = rp(s1, c1, rnd); u2 = rp(s2, c2, rnd)
            p1 = [sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3])) for _ in range(4)]
            kv = [sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3])) for _ in range(4)]
            res = C.evaluate_scaled4(dict(sp=s1, u=u1), s2, u2, p1, 0, kvec=kv, om=True)
            rec = dict(w1=[str(x) for x in w1], w2=[str(x) for x in w2], u1=[str(x) for x in u1], u2=[str(x) for x in u2],
                       p1=[str(x) for x in p1], k=[str(x) for x in kv],
                       res={str(n): dict(alpha={kk: str(v) for kk, v in res[n]['alpha'].items()},
                                         beta=[{kk: str(v) for kk, v in b.items()} for b in res[n]['beta']]) for n in res},
                       time=round(time.time() - t0, 1))
            out[sec].append(rec); print(sec, k, json.dumps(rec), flush=True)
            json.dump(out, open(fn, 'w'), indent=1)
