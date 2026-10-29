"""Translation part of K5 (beta^mu, mu given): split into the subset predicted by the closed form
LHS_P - (D_P O)[infty] = lim (hbar/2) d_{K''} (Q O_W)_tree  [T between two P edges (an A-c line), Khat adjacent to O]
and the complement (must sum to zero if both exact decompositions hold)."""
import json, sys, time
import sympy as sp
import delta_canon as C
D = C.D
I = sp.I
def run(sample, mu):
    w1 = [sp.sympify(x) for x in sample['w1']]; w2 = [sp.sympify(x) for x in sample['w2']]
    D.set_O(w1, w2)
    u1 = [sp.sympify(x) for x in sample['u1']]; u2 = [sp.sympify(x) for x in sample['u2']]
    p1 = [sp.sympify(x) for x in sample['p1']]; kv = [sp.sympify(x) for x in sample['k']]
    p1s = [C.LAMB * x for x in p1]; ks_ = [C.LAMB * x for x in kv]
    legs = C.build_legs(dict(sp='A', u=u1), 'A', u2, p1s, ks_, True, [0] * 4)
    tot = {'subset': 0, 'complement': 0}
    detail = {}
    for part in C.partitions(legs, D.T_TERMS(mu), True):
        for r in D.enumerate_terms(legs, part, 'K', [0, 1]):
            n = len(r['word'])
            v = C.special_scaled(C.trunc_lam(r['Num'], 1), r['curs'], r['s'], {'b': lambda e: sp.expand(-I * e)}, 1)
            val = r['fac'] * (v['b'][0]['1'] + 0 * v['b'][0]['L'])
            Lval = r['fac'] * v['b'][0]['L']
            kinds = [w[0] for w in r['word']]
            iT = kinds.index('SM'); iO = kinds.index('O')
            s = r['s']
            T_edges = (r['edges'][iT], r['edges'][(iT + 1) % n])
            T_onP = (iT != s and (iT + 1) % n != s)
            Kadj_O = (s == iO or s == (iO + 1) % n)
            key = 'subset' if (T_onP and Kadj_O) else 'complement'
            tot[key] += val
            cl = f'T_onP={T_onP},Kadj_O={Kadj_O},' + ','.join(kinds) + '|' + ','.join(r['edges']) + f'|s={s}'
            detail[cl] = sp.expand(detail.get(cl, 0) + val)
            if Lval != 0:
                detail[cl + '|L'] = sp.expand(detail.get(cl + '|L', 0) + Lval)
    return {k: str(sp.expand(v)) for k, v in tot.items()}, {k: str(v) for k, v in detail.items() if v != 0}
if __name__ == '__main__':
    d = json.load(open(sys.argv[1]))['ApAp'][int(sys.argv[2])]
    mu = int(sys.argv[3])
    t0 = time.time()
    tot, det = run(d, mu)
    out = dict(sample=int(sys.argv[2]), mu=mu, totals=tot, by_class=det, k=d['k'], time=round(time.time() - t0, 1))
    json.dump(out, open(f'delta_K5_class_s{sys.argv[2]}_mu{mu}.json', 'w'), indent=1)
    print(json.dumps(out, indent=1))
