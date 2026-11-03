"""Gate G-reg0: arity-zero obstruction from the special-edge cycles of I alone (two and three legs)."""
import json, os, random
from fractions import Fraction as Fr
import numpy as np
import eb_core as E, eb_cycle as Y, eb_fit as F, eb_qb as Q

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')
L = '(log3 - log2)'


def mom3():
    p1 = [Fr(random.randint(-3, 3)) for _ in range(4)]; p2 = [Fr(random.randint(-3, 3)) for _ in range(4)]
    return [-(p1[m] + p2[m]) for m in range(4)], p1, p2


def gq(z):
    return (Fr(z.real), Fr(z.imag))


def run(seed=5, log=print):
    random.seed(seed)
    out = {}
    # two legs (c, A_V): basis int(d*dc, d*A_V); pole: Q of int(A,A)
    vals0, vals1, c0, c1 = [], [], [], []
    for _ in range(3):
        p = [Fr(random.randint(-3, 3)) for _ in range(4)]
        a = np.zeros(10, complex); a[:4] = [random.randint(-2, 2) for _ in range(4)]
        legs = [dict(sp='c', u=np.array([1], complex), p=p), dict(sp='A', u=a, p=[-x for x in p])]
        r = Y.evaluate(legs, [], mlist=(0, 1))
        vals0.append(r.get('k', {}).get(0, E.CLog())); vals1.append(r.get('k', {}).get(1, E.CLog()))
        c0.append([Q.to_fr(Q.F_dsdc_dsA(legs))])
        c1.append([Q.to_fr(Q.apply_Q(Q.F_AA)(legs))])
    f0, ok0 = F.fit(vals0, c0, ['int(d*dc,d*A_V)'])
    f1, ok1 = F.fit(vals1, c1, ['Q int(A,A)'])
    out['two_legs'] = dict(eps0={k: {l: str(x) for l, x in d.items()} for k, d in f0.items()}, ok0=ok0,
                           pole={k: {l: str(x) for l, x in d.items()} for k, d in f1.items()}, ok1=ok1)
    log(f'G-reg0 two legs: {out["two_legs"]}')
    # three legs: G1..G5 of Prop c9-basis, values with f = ([X0,X1],X2)
    def dot(a, b):
        return sum(complex(a[m]) * complex(b[m]) for m in range(4))
    R, MP, RB = Q.R, Q.MP, Q.RB
    sectors = {}
    for name, nconf in (('cAVAV', 4), ('cApAp', 3), ('clamlam', 4), ('ccAsV', 3)):
        vals, cands = [], []
        for _ in range(nconf):
            p0, p1, p2 = mom3()
            if name == 'cAVAV':
                a1 = np.zeros(10, complex); a2 = np.zeros(10, complex)
                a1[:4] = [random.randint(-2, 2) for _ in range(4)]; a2[:4] = [random.randint(-2, 2) for _ in range(4)]
                legs = [dict(sp='c', u=np.array([1], complex), p=p0), dict(sp='A', u=a1, p=p1), dict(sp='A', u=a2, p=p2)]
                G1 = (dot(p2, p2) - dot(p1, p1)) * dot(a1, a2) - dot(a1, p2) * dot(a2, p2) + dot(a1, p1) * dot(a2, p1)
                G2 = dot(a1, p2) * dot(a2, p2) - dot(a1, p1) * dot(a2, p1)
                cands.append([gq(G1), gq(G2)]); labs = ['G1', 'G2']
            elif name == 'cApAp':
                w1 = np.zeros(10, complex); w2 = np.zeros(10, complex)
                w1[4:] = [random.randint(-2, 2) for _ in range(6)]; w2[4:] = [random.randint(-2, 2) for _ in range(6)]
                legs = [dict(sp='c', u=np.array([1], complex), p=p0), dict(sp='A', u=w1, p=p1), dict(sp='A', u=w2, p=p2)]
                cands.append([gq(np.dot(w1, w2) * (dot(p2, p2) - dot(p1, p1)))]); labs = ['G3']
            elif name == 'clamlam':
                l1 = np.array([random.randint(-2, 2) for _ in range(16)], complex)
                l2 = np.array([random.randint(-2, 2) for _ in range(16)], complex)
                legs = [dict(sp='c', u=np.array([1], complex), p=p0), dict(sp='lam', u=l1, p=p1), dict(sp='lam', u=l2, p=p2)]
                rho = lambda p: sum(complex(p[m]) * R[m] for m in range(4))
                G5 = -1j * (l1 @ MP @ rho(p0) @ l2)
                odd = 1j * (l1 @ MP @ rho([p2[m] - p1[m] for m in range(4)]) @ RB[0] @ R[1] @ RB[2] @ R[3] @ l2)
                cands.append([gq(G5), gq(odd)]); labs = ['G5', 'parity-odd']
            else:
                b = np.zeros(10, complex); b[:4] = [random.randint(-2, 2) for _ in range(4)]
                legs = [dict(sp='c', u=np.array([1], complex), p=p0), dict(sp='c', u=np.array([1], complex), p=p1),
                        dict(sp='As', u=b, p=p2)]
                cands.append([gq(-1j * dot(p2, b))]); labs = ['G4']
            r = Y.evaluate(legs, [], mlist=(0, 1))
            vals.append(r.get('f', {}).get(0, E.CLog()))
            assert r.get('f', {}).get(1, E.CLog()).is_zero()
        f, ok = F.fit(vals, cands, labs)
        sectors[name] = dict(eps0={k: {l: str(x) for l, x in d.items()} for k, d in f.items()}, consistent=ok)
        log(f'G-reg0 three legs {name}: {sectors[name]}')
    out['three_legs'] = sectors
    return out


if __name__ == '__main__':
    with open(os.path.join(OUT, 'G_reg0.json'), 'w') as f:
        json.dump(run(), f, indent=1)
