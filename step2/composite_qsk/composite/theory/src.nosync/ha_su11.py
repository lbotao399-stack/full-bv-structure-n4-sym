"""su(1|1) sector Z = Phi_34 = c3+c4+, psi = a1+ c4+ (B03c oscillators): two-site harmonic action (4.27)-(4.29)
restricted to the sector, then H = sum_l H_{l,l+1} on graded cyclic single-trace states of length L.
Output: ha_su11.json (exact spectra, L <= 6)."""
import json, os, itertools
import sympy as sp
import ha_predict as H
HERE = os.path.dirname(os.path.abspath(__file__))
Z = H.letter('Phi', 0, cs=('c3', 'c4')); PSI = {'a1': 1, 'c4': 1}
LET = {'Z': Z, 'p': PSI}
PAR = {'Z': 0, 'p': 1}
def key(x): return H.canon_key(x)
two = {}
for a, b in itertools.product('Zp', repeat=2):
    img = H.H12((LET[a], LET[b]))
    for k, v in img.items():
        for c, d in itertools.product('Zp', repeat=2):
            if key((LET[c], LET[d])) == k:
                two[(a + b, c + d)] = v
def ham_word(w):
    """H acting on a word (open chain sites 0..L-1, plus the bond (L-1,0) with the cyclic sign)"""
    L = len(w); out = {}
    for l in range(L):
        m = (l + 1) % L
        if m != 0:
            pre, pair, post = w[:l], w[l:l + 2], w[l + 2:]
            for (src, dst), v in two.items():
                if src == pair:
                    nw = pre + dst + post
                    out[nw] = out.get(nw, 0) + v
        else:
            # bond (L-1, 0): rotate the last letter to the front with the Koszul sign, act on (0,1), rotate back
            last = w[-1]; rest = w[:-1]
            sg = (-1) ** (PAR[last] * sum(PAR[x] for x in rest))
            rw = last + rest
            for (src, dst), v in two.items():
                if src == rw[:2]:
                    nw2 = dst + rw[2:]
                    first = nw2[0]; rest2 = nw2[1:]
                    sg2 = (-1) ** (PAR[first] * sum(PAR[x] for x in rest2))
                    nw = rest2 + first
                    out[nw] = out.get(nw, 0) + sg * sg2 * v
    return out
def spectrum(L, M):
    words = [''.join(w) for w in itertools.product('Zp', repeat=L) if w.count('p') == M]
    idx = {w: i for i, w in enumerate(words)}; N = len(words)
    Hm = sp.zeros(N, N); C = sp.zeros(N, N)
    for j, w in enumerate(words):
        for nw, v in ham_word(w).items():
            Hm[idx[nw], j] += v
        first = w[0]; rest = w[1:]
        sg = (-1) ** (PAR[first] * sum(PAR[x] for x in rest))
        C[idx[rest + first], j] = sg
    # cyclic sector: eigenvalue 1 of the graded shift C
    Pr = sum((C ** k for k in range(L)), sp.zeros(N, N)) / L
    cols = Pr.columnspace()
    if not cols: return []
    B = sp.Matrix.hstack(*cols)
    Hc = (B.T * B).inv() * B.T * Hm * B
    return sorted([(str(sp.nsimplify(k)), m) for k, m in Hc.eigenvals().items()], key=lambda t: float(sp.sympify(t[0])))
res = {'two_site_H12': {'%s->%s' % k: str(v) for k, v in two.items()}}
for L in range(2, 7):
    for M in range(0, L + 1):
        res['L%d_M%d' % (L, M)] = spectrum(L, M)
json.dump(res, open(os.path.join(HERE, 'ha_su11.json'), 'w'), indent=1)
print(json.dumps(res, indent=0))
