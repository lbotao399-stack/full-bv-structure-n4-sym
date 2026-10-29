import pickle, sympy as sp, sys
sys.path.insert(0,'.')
d = pickle.load(open(sys.argv[1],'rb'))
raw = d['raw']['C']
I = sp.I
def ibp(dD, dP, smax):
    bad = []
    cands = set(k[0] for k in list(dD)+list(dP))
    ext = set()
    for s in cands:
        for i in range(len(s)):
            nw = list(s); nw[i] = (nw[i][0], nw[i][1]+1); ext.add(tuple(nw))
    cands |= ext
    for S in range(smax):
        for a in range(S+1):
            b = S-a
            for sites in cands:
                lhs = dD.get((sites,(a+1,b)),0) + dD.get((sites,(a,b+1)),0)
                rhs = 2*dP.get((sites,(a,b)),0)
                for i in range(len(sites)):
                    if sites[i][1] >= 1:
                        nw = list(sites); nw[i] = (nw[i][0], nw[i][1]-1)
                        rhs += I*dD.get((tuple(nw),(a,b)),0)
                if sp.expand(lhs-rhs) != 0:
                    bad.append(((a,b), sites, sp.expand(lhs), sp.expand(rhs)))
    return bad
for pr, v in raw.items():
    N = max([k[1][0]+k[1][1] for k in v['D0']] + [0])
    b = ibp(v['D0'], v['Pb'], N)
    print(pr, 'max letter order', N, 'violations', len(b), b[:2])
