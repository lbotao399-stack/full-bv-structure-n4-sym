import pickle, sympy as sp, sys
d = pickle.load(open(sys.argv[1], 'rb'))
hP, hD, one = d['hP'], d['hD'], d['one']
Smax = int(sys.argv[2])
def get(h, site, ins):
    return h.get((site, ins), 0)
bad = []
for S in range(0, Smax):
    for a in range(S + 1):
        b = S - a
        for c in range(S + 2):
            dd = S + 1 - c
            lhs = get(hD, (c, dd), (a + 1, b)) + get(hD, (c, dd), (a, b + 1))
            rhs = (get(hD, (c - 1, dd), (a, b)) if c >= 1 else 0) + (get(hD, (c, dd - 1), (a, b)) if dd >= 1 else 0) + 2 * get(hP, (c, dd), (a, b))
            if sp.simplify(lhs - rhs) != 0:
                bad.append(((a, b), (c, dd), lhs, rhs))
print('two-output IBP identity violations:', len(bad))
for x in bad[:20]:
    print(x)
# one-output
for tag in ('AZ', 'ZA'):
    P, D = one[tag + 'P'], one[tag + 'D']
    bad1 = []
    for S in range(0, Smax - 1):
        for a in range(S + 1):
            b = S - a
            c = S + 2
            lhs = D.get(((c,), (a + 1, b)), 0) + D.get(((c,), (a, b + 1)), 0)
            rhs = D.get(((c - 1,), (a, b)), 0) + 2 * P.get(((c,), (a, b)), 0)
            if sp.simplify(lhs - rhs) != 0:
                bad1.append(((a, b), lhs, rhs))
    print(tag, 'one-output IBP violations', bad1)
