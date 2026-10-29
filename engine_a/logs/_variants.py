exec(open('verify_oneloop0d.py').read().split('# variants')[0])
import itertools
def trace(X): return sum((X[(a, a)] for a in V), GA())
def hessv(I, order):
    if order == 0: return {(b, c): I.dleft(c).dleft(b) for b in V for c in V}
    return {(b, c): I.dleft(b).dleft(c) for b in V for c in V}
def PSHv(H, sidx, transP):
    M = {}
    for a in V:
        for c in V:
            s = GA()
            for b in V:
                pab = P.get((b, a), 0) if transP else P.get((a, b), 0)
                if pab != 0:
                    sg = (-1) ** PAR[b] if sidx == 0 else ((-1) ** PAR[a] if sidx == 1 else ((-1)**PAR[c] if sidx==2 else 1))
                    s = s + H[(b, c)] * (pab * sg)
            M[(a, c)] = s
    return M
for order, sidx, transP, tr in itertools.product([0, 1], [0, 1, 2, 3], [0, 1], [0, 1]):
    X3 = PSHv(hessv(I3, order), sidx, transP); X4 = PSHv(hessv(I4, order), sidx, transP)
    T = strace if tr else trace
    f = (T(X4) * sp.Rational(1, 2) + T(mm(X3, X3)) * sp.Rational(1, 4)).degree_part(2)
    d = (W1q - f).map(sp.simplify)
    f2 = (T(X4) * sp.Rational(-1, 2) + T(mm(X3, X3)) * sp.Rational(-1, 4)).degree_part(2)
    d2 = (W1q - f2).map(sp.simplify)
    print(order, sidx, transP, tr, len(d.t), len(d2.t))
