exec(open('check0d_special_triangle.py').read().split('# ---- D-class')[0])
def rpoly4(n):
    F = GA(); ms = [m for m in monos(4) if sum(PAR[v] for v in m) % 2 == 0]; random.shuffle(ms)
    for m in ms[:n]:
        t = GA.const(rnd())
        for v in m: t = t * GA.var(v)
        F = F + t
    return F
I4 = rpoly4(12); H4 = hess(I4)
D = brK(I4, J2, 0)
lhs = deg2(dP(D, P, V))
rhs = deg2((strace(mm(KSH(K, H4), KSH(P, H2))) + strace(mm(KSH(K, H2), KSH(P, H4)))) * sp.Rational(1, 2))
# I4 self-contraction with K inside D would need a field-antifield pair in I4 at one vertex; tree-leg graphs absent at deg 2
x = sp.symbols('x')
T0 = GA({k: c for k, c in (dP(I4, P, V)).t.items()})   # P-tadpole of I4 (quadratic)
extra = deg2(brK(T0, J2, 0))
res = lhs - rhs - extra * x
sol = sp.solve([sp.expand(c) for c in res.t.values()], [x], dict=True)
print('quartic-J2 special bubble: d_P {I4,J2}_K = half[str(K S H4 P S H2)+str(K S H2 P S H4)] + x*{d_P I4, J2}_K :', sol, len(lhs.t), len(rhs.t))
fl = sp.symbols('fl')
