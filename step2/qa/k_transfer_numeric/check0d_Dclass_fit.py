exec(open('check0d_special_triangle.py').read().split('# ---- D-class')[0])
D0 = brK(I3, I3, 0) * sp.Rational(1, 2)
LD = deg2(L1(D0, J2)); PRD = deg2(X(dP(D0, P, V), J2)); onepiD = deg2(LD - PRD)
def deg1(F): return GA({k: sp.expand(c) for k, c in F.t.items() if sum(1 for v in k if v != GAM) == 1 and sp.expand(c) != 0})
KS3, PS3, PS2, KS2 = KSH(K, H3), KSH(P, H3), KSH(P, H2), KSH(K, H2)
cands = {
 'str(K3 P2 P3)': strace(mm(mm(KS3, PS2), PS3)),
 'str(P2 P3 K3)': strace(mm(mm(PS2, PS3), KS3)),
 'str(P3 K3 P2)': strace(mm(mm(PS3, KS3), PS2)),
 'str(K3 P3 P2)': strace(mm(mm(KS3, PS3), PS2)),
 'str(P3 P2 K3)': strace(mm(mm(PS3, PS2), KS3)),
 'str(P2 K3 P3)': strace(mm(mm(PS2, KS3), PS3)),
}
T1 = deg1(strace(mm(PS3, PS2)) * sp.Rational(1, 2))
cands['{T1,I3}_K'] = brK(T1, I3, 1)
cands['{I3,T1}_K'] = brK(I3, T1, 0)
syms = sp.symbols('x0:%d' % len(cands))
resid = onepiD
for s_, (nm, c) in zip(syms, cands.items()):
    resid = resid - deg2(c) * s_
sol = sp.solve([sp.expand(c) for c in resid.t.values()], list(syms), dict=True)
print({nm: (sol[0].get(s_, s_) if sol else None) for s_, nm in zip(syms, cands)})
# also: is the special tadpole alone equal to the difference?
d1 = deg2(onepiD - deg2(cands['str(K3 P2 P3)']) * sp.Rational(1, 2))
print('residual after half str(K3 P2 P3):', len(d1.t))
d2 = deg2(onepiD + deg2(cands['str(K3 P2 P3)']) * sp.Rational(1, 2))
print('residual after -half str(K3 P2 P3):', len(d2.t))
