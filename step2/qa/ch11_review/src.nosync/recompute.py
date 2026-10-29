"""Independent cheap recomputations for the Chapter 11 review (no source code reused).
1. k=1 weights W_0, W_m on [0,1/3] (Prop c11-k1) from the region beta smallest in the 2-simplex.
2. Bilinear composite pair kernel arithmetic (Ex c11-bilinear, Prop c11-scalar, Lemma c11-onebody, Thm c11-P7 (i)).
3. Twist-two two-site matrix of [B03c,(3.14)] (Prop c10-twist): eigenvalues of H'_12+H'_21 on cyclic states, D_1=-8 E.
4. su(1|1) spectra: sum_l (1 - graded exchange) on graded cyclic states of tr(Z^{L-M} psi^M), L<=6.
5. Harmonic numbers for the D_1 tables."""
import json, itertools
from fractions import Fraction as Fr
import sympy as sp
from pathlib import Path
OUT = Path(__file__).resolve().parent.parent / 'recompute.json'
res = {}
b, g = sp.symbols('beta gamma', positive=True)
# 1. W_0: region beta<=gamma', beta<=gamma''=1-beta-gamma'  -> gamma' in [beta, 1-2beta], beta<=1/3
W0_sym = 2*sp.integrate(sp.log(g/b), (g, b, (1-b)/2))
W0_claim = (1-b)*sp.log((1-b)/(2*b)) - (1-3*b)
res['W0_symmetric_form_equals_claim'] = sp.simplify(sp.expand_log(W0_sym - W0_claim, force=True)) == 0
# numeric check of the min form at a few beta
import mpmath as mp
def W0_num(bb):
    bb = mp.mpf(bb.numerator)/bb.denominator
    return mp.quad(lambda x: mp.log(min(x, 1-bb-x)/bb), [bb, (1-bb)/2, 1-2*bb])
res['W0_numeric'] = [[float(bb), float(W0_num(bb)), float(W0_claim.subs(b, bb))] for bb in (Fr(1,10), Fr(1,5), Fr(3,10))]
# beta smallest requires beta<=1/3: at beta=1/3 interval [1/3,1/3] degenerates
res['W0_at_one_third'] = sp.nsimplify(W0_claim.subs(b, sp.Rational(1,3)))
# W_m
m = sp.symbols('m', positive=True, integer=True)
for mm in (1, 2, 3):
    lhs = sp.Rational(1, mm)*2*(sp.integrate(g**mm, (g, b, (1-b)/2)) - b**mm*((1-b)/2 - b))
    claim = sp.Rational(2, mm)*((((1-b)/2)**(mm+1) - b**(mm+1))/(mm+1) - b**mm*((1-b)/2 - b))
    res[f'W{mm}_matches'] = sp.simplify(lhs - claim) == 0
    # direct min form numerically
    bb = mp.mpf(1)/7
    num = mp.quad(lambda x: (min(x, 1-bb-x)**mm - bb**mm)/mm, [bb, (1-bb)/2, 1-2*bb])
    res[f'W{mm}_numeric_vs_claim_at_1/7'] = [float(num), float(claim.subs(b, sp.Rational(1, 7)))]
# moments of W0 over [0,1/3] lie in Q + Q log2 + Q log3
mom = {}
for a_ in range(3):
    for c_ in range(3):
        v = sp.integrate(b**a_*(1-b)**c_*W0_claim, (b, 0, sp.Rational(1, 3)))
        v = sp.nsimplify(sp.expand(sp.expand_log(sp.simplify(v), force=True)))
        mom[f'b^{a_}(1-b)^{c_}'] = str(v)
res['W0_moments'] = mom
ok = True
for k, v in mom.items():
    e = sp.sympify(v)
    rest = sp.expand(e).subs({sp.log(2): 0, sp.log(3): 0})
    coeffs_ok = all(isinstance(t, (sp.Rational, sp.Integer)) or t.is_rational for t in [sp.expand(e).coeff(sp.log(2)), sp.expand(e).coeff(sp.log(3)), rest])
    ok = ok and coeffs_ok
res['W0_moments_in_Q_log2_log3'] = ok
# 2. bilinear arithmetic (units)
wO = 2  # two letters A_perp of weight 1
C_Delta_trinv = Fr(-1, 3)*2*wO - Fr(4, 3)
res['bilinear_translation_invariant_part'] = str(C_Delta_trinv)  # expect -8/3
res['bilinear_letter_count_identity_c'] = str(Fr(-2, 3)*(2+2))       # expect -8/3
res['bilinear_route_value_traceless'] = str(16 - Fr(4, 3))           # Ch9 table first row 44/3
res['bilinear_route_value_trace'] = str(-32 - Fr(4, 3))              # -100/3
res['bilinear_Pb'] = str(Fr(-1, 3)*2/2)                               # position part -2/3 d_b, half -> -1/3
res['Pb_composite_triangle_plus_two_point'] = str(-2 + Fr(5, 3))
# per adjacent pair from Ch10 Prop c10-scalarsector and Rem c9-Pi: -v. + (1/2)(5/3) v. per letter share; -2/3 X per letter halved
res['scalar_pair_v_coefficient'] = str(-1 + Fr(5, 6))   # -1/6
res['scalar_pair_constant'] = str(2*Fr(-2, 3)/2)        # -2/3
res['scalar_pair_trinv'] = str(Fr(-1, 6)*2*2 + Fr(-2, 3))  # -4/3
# Lemma c11-onebody at theta=-1/3: c_phi + (2/3) w = (2/3)(2+|alpha|)
chk = []
for cphi, w0 in ((Fr(2, 3), 1), (Fr(1, 3), Fr(3, 2)), (Fr(4, 3), 0)):
    for al in range(4):
        chk.append(cphi + Fr(2, 3)*(w0+al) == Fr(2, 3)*(2+al))
res['onebody_universal'] = all(chk)
# c_phi from Ch10 table: -4*weight + c_phi = -22/3, -47/3, -8/3 with field-counting weights 2,4,1
res['onebody_cphi_from_ch10_table'] = [str(Fr(-22, 3)+8), str(Fr(-47, 3)+16), str(Fr(-8, 3)+4)]
# 3. twist two, Beisert H'_12 [B03c (3.14)] as quoted in Prop c10-twist
def h(n):
    return sum(Fr(1, j) for j in range(1, n+1))
tw = {}
for Sp in range(0, 9):
    n = Sp+1
    H12 = sp.zeros(n, n)
    for k in range(n):
        for kp in range(n):
            H12[kp, k] = (h(k)+h(Sp-k)) if k == kp else -Fr(1, abs(k-kp))
    # basis |k> = D^k Z D^{S'-k} Z /(k!(S'-k)!); H21 = swap H12 swap; swap: k -> S'-k
    Pm = sp.zeros(n, n)
    for k in range(n):
        Pm[Sp-k, k] = 1
    H = H12 + Pm*H12*Pm
    # cyclic (symmetric) states: v with Pm v = v
    sym = [sp.Matrix([1 if (i == k or i == Sp-k) else 0 for i in range(n)]) for k in range(n//2+1) if k <= Sp-k]
    B = sp.Matrix.hstack(*sym)
    # restrict: H B = B M
    M = (B.T*B).inv()*B.T*H*B
    assert sp.simplify(H*B - B*M) == sp.zeros(n, len(sym))
    ev = sorted([sp.nsimplify(e) for e, mlt in M.eigenvals().items() for _ in range(mlt)])
    tw[Sp] = {'E': [str(e) for e in ev], 'D1': [str(-8*e) for e in ev],
              'target_4h(S)_S_even': [str(4*h(S)) for S in range(0, Sp+1, 2)]}
res['twist_two'] = tw
res['minus32h'] = {S: str(-32*h(S)) for S in (2, 4, 6, 8)}
# 4. su(1|1): letters Z (even), psi (odd); graded cyclic states; H = sum_l (1 - Pi_{l,l+1}) graded exchange
def su11(L, M):
    words = [w for w in itertools.product((0, 1), repeat=L) if sum(w) == M]
    idx = {w: i for i, w in enumerate(words)}
    N = len(words)
    def exch(w, l):
        l2 = (l+1) % L
        if l2 == 0:
            # bond between site L-1 and 0 across the trace: handle via cyclic rotation sign
            return None
        w2 = list(w); w2[l], w2[l2] = w2[l2], w2[l]
        s = -1 if (w[l] == 1 and w[l2] == 1) else 1
        return tuple(w2), s
    # cyclic shift C: tr(X1...XL) -> tr(XL X1...X_{L-1}) with sign (-1)^{|XL|(|X1|+...+|X_{L-1}|)}
    def shift(w):
        s = -1 if (w[-1] == 1 and sum(w[:-1]) % 2 == 1) else 1
        return (w[-1],)+w[:-1], s
    Cm = sp.zeros(N, N)
    for w in words:
        w2, s = shift(w); Cm[idx[w2], idx[w]] += s
    Hm = sp.zeros(N, N)
    for w in words:
        for l in range(L-1):
            w2, s = exch(w, l)
            Hm[idx[w], idx[w]] += 1
            Hm[idx[w2], idx[w]] -= s
    # bond (L-1,0): conjugate the bond (0,1) by the cyclic shift: H_{L,1} = C^{-1} H_{1,2} C? use H_{l,l+1} = C^l H_{0,1} C^{-l}
    H01 = sp.zeros(N, N)
    for w in words:
        w2, s = exch(w, 0)
        H01[idx[w], idx[w]] += 1
        H01[idx[w2], idx[w]] -= s
    Hfull = sp.zeros(N, N)
    Ck = sp.eye(N)
    Cinv = Cm.inv()
    for l in range(L):
        Hfull += Ck*H01*Ck.inv()
        Ck = Ck*Cinv
    # projector onto graded cyclic states: average of C^k
    P = sp.zeros(N, N); Ck = sp.eye(N)
    for k in range(L):
        P += Ck; Ck = Ck*Cm
    P = P/L
    cols = P.columnspace()
    if not cols:
        return None
    B = sp.Matrix.hstack(*cols)
    Mr = (B.T*B).inv()*B.T*Hfull*B
    assert sp.simplify(Hfull*B - B*Mr) == sp.zeros(N, B.shape[1])
    return sorted([sp.nsimplify(sp.simplify(e)) for e, mlt in Mr.eigenvals().items() for _ in range(mlt)], key=lambda e: float(e))
s11 = {}
for L in range(2, 7):
    for M in range(0, L+1):
        r = su11(L, M)
        s11[f'({L},{M})'] = None if r is None else [str(e) for e in r]
res['su11'] = s11
# Heisenberg check of the same code path with bosonic letters (Ch10 Prop c10-su2 (4,2): 0,6)
OUT.write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str))
