"""QA (Part Q): project engine B's TQ3 beta-beta output onto the twisted letters.
Reads engine B's out.nosync/tq3_twist_n0.json (read-only) and B's realization (read-only import).
Twisted outputs: A_perp output on gamma^K = A_{K4} (dual vectors w.r.t. the letter basis {A_{K4}, A_{IJ}}),
lam_- output on lambar_{alphadot 4} (components 11, 15, already restricted by B), output momenta restricted to the
exact span E = span{Gamma(xi, xi')|_V} (this kills every exact derivative u^alpha sigma^m d_m and keeps the
holomorphic ones).  Prints the projected values."""
import os, sys, json, itertools, random
import sympy as sp

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
sys.path.insert(0, os.path.join(ROOT, 'step2/composite_qsk/composite/engine_b/src.nosync/vendor'))
import realization as RZ

Mp = RZ.pairing()
xi = sp.zeros(16, 1); xi[3] = 1
# gamma letters: w_K = Sig[K][J-1,3]; check rho(w) xi = 0
gam = []
for J in (1, 2, 3):
    w = [RZ.Sig[K][J - 1, 3] for K in range(6)]
    v = sp.zeros(16, 1)
    for K in range(6):
        v += w[K] * RZ.R[4 + K] * xi
    assert v == sp.zeros(16, 1), ('gamma not annihilating', J)
    gam.append(w)
# complement letters A_{IJ}, I<J in 1..3
comp = []
for (I_, J) in ((1, 2), (1, 3), (2, 3)):
    comp.append([RZ.Sig[K][I_ - 1, J - 1] for K in range(6)])
# check the complement does not annihilate xi
for w in comp:
    v = sp.zeros(16, 1)
    for K in range(6):
        v += w[K] * RZ.R[4 + K] * xi
    assert v != sp.zeros(16, 1)
L = sp.Matrix(gam + comp)          # rows: letter functionals on the 6 components
Linv = L.inv()                      # columns: dual vectors
dual_gamma = [Linv[:, k] for k in range(3)]
# exact span E
E = []
for c in range(16):
    xp = sp.zeros(16, 1); xp[c] = 1
    vec = [(xi.T * Mp * RZ.R[M] * xp)[0] for M in range(4)]
    E.append(vec)
Em = sp.Matrix(E)
rk = Em.rank()
Ebasis = [list(r) for r in Em.T.columnspace()]
print('exact span rank', rk, 'basis', Ebasis)
# isotropy check
for a in Ebasis:
    for b in Ebasis:
        assert sp.simplify(sum(x * y for x, y in zip(a, b))) == 0
data = json.load(open(os.path.join(ROOT, 'step2/composite_qsk/qsk/engine_b/out.nosync/tq3_twist_n0.json')))
bb = data['beta beta (lam_+, lam_+) -> A_perp lamM']
random.seed(7)
out = {}
for pair, ent in bb.items():
    if isinstance(ent, str):
        out[pair] = ent
        continue
    res = {}
    for key, s in ent.items():
        e = sp.sympify(s)
        vals = []
        for K in range(3):
            for trial in range(2):
                sub = {}
                for c in range(4, 10):
                    sub[sp.Symbol('o0_%d' % c)] = dual_gamma[K][c - 4]
                # lam_- output: keep 11, 15 with random weights
                sub[sp.Symbol('o1_11')] = random.randint(1, 5)
                sub[sp.Symbol('o1_15')] = random.randint(1, 5)
                for j in (0, 1):
                    a, b = random.randint(-4, 4), random.randint(-4, 4)
                    r = [a * Ebasis[0][m] + b * Ebasis[1][m] for m in range(4)]
                    for m in range(4):
                        sub[sp.Symbol('r%d_%d' % (j, m))] = r[m]
                vals.append(str(sp.nsimplify(sp.expand(e.subs(sub)))))
        # also: values on the complement (non-twisted) A_perp outputs, holomorphic momenta
        valsc = []
        for K in range(3):
            sub = {}
            for c in range(4, 10):
                sub[sp.Symbol('o0_%d' % c)] = Linv[:, 3 + K][c - 4]
            sub[sp.Symbol('o1_11')] = 1; sub[sp.Symbol('o1_15')] = 2
            for j in (0, 1):
                r = [Ebasis[0][m] + 2 * Ebasis[1][m] for m in range(4)]
                for m in range(4):
                    sub[sp.Symbol('r%d_%d' % (j, m))] = r[m]
            valsc.append(str(sp.nsimplify(sp.expand(e.subs(sub)))))
        res[key] = {'twisted_gamma_holomorphic': vals, 'complement_A_IJ_holomorphic': valsc}
    out[pair] = res
json.dump(out, open(os.path.join(ROOT, 'step2/composite_qsk/qsk/qa/out.nosync/bb_projection.json'), 'w'), indent=1)
for p, r in out.items():
    if isinstance(r, str):
        print(p, r); continue
    for k, v in r.items():
        print(p, k, 'gamma:', v['twisted_gamma_holomorphic'][:6], '| A_IJ:', v['complement_A_IJ_holomorphic'])
