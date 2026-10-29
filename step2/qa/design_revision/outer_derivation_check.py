# Design revision (2026-10-28): checks for the outer derivation B of psl(4|4) used in
# step2/theory/odd_components_design.tex, Prop prop:o-laws (e) and Prop prop:o-GB.
#
# (1) In the exact realisation of Chapter 2 (superconformal/algebra.py, read-only, cR = 1 as in verify.py):
#     split the odd part Sigma (+) Sigma^* into the four chirality blocks
#       Sigma_+ = Omega-eigenspace +1 of Sigma  (S_+ (x) W),  Sigma_- (S_- (x) W^*),
#       Sigma^*_a, Sigma^*_b (eigenspaces of the chirality operator on Sigma^*),
#     and verify the bracket pattern that makes B := (+1 on one block pair, -1 on the other, 0 on the even part)
#     an even derivation:  [x,y] = 0 whenever x,y lie in blocks with the same B-value; [even, block] stays in the
#     B-eigenspace.  The block pair with B = +1 is the one containing Sigma_+ (S_+ (x) W); its Sigma^* partner is
#     the block whose bracket with Sigma_+ vanishes (expected: S_- (x) W).
# (2) Non-exactness of G_B: the linear system for a field-type rescaling N (weights a_V, a_W, a_lam, a_lamb, a_c)
#     with {S^(1)(xi), N} = kappa S^(1)(xi) term by term (eq:o-S1xi field content); the cyclic sum of the four
#     linear-term conditions forces kappa = 0.
import os, sys, json, random
import sympy as sp
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
sys.path.insert(0, os.path.join(ROOT, 'superconformal'))
import model as Mo
import algebra as Al
Al.set_cR(sp.Integer(1))

out = {}
SIG, SIGB = Al.SIGM, Al.SIGBM
Om = SIGB[0] * SIG[1] * SIGB[2] * SIG[3]          # chirality on Sigma (twist_coeff_check.py)
Oms = SIG[0] * SIGB[1] * SIG[2] * SIGB[3]         # candidate chirality on Sigma^*
I16 = sp.eye(16)
def eigbasis(M, ev):
    ns = (M - ev * I16).nullspace()
    return [sp.Matrix(v) for v in ns]
blocks = {}
for name, M in (('Sig', Om), ('Sigs', Oms)):
    evs = sorted({sp.nsimplify(e) for e in M.eigenvals().keys()}, key=lambda z: (sp.re(z), sp.im(z)))
    for e in evs:
        blocks['%s[%s]' % (name, e)] = eigbasis(M, e)
out['block_dimensions'] = {k: len(v) for k, v in blocks.items()}
rng = random.Random(20261028)
def rvec(basis):
    v = sp.zeros(16, 1)
    for b in basis:
        v += sp.Integer(rng.randint(-3, 3) or 1) * b
    return v
def mk(name, v):
    return Al.Qgen(list(v)) if name.startswith('Sig[') else Al.Sgen(list(v))
names = sorted(blocks)
samples = {n: [mk(n, rvec(blocks[n])) for _ in range(2)] for n in names}
pattern = {}
for i, n1 in enumerate(names):
    for n2 in names[i:]:
        zero = all(Al.bracket(a, b).is_zero() for a in samples[n1] for b in samples[n2])
        pattern['[%s,%s]' % (n1, n2)] = 'zero' if zero else 'nonzero'
out['odd_odd_bracket_pattern'] = pattern
# identify the B=+1 pair: Sigma_+ is Om-eigenvalue +1 (S_+ (x) W, ch01 Prop c1-chirality)
plus_Sig = [n for n in names if n.startswith('Sig[') and n.endswith('[1]')]
assert len(plus_Sig) == 1
plus_Sig = plus_Sig[0]
partners = [n for n in names if n.startswith('Sigs') and pattern.get('[%s,%s]' % tuple(sorted((plus_Sig, n))), pattern.get('[%s,%s]' % (n, plus_Sig))) == 'zero']
out['Sigma_+_commuting_Sigma*_block'] = partners
Bplus = [plus_Sig] + partners
Bminus = [n for n in names if n not in Bplus]
def key(a, b):
    k = '[%s,%s]' % (a, b)
    return k if k in pattern else '[%s,%s]' % (b, a)
same_block_zero = all(pattern[key(a, b)] == 'zero' for blk in (Bplus, Bminus) for a in blk for b in blk)
out['B_blocks'] = {'+1': Bplus, '-1': Bminus}
out['brackets_within_each_B_block_vanish'] = bool(same_block_zero)
# even generators preserve the B-eigenspaces
def which_block(o):
    res = []
    for n in names:
        vec = o.Q if n.startswith('Sig[') else o.eta
        other = o.eta if n.startswith('Sig[') else o.Q
        if all(sp.expand(z) == 0 for z in other):
            if all(sp.expand(z) == 0 for z in vec):
                continue
            Mb = sp.Matrix.hstack(*blocks[n])
            try:
                Mb.gauss_jordan_solve(vec)
                res.append(n)
            except ValueError:
                pass
    return res
evens = {'P': Al.P([1, 2, -1, 3]), 'K': Al.Kb([2, -1, 1, 1]), 'D': Al.Dil(),
         'M': Al.Mrot([[0, 1, 2, 0], [-1, 0, 0, 1], [-2, 0, 0, 3], [0, -1, -3, 0]]),
         'R': Al.Rsym({(4, 5): 1, (5, 4): -1, (6, 9): 2, (9, 6): -2})}
preserve = {}
for en, ev in evens.items():
    for n in names:
        img = Al.bracket(ev, samples[n][0])
        if img.is_zero():
            preserve['[%s,%s]' % (en, n)] = 'zero'
            continue
        # split image into Q-part and eta-part and locate the blocks
        parts = []
        if not all(sp.expand(z) == 0 for z in img.Q):
            parts += which_block(Al.Qgen(list(img.Q)))
        if not all(sp.expand(z) == 0 for z in img.eta):
            parts += which_block(Al.Sgen(list(img.eta)))
        preserve['[%s,%s]' % (en, n)] = parts
out['even_on_blocks'] = preserve
def bval(n): return 1 if n in Bplus else -1
ok_even = all(v == 'zero' or (len(v) >= 1 and all(bval(m) == bval(k.split(',')[1][:-1]) for m in v))
              for k, v in preserve.items())
out['even_part_preserves_B_eigenspaces'] = bool(ok_even)
out['conclusion_B_is_even_derivation'] = bool(same_block_zero and ok_even)
# (2) rescaling system
aV, aW, al, alb, kap = sp.symbols('a_V a_W a_lam a_lamb kappa')
# terms of S^(1)(xi), xi in S_+ (x) W (eq:o-S1xi): (field content, N-weight difference)
terms = {'(Gamma(xi,lambar)|_V, A*_V)': alb - aV, '(Gamma(xi,lam)|_W, A*_W)': al - aW,
         '(rho(dA_V) xi, lam*|_{S+W})': aV - al, '(rho(dA_W) xi, lam*|_{S-W*})': aW - alb,
         '(rho([A_W,A_W]) xi, lam*|_{S+W})': 2*aW - al, '(rho([A_V,A_W]) xi, lam*|_{S-W*})': aV + aW - alb,
         '(rho([A_V,A_V]) xi, lam*|_{S+W})': 2*aV - al}
lin = ['(Gamma(xi,lambar)|_V, A*_V)', '(Gamma(xi,lam)|_W, A*_W)', '(rho(dA_V) xi, lam*|_{S+W})', '(rho(dA_W) xi, lam*|_{S-W*})']
out['rescaling_system'] = {k: str(v) for k, v in terms.items()}
out['cyclic_sum_of_linear_terms'] = str(sp.simplify(sum(terms[k] for k in lin)))
sol = sp.solve([terms[k] - kap for k in lin], [aV, aW, al, alb, kap], dict=True)
out['linear_terms_equal_kappa_solutions'] = [{str(k): str(v) for k, v in s.items()} for s in sol]
out['kappa_forced_zero'] = all(sp.simplify(s.get(kap, kap)) == 0 for s in sol) if sol else False
out['weight_4.5_degree_-2_local_functionals'] = 'none: c* has weight 4, lam* lam* weight 5, lam* A* weight 11/2; no field of weight 1/2'
print(json.dumps(out, indent=1))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'outer_derivation_check.json'), 'w'), indent=1)
