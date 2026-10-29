"""Audit script a2: sign identities of the odd bracket of ch01 def:c1-bracket,
{f,g} = Delta(fg) - Delta(f)g - (-1)^{|f|} f Delta(g),
in a finite-dimensional bigraded model with bosonic and fermionic fields, antifields, ghost.

Checks the identities used by the revised design (prop:o-laws (a), eq:o-twistlaws, thm:o-twist (ii))
and by ch08 Prop c8-laws (a):
  [s, ad_phi] = - ad_{s phi}       for s = ad_S, S even, phi odd   (claimed)
  [s, ad_G]   = + ad_{s G} ?       for G even
  [N_gh, ad_f] = ad_{N_gh f} + ad_f
  ad_f ad_f = (1/2) ad_{{f,f}}     for f even  ((Q_0^xi)^2 = 1/2 {{S^xi,S^xi},-})
  graded symmetry {g,f} = (-1)^{|f||g|}{f,g}
"""
import json, random, sys, os
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grass import Alg

names = ['x1', 'x2', 'X1s', 'X2s', 't1', 't2', 'T1s', 'T2s', 'c', 'cs']
par = [0, 0, 1, 1, 1, 1, 0, 0, 1, 0]
deg = [0, 0, -1, -1, 0, 0, -1, -1, 1, -2]
pairs = [(0, 2), (1, 3), (4, 6), (5, 7), (8, 9)]
A = Alg(names, par, deg, pairs)
rng = random.Random(20261029)

res = {}


def rnd(p, d, n=4, md=3):
    while True:
        e = A.random_hom(rng, p, d, n, md)
        if e:
            return e


# 1. Delta^2 = 0
ok = True
for _ in range(30):
    p = rng.randint(0, 1); d = rng.choice([-2, -1, 0, 1])
    try:
        f = rnd(p, d)
    except Exception:
        continue
    if not A.is_zero(A.Delta(A.Delta(f))):
        ok = False
res['Delta_squared_zero'] = ok

# 2. graded symmetry
cnt = {'(-1)^{|f||g|}': 0, 'other': 0}
for _ in range(60):
    pf, pg = rng.randint(0, 1), rng.randint(0, 1)
    f = rnd(pf, rng.choice([-1, 0, 1])); g = rnd(pg, rng.choice([-1, 0, 1]))
    a = A.br(f, g); b = A.br(g, f)
    if A.is_zero(A.add(b, a, -((-1) ** (pf * pg)))):
        cnt['(-1)^{|f||g|}'] += 1
    else:
        cnt['other'] += 1
res['symmetry_{g,f}=(-1)^{|f||g|}{f,g}'] = cnt


# 3. Jacobi: fit {f,{g,h}} = alpha {{f,g},h} + beta {g,{f,h}} per parity pattern
def fit_jacobi(pf, pg, ph, trials=12):
    found = None
    for _ in range(trials):
        f = rnd(pf, rng.choice([-1, 0])); g = rnd(pg, rng.choice([-1, 0])); h = rnd(ph, rng.choice([-1, 0, 1]))
        L = A.br(f, A.br(g, h)); T1 = A.br(A.br(f, g), h); T2 = A.br(g, A.br(f, h))
        sols = []
        for al in (1, -1):
            for be in (1, -1):
                if A.is_zero(A.add(A.add(L, T1, -al), T2, -be)):
                    sols.append((al, be))
        if not (A.is_zero(T1) and A.is_zero(T2)):
            if found is None:
                found = set(sols)
            else:
                found &= set(sols)
    return sorted(found) if found is not None else None


jac = {}
for pf in (0, 1):
    for pg in (0, 1):
        for ph in (0, 1):
            jac[f'pf={pf},pg={pg},ph={ph}'] = fit_jacobi(pf, pg, ph)
res['jacobi_fits_(alpha,beta)'] = jac

# 4. [ad_S, ad_phi] = - ad_{S phi}?  S even (degree 0), phi odd (degree -1)
ok_minus, ok_plus, n = 0, 0, 0
for _ in range(25):
    S = rnd(0, 0); phi = rnd(1, -1); O = rnd(rng.randint(0, 1), rng.choice([-1, 0, 1]))
    lhs = A.add(A.br(S, A.br(phi, O)), A.br(phi, A.br(S, O)), -1)
    sphi = A.br(S, phi)
    rhs = A.br(sphi, O)
    if A.is_zero(lhs) and A.is_zero(rhs):
        continue
    n += 1
    if A.is_zero(A.add(lhs, rhs, 1)):
        ok_minus += 1
    if A.is_zero(A.add(lhs, rhs, -1)):
        ok_plus += 1
res['[ad_S,ad_phi] (S even, phi odd): s{phi,O}-{phi,sO}'] = {'= -{s phi,O}': ok_minus, '= +{s phi,O}': ok_plus, 'samples': n}

# 5. G even: s{G,O} - {G,sO} vs {sG,O}
ok_minus, ok_plus, n = 0, 0, 0
for _ in range(25):
    S = rnd(0, 0); G = rnd(0, 0); O = rnd(rng.randint(0, 1), rng.choice([-1, 0, 1]))
    pO = A.parity(O)
    # graded commutator: ad_S odd, ad_G odd => [ad_S, ad_G] = ad_S ad_G + ad_G ad_S
    lhs = A.add(A.br(S, A.br(G, O)), A.br(G, A.br(S, O)), 1)
    rhs = A.br(A.br(S, G), O)
    if A.is_zero(lhs) and A.is_zero(rhs):
        continue
    n += 1
    if A.is_zero(A.add(lhs, rhs, 1)):
        ok_minus += 1
    if A.is_zero(A.add(lhs, rhs, -1)):
        ok_plus += 1
res['G even: s{G,O}+{G,sO} (graded commutator of two odd maps)'] = {'= -{sG,O}': ok_minus, '= +{sG,O}': ok_plus, 'samples': n}


# 6. N_gh
def Ngh(a):
    r = {}
    for (e, o), c in a.items():
        d = sum(A.deg[A.even[k]] * x for k, x in enumerate(e)) + sum(A.deg[i] for i in o)
        if d:
            r[(e, o)] = c * d
    return r


ok = 0; n = 0
for _ in range(25):
    f = rnd(rng.randint(0, 1), rng.choice([-1, 0])); g = rnd(rng.randint(0, 1), rng.choice([-1, 0, 1]))
    lhs = Ngh(A.br(f, g)); rhs = A.add(A.add(A.br(Ngh(f), g), A.br(f, Ngh(g))), A.br(f, g))
    n += 1
    ok += A.is_zero(A.add(lhs, rhs, -1))
res['N_gh{f,g} = {N_gh f,g}+{f,N_gh g}+{f,g}'] = f'{ok}/{n}'

# 7. N = ad_psi, psi odd (like S^(1)(R_xi): degree -1, odd parity): [N, ad_f] = ad_{N f}?
ok = 0; n = 0
for _ in range(25):
    psi = rnd(1, -1); f = rnd(0, rng.choice([-1, 0])); O = rnd(rng.randint(0, 1), rng.choice([-1, 0, 1]))
    # N even (odd psi + odd bracket), ad_f odd for f even
    lhs = A.add(A.br(psi, A.br(f, O)), A.br(f, A.br(psi, O)), -1)
    rhs = A.br(A.br(psi, f), O)
    if A.is_zero(lhs) and A.is_zero(rhs):
        continue
    n += 1
    ok += A.is_zero(A.add(lhs, rhs, -1))
res['[ad_psi, ad_f] = ad_{ {psi,f} } (psi odd, f even)'] = f'{ok}/{n}'

# 8. ad_f ad_f = 1/2 ad_{{f,f}} for f even
ok = 0; n = 0; okm = 0
for _ in range(25):
    f = rnd(0, rng.choice([-1, 0])); O = rnd(rng.randint(0, 1), rng.choice([-1, 0, 1]))
    lhs = A.br(f, A.br(f, O)); rhs = A.scale(A.br(A.br(f, f), O), Fraction(1, 2))
    if A.is_zero(lhs) and A.is_zero(rhs):
        continue
    n += 1
    ok += A.is_zero(A.add(lhs, rhs, -1))
    okm += A.is_zero(A.add(lhs, rhs, 1))
res['ad_f ad_f = +1/2 ad_{f,f} (f even)'] = f'{ok}/{n}'
res['ad_f ad_f = -1/2 ad_{f,f} (f even)'] = f'{okm}/{n}'

# 9. symmetric bracket of S^(1)(R)-like odd psi with even S^xi-like: {psi, Sx} = (-1)^{1*0}{Sx, psi}
ok = 0; n = 0
for _ in range(20):
    psi = rnd(1, -1); Sx = rnd(0, rng.choice([0, -1]))
    n += 1
    ok += A.is_zero(A.add(A.br(psi, Sx), A.br(Sx, psi), -1))
res['{psi,Sxi} = {Sxi,psi} (psi odd, Sxi even): N_R S^xi = Q_0^xi(S^(1)(R_xi))'] = f'{ok}/{n}'

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a2_bv_signs.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
