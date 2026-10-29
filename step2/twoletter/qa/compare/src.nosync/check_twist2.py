"""Task 2, twist two (predictions P5, P6b):
 (i)  low-order entries (a+b <= 2) of engine B's exported sl(2) data (two_site.json: sl2_null density, AZ gauge-letter
      kernels) re-derived from the derivative tables kernels_W0/kernels_deriv (equal to engine A, cmp_W0/cmp_deriv);
 (ii) comparator's own assembly of lc^natural on O_k = tr(D^k Z D^{N-k} Z), D = n.d + [n.A, -], to O(A_n), two-field
      (ZZ) part, from the sl2_null density and the AZ/ZA kernels; spectrum vs -8 h(S) u (P5) for N <= 8;
 (iii) eigenvalues of engine A's exported twist-two matrices (kernels.json) recomputed exactly."""
import os, json, itertools, math
from fractions import Fraction
import sympy as sp
import cmp_lib as L
import cmp_deriv_parse as CP

TS = json.load(open(os.path.join(L.EB, 'two_site.json')))
BW = json.load(open(os.path.join(L.EB, 'kernels_W0.json')))
BD = json.load(open(os.path.join(L.EB, 'kernels_deriv.json')))
KA = json.load(open(os.path.join(L.EA, 'out.nosync', 'kernels.json')))
P = json.load(open(os.path.join(L.TL, 'theory', 'predictions.json')))
I = sp.I
res = {}

# ------------------------------------------------------------------ (i) re-derive low-order sl(2) data (engine-B basis)
n_vec = [1, I, 0, 0]                    # n in engine-B directions (A_n = A_0 + i A_1)
nbar_half = [sp.Rational(1, 2), -I / 2, 0, 0]   # p = t nbar/2  ->  n.p = t
ZETA = {'A_I1': 1, 'A_I2': I}           # Z = A_I1 + i A_I2
ZPROJ = {'A_I1': sp.Rational(1, 2), 'A_I2': -I / 2}   # coefficient of Z in A_I
AN = {'A_0': 1, 'A_1': I}


def evaluate(parsed, outs_n, zeta_a, zeta_b, a, b, nder):
    """sum over terms: letters contracted with zeta, outputs projected on Z, momenta p_j = t_j nbar/2, s = u = n,
    s-degree a, u-degree b; returns polynomial in t_j (dict exps -> coeff)"""
    t = sp.symbols('t0:%d' % outs_n)
    tot = 0
    for key, poly in parsed.items():
        za = zeta_a.get(key[0], 0); zb = zeta_b.get(key[1], 0)
        if za == 0 or zb == 0:
            continue
        f = za * zb
        for j in range(outs_n):
            f *= ZPROJ.get(key[2 + j], 0)
        if f == 0:
            continue
        for mon, c in poly.items():
            sd = sum(mon[4 * outs_n:4 * outs_n + 4]) if nder else 0
            ud = sum(mon[4 * outs_n + 4:4 * outs_n + 8]) if nder else 0
            if sd != a or ud != b:
                continue
            v = sp.nsimplify(complex(c).real) + I * sp.nsimplify(complex(c).imag)
            for j in range(outs_n):
                for m in range(4):
                    e = mon[4 * j + m]
                    if e:
                        v *= (t[j] * nbar_half[m]) ** e
            if nder:
                for m in range(4):
                    if mon[4 * outs_n + m]:
                        v *= n_vec[m] ** mon[4 * outs_n + m]
                    if mon[4 * outs_n + 4 + m]:
                        v *= n_vec[m] ** mon[4 * outs_n + 4 + m]
            tot += f * v
    return sp.Poly(sp.expand(tot * math.factorial(a) * math.factorial(b)), *t) if tot != 0 else None, t


def table(pair, nder):
    if nder == 0:
        return BW[pair]
    return BD.get(f'{pair} @ W0+{nder}')


low = {}
ok_low = True
# Z-Z density: pair AW|AW, leg set 'AW AW', B word '1,0' (A order: Y1 = label 0, Y2 = label 1), r = 2
for nder in (0, 1, 2):
    tb = table('AW|AW', nder)
    s = tb['leg_sets']['AW AW']['total'].get('1,0', {}).get('log', '')
    parsed = CP.parse_Bd(s, 2) if nder else {k: {tuple(list(m) + [0] * 8): c for m, c in v.items()}
                                               for k, v in L.parse_B(s, ['AW', 'AW'], 2).items()}
    for a in range(nder + 1):
        b = nder - a
        pol, t = evaluate(parsed, 2, ZETA, ZETA, a, b, 1)
        for c in range(nder + 1):
            d = nder - c
            coeff = pol.coeff_monomial(t[0] ** c * t[1] ** d) if pol is not None else 0
            # value of K on outputs (n.d)^c Z (n.d)^d Z: monomial t0^c t1^d = (n.p0)^c (n.p1)^d = (-i)^{c+d} (n.d)^c (n.d)^d
            val = sp.nsimplify(coeff * (-I) ** (c + d))
            dens = sp.nsimplify(val / 2)
            exp = sp.nsimplify(TS['sl2_null']['entries'].get(f'{c},{d}<-{a},{b}', '0'))
            okk = sp.simplify(dens - exp) == 0
            ok_low &= okk
            low[f'ZZ {c},{d}<-{a},{b}'] = dict(from_tables=str(dens), two_site_json=str(exp), match=bool(okk))
# A_n - Z gauge-letter kernel: pair AV|AW, leg set 'AW', r = 1: planar coefficient = -(value/2)
for nder in (0, 1, 2):
    tb = table('AV|AW', nder)
    s = tb['leg_sets']['AW']['total'].get('0', {}).get('log', '')
    parsed = CP.parse_Bd(s, 1) if nder else {k: {tuple(list(m) + [0] * 8): c for m, c in v.items()}
                                               for k, v in L.parse_B(s, ['AW'], 1).items()}
    for r in range(nder + 1):
        s_ = nder - r
        pol, t = evaluate(parsed, 1, AN, ZETA, r, s_, 1)
        coeff = pol.coeff_monomial(t[0] ** (nder + 1)) if pol is not None else 0
        val = sp.nsimplify(coeff * (-I) ** (nder + 1))
        plan = sp.nsimplify(-val / 2)
        exp = sp.nsimplify(TS['sl2_null']['gauge_letter_kernels'].get(f'AZ,{r},{s_}', '0').replace('+0*I', ''))
        okk = sp.simplify(plan - exp) == 0
        ok_low &= okk
        low[f'AZ {r},{s_}'] = dict(from_tables=str(plan), two_site_json=str(exp), match=bool(okk))
res['low_order_two_site_vs_tables'] = dict(entries=low, pass_=bool(ok_low))
print('low-order two_site entries vs agreed tables:', ok_low)
for k, v in low.items():
    if not v['match']:
        print('  MISMATCH', k, v)

# ------------------------------------------------------------------ (ii) own assembly of the twist-two spectrum
h = {}
for k, v in TS['sl2_null']['entries'].items():
    cd, ab = k.split('<-')
    h[(tuple(map(int, ab.split(','))), tuple(map(int, cd.split(','))))] = sp.nsimplify(v)
AZ = {}; ZA = {}
for k, v in TS['sl2_null']['gauge_letter_kernels'].items():
    typ, r, s_ = k.split(',')
    (AZ if typ == 'AZ' else ZA)[(int(r), int(s_))] = sp.nsimplify(v.replace('+0*I', ''))


def lc_nat_ZZ(k, N):
    """ZZ part of lc^natural(O_k), O_k = tr(D^k Z D^{N-k} Z), as dict {(a,b): coeff} (cyclic words tr(Z_a Z_b))"""
    out = {}
    def add(a, b, c):
        key = (min(a, b), max(a, b))
        out[key] = out.get(key, 0) + c
    # two-letter ZZ term tr(Z_k Z_{N-k}): links (1,2) and (2,1)
    for (x, y) in ((k, N - k), (N - k, k)):
        for (ab, cd), v in h.items():
            if ab == (x, y):
                add(cd[0], cd[1], v)
    add(k, N - k, -4)                       # - u {F, O}: n_F = 2 per scalar letter
    # O(A_n) terms: D^m Z = d^m Z + sum_{j<m} d^j [A, d^{m-1-j} Z] + ...;  d^j [A, Z_q] = sum_i C(j,i) [A_i, Z_{j-i+q}]
    def gauge_terms(m):
        terms = []
        for j in range(m):
            q = m - 1 - j
            for i in range(j + 1):
                terms.append((i, j - i + q, math.comb(j, i)))
        return terms
    # first factor covariant: tr([A_i, Z_m] Z_l) = tr(A_i Z_m Z_l) - tr(Z_m A_i Z_l)
    for (i, m, c) in gauge_terms(k):
        l = N - k
        add(i + m + 1, l, c * AZ.get((i, m), 0))      # tr(A_i Z_m Z_l): pair (A_i, Z_m) -> Z_{i+m+1}
        add(m, i + l + 1, c * ZA.get((i, l), 0))      # cyclic pair (Z_l, A_i) -> Z_{i+l+1}: tr(Z_m Z_{i+l+1})
        add(i + m + 1, l, -c * ZA.get((i, m), 0))     # -tr(Z_m A_i Z_l): pair (Z_m, A_i)
        add(m, i + l + 1, -c * AZ.get((i, l), 0))     # -tr(Z_m A_i Z_l): pair (A_i, Z_l)
    # second factor covariant: tr(Z_k [A_i, Z_m]) = tr(Z_k A_i Z_m) - tr(Z_k Z_m A_i)
    for (i, m, c) in gauge_terms(N - k):
        kk = k
        add(kk, i + m + 1, c * AZ.get((i, m), 0))     # tr(Z_k A_i Z_m): pair (A_i, Z_m)
        add(i + kk + 1, m, c * ZA.get((i, kk), 0))    # pair (Z_k, A_i)
        add(kk, i + m + 1, -c * ZA.get((i, m), 0))    # -tr(Z_k Z_m A_i): pair (Z_m, A_i)
        add(i + kk + 1, m, -c * AZ.get((i, kk), 0))   # cyclic pair (A_i, Z_k): tr(Z_{i+k+1} Z_m)
    return out


def hsum(S):
    return sum(sp.Rational(1, j) for j in range(1, S + 1))


tw = {}
ok_tw = True
for N in range(0, 9):
    basis = [(a, N - a) for a in range(0, N // 2 + 1)]
    bi = {b: i for i, b in enumerate(basis)}
    M = sp.zeros(len(basis), len(basis))
    for j, (a, b) in enumerate(basis):
        img = lc_nat_ZZ(a, N)
        for key, c in img.items():
            M[bi[key], j] += c
    ev = M.eigenvals()
    evs = sorted([sp.nsimplify(e) for e, mlt in ev.items() for _ in range(mlt)])
    target = sorted([-8 * hsum(S) for S in range(0, N + 1, 2)])
    okN = bool(len(evs) == len(target) and all(sp.simplify(x - y) == 0 for x, y in zip(evs, target)))
    ok_tw &= okN
    tw[str(N)] = dict(matrix=[[str(M[i, j]) for j in range(M.cols)] for i in range(M.rows)],
                      eigenvalues=[str(e) for e in evs], target=[str(t) for t in target], match=okN)
    print('N', N, 'eig', [str(e) for e in evs], 'target', [str(t) for t in target], okN)
# Z-Z only (no gauge-letter kernels): eigenvalues, cf. engine B note 2u/(j+1)-type values and engine A G8 diagnostic
zz_only = {}
for N in (2, 4):
    basis = [(a, N - a) for a in range(0, N // 2 + 1)]
    bi = {b: i for i, b in enumerate(basis)}
    M = sp.zeros(len(basis), len(basis))
    for j, (a, b) in enumerate(basis):
        for (x, y) in ((a, b), (b, a)):
            for (ab, cd), v in h.items():
                if ab == (x, y):
                    M[bi[(min(cd), max(cd))], j] += v
        M[j, j] -= 4
    zz_only[str(N)] = sorted(str(sp.nsimplify(e)) for e in M.eigenvals())
res['twist_two_own_assembly'] = dict(per_N=tw, pass_=bool(ok_tw), ZZ_only_without_gauge_letters=zz_only)
print('ZZ only', zz_only)

# ------------------------------------------------------------------ (iii) engine A's exported matrices
ea = {}
ok_a = True
for N, d in KA['twist_two'].items():
    M = sp.Matrix([[sp.nsimplify(x) for x in row] for row in d['matrix_units']])
    evs = sorted(sp.nsimplify(e) for e, mlt in M.eigenvals().items() for _ in range(mlt))
    target = sorted([-8 * hsum(S) for S in range(0, int(N) + 1, 2)])
    okN = len(evs) == len(target) and all(sp.simplify(x - y) == 0 for x, y in zip(evs, target))
    ok_a &= okN
    ea[N] = dict(eigenvalues=[str(e) for e in evs], match=bool(okN))
res['engine_A_matrices_eigenvalues'] = dict(per_N=ea, pass_=bool(ok_a))
print('engine A matrices', ok_a)
# engine B reported eigenvalues
eb = {N: d['eig_lc_natural_units'] for N, d in TS['twist_two'].items()}
res['engine_B_reported'] = eb
res['predictions_P5_twist_two'] = {S: v['lc_natural_units'] for S, v in P['P5_spectrum_Beisert']['twist_two']['values'].items()}
json.dump(res, open(os.path.join(L.QA, 'out.nosync', 'twist2_check.json'), 'w'), indent=1)
