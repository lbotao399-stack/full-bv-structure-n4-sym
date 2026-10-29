"""Scheme-independent spectral targets from the harmonic action of [B03c] (N. Beisert, hep-th/0307015):
  H_12 |s_1..s_n; A> = sum_{s'} c_{n,n12,n21} delta_{C_1,0} delta_{C_2,0} |s'_1..s'_n; A>          (4.27)
  c_{n,n12,n21} = (-1)^{1+n12 n21} Gamma((n12+n21)/2) Gamma(1+(n-n12-n21)/2) / Gamma(1+n/2)       (4.28)
  c_{n,0,0} = h(n/2)                                                                                (4.29)
letters (4.9): D^k F = (a+)^{k+2}(b+)^k, D^k Psi = (a+)^{k+1}(b+)^k c+, D^k Phi = (a+)^k(b+)^k (c+)^2, ...;
central charge per site C = n_a - n_b + n_c - 2; delta D = (g_YM^2 N / 8 pi^2) H (1.2), H = sum_l H_{l,l+1}.
Oscillators a_alpha, b_alphadot bosonic, c_a fermionic; sites are graded tensor factors, the state
|s;A> = A_1^+(s_1) ... A_n^+(s_n)|00> is brought to the canonical order (site 1 block first, inside a site the
order of the oscillator list TYPES) with the Koszul sign of the fermionic oscillators.
Prediction used in composite_theory.tex Sec. 4: on a planar family, D_1 = -8 u E_H on eigenstates (u the unit),
E_H an eigenvalue of H (D_1 = 4 lc^natural under the hypothesis tested there, lc^natural = -2 u E_H, Cor c10-beisert).
Exact arithmetic (fractions); spectra by sympy on the cyclic sector.  Output: ha_predictions.json."""
import json, os, itertools
from fractions import Fraction as Fr
from math import comb
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
TYPES = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'c3', 'c4']
FERM = {t: t.startswith('c') for t in TYPES}
TIDX = {t: i for i, t in enumerate(TYPES)}


def h(n):
    return sum((Fr(1, k) for k in range(1, int(n) + 1)), Fr(0))


def gamma_half(m2):
    """Gamma(m2/2) for a positive integer m2, as sympy exact"""
    return sp.gamma(sp.Rational(m2, 2))


def coef(n, n12, n21):
    if n12 == 0 and n21 == 0:
        return sp.Rational(h(Fr(n, 2)).numerator, h(Fr(n, 2)).denominator) if n % 2 == 0 else sp.harmonic(sp.Rational(n, 2))
    val = (-1) ** (1 + n12 * n21) * gamma_half(n12 + n21) * sp.gamma(1 + sp.Rational(n - n12 - n21, 2)) / sp.gamma(1 + sp.Rational(n, 2))
    return sp.nsimplify(sp.simplify(val))


def central(site):
    na = site.get('a1', 0) + site.get('a2', 0)
    nb = site.get('b1', 0) + site.get('b2', 0)
    nc = sum(site.get(t, 0) for t in ('c1', 'c2', 'c3', 'c4'))
    return na - nb + nc - 2


def canon_key(sites):
    return tuple(tuple(sorted((t, n) for t, n in s.items() if n)) for s in sites)


def H12(state):
    """state: tuple(site1 dict, site2 dict) of oscillator occupation numbers.  Returns dict key -> coefficient."""
    s1, s2 = state
    olist = []        # (type, original site) in canonical order: site 1 block, then site 2 block, type order
    for site_idx, s in ((1, s1), (2, s2)):
        for t in TYPES:
            for _ in range(s.get(t, 0)):
                olist.append((t, site_idx))
    n = len(olist)
    out = {}
    # bosons: group identical (type, site) to use binomials; fermions: explicit
    groups = {}
    for i, (t, si) in enumerate(olist):
        groups.setdefault((t, si), []).append(i)
    bos_groups = [(k, v) for k, v in groups.items() if not FERM[k[0]]]
    fer_idx = [i for i, (t, si) in enumerate(olist) if FERM[t]]
    for bos_choice in itertools.product(*[range(len(v) + 1) for k, v in bos_groups]):
        mult = 1
        for (k, v), m in zip(bos_groups, bos_choice):
            mult *= comb(len(v), m)
        for fer_choice in itertools.product((0, 1), repeat=len(fer_idx)):
            # m hopping bosons per group, fermion i hops if fer_choice = 1
            new = [dict(), dict()]
            n12 = n21 = 0
            for (k, v), m in zip(bos_groups, bos_choice):
                t, si = k
                stay = len(v) - m
                new[si - 1][t] = new[si - 1].get(t, 0) + stay
                new[2 - si][t] = new[2 - si].get(t, 0) + m
                if si == 1:
                    n12 += m
                else:
                    n21 += m
            ferm_new = []
            ok = True
            for i, hop in zip(fer_idx, fer_choice):
                t, si = olist[i]
                ns = si if not hop else 3 - si
                if hop:
                    if si == 1:
                        n12 += 1
                    else:
                        n21 += 1
                if new[ns - 1].get(t, 0) >= 1:
                    ok = False
                    break
                new[ns - 1][t] = 1
                ferm_new.append((ns, TIDX[t]))
            if not ok:
                continue
            if central(new[0]) != 0 or central(new[1]) != 0:
                continue
            # Koszul sign: sort fermions (in list order) to canonical (site, type)
            perm = sorted(range(len(ferm_new)), key=lambda j: ferm_new[j])
            sign = 1
            seen = [False] * len(perm)
            for j in range(len(perm)):
                if not seen[j]:
                    L = 0
                    x = j
                    while not seen[x]:
                        seen[x] = True
                        x = perm[x]
                        L += 1
                    if L % 2 == 0:
                        sign = -sign
            c = coef(n, n12, n21)
            key = canon_key(new)
            out[key] = out.get(key, 0) + sign * mult * c
    return {k: v for k, v in out.items() if v != 0}


def site_parity(site):
    return sum(site.get(t, 0) for t in TYPES if FERM[t]) % 2


def two_site_spectrum(basis_sites):
    """basis_sites: list of (site1, site2).  Graded cyclic sector of H_12 + H_21."""
    keys = [canon_key(b) for b in basis_sites]
    idx = {k: i for i, k in enumerate(keys)}
    N = len(keys)
    M = sp.zeros(N, N)
    for j, b in enumerate(basis_sites):
        for k, v in H12(b).items():
            if k not in idx:
                raise RuntimeError('state left the family: %s' % (k,))
            M[idx[k], j] += v
    # graded swap P: |X1 X2> -> (-1)^{|X1||X2|} |X2 X1>
    P = sp.zeros(N, N)
    for j, b in enumerate(basis_sites):
        s1, s2 = b
        sg = (-1) ** (site_parity(s1) * site_parity(s2))
        P[idx[canon_key((s2, s1))], j] = sg
    H = M + P * M * P
    # cyclic sector: image of (1+P)/2
    Pr = (sp.eye(N) + P) / 2
    cols = Pr.columnspace()
    if not cols:
        return []
    B = sp.Matrix.hstack(*cols)
    Hc = (B.T * B).inv() * B.T * H * B
    ev = Hc.eigenvals()
    return sorted([(sp.nsimplify(k), m) for k, m in ev.items()], key=lambda t: float(t[0]))


def letter(kind, k, alpha=1, adot=1, cs=()):
    s = {}
    a = 'a%d' % alpha
    b = 'b%d' % adot
    if kind == 'F':
        s[a] = k + 2; s[b] = k
    elif kind == 'Psi':
        s[a] = k + 1; s[b] = k
    elif kind == 'Phi':
        s[a] = k; s[b] = k
    for c in cs:
        s[c] = 1
    return {t: n for t, n in s.items() if n}


def family(kind1, cs1, kind2, cs2, Ntot):
    basis = []
    for k1 in range(Ntot + 1):
        basis.append((letter(kind1, k1, cs=cs1), letter(kind2, Ntot - k1, cs=cs2)))
    if (kind1, cs1) != (kind2, cs2):
        for k1 in range(Ntot + 1):
            basis.append((letter(kind2, k1, cs=cs2), letter(kind1, Ntot - k1, cs=cs1)))
    # dedupe
    seen = {}
    for b in basis:
        seen[canon_key(b)] = b
    return list(seen.values())


def closed_family(seed, Ntot_states=None):
    """close a set of two-site states under H12 and the swap"""
    todo = list(seed)
    found = {canon_key(b): b for b in seed}
    while todo:
        b = todo.pop()
        for cand in [(b[1], b[0])]:
            k = canon_key(cand)
            if k not in found:
                found[k] = cand; todo.append(cand)
        for k in H12(b):
            if k not in found:
                st = (dict(k[0]), dict(k[1]))
                found[k] = st; todo.append(st)
    return list(found.values())


def main():
    res = {'source': '[B03c] (4.9), (4.27)-(4.29), (1.2)', 'note': 'E_H = eigenvalue of H on the graded cyclic two-site sector; D_1 = -8 u E_H predicted'}
    hh = lambda n: sp.Rational(h(n).numerator, h(n).denominator)
    # 1. scalar sl(2): Z = Phi_34 = c3 c4, derivatives a1 b1
    sc = {}
    for N in range(0, 7):
        fam = closed_family(family('Phi', ('c3', 'c4'), 'Phi', ('c3', 'c4'), N))
        sc[N] = [str(e) + ('^%d' % m if m > 1 else '') for e, m in two_site_spectrum(fam)]
    res['scalar_sl2_twist_two'] = dict(letters='D^k Phi_34 = (a1+ b1+)^k c3+ c4+', spectra_by_total_derivatives=sc,
                                       expected='4 h(S), S even <= N (B03c (2.31))',
                                       check=[str(4 * hh(S)) for S in range(0, 7, 2)])
    # 2. gluon sector F_11 = (a1+)^2 with D_11 = a1+ b1+
    gl = {}
    for N in range(0, 6):
        fam = closed_family(family('F', (), 'F', (), N))
        gl[N] = [str(e) + ('^%d' % m if m > 1 else '') for e, m in two_site_spectrum(fam)]
    res['gluon_F11'] = dict(letters='D^k F_11 = (a1+)^{k+2}(b1+)^k', spectra_by_total_derivatives=gl,
                            compare=[str(4 * hh(S + 2)) for S in range(0, 6, 2)])
    # 3. fermion sector Psi_1 with c4, closed under H
    fe = {}
    for N in range(0, 6):
        fam = closed_family(family('Psi', ('c4',), 'Psi', ('c4',), N))
        fe[N] = dict(dim=len(fam), spectrum=[str(e) + ('^%d' % m if m > 1 else '') for e, m in two_site_spectrum(fam)])
    res['fermion_Psi1_c4'] = dict(letters='D^k Psi = (a1+)^{k+1}(b1+)^k c4+ (one R component)', by_N=fe)
    fe2 = {}
    for N in range(0, 5):
        fam = closed_family(family('Psi', ('c3',), 'Psi', ('c4',), N))
        fe2[N] = dict(dim=len(fam), spectrum=[str(e) + ('^%d' % m if m > 1 else '') for e, m in two_site_spectrum(fam)])
    res['fermion_Psi1_c3_c4'] = dict(letters='D^k Psi with c3+ and with c4+ (two R components), closed under H', by_N=fe2)
    # 4. mixed (F_11, Psi_1 c4)
    mx = {}
    for N in range(0, 5):
        fam = closed_family(family('F', (), 'Psi', ('c4',), N))
        mx[N] = dict(dim=len(fam), spectrum=[str(e) + ('^%d' % m if m > 1 else '') for e, m in two_site_spectrum(fam)])
    res['mixed_F11_Psi1'] = dict(letters='D^k F_11 and D^k Psi_1 c4 (closed under H)', by_N=mx)
    # 5. Konishi / 20': two scalars without derivatives, all six Phi_ab
    phis = [('c%d' % i, 'c%d' % j) for i in range(1, 5) for j in range(i + 1, 5)]
    basis = [(letter('Phi', 0, cs=p), letter('Phi', 0, cs=q)) for p in phis for q in phis]
    res['scalars_no_derivatives'] = [str(e) + ('^%d' % m if m > 1 else '') for e, m in two_site_spectrum(basis)]
    # 6. su(1|1): Z = c3 c4, psi = a1 c4: two sites
    Z = letter('Phi', 0, cs=('c3', 'c4'))
    psi = {'a1': 1, 'c4': 1}
    res['su11_two_site'] = dict(ZZ=[str(e) for e, m in two_site_spectrum([(Z, Z)])],
                                Zpsi=[str(e) + ('^%d' % m if m > 1 else '') for e, m in two_site_spectrum(closed_family([(Z, psi)]))])
    json.dump(res, open(os.path.join(HERE, 'ha_predictions.json'), 'w'), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == '__main__':
    main()
