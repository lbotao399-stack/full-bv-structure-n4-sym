"""Closed forms, comparator's own implementations:
 (a) k = 0 (Lemma lem:tl-k0): point-split self-contraction 2 (Xa.Xb) 16pi^2 S_eps(u - s),
     16pi^2 S_eps(v) = sum_m eps^-m (-1)^{m+1} (|v|^2/4)^{m-1}/m!   (derived from int_eps dt t^-2 e^{-|v|^2/4t});
     compared with engine B's '(none)' leg sets at W0+1, W0+2 (AV|AV, AW|AW) and W0.
 (b) k = 1 quartic bubble (eq:tl-bubbleclosed), derived independently (Gaussian loop integral, total-time cutoff;
     log coefficient shape independent):  [log] = -N_V sum_{a,b,j} (i s.R)^a (i u.R)^b |u-s|^{2j} (|R|^2/4)^j
     (a+j)!(b+j)!/(a! b! (j!)^2 (a+b+2j+1)!),  iR = d on the product of the outputs (R = sum of the output
     momenta, as stated after Prop prop:tl-genfun);
     compared with engine A's Q-block evaluated on demand (AW|AW, AV|AV, AV|AW -> two outputs), derivatives n <= 3."""
import os, sys, json, math, itertools
from fractions import Fraction
import numpy as np
import sympy as sp
import cmp_lib as L
import cmp_deriv_parse as CP

sys.path.insert(0, os.path.join(L.EA, 'src.nosync'))
import ea_core as E, ea_kern as K, ea_tables as TB

BD = json.load(open(os.path.join(L.EB, 'kernels_deriv.json')))
BW = json.load(open(os.path.join(L.EB, 'kernels_W0.json')))
res = {}

# ---------------------------------------------------------------- (a) k = 0
s_ = sp.symbols('s0:4'); u_ = sp.symbols('u0:4')
v2 = sum((u_[m] - s_[m]) ** 2 for m in range(4))
ok_a = True
rows = []
for pair, comps in (('AV|AV', ['A_0', 'A_1', 'A_2', 'A_3']), ('AW|AW', ['A_I%d' % i for i in range(1, 7)])):
    for nder in (0, 1, 2):
        if nder == 0:
            d = BW[pair]['leg_sets']['(none)']['total'].get('', {})
        else:
            d = BD[f'{pair} @ W0+{nder}']['leg_sets'].get('(none)', {}).get('total', {}).get('', {})
        W = 2 + nder
        for m in range(1, 4):
            key = f'eps^-{m}'
            s = d.get(key, '')
            parsed = CP.parse_Bd(s, 0) if s else {}
            # B polynomial: sum over Xa[c] Xb[c'] of poly in s,u
            got = {}
            for k_, poly in parsed.items():
                ex = 0
                for mon, c in poly.items():
                    t = sp.nsimplify(complex(c).real)
                    for i in range(4):
                        t *= s_[i] ** mon[i] * u_[i] ** mon[4 + i]
                    ex += t
                got[(k_[0], k_[1])] = sp.expand(ex)
            # prediction: total order exactly nder in s,u; only m with 2(m-1) = nder contributes
            pred_poly = sp.expand(2 * (-1) ** (m + 1) * (v2 / 4) ** (m - 1) / sp.factorial(m)) if 2 * (m - 1) == nder else 0
            okm = True
            for ca, cb in itertools.product(comps, comps):
                exp = pred_poly if ca == cb else 0
                if sp.expand(got.get((ca, cb), 0) - exp) != 0:
                    okm = False
            ok_a &= okm
            rows.append(dict(pair=pair, W=W, kind=key, match=bool(okm), predicted=str(pred_poly)))
res['k0_lemma_vs_engineB_none_legsets'] = dict(rows=rows, pass_=bool(ok_a))
print('k=0 lemma vs engine B (none) leg sets:', ok_a)

# ---------------------------------------------------------------- (b) quartic bubble closed form vs engine A Q block
rng = np.random.default_rng(7)


def closed_form(N_V, s, u, R, nmax):
    """{(adeg,bdeg): value} of sum s^al u^be/(al! be!) K(d^al X_a, d^be X_b) restricted to bidegree (a,b) at the
    given s, u vectors (directional), from eq:tl-bubbleclosed"""
    sR = sum(si * Ri for si, Ri in zip(s, R)); uR = sum(ui * Ri for ui, Ri in zip(u, R))
    v2n = sum((ui - si) ** 2 for si, ui in zip(s, u)); R2 = sum(Ri * Ri for Ri in R)
    out = {}
    # expand |u-s|^{2j} in bidegree: |u-s|^2 = s.s - 2 s.u + u.u
    ss = sum(si * si for si in s); su = sum(si * ui for si, ui in zip(s, u)); uu = sum(ui * ui for ui in u)
    for a in range(nmax + 1):
        for b in range(nmax + 1 - a):
            for j in range(0, (nmax - a - b) // 2 + 1):
                coef = Fraction(math.factorial(a + j) * math.factorial(b + j),
                                math.factorial(a) * math.factorial(b) * math.factorial(j) ** 2 * math.factorial(a + b + 2 * j + 1))
                base = (1j * sR) ** a * (1j * uR) ** b * (R2 / 4) ** j
                # |u-s|^{2j} = sum over (p,q,r) p+q+r=j of multinomial ss^p (-2su)^q uu^r: bidegree (2p+q, q+2r)
                for p in range(j + 1):
                    for q in range(j + 1 - p):
                        r = j - p - q
                        mult = math.factorial(j) // (math.factorial(p) * math.factorial(q) * math.factorial(r))
                        key = (a + 2 * p + q, b + q + 2 * r)
                        out[key] = out.get(key, 0) - N_V * float(coef) * base * mult * ss ** p * (-2 * su) ** q * uu ** r
    return out


def dirpoly(c, deg):
    o = {(0, 0, 0, 0): 1}
    for _ in range(deg):
        new = {}
        for e, v in o.items():
            for mu in range(4):
                if c[mu]:
                    ee = list(e); ee[mu] += 1; ee = tuple(ee)
                    new[ee] = new.get(ee, 0) + v * 1j * c[mu]
        o = new
    return o


rows = []
ok_b = True
for (ta, tb, word) in (('AW', 'AW', ('AW', 'AW')), ('AV', 'AV', ('AV', 'AV')), ('AV', 'AW', ('AV', 'AW')),
                       ('AV', 'AV', ('AW', 'AW')), ('AW', 'AW', ('AV', 'AV'))):
    va, vb = TB.letter_vecs(ta), TB.letter_vecs(tb)
    for trial in range(3):
        moms = [rng.integers(-3, 4, 4) for _ in range(2)]
        pols = [rng.integers(-2, 3, L.NA[t]) + 1j * rng.integers(-2, 3, L.NA[t]) for t in word]
        pv = [sum(pols[j][i] * TB.stripped_pol(word[j], i) for i in range(L.NA[word[j]])) for j in range(2)]
        c = rng.integers(-2, 3, 4); d = rng.integers(-2, 3, 4)
        while not c.any():
            c = rng.integers(-2, 3, 4)
        while not d.any():
            d = rng.integers(-2, 3, 4)
        blocks = [('qua', (pv[0], moms[0]), (pv[1], moms[1]))]
        # N_V: the W0 value (no derivatives)
        r0 = K.kernel(blocks, ta, tb, kinds=('log',))
        for nd in range(0, 4):
            for adeg in range(nd + 1):
                bdeg = nd - adeg
                rA = K.kernel(blocks, ta, tb, kinds=('log',), fa=dirpoly(c, adeg) if adeg else None,
                              fb=dirpoly(d, bdeg) if bdeg else None)
                for ia, za_ in enumerate(va[:3]):
                    for ib, zb_ in enumerate(vb[:3]):
                        NV = K.contract(r0['log'], za_, zb_) if r0 and r0.get('log') is not None else (0, 0)
                        NV = -(float(NV[0]) + 1j * float(NV[1]))       # [log] at n = 0 is -N_V
                        vA = K.contract(rA['log'], za_, zb_) if rA and rA.get('log') is not None else (0, 0)
                        vA = float(vA[0]) + 1j * float(vA[1])
                        Rsum = moms[0] + moms[1]             # iR = d acting on the outputs (theory, after eq:tl-genfun)
                        cf = closed_form(NV, c, d, Rsum, nd)
                        vC = cf.get((adeg, bdeg), 0) * math.factorial(adeg) * math.factorial(bdeg)
                        cf2 = closed_form(NV, c, d, -Rsum, nd)
                        vC2 = cf2.get((adeg, bdeg), 0) * math.factorial(adeg) * math.factorial(bdeg)
                        okk = abs(vA - vC) <= 1e-9 * max(1, abs(vA))
                        ok_b &= okk
                        rows.append(dict(pair=f'{ta}|{tb}', word=' '.join(word), n=nd, a=adeg, b=bdeg, engineA=str(vA),
                                         closed_form=str(vC), match=bool(okk),
                                         match_if_R_sign_reversed=bool(abs(vA - vC2) <= 1e-9 * max(1, abs(vA)))))
res['bubble_closed_form_vs_engineA_Qblock'] = dict(n_checks=len(rows), n_nonzero=sum(r['engineA'] not in ('0j', '0') for r in rows),
                                                   pass_=bool(ok_b), rows=rows[:60])
print('bubble closed form vs engine A Q block:', ok_b, len(rows), 'checks;',
      'R sign reversed would match:', all(r['match_if_R_sign_reversed'] for r in rows))
bad = [r for r in rows if not r['match']]
for r in bad[:10]:
    print('  MISMATCH', r)
json.dump(res, open(os.path.join(L.QA, 'out.nosync', 'closed_forms_check.json'), 'w'), indent=1)
