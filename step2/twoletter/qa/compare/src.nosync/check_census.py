"""Task 3: census and output bound of theory/census.json against what the engines found.
For every nonzero kernel (engine A W0 tables, engine B W0 tables, engine B derivative tables): the leg set must be a
predicted log leg set (log) or a predicted eps^{-m} leg set with m <= eps_max_m, and every output-momentum monomial
must have degree W - w_out - 2m (m = 0 for log eps; letter derivatives counted in W).  Predicted leg sets with
identically zero kernels are listed (expected zeros vs others)."""
import os, json
from fractions import Fraction
import cmp_lib as L

WEYL = {'c': Fraction(0), 'AV': Fraction(1), 'AW': Fraction(1), 'lamP': Fraction(3, 2), 'lamM': Fraction(3, 2),
        'AsV': Fraction(3), 'AsW': Fraction(3), 'lamsP': Fraction(5, 2), 'lamsM': Fraction(5, 2), 'cs': Fraction(4)}
GH = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 0, 'lamM': 0, 'AsV': -1, 'AsW': -1, 'lamsP': -1, 'lamsM': -1, 'cs': -2}
cen = json.load(open(os.path.join(L.TL, 'theory', 'census.json')))['output_leg_sets_by_pair_weight']
A = json.load(open(os.path.join(L.EA, 'out.nosync', 'kernels_W0_raw.json')))
B = json.load(open(os.path.join(L.EB, 'kernels_W0.json')))
Bd = json.load(open(os.path.join(L.EB, 'kernels_deriv.json')))


def canon(ls):
    return ' '.join(sorted(ls.split())) if ls not in ('', '(none)') else '(none)'


def predicted(pair, W):
    e = cen[pair].get(str(W) if Fraction(W).denominator > 1 else str(int(Fraction(W))))
    if e is None:
        return None
    return dict(log=set(canon(x) for x in e['log']), eps={canon(k): v for k, v in e['eps_max_m'].items()},
                nder=e['letter_derivatives'])


findings, nonzero, checked_monos = [], [], 0
seen = {}


def check(src, pair, W, ls, kind, mons, nvar_out, extra_deg=0):
    """mons: iterable of monomials (tuples; first 4*n entries = output momenta)"""
    global checked_monos
    outs = ls.split() if ls not in ('', '(none)') else []
    pr = predicted(pair, W)
    m = 0 if kind == 'log' else int(str(kind).replace('eps^-', ''))
    key = (src, pair, str(W), canon(ls), str(kind))
    seen[key] = True
    if pr is None:
        findings.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind=str(kind), issue='pair weight not in census'))
        return
    if kind == 'log':
        if canon(ls) not in pr['log']:
            findings.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind='log', issue='nonzero log kernel outside predicted log leg sets'))
    else:
        mm = pr['eps'].get(canon(ls))
        if mm is None or m > mm:
            findings.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind=str(kind), issue='nonzero eps^-m outside predicted sets'))
    wout = sum((WEYL[t] for t in outs), Fraction(0))
    target = Fraction(W) - wout - 2 * m
    gh = GH[pair.split('|')[0]] + GH[pair.split('|')[1]]
    if sum(GH[t] for t in outs) != gh:
        findings.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind=str(kind), issue='ghost number not conserved'))
    if len(outs) > Fraction(W) + gh:
        findings.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind=str(kind), issue='n_out > W + d_ab'))
    for mon in mons:
        checked_monos += 1
        if sum(mon[:nvar_out]) != target:
            findings.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind=str(kind),
                                 issue=f'momentum degree {sum(mon[:nvar_out])} != W - w_out - 2m = {target}'))
            break
    nonzero.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind=str(kind)))


for pair, v in A.items():
    W = v['W']
    for ls, x in v['leg_sets'].items():
        n = len(ls.split()) if ls != '(none)' else 0
        tot = {}
        for key, d in x.items():
            for kind, entry in d.items():
                for idx, s in entry.items():
                    for mon, c in L.parse_A(s, n).items():
                        if abs(c) > 0:
                            tot.setdefault(kind, set()).add(mon)
        for kind, mons in tot.items():
            check('A_W0', pair, W, ls, kind, mons, 4 * n)
for pair, v in B.items():
    W = v['W']
    for ls, x in v['leg_sets'].items():
        n = len(x['outputs'])
        tot = {}
        for w, d in x['total'].items():
            for kind, s in d.items():
                for key, poly in L.parse_B(s, x['outputs'], n).items():
                    tot.setdefault(kind, set()).update(m for m, c in poly.items() if abs(c) > 0)
        for kind, mons in tot.items():
            if mons:
                check('B_W0', pair, W, ls, kind, mons, 4 * n)
import cmp_deriv_parse as CP
for key, v in Bd.items():
    pair = key.split(' @ ')[0]; W = v['W']
    for ls, x in v['leg_sets'].items():
        n = len(x['outputs'])
        tot = {}
        for w, d in x['total'].items():
            for kind, s in d.items():
                for kk, poly in CP.parse_Bd(s, n).items():
                    tot.setdefault(kind, set()).update(m for m, c in poly.items() if abs(c) > 0)
        for kind, mons in tot.items():
            if mons:
                check('B_deriv', pair, W, ls, kind, mons, 4 * n)

# predicted but zero
zeros = []
for src, data in (('A_W0', A), ('B_W0', B)):
    for pair, v in data.items():
        W = v['W']
        pr = predicted(pair, W)
        for ls in sorted(pr['log']):
            if (src, pair, str(W), ls, 'log') not in seen:
                zeros.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind='log'))
        for ls, mm in pr['eps'].items():
            for m in range(1, mm + 1):
                if (src, pair, str(W), ls, f'eps^-{m}') not in seen:
                    zeros.append(dict(src=src, pair=pair, W=str(W), legset=ls, kind=f'eps^-{m}'))
for key, v in Bd.items():
    pair = key.split(' @ ')[0]; W = v['W']
    pr = predicted(pair, W)
    for ls in sorted(pr['log']):
        if len(ls.split()) <= 2 and ls != '(none)' and ('B_deriv', pair, str(W), ls, 'log') not in seen:
            zeros.append(dict(src='B_deriv', pair=pair, W=str(W), legset=ls, kind='log'))
    for ls, mm in pr['eps'].items():
        if len(ls.split()) <= 2 or ls == '(none)':
            for m in range(1, mm + 1):
                if ('B_deriv', pair, str(W), ls, f'eps^-{m}') not in seen:
                    zeros.append(dict(src='B_deriv', pair=pair, W=str(W), legset=ls, kind=f'eps^-{m}'))
res = dict(n_nonzero_kernels=len(nonzero), n_monomials_checked=checked_monos, findings=findings,
           predicted_but_zero=zeros)
json.dump(res, open(os.path.join(L.QA, 'out.nosync', 'census_check.json'), 'w'), indent=1)
print('nonzero kernels', len(nonzero), 'monomials', checked_monos, 'findings', len(findings))
for f in findings:
    print('FINDING', f)
for z in zeros:
    print('ZERO', z)
