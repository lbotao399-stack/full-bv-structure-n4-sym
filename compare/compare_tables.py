"""Engine A vs engine B vs classical census vs FeynArts manifest: exact comparison of every reported coefficient.
Reads only the project files named in the task; writes compare_tables.json."""
import json, os
import sympy as sp

R = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
A = json.load(open(f'{R}/engine_a/results.json'))
B = json.load(open(f'{R}/engine_b/results.json'))
M = json.load(open(f'{R}/diagrams/manifest.json'))
C = json.load(open(f'{R}/classical/census.json'))
CC = json.load(open(f'{R}/diagrams/data/census_comparison.json'))
SA2 = json.load(open(f'{R}/engine_a/logs/sector_a_su2.json'))
SA3 = json.load(open(f'{R}/engine_a/logs/sector_a_su3.json'))

def q(x):
    return sp.nsimplify(sp.sympify(x)) if x not in (None, '') else sp.Integer(0)

Ag = {g['id']: g for g in A['graphs']}
Bg = {g['id']: g for g in B['graphs']}

# ---- matched list: census class -> FeynArts id -> engine A id -> engine B id
TWO = [
 (1, 'AV,AV:g3[VVV],g3[VVV]:P_AA(V),P_AA(V)', 'AV-AV:AAA+AAA:AV,AV', 'AV-AV:AAA+AAA:AV,AV'),
 (2, 'AV,AV:g3[VWW],g3[VWW]:P_AA(W),P_AA(W)', 'AV-AV:AAA+AAA:Aperp,Aperp', 'AV-AV:AAA+AAA:Ap,Ap'),
 (3, 'AV,AV:g4[VVVV]:P_AA(V)', 'AV-AV:AAAA:AV', 'AV-AV:AAAA:AV'),
 (4, 'AV,AV:g4[VVWW]:P_AA(W)', 'AV-AV:AAAA:Aperp', 'AV-AV:AAAA:Ap'),
 (5, 'AV,AV:gh[V],gh[V]:P_cA*,P_cA*', 'AV-AV:AcAs+AcAs:c,As', 'AV-AV:AcAs+AcAs:c-As,c-As'),
 (6, 'AV,AV:y[V],y[V]:P_lamlam,P_lamlam', 'AV-AV:llA+llA:l,l', 'AV-AV:lamlamA+lamlamA:lam,lam'),
 (7, 'AV,AW:y[V],y[W]:P_lamlam,P_lamlam', None, 'AV-Ap:AAA+AAA:all'),
 (8, 'AW,AW:g3[VWW],g3[VWW]:P_AA(V),P_AA(W)', 'Aperp-Aperp:AAA+AAA:AV,Aperp', 'Ap-Ap:AAA+AAA:Ap,AV'),
 (9, 'AW,AW:g4[VVWW]:P_AA(V)', 'Aperp-Aperp:AAAA:AV', 'Ap-Ap:AAAA:AV'),
 (10, 'AW,AW:g4[WWWW]:P_AA(W)', 'Aperp-Aperp:AAAA:Aperp', 'Ap-Ap:AAAA:Ap'),
 (11, 'AW,AW:y[W],y[W]:P_lamlam,P_lamlam', 'Aperp-Aperp:llA+llA:l,l', 'Ap-Ap:lamlamA+lamlamA:lam,lam'),
 (12, 'c,AsV:gh[V],gh[V]:P_AA(V),P_cA*', 'As-c:AcAs+AcAs:AV,c-As', 'As-c:AcAs+AcAs:AV,c-As'),
 (13, 'lam,lam:y[V],y[V]:P_AA(V),P_lamlam', 'l-l:llA+llA:l,AV', 'lam-lam:lamlamA+lamlamA:lam,AV'),
 (14, 'lam,lam:y[W],y[W]:P_AA(W),P_lamlam', 'l-l:llA+llA:l,Aperp', 'lam-lam:lamlamA+lamlamA:lam,Ap'),
]
EXTRA_ZERO = [  # entries carried by an engine that are not census graphs (a propagator component they need is zero)
 ('Aperp-Aperp:AcAs+AcAs:c,As', 'Ap-Ap:AcAs+AcAs:c-As,c-As', 'G^W vertex needs c (x) A*|_W-perp: no such component of P'),
 ('As-c:AcAs+AcAs:Aperp,c-As', None, 'same reason (FeynArts job tp_cAsW: 0 graphs)'),
 (None, 'AV-Ap:AAAA:all', 'no quartic two-point graph with external (A|_V, A|_W-perp): P has no A|_V (x) A|_W-perp component'),
]
INS = [
 ('O1', 'O:O:P_AA(W)', 'vacuum:O:Aperp', 'vac:O:Ap'),
 ('O2', 'O,AV:O,g3[VWW]:P_AA(W),P_AA(W)', None, None),
 ('O3', 'O,AV,AV:O,g4[VVWW]:P_AA(W),P_AA(W)', 'AV-AV:O+AAAA:Aperp,Aperp', 'AV-AV+O:AAAA+O:Ap,Ap'),
 ('O4', 'O,AV,AV:O,g3[VWW],g3[VWW]:P_AA(W),P_AA(W),P_AA(W)', 'AV-AV:O+AAA+AAA:Aperp,Aperp,Aperp', 'AV-AV+O:AAA+AAA+O:Ap,Ap,Ap'),
 ('O5', 'O,AW,AW:O,g4[WWWW]:P_AA(W),P_AA(W)', 'Aperp-Aperp:O+AAAA:Aperp,Aperp', 'Ap-Ap+O:AAAA+O:Ap,Ap'),
 ('O6', 'O,AW,AW:O,g3[VWW],g3[VWW]:P_AA(V),P_AA(W),P_AA(W)', 'Aperp-Aperp:O+AAA+AAA:Aperp,AV,Aperp', 'Ap-Ap+O:AAA+AAA+O:Ap,AV,Ap'),
]
EXTRA_INS = [('Aperp-Aperp:O+AAA+AAA:Aperp,Aperp,Aperp', 'needs an (A|_W-perp)^3 vertex: none'),
             ('AV-AV:O+AAA+AAA:Aperp,AV,Aperp', 'needs an A|_V A|_V A|_W-perp vertex: none')]
KEYMAP = {'O': 'O', 'trace': 'trace_mix', 'AV_trace': 'AV_mix'}

fa_ids = {d['id'] for d in M['diagrams']}
cen2 = C['two_point']
out = {'two_point': [], 'insertion': [], 'extra_engine_entries': [], 'mismatches': []}

def coeffs(g, which, keymap=None):
    d = g.get(which, {}) if g else {}
    if not isinstance(d, dict):
        return {}
    r = {}
    for kk, v in d.items():
        kk2 = keymap.get(kk, kk) if keymap else kk
        if isinstance(v, str) and ('[' in v or 'dim' in v):
            v = v.split()[0]
        r[kk2] = q(v)
    return {kk: v for kk, v in r.items() if v != 0}

for (ci, fa, ia, ib) in TWO:
    ga, gb = Ag.get(ia), Bg.get(ib)
    cen = cen2[ci - 1]
    ent = dict(census_index=ci, census_Aut=cen['Aut'], census_omega=cen['omega'], feynarts_id=fa,
               feynarts_present=fa in fa_ids, engine_A_id=ia, engine_B_id=ib,
               A_sym=ga['symmetry_factor'] if ga else None, B_sym=gb['symmetry_factor'] if gb else None)
    for which in ('eps_inv', 'log_eps'):
        ca = coeffs(ga, which); cb = coeffs(gb, which)
        if ci == 7:  # engine A: sector total of logs/sector_a (su2, su3); engine B: masked-off-nothing run, raw value exactly 0
            ca = {nm: q(v) for nm, v in SA2['sectors']['AV-Aperp'][which + '_rational'].items() if q(v) != 0}
            cb = {}
        ent[which] = {'A': {kk: str(v) for kk, v in ca.items()}, 'B': {kk: str(v) for kk, v in cb.items()}, 'equal': ca == cb}
        if ca != cb:
            out['mismatches'].append((fa, which, ca, cb))
    out['two_point'].append(ent)

for (ia, ib, why) in EXTRA_ZERO:
    za = (not coeffs(Ag.get(ia), 'eps_inv') and not coeffs(Ag.get(ia), 'log_eps')) if ia else None
    zb = (not coeffs(Bg.get(ib), 'eps_inv') and not coeffs(Bg.get(ib), 'log_eps')) if ib else None
    out['extra_engine_entries'].append(dict(engine_A_id=ia, engine_B_id=ib, A_zero=za, B_zero=zb, reason=why))

for (ci, fa, ia, ib) in INS:
    ga, gb = Ag.get(ia), Bg.get(ib)
    ent = dict(census_index=ci, feynarts_id=fa, feynarts_present=fa in fa_ids, engine_A_id=ia, engine_B_id=ib)
    for which in ('eps_inv', 'log_eps'):
        if ci == 'O1':
            ca = {'vacuum': q(list(ga[which].values())[0])} if ga and ga[which] else {}
            cb = {'vacuum': q(list(gb[which].values())[0].split()[0])} if gb and gb[which] else {}
        else:
            ca = coeffs(ga, which, KEYMAP); cb = coeffs(gb, which)
        ent[which] = {'A': {kk: str(v) for kk, v in ca.items()}, 'B': {kk: str(v) for kk, v in cb.items()}, 'equal': ca == cb}
        if ca != cb:
            out['mismatches'].append((fa, which, ca, cb))
    out['insertion'].append(ent)
for (ia, why) in EXTRA_INS:
    g = Ag[ia]
    out['extra_engine_entries'].append(dict(engine_A_id=ia, engine_B_id=None,
                                           A_zero=not coeffs(g, 'eps_inv', KEYMAP) and not coeffs(g, 'log_eps', KEYMAP), reason=why))

# ---- totals
tot = {}
for kk in ('B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7'):
    a = A['ICT_quadratic'][kk]; b = B['ICT_quadratic'][kk]
    s = {'eps_inv': 0, 'log_eps': 0}
    for ent in out['two_point']:
        for which in s:
            s[which] += q(ent[which]['A'].get(kk, 0))
    tot[kk] = dict(A={w: str(q(a[w])) for w in a}, B={w: str(q(b[w])) for w in b},
                   sum_of_graphs={w: str(v) for w, v in s.items()},
                   equal=all(q(a[w]) == q(b[w]) == s[w] for w in s))
out['ICT_quadratic'] = tot
out['su2_su3_sector_totals_engineA_equal'] = all(SA2['sectors'][s][w + '_rational'] == SA3['sectors'][s][w + '_rational']
                                                for s in SA2['sectors'] for w in ('eps_inv', 'log_eps'))

# ---- intermediate polynomials: engine B C0(x) = 2 g(u), C1 = 4 c  (engine A integrand_and_faces vs engine B formula_chain_data)
fcd = B['formula_chain_data']['two_point_words']
STRUCT = {'struct0': "p^2 (A(-p),A(p))", 'struct1': "(p.A(-p))(p.A(p))"}
inter = []
x, u = sp.symbols('x u')
for (ci, fa, ia, ib) in TWO:
    if ia is None or ib not in fcd or 'integrand_and_faces' not in Ag[ia]:
        continue
    ga = Ag[ia]['integrand_and_faces']; w = fcd[ib]
    C0 = eval(w['C0'], {'I': sp.I, 'x': x}) if w.get('C0') else {}
    C0 = {kk: sp.sympify(v).subs(x, u) for kk, v in C0.items()}
    okg = True
    for s, expr in ga['g(u)'].items():
        e = sp.sympify(expr.replace('I', 'I'))
        name = STRUCT.get(s)
        if ib.startswith('lam') or ib.startswith('As'):
            name = list(C0.keys())[0] if C0 else None
        bval = C0.get(name, 0) if name else 0
        okg &= sp.simplify(bval - 2 * e) == 0
    C1 = eval(w['C1'], {'I': sp.I, 'x': x}) if w.get('C1') else {}
    cA = sp.sympify(list(ga['c'].values())[0]) if ga['c'] else 0
    cB = sp.sympify(list(C1.values())[0]) if C1 else 0
    inter.append(dict(graph=fa, C0_equals_2g=bool(okg), C1_equals_4c=bool(sp.simplify(cB - 4 * cA) == 0)))
out['intermediate_conversion'] = inter

# ---- O(eps), gamma
oa = A['O_eps']['scalar_bilinear']; ob = B['O_eps']['scalar_bilinear']['partial_delta_ICT']
out['partial_delta_ICT'] = {
 'O log': (oa['partial_delta_ICT_log_eps']['O'], ob['O']['log_eps']),
 'trace-mixing log': (oa['partial_delta_ICT_log_eps']['(w1,w2) int omega sum_w (iota_w A, iota_w A)'],
                      ob['(w1,w2) int omega sum_w (iota_w A, iota_w A)  [orthonormal basis of W_R^perp]']['log_eps']),
 'A_V-mixing log': (oa['partial_delta_ICT_log_eps']['(w1,w2) int omega (A_V,A_V)'], ob['(w1,w2) int omega (A_V, A_V)']['log_eps'].split()[0]),
 'vacuum eps^-1 [units hbar/16pi^2, x dim g (w1,w2)]': (oa['partial_delta_ICT_eps_inv']['vacuum: (w1,w2) dim(g) int omega, units hbar/(16 pi^2)'],
                                                        ob['constant (vacuum)']['eps_inv'].split()[0]),
 'field-dependent eps^-1': (oa['partial_delta_ICT_eps_inv']['field-dependent'], ob['O']['eps_inv'])}
out['partial_delta_ICT_equal'] = all(q(a) == q(b) for a, b in out['partial_delta_ICT'].values())
da = A['dilatation']; db = B['dilatation']
out['gamma'] = {
 'traceless': (da['traceless'], db['traceless']['gamma']),
 'trace': (da['trace'], db['trace']['gamma']),
 'letter iota_w A': (da['letters']['iota_w A (W^perp)'], db['letters']['int omega (X, iota_w A)']['gamma']),
 'letter A|_V transverse': (da['letters']['A_V'].split()[0], db['letters']['int omega (X, A_V(v))']['gamma_transverse']),
 'letter lambda': (da['letters']['lambda'], db['letters']['int omega (X, (zeta, lambda))']['gamma']),
 'letter c': (da['letters']['c'], None)}
out['gamma_equal'] = {kk: (b is not None and q(a) == q(b)) for kk, (a, b) in out['gamma'].items()}

# ---- FeynArts boundary sets
bd = [d for d in M['diagrams'] if d['kind'] == 'boundary']
out['feynarts_boundary'] = {
 'route_i_L_eps_pinned': len([d for d in bd if d['route'] == 'i']),
 'route_ii_K_eps_1PI': len([d for d in bd if d['route'] == 'ii' and d['one_PI']]),
 'route_ii_counterterm_trees': len([d for d in bd if d['route'] == 'ii' and not d['one_PI']]),
 'reducible_insertion_graphs': len([d for d in M['diagrams'] if d['kind'] == 'insertion' and not d['one_PI']]),
 'feynarts_vs_census_all_match': all(e['match'] for e in CC['two_point'] + CC['one_O_insertion_all']),
 'feynarts_only': CC['feynarts_only']}
json.dump(out, open(f'{R}/compare/compare_tables.json', 'w'), indent=1, default=str)
print('mismatches:', out['mismatches'])
print('ICT equal:', {kk: v['equal'] for kk, v in tot.items()})
print('intermediate:', inter)
print('partial_delta_ICT equal:', out['partial_delta_ICT_equal'], out['partial_delta_ICT'])
print('gamma:', out['gamma'], out['gamma_equal'])
print('boundary:', out['feynarts_boundary'])
print('extra:', out['extra_engine_entries'])
print('su2=su3 (A):', out['su2_su3_sector_totals_engineA_equal'])

# ---- engine B intermediate integrals -> final coefficients (structure basis of formula_chain_data)
conv = []
for (ci, fa, ia, ib) in TWO:
    w = fcd.get(ib)
    if not w or 'pole_integral' not in w:
        continue
    g = Bg[ib]
    pole = eval(w['pole_integral'], {'I': sp.I}); lg = eval(w['log_integral_A0'], {'I': sp.I})
    rec = {'graph': fa}
    if pole:
        key = 'B3' if ib.startswith('AV-AV') else 'B5'
        rec['pole_integral/4 == eps_inv'] = bool(sp.nsimplify(list(pole.values())[0]) / 4 == q(g['eps_inv'][key]))
    if ib.startswith('AV-AV') or ib.startswith('Ap-Ap'):
        kin = 'B1' if ib.startswith('AV-AV') else 'B4'
        b1 = q(g['log_eps'].get(kin, 0)); b2 = q(g['log_eps'].get('B2', 0))
        ok = True
        if "p^2 (A(-p),A(p))" in lg:
            ok &= -sp.nsimplify(lg["p^2 (A(-p),A(p))"]) / 4 == b1
        if "(p.A(-p))(p.A(p))" in lg:
            ok &= -sp.nsimplify(lg["(p.A(-p))(p.A(p))"]) / 4 == -b1 + b2
        rec['-log_integral/4 == log_eps (structures p^2(u.v) -> B_kin, (p.u)(p.v) -> -B_kin+B2)'] = bool(ok)
    else:
        val = sp.nsimplify(list(lg.values())[0])
        kk = 'B6' if ib.startswith('lam') else 'B7'
        rec['log_integral / log_eps'] = str(sp.simplify(val / q(g['log_eps'][kk])))
    conv.append(rec)
out['engineB_integral_conversion'] = conv
json.dump(out, open(f'{R}/compare/compare_tables.json', 'w'), indent=1, default=str)
for r in conv:
    print(r)
