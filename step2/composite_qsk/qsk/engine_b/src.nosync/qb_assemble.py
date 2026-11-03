"""Assemble results.json and gates.json of route B, PART Q, from the raw outputs in out.nosync/."""
import os, sys, json, glob, time

HERE = os.path.dirname(os.path.abspath(__file__))
EB = os.path.abspath(os.path.join(HERE, '..'))
OUT = os.path.join(EB, 'out.nosync')


def load(name):
    fn = os.path.join(OUT, name)
    if not os.path.exists(fn):
        return None
    return json.load(open(fn))


import re
from fractions import Fraction as _F
_R1 = re.compile(r'^\((-?\d+(?:/\d+)?)\)(\*I)?(.*)$')
_R2 = re.compile(r'^\((-?\d+(?:/\d+)?)\+\((-?\d+(?:/\d+)?)\)\*I\)(.*)$')


def neg_poly_str(s):
    if s == '0':
        return s
    out = []
    for t in s.split(' + '):
        m = _R2.match(t)
        if m:
            out.append(f'({-_F(m.group(1))}+({-_F(m.group(2))})*I){m.group(3)}')
            continue
        m = _R1.match(t)
        assert m, t[:80]
        out.append(f'({-_F(m.group(1))}){m.group(2) or ""}{m.group(3)}')
    return ' + '.join(out)


def neg_tree(x):
    if isinstance(x, dict):
        return {k: neg_tree(v) for k, v in x.items()}
    if isinstance(x, str):
        return neg_poly_str(x)
    return x


def fix_sign(d):
    """TQ2 files written before the sign fix (no sign_version 2) carry -D for odd x: negate them."""
    if d.get('sign_version') == 2:
        return d
    d = dict(d)
    for key in ('kernel', 'planar_eps0'):
        if d.get(key):
            d[key] = neg_tree(d[key])
    d['sign_version'] = '2 (negated at assembly from v1)'
    return d


def gate_rows(src, rows):
    out = []
    for g in rows or []:
        out.append(dict(id=g['id'], value=g.get('value'), target=g.get('target'), pass_or_fail=('pass' if g.get('pass_') else 'FAIL'),
                        note=g.get('note', ''), source=src))
    return out


def run():
    res = {}
    gates = []
    res['meta'] = dict(
        title='Route B (engine B), PART Q: defect cycles of Q_1 (xi), S_1 (eta), K_1 (K_b); momentum-space numerator algebra, Feynman parameters, numerical eps-extrapolation',
        date=time.strftime('%Y-%m-%d %H:%M'),
        theory='step2/composite_qsk/qsk/theory/qsk_theory.tex, engine_spec.json',
        unit='coefficients in units hbar/(16 pi^2) times the colour words C_w = (1 (x) ad_{Y_w1} ... ad_{Y_wr})(-,-)^{-1} (X_a on the first factor); planar adjacent pair: (-1)^r a/2 units u = -g^2 N/(32 pi^2), outputs in trace order Y_wr ... Y_w1',
        names='Xa[comp], Xb[comp]: letter components; Y<j>[comp]: output polarisations; p<j>_m: output momenta (i p_m = d_m on that output); xi<j>: components of xi in Sigma (realisation of engine_b/realization.py: j = 4a + A, a in {0,1} the S_+ index, A in 0..3 the W index; xi in S_+ (x) W has j = 0..7); eta<j>: Sigma^* components; b<m>: vector; engine vector index 0 <-> two-component m = 4, 1,2,3 <-> 1,2,3 (step2/notation/dictionary.json); comp names of tl_engine: A_m, A_I1..A_I6, lam_a<a>^<i> (lambda_+ = lambda_alpha^i), lambar_ad<a>_<i> (lambda_- = lambdabar_alphadot i), lams^a<a>_<i> (lambda*_+), lambars^ad<a><i> (lambda*_-), As_*, cs',
        regulator='Schwinger cutoff t >= eps on every edge P(eps, oo) (T = infinity at generic momenta), special edge K_eps - K_0 with t in [0, eps] and numerator -|k|^2 k0; no dimensional regularisation; no IR regulator other than T (= infinity here)',
        machine='pure Python (sympy rings over Q(i), fractions, mpmath); one heavy process at a time; code src.nosync/, raw outputs out.nosync/')
    e0 = load('e0_conventions.json')
    if e0:
        res['E0'] = e0['results']
        gates += gate_rows('E0 qb_e0.py', e0['gates'])
    tq0 = load('tq0_directions.json')
    if tq0:
        res['TQ0'] = tq0['results']
        gates += gate_rows('TQ0 qb_tq0.py', tq0['gates'])
    gr = load('greg_Pb_AW.json')
    if gr:
        target = {"word [0, 1] eps^-0 1": {"i(b.p0) (a.Y1)(b.Y0)": "-1/3", "i(b.p1) (a.Y1)(b.Y0)": "-1/3"},
                  "word [1, 0] eps^-0 1": {"i(b.p0) (a.Y0)(b.Y1)": "-1/3", "i(b.p1) (a.Y0)(b.Y1)": "-1/3"}}
        ok = True
        for k, d in target.items():
            got = gr['decomposition'].get(k, {})
            for kk, v in d.items():
                if got.get(kk) != v:
                    ok = False
            for kk, v in got.items():
                if kk not in d and v != '0':
                    ok = False
        extra = [k for k in gr['decomposition'] if k not in target]
        ok = ok and not extra
        res['TQ7'] = dict(G_reg_Pb=dict(pair='(A_perp, A_perp) -> A_perp A_perp', decomposition=gr['decomposition'],
                                        composite_engine_B_target=target, seconds=gr['seconds'],
                                        statement='sum over pairs of Upsilon_{P_b} on the scalar bilinear = C^(2)_{P_b} = -(1/6) u b.d(A A) per adjacent pair = -(1/3) u d_b O (Cor qs-DP); no eps^-1 term'),
                      Delta='not run in route B (needs the position part of S^(1)(Delta) with the injected momentum; open)')
        gates.append(dict(id='G-reg (P_b)', value=gr['decomposition'], target='composite engine B so6_AWAW.json word by word: -1/3 i(b.p0) - 1/3 i(b.p1) on the exchanged structure', pass_or_fail='pass' if ok else 'FAIL',
                          note='calibrates the Koszul ordering for even x: CE coordinate of x first, then the parameter of the special edge', source='TQ7 qb_greg.py'))
    # TQ2
    files = sorted(glob.glob(os.path.join(OUT, 'tq2', '*.json')))
    tq2 = []
    comp = {}
    for fn in files:
        d = fix_sign(json.load(open(fn)))
        if d.get('xi', '').startswith('xiP_e'):
            key = (d['pair'], d['legset'])
            comp.setdefault(key, []).append(d)
            continue
        tq2.append(d)
    have_sym = set((d['pair'], d['legset']) for d in tq2)
    for key, lst in comp.items():
        if key in have_sym:
            continue
        lst = sorted(lst, key=lambda x: x['xi'])
        js = [x['xi'][5:] for x in lst]
        lab = ('by components: Upsilon = sum_j xi_j (entry for xi = e_j), j = 0..7' if len(js) == 8 else
               'at xi = e_j for j in ' + ','.join(js) + ' (xi = e_3 is the rank-one u (x) e_4 of the twist; with symbolic letters and outputs it determines the kernel for every xi by Spin(4) x SU(4) covariance and linearity, S_+ (x) W being irreducible)')
        tq2.append(dict(lst[0], xi=lab,
                        kernel_by_xi_component={x['xi']: x['kernel'] for x in lst},
                        planar_eps0_by_xi_component={x['xi']: x['planar_eps0'] for x in lst},
                        kernel=None, planar_eps0=None, seconds=sum(x['seconds'] for x in lst),
                        nterms=sum(x['nterms'] for x in lst), components_done=len(lst)))
    fsp = load('fsp_gpole.json')
    fs = fsp['results'] if fsp else {}
    for d in tq2:
        key = None
        if d['legset'] in ('lamP', 'lamM') and d['pair'] in ('lamP|lamP', 'lamP|lamM', 'lamM|lamM'):
            key = f"{d['pair']} -> {d['legset']}"
        d['finite_special_vertex'] = fs.get(key, {}).get('F_xi') if key else '0 (no two-letter pole ct_1 with this leg set, no ct1_xi at n = 0: weight count)'
        d['ring_names_note'] = 'kernel strings use the renamed generators of meta.names; ring generator list kept in ring'
    res['TQ2'] = dict(n=0, xi='S_+ (x) W, symbolic xi0..xi7 (or by components)', legsets=tq2,
                      summary=dict(n_legsets=len(tq2), zero_legsets=[f"{d['pair']} [{d['legset']}]" for d in tq2 if d.get('nterms') == 0],
                                   bases=sorted(set(b for d in tq2 for b in d.get('bases', []))),
                                   poles=sorted(set(f"{d['pair']} [{d['legset']}]" for d in tq2 if d.get('poles_present')))),
                      Upsilon='Upsilon_xi = eps^0 coefficient (kernel, key eps^-0) + finite_special_vertex',
                      finite_special_vertices={k: v for k, v in fs.items() if k != 'seconds'})
    if fsp:
        gates += gate_rows('G-pole qb_fsp.py', fsp['gates'])
    nm = load('numeps.json')
    if nm:
        res['numerical_eps'] = nm['results']
        for name, v in nm['results'].items():
            if not isinstance(v, dict):
                continue
            rows = v['rows']
            ok = bool(rows) and all(float(r['diff_eps0']) < 1e-8 and float(r['diff_eps_m1']) < 1e-12 and
                                    abs(complex(r['fit_log_eps'].replace(' ', '').replace('(', '').replace(')', ''))) < 1e-8 for r in rows)
            gates.append(dict(id=f'G-num {name}', value=[(r['seq'], r['sigma'], r['diff_eps0'], r['diff_eps_m1'], r['fit_log_eps']) for r in rows],
                              target='finite-eps quadrature (32 digits), nine points eps = 2^-10..2^-18, fit: eps^0 within 1e-8, eps^-1 within 1e-12 of the exact coefficients, log eps coefficient below 1e-8 (G-noLog numeric)',
                              pass_or_fail='pass' if ok else 'FAIL', note='', source='qb_numeps.py'))
    bil = load('bilinears_xi.json')
    if bil:
        res['TQ4'] = bil['results']
        gates += gate_rows('G-plane qb_bil.py', bil['gates'])
    for fn in sorted(glob.glob(os.path.join(OUT, 'etak_eta_*.json')) + glob.glob(os.path.join(OUT, 'etak_K*.json'))):
        ek = json.load(open(fn))
        res.setdefault('TQ5_TQ6', {}).update(ek['results'])
        gates += gate_rows('qb_etak.py', ek['gates'])
    nt = load('numT.json')
    if nt:
        res['numerical_endpoint_T'] = nt['results']
        okT = True; vals = {}
        for k, v in nt['results'].items():
            if not (isinstance(v, dict) and 'by_T' in v):
                continue
            for T, r in v['by_T'].items():
                d0 = float(r['diff_eps0']); d1 = float(r['diff_eps_m1'])
                conv = True
                if r.get('T_minus_infinity_over_eps_at_eps_2m12_2m16_2m20'):
                    zs = [complex(x.replace(' ', '').replace('(', '').replace(')', '')) for x in r['T_minus_infinity_over_eps_at_eps_2m12_2m16_2m20']]
                    conv = abs(zs[2] - zs[1]) < 1e-3 * abs(zs[2])
                okT = okT and d0 < 1e-8 and d1 < 1e-12 and conv
                vals[f'{k} T={T}'] = dict(diff_eps0=r['diff_eps0'], diff_eps_m1=r['diff_eps_m1'], O_eps_T_dependence=conv)
        gates.append(dict(id='G-KT analogue (endpoint independence, route B numerics)', value=vals,
                          target='at T|p|^2 = 10.7, 42.8, 171 (T = 1, 4, 16) and T = infinity the fitted eps^0 (to < 1e-8, about 10 significant digits) and eps^-1 (< 1e-12) coefficients of the T5 cycle equal the exact T = infinity values; (value_T - value_infinity)/eps converges (T-dependence O(eps), Thm qs-simplex (iv))',
                          pass_or_fail='pass' if okT else 'FAIL', note='quadrature split at the kink eps/alpha_sigma = T/(1 - alpha_sigma)', source='qb_numT.py'))
    gm = load('gmisc.json')
    if gm:
        res['gates_misc'] = gm['results']
        gates += gate_rows('qb_gmisc.py', gm['gates'])
    res['TQ1'] = dict(status='not computed in route B (the action side along xi: quadratic Ob^(1)(xi), {S^(1)(xi), hbar F_1^(0)} and the solution hbar F_1^(1)(xi) are open); consequence: the one-body terms -{hbar F_1^(1)(xi), -} of Q_1 are known only up to the particular solution, the free directions theta^xi_V S^V_2 + theta^xi_perp S^perp_2 (TQ0) being the scheme family')
    res['TQ8'] = dict(status='later (needs Gram_1 from the two-point graphs of Def c8-2pt)')
    t3 = load('tq3_twist_n0.json')
    if t3:
        res['TQ3'] = t3
    # closed forms of the constant-numerator integrals (exact, ce_int.J) and of the n = 1 bubble (Prop qs-small (i))
    sys.path.insert(0, os.path.join(EB, '..', '..', 'composite', 'engine_b', 'src.nosync'))
    import ce_int as CI
    from math import factorial
    from fractions import Fraction as Fr
    cf = {}
    okcf = True
    for n in range(2, 7):
        v0 = CI.J(0, 0, tuple([0] * (n - 1)), tuple([1] * (n - 1)))
        v1 = CI.J(1, 0, tuple([0] * (n - 1)), tuple([1] * (n - 1)))
        target0 = CI.LogNum.log(Fr(n, n - 1)) * Fr(1, factorial(n - 1))
        target1 = Fr(1, n * n * (n - 1) * factorial(n - 1))
        ok0 = (v0 - target0).is_zero(); ok1 = (v1 - CI.LogNum(target1)).is_zero()
        okcf = okcf and ok0 and ok1
        cf[f'n_t = {n}'] = dict(eps0=str(v0), eps_m1=str(v1), eps0_closed='log(n_t/(n_t-1))/(n_t-1)!', eps_m1_closed='1/(n_t^2 (n_t-1) (n_t-1)!)', match=ok0 and ok1)
    res['closed_forms'] = dict(
        constant_numerator=dict(statement='a cycle with n_t edges (one special) whose Gaussian coefficient F_m is a constant F contributes 16 pi^2 [eps^0] = F log(n_t/(n_t-1))/(n_t-1)! and 16 pi^2 [eps^-1] = F/(n_t^2 (n_t-1) (n_t-1)!) (Thm qs-simplex with F constant); this is the whole transcendental content of every leg set with superficial degree D = 0',
                                values=cf),
        bubble='n_t = 2 (Prop qs-small (i)): int_0^{1/2} log((1-alpha)/alpha) dalpha = log 2; the G-cubic chain gives 16 pi^2 x value = (-1/(2 eps) - (3/8)|q|^2) x factor (gate G-cubic)',
        D0_legsets='for D = 0 leg sets the kernel is a rational multiple of log(3/2)/2 (n_t = 3) or of log(4/3)/6 (n_t = 4) chain by chain (only the bases log2, log3 with opposite coefficients occur)')
    gates.append(dict(id='G-closed (constant numerator)', value=cf, target='J_0(n_t) = log(n_t/(n_t-1))/(n_t-1)!, J_1(n_t) = 1/(n_t^2 (n_t-1)(n_t-1)!), n_t = 2..6 exact', pass_or_fail='pass' if okcf else 'FAIL', note='', source='qb_assemble.py (ce_int.J)'))
    st = load('status_extra.json')
    if st:
        res.update(st.get('results', {}))
        gates += st.get('gates', [])
    json.dump(res, open(os.path.join(EB, 'results.json'), 'w'), indent=1)
    json.dump(dict(gates=gates, n_pass=sum(1 for g in gates if g['pass_or_fail'] == 'pass'),
                   n_fail=sum(1 for g in gates if g['pass_or_fail'] == 'FAIL')), open(os.path.join(EB, 'gates.json'), 'w'), indent=1)
    return res, gates


if __name__ == '__main__':
    res, gates = run()
    print('gates', sum(1 for g in gates if g['pass_or_fail'] == 'pass'), 'pass,', sum(1 for g in gates if g['pass_or_fail'] == 'FAIL'), 'fail')
    print('results.json', os.path.getsize(os.path.join(EB, 'results.json')) // 1024, 'KB')
