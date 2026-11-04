"""OP1/OP3 tables for the even generators and the pair conjecture (G-D1pair), route B, from the stored rows.

For every field/ghost pair and leg set (n = 0, and n = 1 where the rows exist):
  Upsilon_{P_b}  (b symbolic: b0..b3), Upsilon^0_Delta, Upsilon^0_{M_w} (w symbolic, w_mn = -w_nm, m < n),
  Upsilon_R (R symbolic, R_MN, 4 <= M < N <= 9), Upsilon^0_{K_b} (b symbolic), each split into the placement classes
  'line' (the composite pair kernel C_x = [eps^0] of the special-edge composite cycles) and 'special', the poles
  (eps^-1, eps^-2), the finite special vertex F_x (S-channel, ob_fsv), and lc (log eps coefficient of the arity-zero
  counterterm of the pair, LC rows) with ct_m.
Gates per row:  G-noLog (no log eps in any defect row: structural, checked on the stored keys);
  G-D1pair: Upsilon^0_Delta = 4 lc + C^0_Delta, Upsilon_P = C_P, Upsilon^0_M = C^0_M, Upsilon_R = 0, i.e.
            [eps^0] special + F = 4 lc (Delta) and 0 (P_b, M_w, R), per colour word and per transcendental base;
  the same statement for K^0 is recorded (no theorem).
Position parts: Upsilon_Delta = 2 x_O^mu Upsilon_{P_mu} + Upsilon^0_Delta etc. (Thm op-b (i)) hold by construction in
route B (the position factor x_O^mu is the omega-part of the injected momentum)."""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ob_asm as A
import ob_fsv as FS
import ob_rows as RW
import ob_gen as G
import sympy as sp

OUTF = os.path.join(RW.OUT, 'tables')
os.makedirs(OUTF, exist_ok=True)
BS = list(sp.symbols('b0:4'))
WS = {(m, n): sp.Symbol(f'w{m}{n}') for m in range(4) for n in range(m + 1, 4)}
RS = {(M, N): sp.Symbol(f'R{M}{N}') for M in range(4, 10) for N in range(M + 1, 10)}


def wmat(m, n):
    w = [[0] * 4 for _ in range(4)]
    w[m][n] = 1; w[n][m] = -1
    return w


def rmat(M, N):
    r = [[0] * 6 for _ in range(6)]
    r[M - 4][N - 4] = 1; r[N - 4][M - 4] = -1
    return r


def to_sym(row, dd):
    """{(cls, word, m, base): ring dict} -> {cls: {word: {m: {base: expr}}}}"""
    out = {}
    for (cls, w, m, b), d in dd.items():
        e = A.to_expr(row, d)
        if e != 0:
            out.setdefault(cls, {}).setdefault(w, {}).setdefault(m, {})[b] = e
    return out


def combine(parts):
    """sum of [(coef, {cls: {w: {m: {b: expr}}}})]."""
    out = {}
    for c, P in parts:
        for cls, tw in P.items():
            for w, tm in tw.items():
                for m, tb in tm.items():
                    for b, e in tb.items():
                        dd = out.setdefault(cls, {}).setdefault(w, {}).setdefault(m, {})
                        dd[b] = sp.expand(dd.get(b, 0) + c * e)
    return prune(out)


def prune(T):
    out = {}
    for cls, tw in T.items():
        for w, tm in tw.items():
            for m, tb in tm.items():
                for b, e in tb.items():
                    if e != 0:
                        out.setdefault(cls, {}).setdefault(w, {}).setdefault(m, {})[b] = e
    return out


def gen_tables(pair, ls, n, gens=('P', 'D0', 'M0', 'R', 'K0')):
    ep = A.Row.load('EP', pair, ls, n)
    esw = A.Row.load('ESW', pair, ls, n)
    lc = A.Row.load('LC', pair, ls, n)
    T = {}
    rowref = ep or esw
    if rowref is None:
        return None
    if ep is not None and 'P' in gens:
        T['P'] = combine([(BS[r], to_sym(ep, A.ups_P(ep, [1 if s == r else 0 for s in range(4)], by_class=True))) for r in range(4)])
    T['D0'] = to_sym(rowref, A.ups0_Delta(ep, esw, by_class=True))
    if 'M0' in gens:
        T['M0'] = combine([(s, to_sym(rowref, A.ups0_M(ep, esw, wmat(m, n_), by_class=True))) for (m, n_), s in WS.items()])
    if esw is not None and 'R' in gens:
        T['R'] = combine([(s, to_sym(esw, A.ups_R(esw, rmat(M, N), by_class=True))) for (M, N), s in RS.items()])
    if 'K0' in gens:
        T['K0'] = combine([(BS[r], to_sym(rowref, A.ups0_K(ep, esw, [1 if s == r else 0 for s in range(4)], by_class=True)))
                       for r in range(4)])
    F = {}
    LC = {}
    if lc is not None:
        for (w, m, b), d in A.plain(lc).items():
            LC.setdefault(w, {})[m] = A.to_expr(lc, d)
        F['P'] = {}
        for r in range(4):
            for w, e in FS.fsv(lc, 'P', dict(b=[1 if s == r else 0 for s in range(4)])).items():
                F['P'][w] = sp.expand(F['P'].get(w, 0) + BS[r] * e)
        F['D0'] = FS.fsv(lc, 'D0', {})
        F['M0'] = {}
        for (m, n_), s in (WS.items() if 'M0' in gens else []):
            for w, e in FS.fsv(lc, 'M0', dict(w=wmat(m, n_))).items():
                F['M0'][w] = sp.expand(F['M0'].get(w, 0) + s * e)
        F['R'] = {}
        for (M, N), s in (RS.items() if 'R' in gens else []):
            for w, e in FS.fsv(lc, 'R', dict(R=rmat(M, N))).items():
                F['R'][w] = sp.expand(F['R'].get(w, 0) + s * e)
        F['K0'] = {}
        for r in (range(4) if 'K0' in gens else []):
            for w, e in FS.fsv(lc, 'K0', dict(b=[1 if s == r else 0 for s in range(4)])).items():
                F['K0'][w] = sp.expand(F['K0'].get(w, 0) + BS[r] * e)
        F = {g: {w: e for w, e in d.items() if e != 0} for g, d in F.items()}
    # evaluation point: if any row of the leg set used fixed numeric polarisations / letters, evaluate all at that point
    subs = {}
    for row in (ep, esw, lc):
        if row is None:
            continue
        if row.d.get('pol', 'sym') != 'sym':
            for j, d in enumerate(row.d['pol']):
                for c, v in d.items():
                    subs[sp.Symbol(f'o{j}_{c}')] = sp.sympify(v)
        lt = row.d.get('letters', 'sym')
        if lt != 'sym':
            for key, d in lt.items():
                for c, v in d.items():
                    subs[sp.Symbol(f'{key}_{c}')] = sp.sympify(v)
    if subs:
        def ev(x):
            if isinstance(x, dict):
                return {k: ev(v) for k, v in x.items()}
            return sp.expand(sp.sympify(x).subs(subs))
        T = prune(ev(T)); F = {g: {w: e for w, e in ev(d).items() if e != 0} for g, d in F.items()}; LC = ev(LC)
    numeric = bool(subs)
    return dict(T=T, F=F, LC=LC, has=dict(EP=ep is not None, ESW=esw is not None, LC=lc is not None),
                D=dict(EP=ep.d['D'] if ep else None, ESW=esw.d['D'] if esw else None),
                pol=('numeric evaluation point: ' + str({str(k): str(v) for k, v in subs.items()}) if numeric else 'sym'), outs=rowref.d['outs'],
                subs=subs)


def d1pair(tab, withF):
    """{gen: (pass, residual dict)}: [eps^0] special (+ F if withF) - target, per word and base."""
    res = {}
    for g, target in (('D0', 'lc'), ('P', None), ('M0', None), ('R', None)):
        if g not in tab['T']:
            continue
        sp_ = tab['T'][g].get('special', {})
        words = set(sp_.keys()) | (set(tab['F'].get(g, {}).keys()) if withF else set())
        if target == 'lc':
            words |= set(w for w, d in tab['LC'].items() if 'log' in d)
        resid = {}
        for w in words:
            bases = set(sp_.get(w, {}).get('0', {}).keys()) | {'1'}
            for b in bases:
                e = sp_.get(w, {}).get('0', {}).get(b, 0)
                if b == '1':
                    if withF:
                        e += tab['F'].get(g, {}).get(w, 0)
                    if target == 'lc':
                        e -= 4 * tab['LC'].get(w, {}).get('log', 0)
                e = sp.expand(e)
                if e != 0:
                    resid[f'{w} | {b}'] = str(e)[:400]
        res[g] = (len(resid) == 0, resid)
    return res


def gpole(tab):
    """eq c11-pole for the pair rows (s ct^(1)_{x,m} = 0 by weight on these rows): [eps^-m] D^0_Delta = -4m ct_m
    (scaling: {S(Delta), ct_m(O)} - ct_m({S(Delta), O}) = 2(W - 2m) - 2W), [eps^-m] D_{P_b}, D^0_{M_w}, D_R = 0
    (covariance of the counterterm map)."""
    res = {}
    for g in ('D0', 'P', 'M0', 'R'):
        if g not in tab['T']:
            continue
        resid = {}
        words = set()
        for cls, tw in tab['T'][g].items():
            words |= set(tw.keys())
        words |= set(tab['LC'].keys())
        for w in words:
            for m in ('1', '2'):
                e = 0
                for cls, tw in tab['T'][g].items():
                    for b, ex in tw.get(w, {}).get(m, {}).items():
                        e += ex if b == '1' else ex * sp.Symbol(b)
                if g == 'D0':
                    e += 4 * int(m) * tab['LC'].get(w, {}).get(m, 0)
                e = sp.expand(e)
                if e != 0:
                    resid[f'{w} | eps^-{m}'] = str(e)[:300]
        res[g] = (len(resid) == 0, resid)
    return res


ENTRY_KEYS = ('pair', 'legset', 'n', 'outs', 'has', 'D', 'pol', 'G_noLog', 'G_D1pair', 'D1pair_residuals',
              'finite_special_vertex_nonzero', 'G_D1pair_with_theory_F', 'G_pole', 'pole_residuals')
RESUME_AFTER = time.mktime((2026, 11, 3, 0, 0, 0, 0, 0, -1))     # tables of the earlier production pass are rebuilt


def run(n=0, pairs=None):
    t0 = time.time()
    pairs = pairs or ['AV|AV', 'AV|AW', 'AV|lamP', 'AV|lamM', 'AV|c', 'AW|AW', 'AW|lamP', 'AW|lamM', 'AW|c', 'lamP|lamP',
                      'lamP|lamM', 'lamP|c', 'lamM|lamM', 'lamM|c', 'c|c']
    summary = []
    for pair in pairs:
        lsets = sorted(set(x['legset'] for x in RW.legsets('EP', pair, n)) | set(x['legset'] for x in RW.legsets('ESW', pair, n)))
        for ls in lsets:
            fn = os.path.join(OUTF, f"{pair.replace('|', '_')}__{ls.replace(' ', '_').replace('(', '').replace(')', '')}__n{n}.json")
            # resume: reuse a table written after its rows (the process can be restarted after the RSS watchdog)
            rowf = [RW.fname(k, pair, ls, n) for k in ('EP', 'ESW', 'LC')]
            rmt = max([os.path.getmtime(f) for f in rowf if os.path.exists(f)] or [0])
            if os.path.exists(fn) and os.path.getmtime(fn) > max(rmt, RESUME_AFTER):
                try:
                    d = json.load(open(fn))
                    if 'G_pole' in d and 'G_D1pair_with_theory_F' in d:
                        summary.append({k: d[k] for k in ENTRY_KEYS})
                        continue
                except Exception:
                    pass
            marker = fn + '.inprogress'
            if os.path.exists(marker):
                # this table exceeded the RSS watchdog in an earlier attempt: recorded as too large, skipped
                summary.append(dict(pair=pair, legset=ls, n=n, outs=[] if ls == '(none)' else ls.split(),
                                    has=dict(EP=False, ESW=False, LC=False), D={}, pol='table exceeded the memory cap',
                                    G_noLog=True, G_D1pair={}, D1pair_residuals={}, finite_special_vertex_nonzero=False,
                                    G_D1pair_with_theory_F='not built', G_pole={}, pole_residuals={}, table='exceeded the 1400 MB watchdog'))
                continue
            open(marker, 'w').write('started ' + time.strftime('%Y-%m-%d %H:%M:%S'))
            tab = gen_tables(pair, ls, n, gens=(('P', 'D0', 'M0', 'R', 'K0') if n == 0 else ('P', 'D0')))
            if tab is None:
                os.remove(marker)
                continue
            g = d1pair(tab, False)
            gF = d1pair(tab, True)
            gp = gpole(tab)
            nolog = all('log' not in m for gg in tab['T'].values() for tw in gg.values() for tm in tw.values() for m in tm)
            hasF = any(bool(d) for d in tab['F'].values())
            entry = dict(pair=pair, legset=ls, n=n, outs=tab['outs'], has=tab['has'], D=tab['D'], pol=tab['pol'],
                         G_noLog=nolog, G_D1pair={k: v[0] for k, v in g.items()},
                         D1pair_residuals={k: v[1] for k, v in g.items() if not v[0]},
                         finite_special_vertex_nonzero=hasF,
                         G_D1pair_with_theory_F={k: v[0] for k, v in gF.items()} if hasF else 'same (F = 0)',
                         G_pole={k: v[0] for k, v in gp.items()},
                         pole_residuals={k: v[1] for k, v in gp.items() if not v[0]})
            ser = {gname: {cls: {w: {m: {b: str(e) for b, e in tm.items()} for m, tm in tw.items()} for w, tw in tc.items()}
                           for cls, tc in gt.items()} for gname, gt in tab['T'].items()}
            json.dump(dict(entry, kernels=ser, F={gname: {w: str(e) for w, e in d.items()} for gname, d in tab['F'].items()},
                           LC={w: {m: str(e) for m, e in d.items()} for w, d in tab['LC'].items()}), open(fn, 'w'))
            summary.append(entry)
            os.remove(marker)
    return summary, round(time.time() - t0, 1)


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    s, sec = run(n)
    json.dump(dict(summary=s, seconds=sec), open(os.path.join(RW.OUT, f'tables_summary_n{n}.json'), 'w'), indent=1)
    for e in s:
        print(e['pair'], '[' + e['legset'] + ']', e['has'], 'noLog', e['G_noLog'], e['G_D1pair'], 'withF', e['G_D1pair_with_theory_F'], 'pole', e['G_pole'])
    print('seconds', sec)
