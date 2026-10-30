"""Compare engine A (exact) and engine B (position space, numerical extrapolation) per class.
Reads ../engine_a/results.json and ../engine_b/results.json; writes comparison.json and comparison.txt.
Engine A's entries are exact sympy strings; engine B's may be exact strings or numeric values with digits."""
import os, json, re
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
A = json.load(open(os.path.join(ROOT, 'engine_a', 'results.json')))
B = json.load(open(os.path.join(ROOT, 'engine_b', 'results.json')))
beta2 = sp.Symbol('beta_2')
NS = {'EulerGamma': sp.EulerGamma, 'log': sp.log, 'pi': sp.pi, 'beta_2': beta2, 'gamma_E': sp.EulerGamma, 'gammaE': sp.EulerGamma}

def parse(v):
    """Return (sympy expr or None, float or None, digits or None)."""
    if v is None:
        return None, None, None
    if isinstance(v, dict):
        ex = v.get('exact'); nu = v.get('numeric') or v.get('value'); dg = v.get('digits')
        e = parse(ex)[0] if ex is not None else None
        f = None
        if nu is not None:
            try:
                f = float(str(nu).replace('*beta_2', '*0'))
            except Exception:
                f = None
        if e is not None and f is None:
            f = float(sp.N(e.subs(beta2, 0)))
        return e, f, dg
    if isinstance(v, (int, float)):
        return sp.nsimplify(v) if float(v).is_integer() else None, float(v), None
    s = str(v).strip()
    try:
        e = sp.sympify(s, locals=NS)
        return e, float(sp.N(e.subs(beta2, 0))), None
    except Exception:
        try:
            return None, float(s), None
        except Exception:
            return None, None, None

def getA(cls, struct, key):
    return A['classes'][cls][struct][key]

def getB(cls, struct, key):
    """Flexible lookup in engine B's file: classes[cls] may carry c_S/c_T/d_T/log_eps_S/log_eps_T, or S12/T12 dicts."""
    c = B['classes'].get(cls) or B['classes'].get(cls.lower()) or {}
    if struct in c and isinstance(c[struct], dict) and key in c[struct]:
        return c[struct][key]
    mp = {('S12', 'const'): 'c_S', ('T12', 'const'): 'c_T', ('T12', 'log_x2'): 'd_T', ('S12', 'log_x2'): 'd_S',
          ('S12', 'log_eps'): 'log_eps_S', ('T12', 'log_eps'): 'log_eps_T'}
    k = mp.get((struct, key))
    if k in c:
        return c[k]
    return None

header = ['class', 'structure', 'coefficient', 'engine A (exact)', 'engine A (num)', 'engine B', 'B digits', '|A-B|', 'agree']
rows = []; maxdiff = 0.0; ncomp = 0; nfail = 0
for cls in ['T', 'SE', 'B', 'EX', 'QU', 'OC', 'F']:
    if cls not in A['classes'] and cls != 'T':
        continue
    for struct in ['S12', 'T12']:
        for key in ['const', 'log_x2', 'log_eps']:
            if cls == 'T':
                # both engines normalise to the tree; the tree class carries no one-loop coefficient (0 in both files)
                ea, fa = sp.Integer(0), 0.0
            else:
                ea, fa, _ = parse(getA(cls, struct, key))
            vb = getB(cls, struct, key)
            eb, fb, db = parse(vb)
            if fb is None:
                rows.append([cls, struct, key, str(ea), '%.12g' % fa if fa is not None else '', 'missing', '', '', 'n/a'])
                continue
            diff = abs(fa - fb) if fa is not None else None
            tol = 10 ** (-(db if db else 6)) * max(1.0, abs(fa)) * 5
            ok = diff is not None and diff <= tol
            ncomp += 1; nfail += (0 if ok else 1); maxdiff = max(maxdiff, diff or 0)
            rows.append([cls, struct, key, str(ea), '%.12g' % fa, ('%.12g' % fb) if eb is None else str(eb), str(db or ''), '%.2e' % diff, 'yes' if ok else 'NO'])
# totals
for struct in ['S12', 'T12']:
    for key in ['const', 'log_x2', 'log_eps']:
        ea, fa, _ = parse(A['total'][struct][key])
        tb = B.get('total', {})
        vb = tb.get(struct, {}).get(key) if isinstance(tb.get(struct), dict) else None
        if vb is None:
            mp = {('S12', 'const'): 'c_S', ('T12', 'const'): 'c_T', ('T12', 'log_x2'): 'd_T', ('S12', 'log_x2'): 'd_S', ('S12', 'log_eps'): 'log_eps_S', ('T12', 'log_eps'): 'log_eps_T'}
            vb = tb.get(mp[(struct, key)]) or B.get(mp[(struct, key)])
        eb, fb, db = parse(vb)
        if fb is None:
            rows.append(['TOTAL', struct, key, str(ea), '%.12g' % fa, 'missing', '', '', 'n/a']); continue
        diff = abs(fa - fb); tol = 10 ** (-(db if db else 6)) * max(1.0, abs(fa)) * 5; ok = diff <= tol
        ncomp += 1; nfail += (0 if ok else 1); maxdiff = max(maxdiff, diff)
        rows.append(['TOTAL', struct, key, str(ea), '%.12g' % fa, ('%.12g' % fb) if eb is None else str(eb), str(db or ''), '%.2e' % diff, 'yes' if ok else 'NO'])
exA = parse(getA('EX','S12','const'))[1]; exB = B['classes']['EX'].get('numeric', {})
res_extra = {'EX_S12_const_engineB_raw_numeric': exB.get('value'), 'EX_S12_const_engineA_exact_float': exA,
             'abs_difference': (abs(exA - exB['value']) if exB.get('value') is not None else None), 'engineB_digits': exB.get('digits'),
             'note': 'engine B identifies the EX constant by PSLQ on {1, gamma_E, log 2, log 3} at 12 digits; engine A derives it in closed form'}
res = {'EX_numeric_check': res_extra, 'header': header, 'rows': rows, 'entries_compared': ncomp, 'entries_disagreeing': nfail, 'max_abs_difference': maxdiff,
       'note': 'engine A exact (beta_2 set to 0 for the numeric column); engine B numeric with its stated digits, tolerance 5 x 10^-digits x max(1,|A|)'}
json.dump(res, open(os.path.join(HERE, 'comparison.json'), 'w'), indent=1)
with open(os.path.join(HERE, 'comparison.txt'), 'w') as f:
    f.write('\t'.join(header) + '\n')
    for r in rows:
        f.write('\t'.join(r) + '\n')
    f.write('compared %d, disagreeing %d, max |A-B| = %.3e\n' % (ncomp, nfail, maxdiff))
print(open(os.path.join(HERE, 'comparison.txt')).read())
