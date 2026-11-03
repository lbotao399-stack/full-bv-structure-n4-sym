"""QA: intrinsic K_b kernels, engine A (TQ6 rows, A basis, b = A's vector) versus engine B (symbolic polynomials in
results.json TQ5_TQ6, B basis), evaluated at A's data mapped by the intertwiner S (A's lamP = B's lamM).
B's polynomials are in tl_engine component names; the inverse map name -> slot component is built from
ce_job.T.comp_name (read-only import)."""
import os, sys, json, re
import sympy as sp
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'qa/src.nosync'))
import qq_ab_compare as AB          # noqa: E402 (engine B read-only, intertwiner)
CJ = AB.QJ.CJ
SLOT = AB.QJ.SLOT
B = json.load(open(os.path.join(QSK, 'engine_b/results.json')))['TQ5_TQ6']
A = json.load(open(os.path.join(QSK, 'engine_a/results.json')))['tasks']['TQ6']
bvec = A['b']


def inv_names(slot):
    spc, comps = SLOT[slot]
    return {CJ.T.comp_name(slot, c): c for c in comps}


def evalB(expr_str, pre, slotmap, data):
    """data: {prefix: {comp: value}}; pre: list of (prefix, slot)."""
    s = expr_str
    syms = {}
    def rep(m):
        p, nm = m.group(1), m.group(2)
        key = '%s__%s' % (p, re.sub(r'[^A-Za-z0-9]', '_', nm))
        syms[key] = (p, nm)
        return key
    s2 = re.sub(r'(Y0|Y1|Xa|Xb)\[([^\]]+)\]', rep, s)
    e = sp.sympify(s2)
    sub = {}
    for key, (p, nm) in syms.items():
        comp = slotmap[p][nm]
        sub[sp.Symbol(key)] = data[p].get(comp, 0)
    for m in range(4):
        sub[sp.Symbol('b%d' % m)] = bvec[m]
    for j in range(2):
        for m in range(4):
            sub[sp.Symbol('p%d_%d' % (j, m))] = data.get('mom%d' % j, [0, 0, 0, 0])[m]
    return sp.nsimplify(sp.expand(e.subs(sub)))


out = []
for r in A['intrinsic_rows']:
    if r['pair'] not in (['AW', 'lamP'], ['AW', 'lamM']) or len(r['ordering']) != 1:
        continue
    ta, tb = r['pair']
    d = r['data']
    tBa, za = AB.to_B(ta, d['zeta_a'], 'letter'); tBb, zb = AB.to_B(tb, d['zeta_b'], 'letter')
    to, po = AB.to_B(d['legs'][0][0], d['legs'][0][1], 'out')
    key = 'K %s|%s -> %s' % (tBa, tBb, to)
    ent = B[key]['intrinsic']
    slotmap = {'Xa': inv_names(tBa), 'Xb': inv_names(tBb), 'Y0': inv_names(to)}
    data = {'Xa': za, 'Xb': zb, 'Y0': po, 'mom0': d['legs'][0][2]}
    # c-power: one lam letter (covector) and one lam output (vector): e = 0
    tot = 0
    parts = {}
    for k, v in ent.items():
        basis = k.split('|')[-1].strip()
        val = evalB(v, None, slotmap, data)
        parts[basis] = str(val)
        tot += val * {'1': 1, 'log2': sp.log(2), 'log3': sp.log(3)}[basis]
    aval = AB.Aval(r['value'])
    rec = {'A_pair': r['pair'], 'A_leg': r['leg_set'], 'B_key': key, 'A_value': str(aval), 'B_value_at_A_data': str(sp.expand(tot)),
           'B_parts': parts, 'agree': sp.expand(tot - aval) == 0, 'agree_with_minus': sp.expand(tot + aval) == 0}
    print(rec, flush=True)
    out.append(rec)
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/K_compare.json'), 'w'), indent=1)
