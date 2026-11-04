"""Reviser: revised gate verdicts after the dispositions (no number changes; relabels and re-evaluations).
Writes gates_revised.json next to each engine's gates.json; originals untouched."""
import json, os, re
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
def load(p): return json.load(open(os.path.join(ROOT, p)))
def dump(p, d): json.dump(d, open(os.path.join(ROOT, p), 'w'), indent=1)

# ---------- operator engine B ----------
gB = load('operator/engine_b/gates.json')
cmp = load('operator/qa/out.nosync/cmp_tables.json')
incomplete = {(c['gen'], c['key'][0], c['key'][1], c['key'][2]) for c in cmp['conj_disagree']}
inc_labels = set()
for (gen, pair, outs, n) in incomplete:
    a, b = eval(pair); o = eval(outs)
    inc_labels.add((gen, n, f"{a}|{b} [{' '.join(o)}]"))
    inc_labels.add((gen, n, f"{b}|{a} [{' '.join(o)}]"))
rev = []
for x in gB['gates']:
    y = dict(x); gid = x['id']; y['original_pass'] = x['pass_']
    if 'with the theory finite vertex F_x' in gid:
        y['revised'] = 'withdrawn: tests Upsilon = [eps^0]D + F_x, the double-counted definition (QA-F1, accepted); its failure is evidence for the correction'
        y['pass_revised'] = None
    elif gid.startswith('G-pole reduced form') and not x['pass_']:
        m = re.search(r'\((n=\d)\)', gid); n = m.group(1)[2:]; gen = gid.split()[3]
        fails = re.findall(r"'([^']+)'", str(x['value']))
        nov = [f for f in fails if (gen, n, f) in inc_labels or ('lam' in f.split('|')[0] and ' ' in f.split('[')[1])]
        line = [f for f in fails if f not in nov]
        y['no_verdict_rows'] = nov
        y['line_class_rows_B_F3'] = line
        y['revised'] = ('rows with a missing EP/ESW class: no verdict (QA-F6); remaining failures are the (A_V,X)->X line-class rows, '
                        'for which the full eq c11-pole (ct^(1)(sO) and spectator terms) is not implemented: open (QA-F8)')
        y['pass_revised'] = 'open'
    else:
        y['pass_revised'] = x['pass_']
    rev.append(y)
cnt = dict(pass_=sum(1 for y in rev if y['pass_revised'] is True), withdrawn=sum(1 for y in rev if y['pass_revised'] is None),
           open=sum(1 for y in rev if y['pass_revised'] == 'open'), fail=sum(1 for y in rev if y['pass_revised'] is False))
dump('operator/engine_b/gates_revised.json', dict(gates=rev, counts=cnt, incomplete_rows_no_verdict=sorted(map(list, incomplete)),
     note='G-D1pair verdicts on the 11 rows with a missing class are "no verdict" in route B; engine A has them complete and the pair form holds (qa cmp_tables).'))
print('operator B', cnt)

# ---------- action engine A ----------
gA = load('action/engine_a/gates.json')
dxi = load('action/qa/q_dxi_fit.json')
revA = []
for x in gA:
    y = dict(x); y['original_pass'] = x['pass']
    if x['id'].startswith('G-WZ2 (Delta,xi'):
        y['revised'] = ('replaced: the gate demanded Q-closedness of Ob2 and dropped the position part (finding F2). Corrected gate on the full functional: '
                        'Q Ob2(Delta,xi)|_2 = -25/3 u(g1+g3) = -{Ob1(Delta)|_2,S(xi)} (WZ2 with corrected signs), the log(4/3) part = log(4/3) Q(psi1+psi2); '
                        'pass in qa/q_dxi_fit.json and in action/complete routes A and B (Ob2(Delta,xi+-): solvable, theta^Delta_d = -1/3).')
        y['pass_revised'] = True
    elif x['id'].startswith('G-logex'):
        y['revised'] = 'closed combination lc1(x) - {uN,S(x)} (corrected sign, finding F1); the printed + fails; consistent with the revised Prop at-log'
        y['pass_revised'] = True
    elif x['id'].startswith('G-P2') or 'G-pos' in x['id']:
        y['pass_revised'] = x['pass']; y['revised'] = 'direct P_b cycles computed (not structural)'
    else:
        y['pass_revised'] = x['pass']
    revA.append(y)
revA.append(dict(id='G-log (finite-eps extrapolation)', pass_revised='not run', revised='neither action engine ran it; closed forms of Prop at-small checked by Monte Carlo in qa (does not replace it)'))
dump('action/engine_a/gates_revised.json', revA)
print('action A', sum(1 for y in revA if y['pass_revised'] is True), 'pass of', len(revA))

# ---------- action engine B ----------
gB2 = load('action/engine_b/gates.json')
items = gB2 if isinstance(gB2, list) else gB2.get('gates', gB2)
revB = []
for x in items:
    y = dict(x)
    if x['id'].startswith('G-H'):
        y['id_revised'] = 'G-S1xi2 (two-leg cohomology along Sigma, Prop at-xi2)'
        y['revised'] = 'renamed (finding F10): not Hypothesis at-H at three legs'
    if x['id'].startswith('G-sign'):
        y['revised'] = ('convention text corrected (finding F6): the code gives Q int(lam,lam*) = -int(lam,dslash lam), and (lam*,lam) = +(lam,lam*) since lam* is even; '
                        'Q N = [-1,4,-1] (rerun of qa q_lcN.py)')
    if x['id'].startswith('G-WZ1'):
        y['revised'] = 'passes with Lambda Ob1(x) = -{Ob0,S(x)}, the corrected Lemma at-sformulas (ii)'
    revB.append(y)
revB.append(dict(id='G-log (finite-eps extrapolation)', status='not run', revised='structural only'))
dump('action/engine_b/gates_revised.json', revB)
print('action B', len(revB))
