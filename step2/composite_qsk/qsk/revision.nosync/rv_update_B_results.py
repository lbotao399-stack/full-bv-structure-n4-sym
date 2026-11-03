"""Revision (Part Q): record the revision in engine B's results.json (finding F8 restated, T5 placement, reruns, G-KT).
The first-pass tables are kept; the revision entries point to revision.nosync/out."""
import os, json
QSK = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/qsk'
RV = os.path.join(QSK, 'revision.nosync/out')
p = os.path.join(QSK, 'engine_b/results.json')
d = json.load(open(p))
F8_new = ('F8 (closed form, restated in the revision): 16 pi^2 [eps^0] = F_0 log(n_t/(n_t-1))/(n_t-1)! for a cycle with n_t edges whose '
          'eps^0 integrand F_0 is a constant (the leg sets with superficial degree D = 0, which carry no pole); 16 pi^2 [eps^-1] = '
          'F_1/(n_t^2 (n_t-1)(n_t-1)!) for a constant m = 1 integrand F_1 (superficial degree D = 2); verified exactly for n_t = 2..6.')
for i, f in enumerate(d.get('findings', [])):
    if isinstance(f, str) and f.startswith('F8 (closed form)'):
        d['findings'][i] = F8_new
if isinstance(d.get('closed_forms', {}).get('constant_numerator'), dict):
    d['closed_forms']['constant_numerator']['statement'] = F8_new


def load(n):
    try:
        return json.load(open(os.path.join(RV, n)))
    except Exception as e:
        return {'unavailable': str(e)}


d['revision_2026_10_30'] = {
    'T5_lamstar_placement': 'qb_job: T5 also placed with its lam* leg as the output (finding Q4); qb_job.T5_LAMS_OUTPUT = False reproduces the first pass',
    'n0_lamstar_tables': load('b_lamstar_n0.json'),
    'AB_n0_lamstar': load('ab_compare_n0_lamstar.json'),
    'AB_n1_lamstar': load('ab_compare_n1_lamstar.json'),
    'G_KT_Kb': (load('gkt_Kb.json') or {}).get('gate'),
    'note': 'the first-pass n = 0 tables are unchanged by the new placement where b_lamstar_n0.json says new_equals_old_code; n >= 1 needs the fixed code'}
json.dump(d, open(p, 'w'), indent=1)
print('updated')

# revision gates appended to engine B's gates.json
gp = os.path.join(QSK, 'engine_b/gates.json')
G = json.load(open(gp))
G['gates'] = [g for g in G['gates'] if not str(g.get('id', '')).startswith('REV ')]
bl = load('b_lamstar_n0.json'); n0 = load('ab_compare_n0_lamstar.json'); n1 = load('ab_compare_n1_lamstar.json'); gk = load('gkt_Kb.json')
new = [dict(id='REV G-T5lamstar n = 0 lam*-output leg sets unchanged by the T5 lam*-output placement', value='%d/%d' % (sum(r['new_equals_old_code'] for r in bl), len(bl)),
            target='all (the new configurations cancel at n = 0)', pass_or_fail='pass' if all(r['new_equals_old_code'] for r in bl) else 'FAIL', note='revision 2026-10-30', source='revision.nosync/rv_b_lamstar.py'),
       dict(id='REV G-AB n = 0 lam*-output rows, fixed route B vs route A', value='%d/%d' % (sum(r['agree'] for r in n0), len(n0)), target='all agree',
            pass_or_fail='pass' if all(r['agree'] for r in n0) else 'FAIL', note='revision 2026-10-30', source='revision.nosync/rv_ab.py'),
       dict(id='REV G-AB n = 1 lam*-output rows, fixed route B vs route A', value='%d/%d' % (sum(r['agree'] for r in n1), len(n1)), target='all agree (first-pass B: 8/10)',
            pass_or_fail='pass' if all(r['agree'] for r in n1) else 'FAIL', note='revision 2026-10-30', source='revision.nosync/rv_ab.py'),
       dict(id='REV G-KT along K_b', value=(gk.get('gate') or {}).get('n_configs'), target=(gk.get('gate') or {}).get('target'),
            pass_or_fail='pass' if (gk.get('gate') or {}).get('pass_') else 'FAIL', note='revision 2026-10-30', source='revision.nosync/rv_gkt.py')]
for g in G['gates']:
    if g.get('id') == 'G-closed (constant numerator)':
        g['note'] = ('revision: the eps^0 formula holds for a constant eps^0 integrand F_0 (D = 0 leg sets, no pole), the eps^-1 formula '
                     'for a constant m = 1 integrand F_1 (D = 2); each was evaluated with its own constant integrand')
G['gates'] += new
G['n_pass'] = sum(1 for g in G['gates'] if g.get('pass_or_fail') == 'pass')
G['n_fail'] = sum(1 for g in G['gates'] if g.get('pass_or_fail') != 'pass')
json.dump(G, open(gp, 'w'), indent=1)
print('gates', G['n_pass'], G['n_fail'])
