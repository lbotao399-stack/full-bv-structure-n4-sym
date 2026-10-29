"""Engine A, revision 2026-10-30: extended family runs (QA finding 5: families that rested on route B alone).
fermion sl(2) S' = 4, 5; gluon sl(2) S' = 2; twist two S' = 5, 6.  Same code and pass criteria as run_all.py stages 3, 4, 6,
vertex-vertex sign -1 (notes.txt sec. 2b).  Single process; results written incrementally to results_ext.json, gates_ext.json.
Usage: python3 ext_families.py [job ...]   jobs: f4 g2 t5 f5 t6 (default: all, in this order)"""
import os, sys, json, time, resource
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'src.nosync'))
import numpy as np
import sympy as sp
import cq_core as Q, cq_gval as GV, cq_twist as TW, cq_fsl2 as FS, cq_gluon as GL
import ea_spectra as ES
Q.RULE['koszul'] = 'R2'; GV.INTR['on'] = True
OUT = os.path.join(HERE, 'out.nosync')
T0 = time.time()


def res():
    return {'wall_s': round(time.time() - T0, 1), 'max_rss_MB': round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20, 1)}


def log(*a):
    line = '[%8.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(line, flush=True)
    with open(os.path.join(OUT, 'ext_families_log.txt'), 'a') as f:
        f.write(line + '\n')


def progress(msg):
    with open(os.path.join(HERE, 'PROGRESS.txt'), 'a') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S') + ' ext_families: ' + msg + '\n')


R = json.load(open(os.path.join(HERE, 'results_ext.json'))) if os.path.exists(os.path.join(HERE, 'results_ext.json')) else {}
G = json.load(open(os.path.join(HERE, 'gates_ext.json'))) if os.path.exists(os.path.join(HERE, 'gates_ext.json')) else {}


def dump():
    for d in (HERE, OUT):
        json.dump(R, open(os.path.join(d, 'results_ext.json'), 'w'), indent=1, default=str)
        json.dump(G, open(os.path.join(d, 'gates_ext.json'), 'w'), indent=1, default=str)


def fermion(N):
    r = FS.run(N, n=(1, -1j, 0, 0), log=log)
    g10 = sp.Matrix(sp.sympify(r['M0']['matrix'])) == -sp.Rational(2, 3) * (4 + N) * sp.eye(len(r['M0']['matrix']))
    R['F2_fermion_sl2_N%d' % N] = r
    G['F2_fermion_sl2_N%d' % N] = {'P7c_sumC0_eq_-2/3(4+N)': bool(g10), 'P7b_n': bool(r['Pn']['pass']), 'D1_eq_-8uE_H': bool(r['D1']['pass']),
                                   'pass': bool(g10 and r['Pn']['pass'] and r['D1']['pass']), 'resources': res()}


def gluon(N):
    r = GL.run(N, log=log)
    g10 = sp.Matrix(sp.sympify(r['M0']['letter_matrix'])) == -sp.Rational(2, 3) * (6 + N) * sp.eye(len(r['M0']['letter_matrix']))
    vv0 = all(sp.sympify(x) == 0 for row in r['M0']['vv_matrix'] for x in row)
    gd = r['D1']['kappa=-1']['pass']
    gw = all(v['letter_part_zero'] and v['vv_part_zero'] for v in r['ward'].values())
    R['F3_gluon_N%d' % N] = r
    G['F3_gluon_N%d' % N] = {'P7c_sumC0_eq_-2/3(6+N)': bool(g10), 'vv_part_zero_on_shell': bool(vv0), 'P7b_n_letter_zero': bool(r['Pn']['letter_zero']),
                             'D1_eq_-8uE_H_vvsign_-1': bool(gd), 'ward': bool(gw), 'pass': bool(g10 and r['Pn']['letter_zero'] and gd and gw),
                             'resources': res()}


def twist(N):
    r = TW.run(N, kinds=('D0', 'Pn'), log=log)
    M0 = sp.Matrix([[sp.sympify(x) for x in row] for row in r['D0']['matrix_by_logbasis']['1']])
    logs0 = all(all(sp.sympify(x) == 0 for row in r['D0']['matrix_by_logbasis'][b] for x in row) for b in ['log2', 'log3', 'log5', 'log7'])
    g10 = logs0 and M0 == -sp.Rational(2, 3) * (4 + N) * sp.eye(M0.shape[0])
    Pn = sp.Matrix([[sp.sympify(x) for x in row] for row in r['Pn']['matrix_by_logbasis']['1']])
    a1, a2 = sp.symbols('a1 a2')
    g9 = True
    for ik, k in enumerate(range(N // 2 + 1)):
        img = sum(Pn[a, ik] * a1 ** a * a2 ** (N + 1 - a) for a in range(N + 2))
        Bk = 4 * (a1 ** k * a2 ** (N - k) + a1 ** (N - k) * a2 ** k) if 2 * k != N else 8 * a1 ** k * a2 ** k
        g9 &= sp.expand(img - sp.Rational(-1, 3) * (a1 + a2) * Bk) == 0
    lc = ES.twist_two(N)
    Lm = sp.Matrix([[sp.sympify(x) for x in row] for row in lc['matrix_units']])
    D1 = 4 * Lm + M0 + sp.Rational(2, 3) * (4 + N) * sp.eye(M0.shape[0])
    ev = sorted([sp.nsimplify(e) for e, m in D1.eigenvals().items() for _ in range(m)], key=float)
    target = sorted([-32 * sum(sp.Rational(1, j) for j in range(1, S + 1)) for S in range(0, N + 1, 2)], key=float)
    r['D1'] = {'matrix': [[str(D1[i, j]) for j in range(D1.shape[1])] for i in range(D1.shape[0])], 'eigenvalues': [str(x) for x in ev],
               'target_-32h(S)': [str(x) for x in target], 'pass': ev == target}
    R['F1_twist_two_N%d' % N] = r
    G['F1_twist_two_N%d' % N] = {'P7c_sumC0_eq_-2/3(4+N)': bool(g10), 'P7b_n': bool(g9), 'D1_eq_-32h(S)': bool(ev == target),
                                 'pass': bool(g10 and g9 and ev == target), 'resources': res()}


JOBS = {'f4': (fermion, 4), 'g2': (gluon, 2), 't5': (twist, 5), 'f5': (fermion, 5), 't6': (twist, 6)}
if __name__ == '__main__':
    jobs = sys.argv[1:] or ['f4', 'g2', 't5', 'f5', 't6']
    progress('started %s' % jobs)
    for j in jobs:
        f, N = JOBS[j]
        log('job', j, 'start')
        progress('job %s start' % j)
        f(N)
        dump()
        log('job', j, 'done', G[[k for k in G if k.endswith('N%d' % N) and k.startswith({'f': 'F2', 'g': 'F3', 't': 'F1'}[j[0]])][0]], res())
        progress('job %s done: %s' % (j, res()))
    progress('finished')
