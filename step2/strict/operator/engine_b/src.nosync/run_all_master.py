"""Operator engine B (route B), strict Costello route: run everything (resumable).

  python3 run_all.py [--skip-queue] [--force]   (--force re-runs the steps whose output exists)

Steps (each resumable; outputs in out.nosync/, final results.json, gates.json here):
  0. census extraction (out.nosync/census_n01.json) from step2/composite_qsk/qsk/theory/src.nosync/census_qsk.json
  1. OP0: operators on plane waves for every generator (src.nosync/ob_e0.py)
  2. production queue of the pair rows (src.nosync/ob_queue.py --worker p1 p2 p5 p4 p3; one worker process at a time,
     restarted when its peak RSS passes 600 MB; heavy rows retried with fixed numeric output polarisations, then with
     numeric letters)
  3. tables n = 0, 1 with G-noLog, G-D1pair (both readings), reduced G-pole (src.nosync/ob_tables.py)
  4. F0 (ob_f0.py, ob_xif0.py), F1 (ob_fam.py), second evaluation (ob_conj.py), K_b decomposition (ob_pos.py),
     OP7 one-sided relations (ob_r7.py), OP8 finite T (ob_kt.py); F3, F4 (ob_sl2.py Z / psi, null-mode defect cycles),
     F2 (ob_su11.py)
  5. assembly (ob_assemble.py) -> results.json, gates.json (copies in out.nosync/final/)
Memory: in-process cap 1250 MB, RSS watchdog 1400 MB (polled every 0.3 s) in every engine process (env OB_MEMCAP_MB,
OB_WATCHDOG_MB); one engine process at a time.  Prints the wall time and the peak RSS of the child processes.
The master copy of this script is src.nosync/run_all_master.py (iCloud may evict run_all.py)."""
import os, sys, json, time, subprocess, resource
HERE = os.path.dirname(os.path.abspath(__file__))
if os.path.basename(HERE) == 'src.nosync':
    HERE = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'src.nosync')
OUT = os.path.join(HERE, 'out.nosync')
PROG = os.path.join(HERE, 'PROGRESS.txt')
SL2_ZMAX = os.environ.get('OB_SL2_ZMAX', '6')
SL2_PSIMAX = os.environ.get('OB_SL2_PSIMAX', '5')


def log(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S') + ' run_all: ' + msg
    print(line, flush=True)
    with open(PROG, 'a') as f:
        f.write(line + '\n')


def run(args, desc, timeout=None):
    t0 = time.time()
    r = subprocess.run([sys.executable] + args, cwd=SRC, capture_output=True, text=True, timeout=timeout)
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 2 ** 20
    log(f'{desc}: exit {r.returncode}, {round(time.time() - t0)} s, children peak RSS {peak:.0f} MB')
    if r.returncode != 0:
        log('  stderr: ' + (r.stderr.strip().splitlines() or [''])[-1][:300])
    return r


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    os.environ.setdefault('OB_MEMCAP_MB', '1250')
    os.environ.setdefault('OB_WATCHDOG_MB', '1400')
    if not os.path.exists(os.path.join(OUT, 'census_n01.json')):
        src = os.path.join(HERE, '..', '..', '..', 'composite_qsk', 'qsk', 'theory', 'src.nosync', 'census_qsk.json')
        d = json.load(open(src))
        out = {'source': src, 'description': d['description'][:4000], 'generators': {}}
        for g in ['even_P', 'even_0', 'xi', 'xibar', 'eta_intrinsic', 'etabar_intrinsic', 'K_intrinsic']:
            gg = {}
            for k, e in d['generators'][g]['pairs'].items():
                gg[k] = {'W0': e.get('W0'), 'pair_degree': e.get('pair_degree'), 'field_pair': e.get('field_pair'), 'by_n': {}}
                for nd in ['0', '1']:
                    b = e['by_letter_derivatives'].get(nd)
                    if b is not None:
                        gg[k]['by_n'][nd] = {'W': b['W'], 'defect_leg_sets': b['defect_leg_sets'], 'counterterm_leg_sets': b.get('counterterm_leg_sets', {})}
            out['generators'][g] = gg
        json.dump(out, open(os.path.join(OUT, 'census_n01.json'), 'w'))
        log('census extracted')
    if not os.path.exists(os.path.join(OUT, 'op0.json')):
        run(['ob_e0.py'], 'OP0')
    if '--skip-queue' not in sys.argv:
        for ph, sec in (('p1', 1500), ('p2', 1500), ('p5', 1500), ('p4', 2400), ('p3', 2400)):
            run(['ob_queue.py', '--worker', ph, str(sec)], f'queue {ph}')
    for n in (0, 1):
        for attempt in range(6):          # resumable: restarted after the RSS watchdog, finished tables reused
            r = run(['ob_tables.py', str(n)], f'tables n={n} (attempt {attempt + 1})')
            if r.returncode == 0:
                break
    outs = {'ob_f0.py': 'f0.json', 'ob_xif0.py': 'f0_odd.json', 'ob_fam.py': 'fam_f1.json', 'ob_conj.py': 'conj.json',
            'ob_pos.py': 'pos.json', 'ob_r7.py': 'r7.json', 'ob_kt.py': 'op8_kt.json'}
    for script in ('ob_f0.py', 'ob_xif0.py', 'ob_fam.py', 'ob_conj.py', 'ob_pos.py', 'ob_r7.py', 'ob_kt.py'):
        if '--force' in sys.argv or not os.path.exists(os.path.join(OUT, outs[script])):
            run([script], script)
    # strict route on the sl(2) and su(1|1) families (defect cycles in the null mode; rows cached in out.nosync/sl2)
    for args in (['ob_sl2.py', 'Z', SL2_ZMAX], ['ob_sl2.py', 'psi', SL2_PSIMAX], ['ob_su11.py', '6']):
        o = os.path.join(OUT, {'Z': f'sl2_Z_S{SL2_ZMAX}.json', 'psi': f'sl2_psi_S{SL2_PSIMAX}.json', '6': 'su11_L6.json'}[args[1]])
        if '--force' in sys.argv or not os.path.exists(o):
            run(args, ' '.join(args))
    run(['ob_assemble.py'], 'assembly')
    log(f'finished in {round(time.time() - t0)} s')


if __name__ == '__main__':
    main()
