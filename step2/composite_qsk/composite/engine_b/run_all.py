"""Route B (engine B), PART C: run everything, one process at a time, printing wall/cpu/peak RSS per step.
Usage: python3 run_all.py [--force]   (without --force, steps whose output exists in out.nosync/ are skipped)."""
import os, sys, subprocess, time, resource, json
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'src.nosync')
OUT = os.path.join(HERE, 'out.nosync')
STEPS = [
    (['ce_gates_basic.py'], 'basic_gates.json'),       # core identities, exact integrals, GC3 bilinear, Lemma (ii), colour
    (['ce_so6.py'], 'so6_AWAW.json'),                  # F0: general SO(6) scalar pair, planar density
    (['ce_su11.py'], 'su11.json'),                     # F2: su(1|1) L <= 6
    (['ce_sl2.py', 'Z', '4'], 'sl2_Z_S4.json'),         # F1: twist two S' <= 4
    (['ce_sl2.py', 'psi', '5'], 'sl2_psi_S5.json'),     # F2: fermion sl(2) S' <= 5
    (['ce_gluon.py', '2'], 'gluon_S2.json'),            # F3: gluon sl(2) S' <= 2
    (['ce_mixed.py', '4'], 'mixed_S4.json'),            # F3: (F11, psi) S' <= 4
    (['ce_psiZF.py', '2'], 'psiZF_S2.json'),            # F2: psi psi' with Z F11, S' <= 2
    (['ce_numeps.py'], 'numeps.json'),                  # GC5b numerical eps extrapolation
    (['ce_assemble.py'], None),                         # gates.json, results.json
]


def step(args, log):
    t0 = time.time()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with open(os.path.join(OUT, log), 'w') as fh:
        r = subprocess.run([sys.executable] + args, cwd=SRC, stdout=fh, stderr=subprocess.STDOUT)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    peak = after.ru_maxrss / 2 ** 20
    msg = f'{" ".join(args):28s} exit {r.returncode}  wall {time.time() - t0:7.1f} s  cpu {cpu:7.1f} s  peak child RSS {peak:6.0f} MB'
    print(msg, flush=True)
    with open(os.path.join(HERE, 'PROGRESS.txt'), 'a') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M ') + 'run_all: ' + msg + '\n')
    return r.returncode


def main():
    force = '--force' in sys.argv
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    for args, outf in STEPS:
        if outf and not force and os.path.exists(os.path.join(OUT, outf)):
            print(f'{" ".join(args):28s} skipped (out.nosync/{outf} exists; --force to recompute)', flush=True)
            continue
        log = 'log_' + '_'.join(a.replace('.py', '') for a in args) + '.txt'
        step(args, log)
    g = json.load(open(os.path.join(HERE, 'gates.json')))
    npass = sum(1 for v in g.values() if v['pass_'] is True)
    nfail = sum(1 for v in g.values() if v['pass_'] is False)
    nrun = sum(1 for v in g.values() if v['pass_'] is None)
    print(f'gates: {npass} pass, {nfail} fail, {nrun} not run; total wall {time.time() - t0:.0f} s; '
          f'peak RSS of this driver {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20:.0f} MB')


if __name__ == '__main__':
    main()
