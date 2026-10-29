"""Engine B (route B), end to end:
  1. tl_integrals self-test (pinned-edge integrals against quadrature)
  2. kernels_W0.py   : all 21 pair types at W = W0 (resumable; finished pairs are skipped)
  2b. kernels_deriv.py: derivative tables W0+1 (all pairs, <= 2 outputs), W0+2 (scalar/gauge pairs)
  3. gates.py        : G1-G9, closed forms, numerical eps residuals (gates.json, two_site.json)
  4. export.py       : kernels.json
Each step runs in its own process, one at a time; the resources of each are printed.
Usage: python3 run_all.py [--skip-W0] [--skip-deriv] [--sections colour bilinears ...]"""
import os, sys, subprocess, time, resource, json
HERE = os.path.dirname(os.path.abspath(__file__))


def step(args, log):
    t0 = time.time()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with open(os.path.join(HERE, log), 'w') as fh:
        r = subprocess.run([sys.executable] + args, cwd=HERE, stdout=fh, stderr=subprocess.STDOUT)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    peak = after.ru_maxrss / 2**20
    msg = f'{" ".join(args):45s} exit {r.returncode}  wall {time.time() - t0:7.1f} s  cpu {cpu:7.1f} s  peak child RSS {peak:6.0f} MB'
    print(msg, flush=True)
    with open(os.path.join(HERE, 'PROGRESS.md'), 'a') as f:
        f.write(f'  run_all: {msg}\n')
    return r.returncode


def main():
    argv = sys.argv[1:]
    secs = []
    if '--sections' in argv:
        secs = argv[argv.index('--sections') + 1:]
    t0 = time.time()
    step(['tl_integrals.py'], 'log_integrals.txt')
    if '--skip-W0' not in argv:
        step(['kernels_W0.py'], 'log_W0.txt')
    if '--skip-deriv' not in argv:
        step(['kernels_deriv.py'], 'log_deriv.txt')
    os.environ.setdefault('SL2_SMAX', '8')
    step(['gates.py'] + secs, 'log_gates.txt')
    step(['export.py'], 'log_export.txt')
    g = json.load(open(os.path.join(HERE, 'gates.json')))
    npass = sum(1 for v in g.values() if v['pass_'])
    print(f'gates: {npass}/{len(g)} pass; total wall {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
