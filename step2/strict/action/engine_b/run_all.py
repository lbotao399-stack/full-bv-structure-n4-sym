#!/usr/bin/env python3
"""Engine B (route B), action side, strict Costello route: single entry point.

Runs every step as a separate process (one heavy process at a time), skips a step whose output exists (resumable),
enforces the memory cap (RSS polled every 0.5 s, the step is killed above 1400 MB) and prints wall time and peak RSS.
Outputs: out.nosync/*.json (raw), results.json, gates.json (assembled), PROGRESS.txt (log)."""
import os, subprocess, sys, time

ROOT = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(ROOT, 'src.nosync')
OUT = os.path.join(ROOT, 'out.nosync')
CAP_MB = 1400
STEPS = [
    ('G-reg0 (arity zero, two and three legs)', ['eb_greg0.py'], ['G_reg0.json']),
    ('T1 Ob1(Delta) two legs', ['eb_t1.py'], ['T1.json']),
    ('T4 Ob1(xi) two legs', ['eb_t4.py'], ['T4.json']),
    ('G-sign, Q images, G-WZ1 along xi', ['eb_gsign.py'], ['G_sign.json']),
    ('B1 raw cycles (Delta,P_b)', ['eb_b1.py'], ['B1_raw_c_cs.json', 'B1_raw_AV_AsV.json', 'B1_raw_Ap_Asp.json',
                                                 'B1_raw_lam+_lams+.json', 'B1_raw_lam-_lams-.json']),
    ('B1/S1 cohomology weight 5 and theta', ['eb_b1fit.py'], ['B1_fit.json']),
    ('T3 R_r and M_omega', ['eb_t3.py'], ['T3.json']),
    ('S1 two-leg cohomology tables', ['eb_s1.py'], ['S1.json']),
    ('G-cyc, G-log (numerical eps)', ['eb_gnum.py'], ['G_num.json']),
]


def log(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S') + ' run_all: ' + msg
    print(line, flush=True)
    with open(os.path.join(ROOT, 'PROGRESS.txt'), 'a') as f:
        f.write(line + '\n')


def rss_mb(pid):
    try:
        out = subprocess.run(['ps', '-o', 'rss=', '-p', str(pid)], capture_output=True, text=True).stdout.strip()
        return int(out) / 1024 if out else 0.0
    except Exception:
        return 0.0


def run_step(name, script, outs):
    if all(os.path.exists(os.path.join(OUT, o)) for o in outs):
        log(f'skip (outputs exist): {name}')
        return True
    if not os.path.exists(os.path.join(SRC, script[0])):
        log(f'missing script for {name}: {script[0]}')
        return False
    t0 = time.time()
    logf = open(os.path.join(OUT, script[0].replace('.py', '.log')), 'a')
    p = subprocess.Popen([sys.executable, '-u'] + script, cwd=SRC, stdout=logf, stderr=subprocess.STDOUT)
    peak = 0.0
    while p.poll() is None:
        m = rss_mb(p.pid)
        peak = max(peak, m)
        if m > CAP_MB:
            p.kill()
            log(f'KILLED {name}: RSS {m:.0f} MB > {CAP_MB} MB')
            return False
        time.sleep(0.5)
    log(f'{name}: exit {p.returncode}, wall {time.time() - t0:.0f} s, peak RSS {peak:.0f} MB')
    return p.returncode == 0


def main():
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    ok = True
    for name, script, outs in STEPS:
        ok = run_step(name, script, outs) and ok
    r = subprocess.run([sys.executable, 'eb_assemble.py'], cwd=SRC, capture_output=True, text=True)
    log('assemble: ' + (r.stdout.strip() or r.stderr.strip()[-300:]))
    import resource
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    log(f'total wall {time.time() - t0:.0f} s, children user CPU {ru.ru_utime:.0f} s, max child RSS {ru.ru_maxrss / 2 ** 20:.0f} MB')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
