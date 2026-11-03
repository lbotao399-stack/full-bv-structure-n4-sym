"""Route B (engine B), PART Q: run everything (resumable) and assemble results.json, gates.json.
Each step runs in its own process (one heavy process at a time); a step whose output exists is skipped unless --force.
Prints wall time and peak RSS of every step."""
import os, sys, subprocess, time, resource, json

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'src.nosync')
OUT = os.path.join(HERE, 'out.nosync')
FORCE = '--force' in sys.argv
STEPS = [
    ('E0 conventions (commutators, arrows, composite line)', ['qb_e0.py', os.path.join(OUT, 'e0_conventions.json')], 'e0_conventions.json'),
    ('TQ0 quadratic closed directions', ['qb_tq0.py', os.path.join(OUT, 'tq0_directions.json')], 'tq0_directions.json'),
    ('G-reg along P_b (scalar pair)', ['qb_greg.py', 'AW'], 'greg_Pb_AW.json'),
    ('TQ2 kernels Upsilon_xi, n = 0, all leg sets', ['qb_tq2.py', '1500'], None),
    ('finite special vertices and G-pole', ['qb_fsp.py'], 'fsp_gpole.json'),
    ('numerical eps-extrapolation', ['qb_numeps.py'], 'numeps.json'),
    ('TQ4 bilinears and G-plane', ['qb_bil.py'], 'bilinears_xi.json'),
    ('TQ3 twist entries n = 0', ['qb_tq3.py'], 'tq3_twist_n0.json'),
    ('TQ5 intrinsic parts along eta (pair 0)', ['qb_etak.py', 'eta:0'], 'etak_eta_0.json'),
    ('TQ5 intrinsic parts along eta (pair 1)', ['qb_etak.py', 'eta:1'], 'etak_eta_1.json'),
    ('TQ5 intrinsic parts along eta (pair 2)', ['qb_etak.py', 'eta:2'], 'etak_eta_2.json'),
    ('TQ5 intrinsic parts along eta (pair 3)', ['qb_etak.py', 'eta:3'], 'etak_eta_3.json'),
    ('TQ5 intrinsic parts along eta (pair 4)', ['qb_etak.py', 'eta:4'], 'etak_eta_4.json'),
    ('TQ5 intrinsic parts along eta (pair 5)', ['qb_etak.py', 'eta:5'], 'etak_eta_5.json'),
    ('TQ6 intrinsic parts along K_b', ['qb_etak.py', 'K'], 'etak_K.json'),
    ('G-cubic, G-swap, G-chirality', ['qb_gmisc.py'], 'gmisc.json'),
]


def main():
    t00 = time.time()
    for name, cmd, out in STEPS:
        if out and os.path.exists(os.path.join(OUT, out)) and not FORCE:
            print(f'[skip] {name}: {out} exists')
            continue
        t0 = time.time()
        r = subprocess.run([sys.executable] + [cmd[0]] + cmd[1:], cwd=SRC)
        ru = resource.getrusage(resource.RUSAGE_CHILDREN)
        print(f'[done] {name}: exit {r.returncode}, {time.time() - t0:.1f} s, peak RSS of children so far {ru.ru_maxrss / 2 ** 20:.0f} MB', flush=True)
    r = subprocess.run([sys.executable, 'qb_assemble.py'], cwd=SRC)
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    print(f'total {time.time() - t00:.1f} s; peak RSS of any child {ru.ru_maxrss / 2 ** 20:.0f} MB; user CPU {ru.ru_utime:.0f} s')


if __name__ == '__main__':
    main()
