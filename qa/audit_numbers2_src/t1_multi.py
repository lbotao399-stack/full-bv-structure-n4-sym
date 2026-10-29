# Driver: runs t1.py at many random integer configurations (entries in [-RANGE, RANGE]) and checks that
# every per-class ratio (log eps coefficient / tree vertex, in units hbar tr_g(ad ad)/((X,Y)16pi^2)) is
# configuration independent (Schwartz-Zippel: a residual structure of degree d vanishes at a random point
# with probability <= d/(2 RANGE + 1)).
import json, subprocess, sys
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
RANGE = int(sys.argv[2]) if len(sys.argv) > 2 else 9
allres = []
for seed in range(100, 100 + N):
    out = subprocess.run([sys.executable, 't1.py', str(seed), str(RANGE)], capture_output=True, text=True, check=True).stdout
    allres.append(json.loads(out))
summary = {}
for block in ['AVAVAV', 'AV_two_point', 'lambda_two_point', 'lam_lam_AWperp']:
    summary[block] = {}
    for k in allres[0][block]:
        vals = [json.dumps(r[block][k]) for r in allres]
        if 'longitudinal log per loop' in k:
            continue
        summary[block][k] = {'values': sorted(set(vals)), 'configuration_independent': len(set(vals)) == 1}
summary['n_configurations'] = N
summary['entry_range'] = [-RANGE, RANGE]
summary['clifford_ok'] = all(r['clifford_ok'] for r in allres)
print(json.dumps(summary, indent=1))
