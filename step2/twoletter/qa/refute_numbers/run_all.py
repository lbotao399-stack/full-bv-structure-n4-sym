"""Referee (numbers) for R8: runs every check script in sequence (one process at a time, < 0.5 GB RSS each).
Outputs out_*.json in this folder; findings.json is written by hand from them (see findings.json 'checks').
Scripts:
  ref_core.py          referee engine R (own vertices from Ch1, own Wick/Koszul signs, Laplace-exact Schwinger integrals)
  ref_cliff.py         own 10d chiral Clifford realisation
  inv_eval.py          evaluator of engine A index-free strings (data only)
  t_basic.py           integration routine vs theory examples (Prop prop:tl-shape)
  t_bos_W0.py          bosonic W0 kernels vs engine A invariant forms (all components, symbolic momenta)
  t_deriv.py           derivative letters (W up to 6), eps^{-1}, eps^{-2}; k = 0 vs Lemma lem:tl-k0
  t_ghost.py           ghost / antifield pairs (odd pairs included) vs engine A forms
  t_su11.py            su(1|1) sector from scratch (fermion test, new)
  t_su11_engA.py       su(1|1) sector with engine A's kernel function
  t_ferm_engA.py       fermion kernels, component level, R (engine A Clifford) vs engine A, W up to 6
  t_twist2.py N        twist-two spectrum from own kernels, N <= 8
  t_p6b_su2.py         P6b matrix level, twist two without O(A_n), su(2) spectra L <= 8
  t_ward1.py           W1 for lambda, c, A_V (full colour)
  t_ward2.py, t_ward2b.py  W2 instances (run by neither engine)
  t_k3.py              k = 3 chains with derivatives vs engine A
  t_colour.py          planar colour reading with explicit gl(N) matrices
  t_bubble_example.py  theory example after Prop prop:tl-bubble
  t_census.py          independent census of log leg sets vs theory/census.json
  t_numeric.py         finite-eps quadrature residuals (no momentum expansion), slow (~10-20 min)
"""
import subprocess, sys, time
SCRIPTS = ['t_basic.py', 't_bos_W0.py', 't_deriv.py', 't_ghost.py', 't_su11.py', 't_su11_engA.py', 't_ferm_engA.py',
           ['t_twist2.py', '8'], 't_p6b_su2.py', 't_ward1.py', 't_ward2.py', 't_ward2b.py', 't_k3.py', 't_colour.py',
           't_bubble_example.py', 't_census.py']
if '--numeric' in sys.argv:
    SCRIPTS += ['t_numeric.py', 't_numeric_fast.py']
for s in SCRIPTS:
    cmd = [sys.executable] + (s if isinstance(s, list) else [s])
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    print('==', ' '.join(cmd[1:]), 'exit', r.returncode, '%.0fs' % (time.time() - t0))
    print('\n'.join(l for l in r.stdout.splitlines() if 'Warning' not in l)[-3000:])
