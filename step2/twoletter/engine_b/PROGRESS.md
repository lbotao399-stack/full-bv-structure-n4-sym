# Engine B (route B) progress log — step2/twoletter/engine_b

2026-10-29  start. Read DESIGN_STEP2 R8, theory/twoletter_theory.tex, engine_spec.json, predictions.json.
  Reuse (allowed, outside twoletter/engine_a): engine_b/vertices.py (action terms, plane-wave vertices, Hessian blocks),
  engine_b/loops.py (propagator numerators by inversion of the free Hessian), engine_b/realization.py (two-component realisation).
  Plan: tl_integrals.py (pinned-edge u-integrals, Feynman parametrisation, exact 1-d integrals, logs tracked),
        tl_engine.py (templates, chains, supertrace, singular extraction, colour words),
        gates / kernels drivers, run_all.py.
  kernels_W0: AV|c done (13.5 s), maxrss 104 MB
  kernels_W0: AW|AW done (0.9 s), maxrss 104 MB
  kernels_W0: AW|c done (0.0 s), maxrss 104 MB
  kernels_W0: c|c done (0.0 s), maxrss 104 MB
  kernels_W0: AV|AV done (14.4 s), maxrss 101 MB
  kernels_W0: AV|AW done (0.2 s), maxrss 101 MB
  kernels_W0: AV|lamP done (3.9 s), maxrss 101 MB
  kernels_W0: AV|lamM done (3.6 s), maxrss 101 MB
  kernels_W0: AV|c done (0.0 s), maxrss 101 MB
  kernels_W0: AV|AsV done (1.0 s), maxrss 101 MB
  kernels_W0: AW|AW done (0.4 s), maxrss 101 MB
  kernels_W0: AW|lamP done (4.8 s), maxrss 110 MB
  kernels_W0: AW|lamM done (4.8 s), maxrss 111 MB
  kernels_W0: AW|c done (0.0 s), maxrss 111 MB
  kernels_W0: AW|AsV done (1.3 s), maxrss 111 MB
  gates: G1-G9 + closed form pass; sl(2) twist two S<=8 passes only with gauge-letter kernels of covariant derivatives
  kernels_W0: lamP|lamP done (328.9 s), maxrss 1224 MB
  kernels_W0: lamP|lamM INCOMPLETE ['AV AW AW'] (393.6 s)
  kernels_W0: lamP|c done (25.6 s)
  kernels_W0: lamP|AsV done (333.3 s)
  kernels_W0: lamM|lamM done (436.0 s)
  kernels_W0: lamM|c done (25.6 s)
  kernels_W0: lamM|AsV done (328.4 s)
  kernels_W0: c|c done (25.6 s)
  kernels_W0: c|AsV done (135.6 s)
  kernels_W0: AsV|AsV done (252.0 s)
  kernels_W0: lamP|lamM done (314.8 s)
  kernels_deriv: AV|AV @ W0+1 done (17 s), peak RSS 142 MB
  kernels_deriv: AV|AW @ W0+1 done (17 s), peak RSS 168 MB
  kernels_deriv: AV|lamP @ W0+1 done (22 s), peak RSS 353 MB
  kernels_deriv: AV|lamM @ W0+1 done (23 s), peak RSS 345 MB
  kernels_deriv: AV|c @ W0+1 done (13 s), peak RSS 102 MB
  kernels_deriv: AV|AsV @ W0+1 done (16 s), peak RSS 144 MB
  kernels_deriv: AW|AW @ W0+1 done (17 s), peak RSS 213 MB
  kernels_deriv: AW|lamP @ W0+1 done (25 s), peak RSS 451 MB
  kernels_deriv: AW|lamM @ W0+1 done (24 s), peak RSS 457 MB
  kernels_deriv: AW|c @ W0+1 done (13 s), peak RSS 113 MB
  kernels_deriv: AW|AsV @ W0+1 done (17 s), peak RSS 181 MB
  kernels_deriv: lamP|c @ W0+1 done (13 s), peak RSS 114 MB
  kernels_deriv: lamP|AsV @ W0+1 done (29 s), peak RSS 401 MB
  kernels_deriv: lamM|c @ W0+1 done (13 s), peak RSS 128 MB
  kernels_deriv: lamM|AsV @ W0+1 done (30 s), peak RSS 177 MB
  kernels_deriv: c|c @ W0+1 done (13 s), peak RSS 143 MB
  kernels_deriv: c|AsV @ W0+1 done (15 s), peak RSS 146 MB
  kernels_deriv: AsV|AsV @ W0+1 done (20 s), peak RSS 179 MB
  kernels_deriv: AV|AV @ W0+2 done (40 s), peak RSS 247 MB
  kernels_deriv: AV|AW @ W0+2 done (42 s), peak RSS 341 MB
  kernels_deriv: AW|AW @ W0+2 done (48 s), peak RSS 371 MB
  kernels_deriv: lamP|lamP @ W0+1 done (51 s), peak RSS 338 MB
  kernels_deriv: lamP|lamM @ W0+1 done (79 s), peak RSS 445 MB
  kernels_deriv: lamM|lamM @ W0+1 done (51 s), peak RSS 339 MB
  run_all: tl_integrals.py                               exit 0  wall    76.6 s  cpu    76.5 s  peak child RSS     22 MB
  run_all: kernels_W0.py                                 exit 0  wall    12.5 s  cpu    12.4 s  peak child RSS    120 MB
  run_all: kernels_deriv.py                              exit 0  wall    12.6 s  cpu    12.5 s  peak child RSS    347 MB
  run_all: gates.py                                      exit 0  wall   389.2 s  cpu   385.9 s  peak child RSS    347 MB
  run_all: export.py                                     exit 0  wall     0.4 s  cpu     0.3 s  peak child RSS    369 MB
2026-10-29  final: run_all.py end to end (resumed W0/deriv parts), gates 45/45 PASS (gates.json). W0 tables: all 21 pair types,
  every census leg set at W0 (kernels_W0.json, parts_W0/). Derivative tables: W0+1 all pairs (<= 2 outputs), W0+2 for
  AV|AV, AV|AW, AW|AW (<= 2 outputs) (kernels_deriv.json). Spectra: SO(6) density u(2P-K); su(2) cyclic L<=8; twist two
  S<=8 (-8 h(S') u) only after including the gauge-letter pair kernels of the covariant derivatives. kernels.json exported.
