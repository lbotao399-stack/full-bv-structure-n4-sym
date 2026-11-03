"""QA bisection, step 1 (engine A, read-only): A's values of the (lamP, lamP) -> (AW, lamM) defect cycles on mixed
data (twisted / generic xi, letters, outputs), real integer momenta.  Writes qa/out.nosync/bisect_A.json."""
import os, sys, json
import numpy as np
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_twist as TW, qa_core as QC, qa_util as U, ea_core as E
T = TW.Twist()
G5S = np.asarray(E.G5S, dtype=complex)
rng = np.random.default_rng(17)
def rv(n):
    return rng.integers(-2, 3, n) + 1j * rng.integers(-2, 3, n)
def chir(z, s):
    return (np.eye(16) + s * G5S) @ z
resA = json.load(open(os.path.join(QSK, 'engine_a/results.json')))
xi_gen = np.array([complex(x) for x in resA['tasks']['TQ2']['xi']])
cases = []
q = [rng.integers(-3, 4, 4).astype(complex) for _ in range(2)]
gen_la, gen_lb = chir(rv(16), 1), chir(rv(16), 1)
gen_o0, gen_o1 = rv(6), chir(rv(16), -1)
spec = {
    'all twisted': (T.xi, T.beta_z[0], T.beta_z[1], T.gam_out[2], T.dc_out[0]),
    'twisted xi, generic rest': (T.xi, gen_la, gen_lb, gen_o0, gen_o1),
    'generic xi, twisted letters+outputs': (xi_gen, T.beta_z[0], T.beta_z[1], T.gam_out[2], T.dc_out[0]),
    'twisted xi+letters, generic outputs': (T.xi, T.beta_z[0], T.beta_z[1], gen_o0, gen_o1),
    'twisted xi+outputs, generic letters': (T.xi, gen_la, gen_lb, T.gam_out[2], T.dc_out[0]),
    'all generic': (xi_gen, gen_la, gen_lb, gen_o0, gen_o1),
    'scaled-8 all twisted': (T.xi, 8 * T.beta_z[0], 8 * T.beta_z[1], 2 * T.gam_out[2], T.dc_out[0]),
}
out = []
for name, (xi, la, lb, o0, o1) in spec.items():
    gx = QC.odd_gen(xi, name='xi_%s' % name.replace(' ', '_'))
    legs = [('AW', o0, q[0]), ('lamM', o1, q[1])]
    v = QC.pair_value(('lamP', la, None), ('lamP', lb, None), legs, gx, what=('fin',))
    c = lambda a: [str(complex(x)) for x in np.asarray(a).ravel()]
    rec = {'case': name, 'xi': c(xi), 'la': c(la), 'lb': c(lb), 'o0': c(o0), 'o1': c(o1), 'q': [c(x) for x in q], 'A': U.vstr(v['fin'])}
    print(name, 'A =', rec['A'], flush=True)
    out.append(rec)
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/bisect_A.json'), 'w'), indent=1)
