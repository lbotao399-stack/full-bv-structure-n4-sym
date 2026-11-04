"""Route B (engine B code path) for the action-side completion: imports engine B's verified modules read-only
(step2/strict/action/engine_b/src.nosync), template cache: engine B's cache read-only, new templates in
complete/out.nosync/cache."""
import os, sys, pickle, time, json
HERE = os.path.dirname(os.path.abspath(__file__))
COMPLETE = os.path.dirname(HERE)
BSRC = os.path.abspath(os.path.join(COMPLETE, '..', 'engine_b', 'src.nosync'))
if BSRC not in sys.path:
    sys.path.insert(0, BSRC)
import numpy as np
import sympy as sp
import eb_core as E
import eb_cycle as Y
import eb_gens as G
import eb_pos as P
import eb_qb as Q
import ce_core as C
MYCACHE = os.path.join(COMPLETE, 'out.nosync', 'cache')
os.makedirs(MYCACHE, exist_ok=True)
V = E.V
LOGF = os.path.join(COMPLETE, 'PROGRESS.txt')


def log(s, tag='routeB'):
    line = time.strftime('%Y-%m-%d %H:%M:%S') + ' ' + tag + ': ' + s
    print(line, flush=True)
    with open(LOGF, 'a') as f:
        f.write(line + '\n')


def template(name, terms, ext, extra_odd=0, gen_syms=None):
    fn = os.path.join(MYCACHE, f'tmpl_{name}.pkl')
    if os.path.exists(fn):
        with open(fn, 'rb') as f:
            return pickle.load(f)
    T = E.build_template(terms, ext, extra_odd=extra_odd, gen_syms=gen_syms, name=None)
    with open(fn, 'wb') as f:
        pickle.dump(T, f)
    return T


R = [np.array(C.Rm[M].tolist(), dtype=complex) for M in range(10)]     # Sigma -> Sigma^*
RB = [np.array(C.Rbm[M].tolist(), dtype=complex) for M in range(10)]   # Sigma^* -> Sigma
MP = np.array(C.Mp.tolist(), dtype=complex)                             # (s, s*) = s^T MP s*
LAMP = list(range(8)); LAMM = list(range(8, 16)); LAMS_M = list(range(8)); LAMS_P = list(range(8, 16))
