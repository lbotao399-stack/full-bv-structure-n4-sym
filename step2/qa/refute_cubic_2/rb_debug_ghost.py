import sys, os, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
import rb_loop as RL
p = np.array([1., 2., -1., 3.]); X0 = np.array([1., 0., 2.]); X1 = np.array([0., 3., 1.])
a = np.zeros((10, 3)); a[:4, :] = np.outer([2., -1., 1., 1.], X1)   # rank one a x X1
ext = [[(X0, p)], [(a.reshape(-1), -p)]]
for g in RL.enumerate_cycles([('c',), ('A',)]):
    vt = [v[0] for v in g['verts']]
    if vt != ['AcAst', 'AcAst']:
        continue
    for parts in RL.a_parts(g):
        fn = lambda l: RL.graph_N(g, parts, ext, l)
        for l in [np.array([1., 0, 0, 0]), np.array([0., 1, 2, 0])]:
            # hand formula: N = -2i (X0,X1) l^2 (l - p, a_vec) for X on the edge with momentum l
            av = np.array([2., -1., 1., 1.])
            print(g['kinds'], [ (v[4], v[5]) for v in g['verts']], parts, fn(l),
                  -2j * (X0 @ X1) * (l @ l) * ((l - p) @ av), -2j * (X0 @ X1) * ((l-p)@(l-p)) * (l @ av))
