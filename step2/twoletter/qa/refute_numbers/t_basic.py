import sympy as sp
from ref_core import *
p = sp.symbols('p0:4')
l2 = sum(x*x for x in L)
# theory example Prop prop:tl-shape: E=2, N=|l|^4, edges l and l-p: per-edge |p|^2/8 eps^-1 + 1/4 eps^-2, log 0
print(singular_parts(l2**2, [[0]*4, [-x for x in p]]))
# triangle N=|l|^4: eps^-1 coefficient 1/3 (per-edge)
q = sp.symbols('q0:4')
print(singular_parts(l2**2, [[0]*4, list(p), list(q)]))
# bubble constant numerator: log coefficient -1
print(singular_parts(sp.Integer(1), [[0]*4, list(p)]))
# J sanity: E=1
print(J_orthant([0],3), J_orthant([2],5), J_orthant([0,0],3), J_orthant([1,0],4))
