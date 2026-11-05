"""closed forms of the s_x Psi vertices and of [x., QGF] in frame Psi."""
import sympy as sp
import bvjet, models as Mo, ops as O, freecx as Fx, gfsolve as GS
from test_equiv import gens_even, odd, b
X = bvjet.X; p = Fx.p
M = Mo.build_model(True); F = Mo.Fields(M, 'Psi')
fieldnames = set(['c','cb']+['A%d'%m for m in range(4)]+['P%d'%I for I in range(6)]+['L%d'%a for a in range(16)])
def field_part(P):
    out={}
    for (e,o),c in P.items():
        js=[j for j,pw in e]+list(o)
        if all(M.fields[j[0]].name in fieldnames for j in js):
            out[(e,o)]=c
    return out
res={}
# (i) even generators: field-field part = -1/2 int (cbar, iota_{grad div v} A|_V)
for k in ('P','M','D','K'):
    v = gens_even[k]
    ff = field_part(Mo.S1_even(F, v=v))
    gd = [sp.diff(Mo.div(v), t) for t in X]
    pred = {}
    for m in range(4):
        if gd[m]!=0:
            pred = M.plus(pred, M.mul(M.const(-sp.Rational(1,2)*gd[m]), M.mul(F.cb, F.A[m])))
    res['s_vPsi(%s) = -1/2 (cbar, iota_grad(div v) A_V)'%k] = M.is_zero_functional(M.plus(ff, M.scale(pred,-1)))
# K_b explicitly 4 (cbar, iota_b A)
pred = {}
for m in range(4):
    pred = M.plus(pred, M.mul(M.const(4*b[m]), M.mul(F.cb, F.A[m])))
res['s_KPsi = 4 (cbar, iota_b A_V)'] = M.is_zero_functional(M.plus(field_part(Mo.S1_even(F, v=gens_even['K'])), M.scale(pred,-1)))
# (ii) odd generators: field-field part = (Gamma(eps,lambda)|_V, d cbar)
for k in ('xi1','eta1'):
    xi, eta = odd[k]
    ff = field_part(Mo.S1_odd(F, xi, eta))
    eps = Mo.killing(xi, eta)
    G = Mo.Gamma_lam(F, eps)
    pred = {}
    for m in range(4):
        pred = M.plus(pred, M.mul(G[m], F.d(F.cb, m)))
    res['s_epsPsi(%s) = (Gamma(eps,lambda)|_V, d cbar)'%k] = M.is_zero_functional(M.plus(ff, M.scale(pred,-1)))
# (iii) [v., QGF] = -1/2 div v QGF + 1/4 rho(d div v) on lambda^* + [V_v, QGF]
Gop = O.op_from_momentum(M, GS.QGF_Psi_red(p), p, 1)
e_all = O.generic_section(M, range(len(M.fields)))
Rb = Mo.Rm_blocks
for k in ('D','K'):
    v = gens_even[k]; dv = Mo.div(v); gd=[sp.diff(dv,t) for t in X]
    xo = O.op_from_functional(M, Mo.S1_even(F, v=v))
    lhs = O.gcomm(xo, Gop, e_all)
    Vfun = {}
    for m in range(4):
        if gd[m]!=0:
            Vfun = M.plus(Vfun, M.mul(M.const(-sp.Rational(1,2)*gd[m]), M.mul(F.cb, F.A[m])))
    rhs = O.sec_scale(Gop(e_all), -sp.Rational(1,2)*dv)
    if Vfun:
        Vop = O.op_from_functional(M, Vfun)
        rhs = O.sec_sub(rhs, O.gcomm(Vop, Gop, e_all), -1)
    Lidx = [i for i,f in enumerate(M.fields) if f.name.startswith('L') and not f.name.startswith('Ls')]
    Lsidx = [i for i,f in enumerate(M.fields) if f.name.startswith('Ls')]
    extra = {}
    for a_,ia in enumerate(Lidx):
        t = 0
        for m in range(4):
            if gd[m]!=0:
                for b_,ib in enumerate(Lsidx):
                    if Rb[m][1][a_,b_]!=0:
                        t += sp.Rational(1,4)*gd[m]*Rb[m][1][a_,b_]*e_all[ib]
        if t!=0: extra[ia]=sp.expand(t)
    rhs = O.sec_sub(rhs, extra, -1)
    res['[v.,QGF] = -1/2 div v QGF + 1/4 rho(d div v)|_lam* + [V_v,QGF]  (%s)'%k] = O.is_zero_section(O.sec_sub(lhs, rhs))
# explicit V_K operator components
Vfun = {}
for m in range(4):
    Vfun = M.plus(Vfun, M.mul(M.const(4*b[m]), M.mul(F.cb, F.A[m])))
VK = O.op_from_functional(M, Vfun)
img = VK(e_all)
res['V_K_components'] = {M.fields[a].name: str(sp.factor(val)) for a,val in img.items()}
img2 = O.gcomm(VK, Gop, e_all)
res['[V_K,QGF]_components'] = {M.fields[a].name: str(sp.factor(val)) for a,val in img2.items()}
for k,v in res.items(): print(k, ':', v)
import json; json.dump(res, open('vertices_results.json','w'), indent=1)
