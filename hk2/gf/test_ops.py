"""
Operator-level checks in frame Psi: the linear action x. of every generator (from the quadratic part of
S^(1)_Psi(x)), its graded commutator with QGF, the composite kernel QGF x. QGF, and the s_x Psi vertex V_x.
"""
import sympy as sp, json, time
import bvjet, models as Mo, ops as O, freecx as Fx, gfsolve as GS
from test_equiv import gens_even, odd, r

X = bvjet.X
p = Fx.p


def slot_indices(M, prefixes):
    return [i for i, f in enumerate(M.fields) if any(f.name == q or (f.name.startswith(q) and f.name[len(q):].isdigit()) for q in prefixes)]


def run(verbose=True):
    t0 = time.time()
    res = {}
    M = Mo.build_model(True)
    F = Mo.Fields(M, 'Psi')
    G = O.op_from_momentum(M, GS.QGF_Psi_red(p), p, 1)
    Qop = O.op_from_functional(M, Mo.S0(F))
    Qmom = O.op_from_momentum(M, GS.Q_Psi_red(p, -sp.Rational(1, 2), -1, 1), p, 1)
    fields = slot_indices(M, ['c', 'A', 'P', 'L', 'cb'])
    fields = [i for i in fields if M.fields[i].name not in ('cs', 'cbs')]
    antif = [i for i in range(len(M.fields)) if i not in fields]
    e_all = O.generic_section(M, range(len(M.fields)))
    # Q from S0 equals Q_Psi of the momentum construction; H = [Q, QGF] = -Laplacian on every slot
    res['Q_from_S0_equals_Q_Psi'] = O.is_zero_section(O.sec_sub(Qop(e_all), Qmom(e_all)))
    H = O.gcomm(Qop, G, e_all)
    lap = {a: -sum(sp.diff(f, t, 2) for t in X) for a, f in e_all.items()}
    res['H=[Q,QGF]=-d^2_every_slot'] = O.is_zero_section(O.sec_sub(H, lap))
    res['QGF^2=0'] = O.is_zero_section(G(G(e_all)))
    gens = {}
    for k in ('P', 'M', 'D', 'K'):
        gens[k] = Mo.S1_even(F, v=gens_even[k])
    gens['R'] = Mo.S1_even(F, r=r)
    gens['xi'] = Mo.S1_odd(F, *odd['xi1'])
    gens['eta'] = Mo.S1_odd(F, *odd['eta1'])
    out_tables = {}
    for k, S1 in gens.items():
        xo = O.op_from_functional(M, S1)
        c = O.gcomm(xo, G, e_all)
        comp = G(xo(G(e_all)))
        if k in ('K', 'xi', 'eta'):
            res['composite QGF x. QGF nonzero along ' + k] = not O.is_zero_section(comp)
        else:
            res['composite QGF x. QGF = 0 along ' + k] = O.is_zero_section(comp)
        # V_x : the field -> antifield part of x.
        e_f = {a: e_all[a] for a in fields}
        xf = xo(e_f)
        Vx = {a: v for a, v in xf.items() if a in antif}
        comp2 = G({a: v for a, v in xo(G(e_all)).items() if a in antif})
        res['composite = QGF V_x QGF along ' + k] = O.is_zero_section(O.sec_sub(comp, comp2))
        out_tables[k] = {M.fields[a].name: str(v) for a, v in c.items()}
        if k == 'D':
            res['[Delta., QGF] = -4 QGF'] = O.is_zero_section(O.sec_sub(c, O.sec_scale(G(e_all), -4)))
            Hd = O.gcomm(xo, Qop, e_all)
            res['[Delta., Q] = 0'] = O.is_zero_section(Hd)
        if k in ('P', 'M', 'R'):
            res['[x., QGF] = 0 along ' + k] = O.is_zero_section(c)
    if verbose:
        print('ops %.1fs' % (time.time() - t0))
    return res, out_tables


if __name__ == '__main__':
    res, tabs = run()
    for k, v in res.items():
        print(k, v)
    json.dump(tabs, open('commutator_tables.json', 'w'), indent=1)
