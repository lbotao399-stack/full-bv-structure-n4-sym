"""Engine B: one-loop graphs with one insertion of O = int omega (iota_{w1} A, iota_{w2} A)  (omega = plane wave of momentum
q_O, amplitude 1) and partial_delta I^CT[eps].  d_delta W_1 = (hbar/2) sum_n STr[(P Pi H)^{n-1} P Pi H_O] written out over all
words (no cyclicity assumed)."""
import json, time, itertools
import sympy as sp
import vertices as V
import loops as L
from sectors import leg, AVc, APc, CUBW, QUAW, fit

I = sp.I
q = L.q
p1 = sp.symbols('p0:4'); p2 = sp.symbols('k0:4')
w1 = [0] * 4 + list(sp.symbols('w1_0:6')); w2 = [0] * 4 + list(sp.symbols('w2_0:6'))
qO = [-p1[m] - p2[m] for m in range(4)]

def L_O(s):
    a, b = s
    return sum(a['u'][M] * w1[M] for M in range(10)) * sum(b['u'][M] * w2[M] for M in range(10))
V.TERMS['O'] = (sp.Integer(1), ('A', 'A'), L_O, ('k', (0, 1)))

def H_of(kind, ext, k, kp):
    if kind == 'cubic':
        return L.full_hessian(V.CUBIC, ext, k, kp, CUBW)
    if kind == 'quartic':
        return L.full_hessian(V.QUARTIC, ext, k, kp, QUAW)
    if kind == 'O':
        return L.full_hessian(['O'], [], k, kp, {'k': 1})

def word_value(word, mask):
    """word: list of (kind, ext legs, injected momentum).  Returns (a,b) of int STr[P Pi h_1 P Pi h_2 ...] (units 1/16pi^2)."""
    mats = []; lines = []; cur = [sp.Integer(0)] * 4      # line momentum = q - cur
    for (kind, ext, inj) in word:
        lines.append(list(cur))
        Nq = L.prop_numerator([q[m] - cur[m] for m in range(4)], mask)
        nxt = [cur[m] + inj[m] for m in range(4)]
        h = H_of(kind, ext, [-(q[m] - cur[m]) for m in range(4)], [q[m] - nxt[m] for m in range(4)])
        mats += [Nq, L.PiH(h)]
        cur = nxt
    assert all(sp.simplify(c) == 0 for c in cur)
    Num = L.str_product(mats)
    if Num == 0:
        return sp.Integer(0), sp.Integer(0)
    a, b = L.boundary_singular(Num, lines)
    return sp.expand(a * 16 * sp.pi**2), sp.expand(b * 16 * sp.pi**2)

def structures(e1, e2):
    def val(Lf):
        V.TERMS['_s'] = (sp.Integer(1), ('A', 'A'), Lf, ('k', (0, 1)))
        v = V.vertex('_s', [e1, e2]).get('k', 0); del V.TERMS['_s']; return sp.expand(v)
    w12 = sum(w1[M] * w2[M] for M in range(10))
    return {
        'O': val(L_O),
        'trace_mix': val(lambda s: w12 * sum(s[0]['u'][M] * s[1]['u'][M] for M in range(4, 10))),
        'AV_mix': val(lambda s: w12 * sum(s[0]['u'][M] * s[1]['u'][M] for M in range(4))),
    }

def fit_struct(target, S, names):
    cs = sp.symbols(f'C0:{len(names)}')
    expr = sp.expand(target - sum(c * S[n] for c, n in zip(cs, names)))
    syms = sorted(expr.free_symbols - set(cs), key=lambda s: s.name)
    eqs = sp.Poly(expr, *syms).coeffs() if syms else [expr]
    eqs = [e for e in eqs if e != 0]
    if not eqs:
        return {n: '0 (graph vanishes)' if target == 0 else 'undetermined' for n in names}
    sol = sp.solve(eqs, cs, dict=True)
    if not sol:
        return None
    return {n: str(sp.nsimplify(sol[0].get(c, c))) for c, n in zip(cs, names)}

def run():
    out = {}
    for sector, comps, masks in (('Ap-Ap', APc, {'tri': {('AV', 'AV'), ('Ap', 'Ap')}, 'bub': {('Ap', 'Ap')}}),
                                 ('AV-AV', AVc, {'tri': {('Ap', 'Ap')}, 'bub': {('Ap', 'Ap')}})):
        e1 = leg('A', 'u', p1, comps); e2 = leg('A', 'v', p2, comps)
        S = structures(e1, e2)
        names = ['O', 'trace_mix'] if sector == 'Ap-Ap' else ['AV_mix']
        # triangle: (hbar/2)[STr(P h1 P h2 P hO) + STr(P h2 P h1 P hO)] (+ the 4 other words equal to these by cyclicity are
        # accounted for by the 1/3 -> each cyclic class appears 3 times); computed here as the two cyclic classes, weight 1/2.
        t0 = time.time()
        aT = 0; bT = 0
        for word in ([('cubic', [e1], list(p1)), ('cubic', [e2], list(p2)), ('O', [], qO)],
                     [('cubic', [e2], list(p2)), ('cubic', [e1], list(p1)), ('O', [], qO)]):
            a, b = word_value(word, masks['tri']); aT += a / 2; bT += b / 2
        gid = f'{sector}+O:AAA+AAA+O:' + ('Ap,AV,Ap' if sector == 'Ap-Ap' else 'Ap,Ap,Ap')
        out[gid] = dict(sector=sector + ' with O', eps_inv=fit_struct(aT, S, names), log_eps=fit_struct(bT, S, names),
                        raw_a=str(aT), raw_b=str(bT)[:300], seconds=round(time.time() - t0, 1))
        print(gid, out[gid]['eps_inv'], out[gid]['log_eps'], flush=True)
        # quartic + O bubble: (hbar/2)(1/2)[STr(P h12 P hO) + STr(P hO P h12)]
        t0 = time.time()
        aB = 0; bB = 0
        p12 = [p1[m] + p2[m] for m in range(4)]
        for word in ([('quartic', [e1, e2], p12), ('O', [], qO)], [('O', [], qO), ('quartic', [e1, e2], p12)]):
            a, b = word_value(word, masks['bub']); aB += a / 4; bB += b / 4
        gid = f'{sector}+O:AAAA+O:Ap,Ap'
        out[gid] = dict(sector=sector + ' with O', eps_inv=fit_struct(aB, S, names), log_eps=fit_struct(bB, S, names),
                        raw_a=str(aB), raw_b=str(bB)[:300], seconds=round(time.time() - t0, 1))
        print(gid, out[gid]['eps_inv'], out[gid]['log_eps'], flush=True)
    # vacuum: (hbar/2) STr(P Pi h_O) with q_O = 0, colour tr((-,-)^{-1}(-,-)) = dim g
    h = H_of('O', [], [-x for x in q], list(q))
    Num = L.str_product([L.prop_numerator(list(q), {('Ap', 'Ap')}), L.PiH(h)])
    a, b = L.boundary_singular(Num, [[0] * 4])
    w12 = sum(w1[M] * w2[M] for M in range(10))
    out['vac:O:Ap'] = dict(sector='vacuum (units hbar dim(g)/(16 pi^2), times int omega)',
                           eps_inv=str(sp.simplify(a * 16 * sp.pi**2 / 2 / w12)) + ' * (w1,w2)',
                           log_eps=str(sp.simplify(b * 16 * sp.pi**2 / 2)))
    print('vac', out['vac:O:Ap'])
    return out

if __name__ == '__main__':
    res = run()
    with open('o_graphs.json', 'w') as fh:
        json.dump(res, fh, indent=1)
