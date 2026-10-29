"""Gaussian-moment data for the O-insertion graphs (report formula chains)."""
import json
import sympy as sp
import vertices as V
import loops as L
import sectors as S
import o_graphs as OG

q = L.q
def word_data(word, mask, e1, e2):
    mats = []; lines = []; cur = [sp.Integer(0)] * 4
    for (kind, ext, inj) in word:
        lines.append(list(cur))
        Nq = L.prop_numerator([q[m] - cur[m] for m in range(4)], mask)
        nxt = [cur[m] + inj[m] for m in range(4)]
        hh = OG.H_of(kind, ext, [-(q[m] - cur[m]) for m in range(4)], [q[m] - nxt[m] for m in range(4)])
        mats += [Nq, L.PiH(hh)]; cur = nxt
    Num = L.str_product(mats)
    E = len(lines)
    betas = sp.symbols(f'beta0:{E}')
    shift = [sum(betas[e] * lines[e][m] for e in range(E)) for m in range(4)]
    C = L.gaussian_moments(Num, shift)
    S_ = OG.structures(e1, e2)
    names = ['O', 'trace_mix'] if e1['u'][5] != 0 else ['AV_mix']
    out = {}
    for k in sorted(C):
        Ck = C[k]
        # coefficient at zero external momenta (only these reach eps^0 for the degree-0 O graphs)
        Ck0 = sp.expand(Ck.subs({**{x: 0 for x in OG.p1}, **{x: 0 for x in OG.p2}}))
        if Ck0 == 0:
            continue
        out[f'C{k}(zero ext. momenta)'] = OG.fit_struct(Ck0, S_, names)
    return out, E

res = {}
for sector, comps, masks in (('Ap-Ap', S.APc, {'tri': {('AV', 'AV'), ('Ap', 'Ap')}, 'bub': {('Ap', 'Ap')}}),
                             ('AV-AV', S.AVc, {'tri': {('Ap', 'Ap')}, 'bub': {('Ap', 'Ap')}})):
    e1 = S.leg('A', 'u', OG.p1, comps); e2 = S.leg('A', 'v', OG.p2, comps)
    p12 = [OG.p1[m] + OG.p2[m] for m in range(4)]
    for name, word, mask in [('triangle (leg p1, leg p2, O)', [('cubic', [e1], list(OG.p1)), ('cubic', [e2], list(OG.p2)), ('O', [], OG.qO)], masks['tri']),
                             ('triangle (leg p2, leg p1, O)', [('cubic', [e2], list(OG.p2)), ('cubic', [e1], list(OG.p1)), ('O', [], OG.qO)], masks['tri']),
                             ('bubble (quartic, O)', [('quartic', [e1, e2], p12), ('O', [], OG.qO)], masks['bub']),
                             ('bubble (O, quartic)', [('O', [], OG.qO), ('quartic', [e1, e2], p12)], masks['bub'])]:
        d, E = word_data(word, mask, e1, e2)
        res[f'{sector} {name}'] = dict(edges=E, moments=d)
        print(sector, name, E, d, flush=True)
json.dump(res, open('report_data_O.json', 'w'), indent=1)
