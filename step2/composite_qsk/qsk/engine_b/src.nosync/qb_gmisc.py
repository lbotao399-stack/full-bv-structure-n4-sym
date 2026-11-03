"""Gates G-cubic (Example qs-cubic), G-swap (graded symmetry of the pair kernel under X_a <-> X_b), G-chirality
(the kernel along xi in S_- (x) W^* on three pairs has the leg sets obtained by the chirality swap)."""
import os, json, time
import sympy as sp
import qb_core as B
import qb_job as QJ

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')


def ex(job, P, swap=False):
    syms = []
    for n in job.names:
        if swap and n.startswith('za_'):
            n = 'zb_' + n[3:]
        elif swap and n.startswith('zb_'):
            n = 'za_' + n[3:]
        syms.append(sp.Symbol(n))
    return sp.expand(P.as_expr(*syms)) if P else sp.Integer(0)


def run():
    t0 = time.time()
    gates = []; res = {}
    xi = B.XI[:8] + [0] * 8
    # G-cubic: the T5 cycle at the special edge, (A_perp, lam_+) -> A_perp: eps^0 = (3/4)|q|^2 x eps^-1, i.e.
    # 16 pi^2 x value = (-1/(2 eps) - (3/8)|q|^2) x factor
    job = QJ.DJob('AW', 'lamP', ['AW'], 'xi', xi, terms='cubic', tag='gc')
    r = job.kernel(mlist=(0, 1))['total']
    P0 = ex(job, r[(0,)][0]['1']); P1 = ex(job, r[(0,)][1]['1'])
    q2 = sum(sp.Symbol(f'r0_{m}') ** 2 for m in range(4))
    okc = sp.expand(P0 - sp.Rational(3, 4) * q2 * P1) == 0 and P1 != 0 and set(r[(0,)][0].keys()) == {'1'}
    res['G-cubic'] = dict(eps_m1=str(P1)[:400] + ' ...', ratio='eps^0 / eps^-1 = (3/4)|q|^2', bases=sorted(r[(0,)][0].keys()))
    gates.append(dict(id='G-cubic', pass_=okc, value='eps^0 - (3/4)|q|^2 eps^-1 = 0, eps^0 rational' if okc else 'mismatch',
                      target='Example qs-cubic: 16 pi^2 x value = (-1/(2 eps) - (3/8)|q|^2) x (T5 factor x fibre factors)'))
    # G-swap on three pairs
    sw = {}
    oks = True
    for sa, sb, outs in [('AW', 'lamP', ['AW']), ('AV', 'lamP', ['AV']), ('lamP', 'lamM', ['lamM'])]:
        j1 = QJ.DJob(sa, sb, outs, 'xi', xi, tag='sw1'); r1 = j1.kernel(mlist=(0, 1))['total']
        j2 = QJ.DJob(sb, sa, outs, 'xi', xi, tag='sw1'); r2 = j2.kernel(mlist=(0, 1))['total']
        pa = B.V.PAR[B.SLOT[sa][0]]; pb = B.V.PAR[B.SLOT[sb][0]]
        ok = True
        for w, dm in r1.items():
            for m, d in dm.items():
                for bb, P in d.items():
                    Q = r2.get(tuple(reversed(w)), {}).get(m, {}).get(bb)
                    s = (-1) ** (pa * pb + len(w))
                    if sp.expand(ex(j2, Q, True) - s * ex(j1, P)) != 0:
                        ok = False
        sw[f'{sa}|{sb} -> {" ".join(outs)}'] = ok
        oks = oks and ok
    res['G-swap'] = sw
    gates.append(dict(id='G-swap', pass_=oks, value=sw,
                      target='Upsilon(X_b,X_a)_w = (-1)^{|X_a||X_b| + |w|} Upsilon(X_a,X_b)_{reversed w} with the letters relabelled'))
    # G-chirality: along xi in S_- (x) W^* the nonzero leg sets of three pairs are those of S_+ (x) W with lamP <-> lamM
    xim = [0] * 8 + B.XI[8:]
    chir = {}
    okch = True
    for sa, sb, outsP, outsM in [('AW', 'AW', ['lamM'], ['lamP']), ('AV', 'AV', ['lamM'], ['lamP']), ('AW', 'lamM', ['AW'], ['AW'])]:
        jp = QJ.DJob(sa, sb, outsP, 'xi', xi, tag='chp'); rp = jp.kernel(mlist=(0, 1))['total']
        sbm = {'lamP': 'lamM', 'lamM': 'lamP'}.get(sb, sb); sam = {'lamP': 'lamM', 'lamM': 'lamP'}.get(sa, sa)
        jm = QJ.DJob(sam, sbm, outsM, 'xi', xim, tag='chm'); rm = jm.kernel(mlist=(0, 1))['total']
        cnt = lambda r: sorted((str(w), m, b, len(P.terms())) for w, dm in r.items() for m, d in dm.items() for b, P in d.items())
        same = cnt(rp) == cnt(rm)
        chir[f'{sa}|{sb}->{outsP} vs {sam}|{sbm}->{outsM}'] = dict(S_plus=cnt(rp), S_minus=cnt(rm), same_structure=same)
        okch = okch and same
    res['G-chirality'] = chir
    gates.append(dict(id='G-chirality', pass_=okch, value={k: v['same_structure'] for k, v in chir.items()},
                      target='kernels along S_- (x) W^* = chirality-swapped kernels along S_+ (x) W (same words, poles, bases and number of terms) on three pairs'))
    res['seconds'] = round(time.time() - t0, 1)
    return res, gates


if __name__ == '__main__':
    res, gates = run()
    json.dump(dict(results=res, gates=gates), open(os.path.join(OUT, 'gmisc.json'), 'w'), indent=1)
    for g in gates:
        print(g['id'], 'PASS' if g['pass_'] else 'FAIL', g['value'])
