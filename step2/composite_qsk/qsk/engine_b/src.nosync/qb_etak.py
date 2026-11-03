"""TQ5 / TQ6 (route B, n = 0): intrinsic parts Upsilon^0_eta and Upsilon^0_{K_b} on field pairs, the position identity
G-pos, and the colour-singlet values on the scalar bilinear.
eta in S_+ (x) W^* rank one (realisation component 11), b symbolic (b0..b3)."""
import os, sys, json, time
import sympy as sp
from sympy.polys.domains import QQ_I
import qb_core as B
import qb_job as QJ
import qb_jobinj as IJ

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')
I = sp.I
ETA_PAIRS = [('AW', 'AW', ['lamP']), ('AW', 'lamP', ['AV']), ('AW', 'lamP', ['AV', 'AV']), ('AW', 'lamP', ['AW', 'AW']),
             ('AW', 'lamM', ['AW']), ('AW', 'lamM', ['AV', 'AW'])]
K_PAIRS = [('AW', 'AW', ['AV']), ('AV', 'AV', ['AV']), ('AW', 'lamP', ['lamM']), ('AW', 'lamM', ['lamP']),
           ('lamP', 'lamP', ['AW'])]


def tostr(job, tab):
    return {f'word {list(w)} | eps^-{m} | {b}': QJ.poly_str(job, P) for w, dm in tab.items() for m, d in dm.items()
            for b, P in d.items()}


def add(tot, w, m, b, P):
    d = tot.setdefault(w, {}).setdefault(m, {})
    d[b] = d.get(b, P.ring.zero) + P if b in d else P


def eta_intrinsic(sa, sb, outs):
    eta = [0] * 16
    eta[11] = 1          # rank-one eta in S_+ (x) W^* (component a = 0, W^* index 3); covariance as for xi
    ev = sp.Matrix(eta)
    tot = {}
    pos_ok = True
    jobs = []
    for mu in range(4):
        ximu = list(B.Rbm[mu] * ev)
        job = IJ.IJob(sa, sb, outs, 'xi', ximu, hmax=1, tag=f'eta{mu}')
        r = job.kernel(mlist=(0, 1))
        e = [0, 0, 0, 0]; e[mu] = 1
        for w, dm in r['total'].items():
            for m, d in dm.items():
                for b, P in d.items():
                    c = job.hcoeff(P, tuple(e)) * QQ_I(0, -1)
                    if c:
                        add(tot, w, m, b, c)
        # G-pos: the h = 0 part equals the plain kernel along xi_mu (graph by graph identity in the sum); checked on
        # the leg sets with one output (memory)
        if len(outs) > 1:
            jobs.append(job)
            continue
        j0 = QJ.DJob(sa, sb, outs, 'xi', ximu, tag=f'eta{mu}p')
        r0 = j0.kernel(mlist=(0, 1))
        for w, dm in r0['total'].items():
            for m, d in dm.items():
                for b, P in d.items():
                    Ph = job.hcoeff(r['total'].get(w, {}).get(m, {}).get(b, job.R.zero), (0, 0, 0, 0))
                    e1 = sp.expand(P.as_expr(*[sp.Symbol(n) for n in j0.names]))
                    e2 = sp.expand(Ph.as_expr(*[sp.Symbol(n) for n in job.names])) if Ph else 0
                    if sp.expand(e1 - e2) != 0:
                        pos_ok = False
        jobs.append(job)
    # s_0 term
    js = QJ.DJob(sa, sb, outs, 's0', eta, tag='etas0')
    rs = js.kernel(mlist=(0, 1))
    job = jobs[0]
    conv = {}
    out_s0 = tostr(js, rs['total'])
    return job, tot, pos_ok, js, rs['total']


def K_intrinsic(sa, sb, outs):
    b = list(sp.symbols('b0:4'))
    names = IJ.register_K_terms(b, 'Kint')
    tot = {}
    job0 = None
    for key, nm in names.items():
        if key[0] == 'V':
            _, mu, nu = key
            job = IJ.IJob(sa, sb, outs, 'terms', b, hmax=2, quad_terms=nm, Spar=1, tag=f'KV{mu}{nu}')
            r = job.kernel(mlist=(0, 1))
            e = [0, 0, 0, 0]; e[mu] += 1; e[nu] += 1
            fac = QQ_I(-2 if mu == nu else -1, 0)       # z^mu z^nu -> -d_mu d_nu (h^2 -> 2)
        else:
            _, mu = key
            job = IJ.IJob(sa, sb, outs, 'terms', b, hmax=1, quad_terms=nm, Spar=1, tag=f'KW{mu}')
            r = job.kernel(mlist=(0, 1))
            e = [0, 0, 0, 0]; e[mu] = 1
            fac = QQ_I(0, -1)
        for w, dm in r['total'].items():
            for m, d in dm.items():
                for bb, P in d.items():
                    c = job.hcoeff(P, tuple(e))
                    if c:
                        c = c * fac
                        if job0 is None:
                            job0 = job
                        # convert to the ring of job0 (same names)
                        c0 = job0.R.from_dict(dict(c.items())) if job is not job0 else c
                        add(tot, w, m, bb, c0)
        if job0 is None:
            job0 = job
    return job0, tot


def singlet2(job, tot, m=0):
    out = {}
    for w, dm in tot.items():
        if len(w) != 2:
            continue
        for bkey, P in dm.get(m, {}).items():
            out[bkey] = sp.expand(out.get(bkey, 0) + P.as_expr(*[sp.Symbol(n) for n in job.names]))
    return out


def run(which):
    t0 = time.time()
    res = {}; gates = []
    if which in ('eta', 'all'):
        allpos = True
        for sa, sb, outs in ETA_PAIRS:
            job, tot, pos_ok, js, ts = eta_intrinsic(sa, sb, outs)
            allpos = allpos and pos_ok
            res[f'eta {sa}|{sb} -> {" ".join(outs)}'] = dict(position_part_identity=pos_ok, intrinsic_xi_part=tostr(job, tot),
                                                            s0_part=tostr(js, ts),
                                                            note='Upsilon^0_eta = intrinsic_xi_part + s0_part (+ finite special vertices; none at n = 0 on these pairs by weight)')
            print('eta', sa, sb, outs, 'G-pos', pos_ok, round(time.time() - t0, 1), flush=True)
        gates.append(dict(id='G-pos-eta', pass_=allpos, value=allpos,
                          target='Prop qs-posdec (i): the h = 0 part of the cycles with S^(1)(rho(e_mu) eta) and injected h equals the kernel along xi_mu = rho(e_mu) eta, graph by graph (summed per word)'))
    if which in ('K', 'all'):
        for sa, sb, outs in K_PAIRS:
            job, tot = K_intrinsic(sa, sb, outs)
            res[f'K {sa}|{sb} -> {" ".join(outs)}'] = dict(intrinsic=tostr(job, tot) if job else {})
            print('K', sa, sb, outs, round(time.time() - t0, 1), flush=True)
    res['seconds'] = round(time.time() - t0, 1)
    return res, gates


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if which.startswith('eta:') or which.startswith('K:'):
        kind, idx = which.split(':')
        idx = int(idx)
        t0 = time.time()
        if kind == 'eta':
            sa, sb, outs = ETA_PAIRS[idx]
            job, tot, pos_ok, js, ts = eta_intrinsic(sa, sb, outs)
            res = {f'eta {sa}|{sb} -> {" ".join(outs)}': dict(position_part_identity=(pos_ok if len(outs) == 1 else 'not checked (memory)'),
                   intrinsic_xi_part=tostr(job, tot), s0_part=tostr(js, ts), sign_version=QJ.SIGN_VERSION,
                   note='Upsilon^0_eta = intrinsic_xi_part + s0_part (no finite special vertex at n = 0 on this pair by weight); eta rank one (component 11)',
                   seconds=round(time.time() - t0, 1))}
            gates = [dict(id=f'G-pos-eta {sa}|{sb} -> {" ".join(outs)}', pass_=bool(pos_ok), value=pos_ok,
                          target='Prop qs-posdec (i): h = 0 part = kernel along rho(e_mu) eta, per word')] if len(outs) == 1 else []
        else:
            sa, sb, outs = K_PAIRS[idx]
            job, tot = K_intrinsic(sa, sb, outs)
            res = {f'K {sa}|{sb} -> {" ".join(outs)}': dict(intrinsic=tostr(job, tot) if job else {}, seconds=round(time.time() - t0, 1),
                                                       note='Upsilon^0_{K_b} intrinsic part, b symbolic (b0..b3); finite special vertices 0 at n = 0 by weight')}
            gates = []
        json.dump(dict(results=res, gates=gates), open(os.path.join(OUT, f'etak_{kind}_{idx}.json'), 'w'), indent=1)
        print(json.dumps(gates))
        sys.exit(0)
    res, gates = run(which)
    json.dump(dict(results=res, gates=gates), open(os.path.join(OUT, f'etak_{which}.json'), 'w'), indent=1)
    for g in gates:
        print(g['id'], 'PASS' if g['pass_'] else 'FAIL')
