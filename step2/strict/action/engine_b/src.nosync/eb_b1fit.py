"""Task B1 / S1 (weight 5, p <= 1, two legs): cohomology of Q on the degree -1 functionals linear in b, and the
decomposition of Ob^(2)(Delta,P_b)|_2 (theory normalisation = -(raw), eb_b1) onto [S(P_b)_d] modulo Q-exact.
theta_d = o_d / (2u) (Theorem at-theta)."""
import json, os, random, sys, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import eb_qb as Q
import eb_coh as H

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')
LOGK = ['1', 'log2', 'log3']


def parse_clog(d):
    out = {}
    for k in LOGK:
        re = Fr(d['re'].get(k, '0'))
        im = Fr(d['im'].get(k, '0'))
        out[k] = (re, im)
    return out


def load_sector(sec, sign=-1):
    recs = json.load(open(os.path.join(OUT, f'B1_raw_{sec}.json')))
    import eb_b1
    s0, s1, r0, r1 = eb_b1.SECT[sec]
    configs, obvec = [], {k: [] for k in LOGK}
    for rec in recs:
        c = rec['config']
        u0 = np.array([complex(x) for x in c['u0']]); u1 = np.array([complex(x) for x in c['u1']])
        p = [Fr(x) for x in c['p']]; b = np.array([complex(Fr(x)) for x in c['b']])

        def legs_at(r, u0=u0, u1=u1, p=p, s0=s0, s1=s1):
            return [dict(sp=s0, u=u0, p=p), dict(sp=s1, u=u1, p=[-p[m] - r[m] for m in range(4)])]
        configs.append((legs_at, b))
        al = parse_clog(rec['alpha']['k']['0']) if 'k' in rec['alpha'] else {k: (Fr(0), Fr(0)) for k in LOGK}
        pole = rec['alpha'].get('k', {}).get('1')
        assert pole is None or all(v == '0' for part in pole.values() for v in part.values()) or not any(pole['re'].values()) and not any(pole['im'].values()), pole
        bes = []
        for mu in range(4):
            bb = rec['beta'][str(mu)]
            bes.append(parse_clog(bb['k']['0']) if 'k' in bb else {k: (Fr(0), Fr(0)) for k in LOGK})
        for k in LOGK:
            obvec[k] += [(sign * al[k][0], sign * al[k][1])] + [(sign * be[k][0], sign * be[k][1]) for be in bes]
    return configs, obvec


def rand_cfg(species, rngs):
    p = [Fr(random.randint(-3, 3)) for _ in range(4)]
    us = []
    for spc, rg in zip(species, rngs):
        u = np.zeros(Q.DIM[spc], complex)
        for c in (rg if rg is not None else range(Q.DIM[spc])):
            u[c] = random.randint(-2, 2)
        us.append(u)
    b = np.array([random.randint(-2, 2) for _ in range(4)], complex)

    def legs_at(r, us=us, p=p):
        return [dict(sp=species[0], u=us[0], p=p), dict(sp=species[1], u=us[1], p=[-p[m] - r[m] for m in range(4)])]
    return legs_at, b


V4 = range(4); P6 = range(4, 10); LP = range(8); LM = range(8, 16)
BLOCKS = {
    'V': dict(sectors=['c_cs', 'AV_AsV'],
              deg1=[(('c', 'cs'), 'SS', 1, False), (('c', 'cs'), 'SS', 2, True),
                    (('A', 'As'), 'VV', 1, False), (('A', 'As'), 'VV', 2, True)],
              deg2=[(('cs', 'A'), 'SV', 0, False), (('cs', 'A'), 'SV', 1, True), (('As', 'As'), 'VV', 0, True)],
              targets=[(('c', 'As'), (None, V4)), (('As', 'c'), (V4, None)), (('A', 'A'), (V4, V4))]),
    'P': dict(sectors=['Ap_Asp'],
              deg1=[(('A', 'As'), 'PP', 1, False), (('A', 'As'), 'PP', 2, True)],
              deg2=[(('As', 'As'), 'PP', 0, True)],
              targets=[(('A', 'A'), (P6, P6))]),
    'L': dict(sectors=['lam+_lams+', 'lam-_lams-'],
              deg1=[(('lam', 'lams'), 'LLs', 1, False), (('lam', 'lams'), 'LLs', 2, True)],
              deg2=[(('lams', 'lams'), 'LsLs', 0, False), (('lams', 'lams'), 'LsLs', 1, True)],
              targets=[(('lam', 'lam'), (None, None))]),
}


def sector_parts():
    """S(P_b)_d densities (theory Def at-names (ii)) as Fn lists with coefficients."""
    bd = lambda s, ctx: np.dot(1j * Q.pv(s[0]), ctx['b'])
    SV = [H.Fn('S(P_b)_V:A', ('A', 'As'), lambda s, ctx: bd(s, ctx) * np.dot(s[0]['u'][:4], s[1]['u'][:4]), False),
          H.Fn('S(P_b)_V:c', ('c', 'cs'), lambda s, ctx: -bd(s, ctx) * s[0]['u'][0] * s[1]['u'][0], False)]
    SP = [H.Fn('S(P_b)_perp', ('A', 'As'), lambda s, ctx: bd(s, ctx) * np.dot(s[0]['u'][4:], s[1]['u'][4:]), False)]
    SL = [H.Fn('S(P_b)_lam', ('lam', 'lams'), lambda s, ctx: -bd(s, ctx) * (s[0]['u'] @ Q.MP @ s[1]['u']), False)]
    return {'V': SV, 'P': SP, 'L': SL}


def sum_rows(rows):
    out = []
    for col in zip(*rows):
        out.append((sum(c[0] for c in col), sum(c[1] for c in col)))
    return out


def to_sp(z):
    return sp.Rational(z[0].numerator, z[0].denominator) + sp.I * sp.Rational(z[1].numerator, z[1].denominator)


def exact_solve(Arows, b):
    """solve sum_k x_k Arows[k] = b exactly over Q(i); returns list or None."""
    M = sp.Matrix([[to_sp(z) for z in row] for row in Arows]).T
    bv = sp.Matrix([to_sp(z) for z in b])
    try:
        sol, params = M.gauss_jordan_solve(bv)
    except ValueError:
        return None
    if params.shape[0]:
        sol = sol.subs({p: 0 for p in params})
    return [sp.nsimplify(x) for x in sol]


def analyse_block(name, log=print, ntarget=6, seed=17):
    random.seed(seed)
    blk = BLOCKS[name]
    configs, obv = [], {k: [] for k in LOGK}
    for sec in blk['sectors']:
        c, o = load_sector(sec)
        configs += c
        for k in LOGK:
            obv[k] += o[k]
    # degree -1 candidates
    cands = []
    for slots, kind, nm, we in blk['deg1']:
        for lab, L in H.make_candidates(slots, kind, nm, we):
            cands.append(H.Fn(lab, slots, L, we))
    t0 = time.time()
    rows = H.evaluate_fns(cands, configs)
    keep = H.independent(rows)
    basis = [cands[i] for i in keep]
    brows = [rows[i] for i in keep]
    log(f'B1/S1 block {name}: {len(cands)} degree -1 candidates, rank {len(basis)} on {len(configs)} configs ({time.time() - t0:.0f} s)')
    # closedness: Q on target configs
    tcfgs = []
    for (spc, rg) in blk['targets']:
        for _ in range(ntarget):
            tcfgs.append(rand_cfg(spc, rg))
    qrows = H.evaluate_fns(basis, tcfgs, apply_q=True)
    ker = Q.nullspace_mod(H.mod_rows(qrows))
    # exact kernel: solve over sympy
    Mq = sp.Matrix([[to_sp(z) for z in row] for row in qrows]).T
    kerx = Mq.nullspace()
    log(f'  closed subspace dim {len(kerx)} (mod p {len(ker)})')
    # degree -2 candidates and their Q images in the degree -1 basis coordinates
    c2 = []
    for slots, kind, nm, we in blk['deg2']:
        for lab, L in H.make_candidates(slots, kind, nm, we):
            c2.append(H.Fn(lab, slots, L, we))
    img = H.evaluate_fns(c2, configs, apply_q=True)
    ik = H.independent(img) if any(any(z != (0, 0) for z in r) for r in img) else []
    imrows = [img[i] for i in ik]
    log(f'  degree -2 candidates {len(c2)}, image dim {len(imrows)}')
    # express image in basis coordinates
    imcoords = []
    for r in imrows:
        x = exact_solve(brows, r)
        assert x is not None, 'Q image not in the span of the degree -1 basis'
        imcoords.append(sp.Matrix(x))
    # Ob in basis coordinates per log component
    obcoords = {}
    for k in LOGK:
        x = exact_solve(brows, obv[k])
        obcoords[k] = None if x is None else sp.Matrix(x)
    # Ob closed?  (Q Ob = 0): check coordinates are in kernel span
    kerM = sp.Matrix.hstack(*kerx) if kerx else sp.zeros(len(basis), 0)
    closed = {}
    for k, x in obcoords.items():
        if x is None:
            closed[k] = None
            continue
        closed[k] = (sp.Matrix.hstack(kerM, x).rank() == kerM.rank())
    # sector parts coordinates
    sp_parts = sector_parts()
    Sd = {}
    for nm2, fns in sp_parts.items():
        rr = sum_rows(H.evaluate_fns(fns, configs))
        x = exact_solve(brows, rr)
        if x is not None and any(v != 0 for v in x):
            Sd[nm2] = sp.Matrix(x)
    # cohomology: H = ker / im ; decompose Ob = sum o_d S_d + exact + other classes
    imM = sp.Matrix.hstack(*imcoords) if imcoords else sp.zeros(len(basis), 0)
    hdim = kerM.rank() - imM.rank()
    res = dict(n_candidates=len(cands), basis_dim=len(basis), closed_dim=kerM.rank(), exact_dim=imM.rank(),
               H_minus1_dim=hdim, ob_in_span={k: obcoords[k] is not None for k in LOGK}, ob_closed=closed)
    # classes of the sector parts independent mod image?
    Scols = list(Sd.items())
    if Scols:
        SM = sp.Matrix.hstack(*[c for _, c in Scols])
        res['sector_parts_closed'] = [bool(sp.Matrix.hstack(kerM, c).rank() == kerM.rank()) for _, c in Scols]
        res['sector_parts_independent_mod_exact'] = bool(sp.Matrix.hstack(imM, SM).rank() == imM.rank() + SM.shape[1])
        # other classes: complete {S_d} to a basis of H
        reps = []
        cur = sp.Matrix.hstack(imM, SM)
        for j in range(kerM.shape[1]):
            v = kerM[:, j]
            if sp.Matrix.hstack(cur, v).rank() > cur.rank():
                reps.append(v)
                cur = sp.Matrix.hstack(cur, v)
        res['other_classes'] = len(reps)
        # decompose Ob
        dec = {}
        for k, x in obcoords.items():
            if x is None:
                dec[k] = None
                continue
            A = sp.Matrix.hstack(SM, *(reps or [sp.zeros(len(basis), 0)]), imM) if (reps or imM.shape[1]) else SM
            sol, params = A.gauss_jordan_solve(x)
            sol = sol.subs({p: 0 for p in params})
            dec[k] = {Scols[i][0]: str(sp.nsimplify(sol[i])) for i in range(len(Scols))}
            dec[k]['other_classes'] = [str(sp.nsimplify(sol[len(Scols) + i])) for i in range(len(reps))]
        res['Ob_decomposition'] = dec
    log(f'  {res}')
    return res


if __name__ == '__main__':
    out = {}
    for nm in (sys.argv[1:] or ['V', 'P', 'L']):
        out[nm] = analyse_block(nm)
    with open(os.path.join(OUT, 'B1_fit.json'), 'w') as f:
        json.dump(out, f, indent=1, default=str)
