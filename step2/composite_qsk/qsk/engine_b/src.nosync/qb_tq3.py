"""TQ3 (route B, n = 0 entries): the twist pairs at xi = u (x) e_4, restricted to the twisted letters.

Realisation: xi = u (x) e_4 with u the first basis vector of S_+ (xi component 3 = 1); u' the dual basis functional
(u'^alpha u_alpha = 1, u' annihilating the second basis vector); beta_I = u'^alpha lam_alpha^I = the lam_+ component
(a = 0, W index I-1); gamma^I = A_{I4} = sum_K Sig^K_{I-1,3} A_{4+K} (Sig = 't Hooft blocks of the realisation, the
two-component Sigma^I_ij of the dictionary); d_alphadot c <-> lambar_{alphadot 4} = the lam_- components with W^* index 3
(realisation components 11, 15); b <-> (xi, lam*) = xi^T Mp lam* (the lam*_+ component paired with xi).
Entries computed at n = 0:
  beta gamma [BK (3.66)]:  Upsilon_xi(A_perp, lam_+) on gamma^J, beta_I -> lambar_4 lambar_4
  beta beta  [BK (3.67)]:  Upsilon_xi(lam_+, lam_+) on beta_I, beta_J -> A_perp lam_-  (A_perp output on the 6 directions)
  b c (b <-> lam*_+) [BK (3.64)]: Upsilon_xi(c, lam*_+) on c, b -> lambar_4 lambar_4
The finite special vertices vanish on these pairs at n = 0 (weight count, notes.txt).  Colour words as tl_engine; the
BK colour structure (1 (x) ad_Y2 ad_Y1)(-,-)^{-1} is the length-two word."""
import os, sys, json, time
import sympy as sp
import qb_core as B
import qb_job as QJ
import realization as RZ

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')
XI = [0] * 16
XI[3] = 1
LB4 = [11, 15]


def gamma_pol(J):
    """polarisation of gamma^J = A_{J4} on the A_perp components 4..9."""
    return {4 + K: RZ.Sig[K][J - 1, 3] for K in range(6)}


def restrict_outputs(job, P, keep):
    """set the output polarisation symbols outside keep[j] to zero; returns sympy expression."""
    sub = {}
    for nm in job.names:
        if nm.startswith('o') and '_' in nm and nm[1:].split('_')[0].isdigit():
            j, c = nm[1:].split('_')
            j, c = int(j), int(c)
            if keep.get(j) is not None and c not in keep[j]:
                sub[sp.Symbol(nm)] = 0
    e = sp.expand(P.as_expr(*[sp.Symbol(n) for n in job.names]))
    return sp.expand(e.subs(sub))


def run():
    t0 = time.time()
    res = {}
    # beta gamma
    bg = {}
    for I_ in (1, 2, 3):
        for J in (1, 2, 3):
            za = gamma_pol(J)
            zb = {c: (1 if c == I_ - 1 else 0) for c in range(8)}
            job = QJ.DJob('AW', 'lamP', ['lamM', 'lamM'], 'xi', XI, za=za, zb=zb, tag='bg')
            r = job.kernel(mlist=(0, 1))
            ent = {}
            for w, dm in r['total'].items():
                for m, d in dm.items():
                    for b, P in d.items():
                        ent[f'word {list(w)} eps^-{m} {b}'] = str(restrict_outputs(job, P, {0: LB4, 1: LB4}))
            bg[f'gamma^{J} beta_{I_}'] = ent or 'identically 0 (all words, all output components)'
        print('beta gamma', I_, round(time.time() - t0, 1), flush=True)
    res['beta gamma (A_perp, lam_+) -> lamM lamM'] = bg
    # beta beta
    bb = {}
    for I_ in (1, 2, 3):
        for J in (1, 2, 3):
            za = {c: (1 if c == I_ - 1 else 0) for c in range(8)}
            zb = {c: (1 if c == J - 1 else 0) for c in range(8)}
            job = QJ.DJob('lamP', 'lamP', ['AW', 'lamM'], 'xi', XI, za=za, zb=zb, tag='bb')
            r = job.kernel(mlist=(0, 1))
            ent = {}
            for w, dm in r['total'].items():
                for m, d in dm.items():
                    for b, P in d.items():
                        ent[f'word {list(w)} eps^-{m} {b}'] = str(restrict_outputs(job, P, {1: LB4}))
            bb[f'beta_{I_} beta_{J}'] = ent or 'identically 0'
        print('beta beta', I_, round(time.time() - t0, 1), flush=True)
    res['beta beta (lam_+, lam_+) -> A_perp lamM'] = bb
    # b c with b <-> (xi, lam*)
    xv = sp.Matrix(XI)
    zbv = (B.Mp.T * xv)
    zb = {c: zbv[c] for c in range(8, 16)}
    job = QJ.DJob('c', 'lamsP', ['lamM', 'lamM'], 'xi', XI, za={0: 1}, zb=zb, tag='bc')
    r = job.kernel(mlist=(0, 1))
    ent = {}
    for w, dm in r['total'].items():
        for m, d in dm.items():
            for b, P in d.items():
                ent[f'word {list(w)} eps^-{m} {b}'] = str(restrict_outputs(job, P, {0: LB4, 1: LB4}))
                ent[f'word {list(w)} eps^-{m} {b} (all lam_- components)'] = QJ.poly_str(job, P)
    res['b c (c, lam*_+), b = (xi, lam*) -> lamM lamM'] = ent or 'identically 0'
    res['seconds'] = round(time.time() - t0, 1)
    return res


if __name__ == '__main__':
    res = run()
    json.dump(res, open(os.path.join(OUT, 'tq3_twist_n0.json'), 'w'), indent=1)
    print(json.dumps(res, indent=1)[:4000])
