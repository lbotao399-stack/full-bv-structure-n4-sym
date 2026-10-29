"""Fermionic off-shell comparison engine A vs engine B (finding 4 of the QA report).
Engine A (route A, Schwinger/simplex, derived vertex-vertex sign -1) and engine B (route B, first-principles Koszul rule
s(c,sigma)) evaluated at the same off-shell points for every field pair with a gaugino letter or output and every leg set
with at most two outputs (census_sp.json, pair weight without derivatives).  Bases: the exact intertwiner S of
spinor_map.json (u_B = S u_A on Sigma; letters are covectors, zeta_B = S^{-T} zeta_A; output polarisations are vectors,
pol_B = S pol_A).  Under S the +1 eigenspace of rho(e_0 e_1 e_2 e_3) on Sigma (A's label lamP) is B's components 8..15
(B's label lamM): the engines use opposite orientations of V, so the labels are exchanged.
Usage: python3 -B ferm_compare.py [npts]"""
import os, sys, json, time, resource, itertools
import numpy as np
import sympy as sp

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
HERE = os.path.dirname(os.path.abspath(__file__))
ASRC = os.path.join(ROOT, 'step2', 'composite_qsk', 'composite', 'engine_a', 'src.nosync')
BSRC = os.path.join(ROOT, 'step2', 'composite_qsk', 'revise.nosync', 'b_vendor_check', 'src.nosync')   # copy of B's code
sys.path.insert(0, ASRC)
import cq_core as Q, cq_gval as GV, cq_eval as V          # noqa: E402
import ea_core as E                                         # noqa: E402
Q.RULE['koszul'] = 'R2'; GV.INTR['on'] = True
sys.path.insert(0, BSRC)
os.chdir(BSRC)
import ce_job as CJ                                         # noqa: E402

npts = int(sys.argv[1]) if len(sys.argv) > 1 else 1
SM = json.load(open(os.path.join(HERE, 'spinor_map.json')))
C2 = sp.sympify(SM['c2'])                       # S = c S0, c^2 = C2; S0 Gaussian rational
S = sp.Matrix(16, 16, lambda i, j: sp.nsimplify(sp.simplify(sp.sympify(SM['S'][i][j]) / sp.sqrt(C2))))
SinvT = S.inv().T
from realization import pairing as _pairing
_Mp = _pairing()
SST0 = (S.T * _Mp).inv()                          # Sst = c^{-1} SST0
SWAP = {'lamP': 'lamM', 'lamM': 'lamP', 'lamsP': 'lamsM', 'lamsM': 'lamsP'}
ONE = {(0, 0, 0, 0): 1}
RN = {'AV': 4, 'AW': 6, 'lamP': 16, 'lamM': 16, 'c': 1, 'cs': 1, 'AsV': 4, 'AsW': 6, 'lamsP': 16, 'lamsM': 16}
COMPS_B = {'AV': list(range(4)), 'AW': list(range(4, 10)), 'AsV': list(range(4)), 'AsW': list(range(4, 10)), 'c': [0], 'cs': [0]}
G5S = np.asarray(E.G5S); G5SS = np.asarray(E.G5SS)
beta = (1, 2, -1, 3)


def rv(n, rng):
    return np.array([complex(int(x), int(y)) for x, y in zip(rng.integers(-2, 3, n), rng.integers(-2, 3, n))])


def chir(t, z):
    """A: vector of definite chirality (letters and outputs; G5S, G5SS symmetric)"""
    if t in ('lamP', 'lamM'):
        return ((np.eye(16, dtype=int) + (1 if t == 'lamP' else -1) * G5S) @ z)
    if t in ('lamsP', 'lamsM'):
        return ((np.eye(16, dtype=int) + (1 if t == 'lamsP' else -1) * G5SS) @ z)
    return z


def symvec(z):
    return sp.Matrix([sp.nsimplify(complex(x).real) + sp.I * sp.nsimplify(complex(x).imag) for x in z])


def to_B(t, z, kind):
    """A slot vector -> (B type, dict comp -> value)"""
    if t in ('lamP', 'lamM', 'lamsP', 'lamsM'):
        zs = symvec(z)
        if t in ('lamP', 'lamM'):
            zb = (SinvT * zs) if kind == 'letter' else (S * zs)
        else:
            zb = (SST0.inv().T * zs) if kind == 'letter' else (SST0 * zs)
        tb = SWAP[t]
        zb = sp.simplify(zb)
        rng_ = range(8, 16) if tb in ('lamM', 'lamsP') else range(0, 8)
        assert all(sp.simplify(zb[i]) == 0 for i in range(16) if i not in rng_), (t, kind, list(zb))
        return tb, {i: sp.nsimplify(zb[i]) for i in rng_}
    return t, {c: sp.nsimplify(complex(x).real) + sp.I * sp.nsimplify(complex(x).imag) for c, x in zip(COMPS_B[t], z)}


LOG = {'1': 1, 'log2': sp.log(2), 'log3': sp.log(3), 'log5': sp.log(5), 'log7': sp.log(7)}


def Aval(v):
    s = V.vstr(v)
    return sp.expand(sp.sympify(s.replace('i', 'I').replace('+-', '-').replace('log2', 'log(2)').replace('log3', 'log(3)')
                                .replace('log5', 'log(5)').replace('log7', 'log(7)')))


def main():
    cen = json.load(open(os.path.join(ROOT, 'step2', 'composite_qsk', 'composite', 'theory', 'src.nosync', 'census_sp.json')))['pairs']
    pairs = ['AV|lamP', 'AV|lamM', 'AW|lamP', 'AW|lamM', 'lamP|lamP', 'lamP|lamM', 'lamM|lamM', 'AV|AV', 'AV|AW']
    rng = np.random.default_rng(29)
    out = []
    t0 = time.time()
    for pname in pairs:
        ta, tb = pname.split('|')
        ent = cen[pname]['by_letter_derivatives']['0']
        lsP = ent['leg_sets_Pb_and_position']; lsT = ent['leg_sets_translation_invariant']
        for ls in lsP:
            types = ls.split()
            if len(types) > 2:
                continue
            if not any(t.startswith('lam') for t in [ta, tb] + types):
                continue
            for w in sorted(set(itertools.permutations(types))):
                for it in range(npts):
                    za, zb = chir(ta, rv(RN[ta], rng)), chir(tb, rv(RN[tb], rng))
                    moms = [rng.integers(-3, 4, 4) for _ in w]
                    legs = [(t, chir(t, rv(RN[t], rng)), m) for t, m in zip(w, moms)]
                    Xa, Xb = (ta, za, ONE), (tb, zb, ONE)
                    rec = {'pair': pname, 'word_A': list(w)}
                    kinds = [('Pb', ('Pv', beta))] + ([('D0', 'D0')] if ls in lsT else [])
                    vA = {}
                    for nm, kd in kinds:
                        a, b = GV.kernel(Xa, Xb, legs, kd)
                        vA[nm] = {'letter': Aval(a), 'vv': Aval(b), 'value': Aval(GV.combine((a, b), -1))}
                    # engine B
                    tBa, zBa = to_B(ta, za, 'letter'); tBb, zBb = to_B(tb, zb, 'letter')
                    outsB, polB = [], []
                    for (t, p, m) in legs:
                        tt, pp = to_B(t, p, 'out'); outsB.append(tt); polB.append(pp)
                    e = (sum(t in ('lamP', 'lamM') for t in w) - sum(t in ('lamP', 'lamM') for t in (ta, tb))
                         - sum(t in ('lamsP', 'lamsM') for t in w) + sum(t in ('lamsP', 'lamsM') for t in (ta, tb)))
                    assert e % 2 == 0
                    fac = C2 ** (e // 2)
                    job = CJ.CJob(tBa, tBb, outsB, pol=polB, za=zBa, zb=zBb)
                    res = job.kernels(mlist=(0,))
                    subsd = {}
                    for j, m in enumerate(moms):
                        for mu in range(4):
                            subsd[sp.Symbol(f'r{j}_{mu}')] = int(m[mu])
                    for mu in range(4):
                        subsd[sp.Symbol(f'a{mu}')] = beta[mu]
                    rec['types_B'] = [tBa, tBb, outsB]
                    for nm, _ in kinds:
                        bw = {}
                        for wd, dm in res[nm].items():
                            tot = 0
                            for bb, P in dm.get(0, {}).items():
                                tot += sp.expand(P.as_expr().subs(subsd)) * LOG.get(bb, sp.Symbol('L_' + str(bb)))
                            bw[str(wd)] = sp.expand(fac * tot)
                        match = {kap: [wd for wd, v in bw.items() if sp.expand(v - (vA[nm]['letter'] + kap * vA[nm]['vv'])) == 0]
                                 for kap in (-1, 1)}
                        if not bw:
                            match = {kap: (['(B zero)'] if sp.expand(vA[nm]['letter'] + kap * vA[nm]['vv']) == 0 else []) for kap in (-1, 1)}
                        rec[nm] = {'A_letter': str(vA[nm]['letter']), 'A_vv': str(vA[nm]['vv']), 'A_value': str(vA[nm]['value']),
                                   'B_words': {k: str(v) for k, v in bw.items()}, 'match_vvsign_-1': match[-1], 'match_vvsign_+1': match[1]}
                    rec['seconds'] = round(time.time() - t0, 1)
                    rec['rss_MB'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20, 1)
                    ok = all(rec[nm]['match_vvsign_-1'] for nm, _ in kinds)
                    rec['agree'] = bool(ok)
                    print(pname, ls, w, 'agree(-1):', ok, 'agree(+1):', all(rec[nm]['match_vvsign_+1'] for nm, _ in kinds),
                          rec['seconds'], 's', rec['rss_MB'], 'MB', flush=True)
                    out.append(rec)
                    json.dump(out, open(os.path.join(HERE, 'ferm_compare.json'), 'w'), indent=1)
    n = len(out); na = sum(r['agree'] for r in out)
    print('ROWS', n, 'AGREE', na)
    json.dump({'rows': out, 'n_rows': n, 'n_agree': na}, open(os.path.join(HERE, 'ferm_compare.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
