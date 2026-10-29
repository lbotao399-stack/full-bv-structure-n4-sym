"""Comparison of engine B's derivative tables (kernels_deriv.json: W0+1 for every pair type, W0+2 for AV|AV, AV|AW,
AW|AW; leg sets with <= 2 outputs) with engine A evaluated on demand (ea_kern.kernel with letter factors fa, fb; the
engine-A code is imported read-only from engine_a/src.nosync, its block cache is only read).

For each (pair, W, leg set, ordered word, kind) and each split (a, b) of the n letter derivatives, NS random sample
points: integer output momenta, random Gaussian-integer polarisation coefficients in engine A's basis, integer
derivative directions c (letter a), d (letter b).  Engine A: K((c.d)^a X_a, (d.d)^b X_b) summed over the block
structures of the word.  Engine B: a! b! [s^a u^b] of the generating polynomial at s = c, u = d (basis map of
cmp_lib, Koszul sign as in cmp_W0)."""
import os, sys, json, time, math, itertools, resource
from fractions import Fraction
import numpy as np
import cmp_lib as L

sys.path.insert(0, os.path.join(L.EA, 'src.nosync'))
import ea_core as E, ea_kern as K, ea_tables as TB

NS = int(os.environ.get('CMP_NS', '2'))
ONLY = sys.argv[1:]
OUT = os.path.join(L.QA, 'out.nosync')
T = L.transforms('engineA', 'D2B')
vm = L.D2B
rng = np.random.default_rng(20261029)
t0 = time.time()
Bd = json.load(open(os.path.join(L.EB, 'kernels_deriv.json')))
print('loaded', round(time.time() - t0, 1), 's, RSS MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 2 ** 20,
      flush=True)


def parse_Bd(s, n):
    """like cmp_lib.parse_B, plus point-splitting variables s0..s3, u0..u3 (appended after the 4n momenta)"""
    out = {}
    if not s or s.strip() == '0':
        return out
    for term in s.split(' + '):
        m = L._TERM.match(term.strip())
        coef = Fraction(m.group(1)); rest = m.group(2)
        im = False; ia = ib = None; opos = [None] * n; mon = [0] * (4 * n + 8)
        facs = rest[1:].split('*') if rest.startswith('*') else rest.split('*')
        k = 0
        while k < len(facs):
            f = facs[k]; e = 1
            if not f:
                k += 1; continue
            if k + 2 < len(facs) and facs[k + 1] == '':
                e = int(facs[k + 2]); k += 2
            if f == 'I':
                im = not im
            elif f.startswith('Xa['):
                ia = f[3:-1]
            elif f.startswith('Xb['):
                ib = f[3:-1]
            elif f.startswith('Y'):
                lab = int(f[1:f.index('[')]); opos[lab] = f[f.index('[') + 1:-1]
            elif f.startswith('p'):
                lab, mm = f[1:].split('_'); mon[4 * int(lab) + int(mm)] += e
            elif f[0] in 'su' and f[1:].isdigit():
                mon[4 * n + (0 if f[0] == 's' else 4) + int(f[1:])] += e
            else:
                raise ValueError(f)
            k += 1
        val = complex(coef) * (1j if im else 1)
        d = out.setdefault((ia, ib) + tuple(opos), {})
        d[tuple(mon)] = d.get(tuple(mon), 0) + val
    return out


def b_eval(parsed, ta, tb, outs, seq, pols, moms, c, d, adeg, bdeg):
    """engine-B value in engine-A letter basis: matrix [ia, ib]"""
    n = len(outs)
    odd = [l for l in seq if L.PARITY[outs[l]]]
    sgn = 1
    for i in range(len(odd)):
        for j in range(i + 1, len(odd)):
            if odd[i] > odd[j]:
                sgn = -sgn
    Ta, Tb = T[('letter', ta)], T[('letter', tb)]
    # output polarisation in B slot coordinates: pol_A coefficients x -> B field coefficient vector y with
    # sum_io x_io T_A[io] = sum_j y_j B[j]  =>  y = x^T T_out
    ys = {}
    for j, l in enumerate(seq):
        ys[l] = np.asarray(pols[j]) @ T[('output', outs[l])]
    # B momenta / directions in B components: p_{lab, vm[mu]} = mom_A[j][mu]
    pB = {}
    for j, l in enumerate(seq):
        for mu in range(4):
            pB[4 * l + vm[mu]] = moms[j][mu]
    sB = [0] * 4; uB = [0] * 4
    for mu in range(4):
        sB[vm[mu]] = c[mu]; uB[vm[mu]] = d[mu]
    na, nb = len(L.SLOT_B[ta][1]), len(L.SLOT_B[tb][1])
    Mb = np.zeros((na, nb), complex)
    for key, poly in parsed.items():
        pa = L.NAME2[(ta, key[0])]; pb = L.NAME2[(tb, key[1])]
        f = 1
        for l in range(n):
            f = f * ys[l][L.NAME2[(outs[l], key[2 + l])]]
        if f == 0:
            continue
        tot = 0
        for mon, cf in poly.items():
            sd = sum(mon[4 * n:4 * n + 4]); ud = sum(mon[4 * n + 4:])
            if sd != adeg or ud != bdeg:
                continue
            v = cf
            for i in range(4 * n):
                if mon[i]:
                    v *= pB[i] ** mon[i]
            for m in range(4):
                if mon[4 * n + m]:
                    v *= sB[m] ** mon[4 * n + m]
                if mon[4 * n + 4 + m]:
                    v *= uB[m] ** mon[4 * n + 4 + m]
            tot += v
        Mb[pa, pb] += f * tot
    return sgn * math.factorial(adeg) * math.factorial(bdeg) * (Ta @ Mb @ Tb.T)


def dirpoly(c, deg):
    """(i c.p)^deg as {exponent tuple: Gaussian int}"""
    out = {(0, 0, 0, 0): 1}
    for _ in range(deg):
        new = {}
        for e, v in out.items():
            for mu in range(4):
                if c[mu]:
                    ee = list(e); ee[mu] += 1; ee = tuple(ee)
                    new[ee] = new.get(ee, 0) + v * 1j * c[mu]
        out = new
    return out


def a_eval(ta, tb, word, pols, moms, c, d, adeg, bdeg, mmax):
    """engine A, summed over the block structures of the word: {kind: matrix [ia, ib]}"""
    va, vb = TB.letter_vecs(ta), TB.letter_vecs(tb)
    res = {}
    for struct in (TB.block_structures(list(word)) if word else [[]]):
        blocks = []
        for blk in struct:
            if len(blk) == 1:
                j = blk[0]
                pv = sum(pols[j][i] * TB.stripped_pol(word[j], i) for i in range(len(pols[j])))
                blocks.append(('cub', word[j], pv, np.array(moms[j], dtype=int)))
            else:
                j, jj = blk
                p1 = sum(pols[j][i] * TB.stripped_pol(word[j], i) for i in range(len(pols[j])))
                p2 = sum(pols[jj][i] * TB.stripped_pol(word[jj], i) for i in range(len(pols[jj])))
                blocks.append(('qua', (p1, np.array(moms[j], dtype=int)), (p2, np.array(moms[jj], dtype=int))))
        r = K.kernel(blocks, ta, tb, kinds=('log', 'eps') if mmax else ('log',), mmax=mmax or None,
                     fa=dirpoly(c, adeg), fb=dirpoly(d, bdeg))
        if not r:
            continue
        for kind, v in r.items():
            if v is None:
                continue
            M = np.zeros((len(va), len(vb)), complex)
            for ia, za in enumerate(va):
                for ib, zb in enumerate(vb):
                    x = K.contract(v, za, zb)
                    M[ia, ib] = complex(float(x[0]), float(x[1]))
            res[kind] = res.get(kind, 0) + M
    return res


rows = []
for key, entry in Bd.items():
    pair = key.split(' @ ')[0]
    if ONLY and key not in ONLY and pair not in ONLY:
        continue
    ta, tb = entry['letters']
    nder = entry['derivatives_on_letters']
    W = Fraction(entry['W'])
    for ls, x in entry['leg_sets'].items():
        outs = x['outputs']; n = len(outs)
        wout = sum((E.WEYL[t] for t in outs), Fraction(0))
        mmax = int((W - wout) // 2) if W - wout >= 0 else 0
        Bw = x['total']
        # all words: B words plus every ordering of the outputs (zero on B's side if absent)
        allw = set(Bw)
        for perm in itertools.permutations(range(n)):
            allw.add(','.join(map(str, perm)))
        if n == 0:
            allw = {''} | set(Bw)
        for w in sorted(allw):
            wl = [int(v) for v in w.split(',')] if w else []
            seq = list(reversed(wl))
            word = tuple(outs[l] for l in seq)
            parsed = {kd: parse_Bd(s, n) for kd, s in Bw.get(w, {}).items()}
            worst, mx, npts = 0.0, 0.0, 0
            fails = []
            for smp in range(NS):
                moms = [tuple(int(v) for v in rng.integers(-3, 4, 4)) for _ in range(n)]
                pols = [tuple(complex(int(a), int(b)) for a, b in zip(rng.integers(-2, 3, L.NA[t]),
                                                                       rng.integers(-2, 3, L.NA[t]))) for t in word]
                c = tuple(int(v) for v in rng.integers(-2, 3, 4)); d = tuple(int(v) for v in rng.integers(-2, 3, 4))
                # a zero direction gives an empty letter polynomial, which ea_kern.kernel reads as 'no factor': avoid it
                while not any(c):
                    c = tuple(int(v) for v in rng.integers(-2, 3, 4))
                while not any(d):
                    d = tuple(int(v) for v in rng.integers(-2, 3, 4))
                for adeg in range(nder + 1):
                    bdeg = nder - adeg
                    RA = a_eval(ta, tb, word, pols, moms, c, d, adeg, bdeg, mmax)
                    kinds = set(k for k in RA) | set(parsed)
                    for kd in kinds:
                        kA = kd if kd in RA else None
                        if isinstance(kd, str) and kd.startswith('eps^-'):
                            kA = int(kd[5:])
                        elif isinstance(kd, int):
                            kA = kd
                        kB = 'log' if kd == 'log' else ('eps^-%d' % (kd if isinstance(kd, int) else int(kd[5:])))
                        MA = RA.get(kA, 0)
                        MB = b_eval(parsed[kB], ta, tb, outs, seq, pols, moms, c, d, adeg, bdeg) if kB in parsed else 0
                        D = np.abs(np.asarray(MA) - np.asarray(MB))
                        dd = float(np.max(D)) if np.size(D) else 0.0
                        if dd > 1e-8 * max(1.0, float(np.max(np.abs(MA))) if np.size(MA) else 1.0):
                            fails.append(dict(sample=smp, moms=moms, c=c, d=d, adeg=adeg, bdeg=bdeg, kind=str(kd),
                                              diff=dd, ratio=str(complex(np.vdot(np.ravel(MB), np.ravel(MA)) /
                                                                         max(1e-300, np.vdot(np.ravel(MB), np.ravel(MB)).real)))
                                              if np.size(MB) > 1 else None))
                        worst = max(worst, dd)
                        mx = max(mx, float(np.max(np.abs(MA))) if np.size(MA) else 0.0,
                                 float(np.max(np.abs(MB))) if np.size(MB) else 0.0)
                        npts += 1
            ok = worst <= 1e-8 * max(1.0, mx)
            rows.append(dict(table=key, legset=ls or '(none)', Bword=w, Aword=' '.join(word) or '(none)',
                             in_B=w in Bw, samples=npts, max_abs=mx, max_diff=worst, agree=bool(ok), fails=fails))
            if not ok:
                print('DISAGREE', key, ls, w, ' '.join(word), mx, worst, fails, flush=True)
    print(key, 'done', round(time.time() - t0, 1), 's RSS MB',
          resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 2 ** 20, flush=True)

summ = dict(n_rows=len(rows), n_agree=sum(r['agree'] for r in rows), n_nonzero=sum(r['max_abs'] > 0 for r in rows),
            samples_per_word=NS, seconds=round(time.time() - t0, 1),
            peak_rss_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 2 ** 20)
tag = ('_' + '_'.join(x.replace('|', '-').replace(' ', '') for x in ONLY)) if ONLY else ''
json.dump(dict(summary=summ, rows=rows), open(os.path.join(OUT, f'cmp_deriv{tag}.json'), 'w'), indent=1)
print(json.dumps(summ))
