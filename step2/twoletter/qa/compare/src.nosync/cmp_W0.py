"""Entry-by-entry comparison of the pair-weight W_0 tables of engine A (kernels_W0_raw.json) and engine B
(kernels_W0.json), in engine A's component basis.  Levels: (pair, leg set, ordered word, block structure, kind) and
(pair, leg set, ordered word, total over block structures, kind); every letter/output component and every
output-momentum monomial.  Usage: python3 cmp_W0.py [variant] [vecmap] [pair...]"""
import os, sys, json, time, itertools, resource
import numpy as np
import cmp_lib as L

VARIANT = sys.argv[1] if len(sys.argv) > 1 else 'engineA'
VECMAP = sys.argv[2] if len(sys.argv) > 2 else 'D2B'
ONLY = [x for x in sys.argv[3:] if x != 'nokoszul']
KOSZUL = 'nokoszul' not in sys.argv
OUT = os.path.join(L.QA, 'out.nosync')
T = L.transforms(VARIANT, VECMAP)
vm = L.D2B if VECMAP == 'D2B' else list(range(10))

t0 = time.time()
A = json.load(open(os.path.join(L.EA, 'out.nosync', 'kernels_W0_raw.json')))
B = json.load(open(os.path.join(L.EB, 'kernels_W0.json')))


def a_block_tensor(entry, ta, tb, word):
    """engine A entry {idx: polystring} -> {(ia,ib,io..): {mon: c}}"""
    n = len(word)
    out = {}
    for idx, s in entry.items():
        parts = [p for p in idx.split(',')]
        ia, ib = int(parts[0]), int(parts[1])
        io = tuple(int(x) for x in parts[2:] if x != '')
        out[(ia, ib) + io] = L.parse_A(s, n)
    return out


def b_to_A(parsed, ta, tb, outs, seq):
    """B parsed {(name_a, name_b, name_lab0..): {monB}} with A order seq (labels from letter a) -> A-basis dense dict"""
    n = len(outs)
    # dense B array over slot positions
    dims_b = [len(L.SLOT_B[ta][1]), len(L.SLOT_B[tb][1])] + [len(L.SLOT_B[outs[l]][1]) for l in seq]
    monsA = {}
    entries = []
    for key, poly in parsed.items():
        pa = L.NAME2[(ta, key[0])]; pb = L.NAME2[(tb, key[1])]
        po = [L.NAME2[(outs[l], key[2 + l])] for l in seq]
        for monB, c in poly.items():
            monA = [0] * (4 * n)
            for j, l in enumerate(seq):
                for mu in range(4):
                    monA[4 * j + mu] = monB[4 * l + vm[mu]]
            monA = tuple(monA)
            if monA not in monsA:
                monsA[monA] = len(monsA)
            entries.append(((pa, pb, *po), monsA[monA], c))
    if not entries:
        return {}
    arr = np.zeros(dims_b + [len(monsA)], complex)
    for idx, mi, c in entries:
        arr[idx + (mi,)] += c
    mats = [T[('letter', ta)], T[('letter', tb)]] + [T[('output', outs[l])] for l in seq]
    for ax, M in enumerate(mats):
        arr = np.moveaxis(np.tensordot(M, arr, axes=([1], [ax])), 0, ax)
    # Koszul sign: engine A orders the odd output parameters along the chain (seq), engine B by increasing label
    odd = [l for l in seq if L.PARITY[outs[l]]]
    sgn = 1
    for i in range(len(odd)):
        for j in range(i + 1, len(odd)):
            if odd[i] > odd[j]:
                sgn = -sgn
    if KOSZUL:
        arr = sgn * arr
    inv = {v: k for k, v in monsA.items()}
    out = {}
    nz = np.argwhere(np.abs(arr) > 1e-12)
    for ix in nz:
        ix = tuple(int(x) for x in ix)
        out.setdefault(ix[:-1], {})[inv[ix[-1]]] = arr[ix]
    return out


def add_into(D, X, sgn=1):
    for k, p in X.items():
        d = D.setdefault(k, {})
        for m, c in p.items():
            d[m] = d.get(m, 0) + sgn * c


def diff(XA, XB):
    keys = set(XA) | set(XB)
    mx, md, nent = 0.0, 0.0, 0
    worst = None
    for k in keys:
        pa, pb = XA.get(k, {}), XB.get(k, {})
        for m in set(pa) | set(pb):
            a, b = pa.get(m, 0), pb.get(m, 0)
            if abs(a) < 1e-12 and abs(b) < 1e-12:
                continue
            nent += 1
            mx = max(mx, abs(a), abs(b))
            if abs(a - b) > md:
                md = abs(a - b); worst = (k, m, complex(a), complex(b))
    return nent, mx, md, worst


def ratio_fit(XA, XB):
    """best scalar r with XA ~ r XB (detects global sign / factor conventions)"""
    num, den = 0, 0
    for k in set(XA) | set(XB):
        pa, pb = XA.get(k, {}), XB.get(k, {})
        for m in set(pa) | set(pb):
            a, b = pa.get(m, 0), pb.get(m, 0)
            num += np.conj(b) * a; den += abs(b) ** 2
    return complex(num / den) if den else None


def orient(Tp, seq):
    """blocks from letter a for a B topology Tp and A order seq; None if inconsistent"""
    res = set()
    for blocks in (Tp, list(reversed(Tp))):
        pos, ok = 0, True
        for blk in blocks:
            if set(seq[pos:pos + len(blk)]) != set(blk):
                ok = False; break
            pos += len(blk)
        if ok and pos == len(seq):
            res.add(''.join('C' if len(b) == 1 else 'Q' for b in blocks))
    return res


rows = []
TOL = 1e-9
for pt in A:
    if ONLY and pt not in ONLY:
        continue
    ta, tb = A[pt]['letters']
    bpt = B.get(pt)
    legsets = sorted(set(A[pt]['leg_sets']) | set(bpt['leg_sets'] if bpt else []))
    for ls in legsets:
        aL = A[pt]['leg_sets'].get(ls, {})
        bL = bpt['leg_sets'].get(ls) if bpt else None
        outs = bL['outputs'] if bL else (ls.split() if ls != '(none)' else [])
        n = len(outs)
        # engine A: per (word types, block) and totals per word types
        A_blk, A_tot = {}, {}
        for key, d in aL.items():
            wt, blk = [x.strip() for x in key.split('|')]
            word = tuple(wt.split()) if wt != '(none)' else ()
            for kind, entry in d.items():
                X = a_block_tensor(entry, ta, tb, word)
                A_blk[(word, blk, kind)] = X
                add_into(A_tot.setdefault((word, kind), {}), X)
        # engine B
        B_blk, B_tot = {}, {}
        if bL:
            for w, d in bL['total'].items():
                wl = [int(x) for x in w.split(',')] if w else []
                seq = list(reversed(wl))
                word = tuple(outs[l] for l in seq)
                for kind, s in d.items():
                    X = b_to_A(L.parse_B(s, outs, n), ta, tb, outs, seq)
                    B_tot[(w, word, kind)] = X
            for tpo in bL['by_topology']:
                Tp = tpo['topology']
                for w, d in tpo['words'].items():
                    wl = [int(x) for x in w.split(',')] if w else []
                    seq = list(reversed(wl))
                    word = tuple(outs[l] for l in seq)
                    bs = orient(Tp, seq) if Tp else {'0'}
                    assert len(bs) == 1, (pt, ls, Tp, w, bs)
                    blk = bs.pop()
                    for kind, s in d.items():
                        X = b_to_A(L.parse_B(s, outs, n), ta, tb, outs, seq)
                        add_into(B_blk.setdefault((w, word, blk, kind), {}), X)
        # compare totals: every B word against A's word types
        seenA = set()
        for (w, word, kind), XB in B_tot.items():
            XA = A_tot.get((word, kind), {})
            seenA.add((word, kind))
            nent, mx, md, worst = diff(XA, XB)
            rows.append(dict(pair=pt, legset=ls, level='total', Bword=w, Aword=' '.join(word) or '(none)', block='all',
                             kind=kind, n_entries=nent, max_abs=mx, max_diff=md,
                             agree=bool(md <= TOL * max(1, mx)), ratio_A_over_B=None if md <= TOL * max(1, mx) else
                             str(ratio_fit(XA, XB)), worst=None if md <= TOL * max(1, mx) else str(worst)))
        for (word, kind), XA in A_tot.items():
            if (word, kind) not in seenA:
                nent, mx, md, worst = diff(XA, {})
                rows.append(dict(pair=pt, legset=ls, level='total', Bword=None, Aword=' '.join(word) or '(none)',
                                 block='all', kind=kind, n_entries=nent, max_abs=mx, max_diff=md,
                                 agree=bool(md <= TOL * max(1, mx)), ratio_A_over_B=None, worst=str(worst)))
        seenA = set()
        for (w, word, blk, kind), XB in B_blk.items():
            blkA = '0' if blk == '0' else blk
            XA = A_blk.get((word, blkA, kind), {})
            seenA.add((word, blkA, kind))
            nent, mx, md, worst = diff(XA, XB)
            ok = md <= TOL * max(1, mx)
            rows.append(dict(pair=pt, legset=ls, level='block', Bword=w, Aword=' '.join(word) or '(none)', block=blk,
                             kind=kind, n_entries=nent, max_abs=mx, max_diff=md, agree=bool(ok),
                             ratio_A_over_B=None if ok else str(ratio_fit(XA, XB)), worst=None if ok else str(worst)))
        for (word, blk, kind), XA in A_blk.items():
            if (word, blk, kind) not in seenA:
                nent, mx, md, worst = diff(XA, {})
                rows.append(dict(pair=pt, legset=ls, level='block', Bword=None, Aword=' '.join(word) or '(none)',
                                 block=blk, kind=kind, n_entries=nent, max_abs=mx, max_diff=md,
                                 agree=bool(md <= TOL * max(1, mx)), ratio_A_over_B=None, worst=str(worst)))
    print(pt, 'done', round(time.time() - t0, 1), 's', flush=True)

rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20
summ = dict(variant=VARIANT, vecmap=VECMAP, meta=T['_meta'], n_rows=len(rows),
            n_agree=sum(r['agree'] for r in rows), n_disagree=sum(not r['agree'] for r in rows),
            seconds=round(time.time() - t0, 1), peak_rss_MB=round(rss))
tag = f'{VARIANT}_{VECMAP}' + ('' if KOSZUL else '_nokoszul') + ('_' + '_'.join(x.replace('|', '-') for x in ONLY) if ONLY else '')
json.dump(dict(summary=summ, rows=rows), open(os.path.join(OUT, f'cmp_W0_{tag}.json'), 'w'), indent=1)
print(json.dumps(summ))
for r in rows:
    if not r['agree']:
        print('DISAGREE', r['pair'], '|', r['legset'], '|', r['level'], r['Bword'], r['Aword'], r['block'], r['kind'],
              'n', r['n_entries'], 'max', round(r['max_abs'], 4), 'diff', round(r['max_diff'], 4), 'A/B', r['ratio_A_over_B'])
