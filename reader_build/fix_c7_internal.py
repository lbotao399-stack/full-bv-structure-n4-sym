#!/usr/bin/env python3
"""One-off: replace the overstated sentence of Remark rem:c7-internal in a reader document (revision-guarded).
Usage: fix_c7_internal.py PORT [--dry-run]. Changes two sentence blocks of origin b-01282/b-01284 and inserts one
new sentence block after the first; every other block must stay identical."""
import hashlib, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import merge_into_live as mil

UNIT = r'\frac{\hbar\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{16\pi^2(X,Y)}'
S2 = ('<p>By Proposition <a href="#b-01244">7.6</a> every kernel, vertex and counterterm entering them is $x_1$-covariant, '
      'but the internal-insertion graphs still contribute: for a translation $x_1=P_b$, on the part of the scalar bilinears '
      r'with two legs $A|_{W_{\mathbb R}^\perp}$, they give $-\big(\tfrac13+\theta_P\big)' + UNIT + r'\{S^{(1)}(P_b),O\}$ '
      r'($\theta_P$ as in Theorem <a href="#c9-01889">9.33</a>), which the finite part of the operators of Proposition '
      r'<a href="#c9-01895">9.34</a> removes; along $\mathfrak{so}(W_{\mathbb R}^\perp)$ the composite kernel vanishes, '
      'since $v.$ acts trivially on the slot of $c$. </p>')
S3 = r'<p>When these graphs vanish or are removed in this way, for $O$ with $(Q+\{I,-\})O=0$,</p>' + '\n'
S4 = (r'<p>The arity-one counterterms along $x_1$ and the internal-insertion graphs are evaluated in Chapter 9 along the '
      r'dilatation and the translations; along $\mathfrak{so}(V)$ they remain to be evaluated.</p>' + '\n')
A, B = 'sentence-8f1edb9425823d2cf3d04581', 'sentence-4c2ac769d16668706e72c202'


def main():
    port = int(sys.argv[1]); dry = '--dry-run' in sys.argv
    st, state = mil.req(port, 'GET'); assert st == 200
    rev, doc = state['revision'], state['document']
    root = Path(__file__).resolve().parents[1]
    (root / f'reader/qa/backup_port{port}_rev{rev:06d}_before_c7_internal.json').write_text(json.dumps(state, ensure_ascii=False))
    old = doc['blocks']; ids = {b['id']: b for b in old}
    assert 'this is the reason these components' in ids[A]['source'] and 'which realise this choice' in ids[B]['source']
    assert not ids[A].get('formatting') and not ids[B].get('formatting')
    assert ids['c9-01889']['source'].find('translation component') >= 0 or 'P_1' in json.dumps(ids['c9-01889'])
    new_id = 'sentence-' + hashlib.sha256(('b-01282\0' + '3\0' + S3).encode()).hexdigest()[:24]
    blocks = []
    for b in old:
        if b['id'] == A:
            a = dict(b); a['source'] = S2; blocks.append(a)
            blocks.append({'id': new_id, 'type': 'paragraph', 'source': S3, 'parentId': b['parentId'],
                           'sentenceSplit': {'version': 1, 'originId': 'b-01282', 'part': 3}})
        elif b['id'] == B:
            c = dict(b); c['source'] = S4; blocks.append(c)
        else:
            blocks.append(b)
    others = [b for b in blocks if b['id'] not in (A, B, new_id)]
    assert others == [b for b in old if b['id'] not in (A, B)]
    print({'revision': rev, 'blocks_before': len(old), 'blocks_after': len(blocks), 'new_id': new_id})
    if dry:
        return
    doc = dict(doc); doc['blocks'] = blocks
    st, res = mil.req(port, 'PUT', {'document': doc, 'revision': rev}); assert st == 200, (st, res)
    st, back = mil.req(port, 'GET')
    print({'revision_after': res['revision'], 'readback_ok': back['document']['blocks'] == blocks})


if __name__ == '__main__':
    main()
