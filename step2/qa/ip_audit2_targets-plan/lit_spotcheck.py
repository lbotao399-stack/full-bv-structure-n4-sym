#!/usr/bin/env python3
"""Spot check of the literature targets used by inner_product_design.tex against the local arXiv copies
(step2/literature/ip_sources/<id>/paper.txt, downloaded from arXiv; version stamps checked in ip_audit2_checks.py A9).
For each target: the equation labels and verbatim key phrases that the design or ip_targets.tex quote, found / not found,
with the line number.  Read-only; writes lit_spotcheck.json next to this file."""
import json, os, re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.abspath(os.path.join(HERE, '..', '..', 'literature', 'ip_sources'))

CHECKS = {
    'BKRS hep-th/0104016': ('hep-th_0104016', ['(9)', '(10)', '(11)', '(17)', '(19)', '(21)', '(24)', '(27)', '(45)'],
                            ['point', 'Fermi-Feynman gauge', 'has no intrinsic normalization', 'a1 = 3N/8π', 'bonus symmetry',
                             'Konishi anomaly [8]', 'choose to work in the latter']),
    'EJSS hep-th/0501077': ('hep-th_0501077', ['(13)', '(15)', '(17)', '(18)', '(19)', '(20)', '(21)'],
                            ['supersymmetric dimensional reduction', 'This may be a scheme-dependent property']),
    'B04 hep-th/0407277': ('hep-th_0407277', ['(2.1)', '(2.2)', '(2.14)', '(2.17)', '(2.19)', '(2.20)', '(2.25)', '(2.39)'],
                           ['the one-loop correction vanishes', 'dimensional']),
    'BKS hep-th/0303060': ('hep-th_0303060', [], ['self-adjoint with respect to the tree-level scalar', 'dimensional regularization',
                                                 'real eigenvalues']),
    'BMV 2306.04693': ('2306.04693', ['(6)', '(7)', '(9)'], ['self-adjoint']),
    'AFGJ hep-th/9708042': ('hep-th_9708042', ['(2.6)', '(2.7)', '(4.15)'], ['determined by the free field content']),
    'EHW hep-th/9905085': ('hep-th_9905085', ['(45)', '(46)', '(47)', '(48)'], ['holds exactly']),
    'Ryz hep-th/0109064': ('hep-th_0109064', ['(13)', '(29)', '(30)', '(31)', '(32)', '(33)', '(34)'],
                           ['simpler equivalent prescription', 'The same result is obtained in dimensional', 'in dimensional regularization']),
    'Mi97 hep-th/9712074': ('hep-th_9712074', ['(2.20)', '(2.24)', '(2.25)', '(2.26)', '(2.27)', '(2.32)', '(2.33)', '(3.2)', '(3.28)', '(3.29)'], []),
    'B03 hep-th/0310252': ('hep-th_0310252', ['(4.11)', '(4.12)', '(4.13)', '(4.14)'],
                           ['will not require an adjoint operation', 'reality constraint']),
    'HSW hep-th/9808162': ('hep-th_9808162', [], ['anomaly considerations determine the coefficient']),
    'BRvR 1304.1803': ('1304.1803', [], ['4(u − v)', 'crossing']),
}

def norm(s):
    return re.sub(r'\s+', ' ', s)

out = {}
for name, (d, eqs, phrases) in CHECKS.items():
    txt = open(os.path.join(SRC, d, 'paper.txt'), errors='ignore').read()
    lines = txt.split('\n')
    flat = norm(txt)
    res = {'equations': {}, 'phrases': {}}
    for e in eqs:
        ln = next((i + 1 for i, l in enumerate(lines) if e in l), None)
        res['equations'][e] = ln
    for p in phrases:
        res['phrases'][p] = (norm(p) in flat) or (norm(p).replace('-', ' ') in flat.replace('-', ' '))
    out[name] = res
missing = {k: [e for e, v in r['equations'].items() if v is None] + [p for p, v in r['phrases'].items() if not v] for k, r in out.items()}
out['_missing'] = {k: v for k, v in missing.items() if v}
json.dump(out, open(os.path.join(HERE, 'lit_spotcheck.json'), 'w'), indent=1)
print(json.dumps(out['_missing'], indent=1))
