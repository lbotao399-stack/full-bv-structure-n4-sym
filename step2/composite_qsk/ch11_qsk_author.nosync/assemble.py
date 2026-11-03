"""Assemble the new chapter 11: Sections 11.1-11.5 unchanged (lines 1-592 of the chapter before this edit) + new 11.6."""
from pathlib import Path
B = Path('ch11_before.tex').read_text().split('\n')
old = lambda a, b: '\n'.join(B[a - 1:b])          # 1-based inclusive
s1to5 = old(1, 592)
assert s1to5.rstrip().endswith('\\end{remark}') and B[592].startswith('\\section{The two-body kernels')
plane = old(664, 671)
plane_add = (" In the evaluation of Theorem \\ref{thm:c11-bilinears} the two-body part on these bilinears, "
             "$-2L\\,\\mathcal B_1+4L(\\mathcal B_2+\\mathcal B_3)$ units, vanishes term by term, with "
             "$\\mathcal B_1=\\{S^{(1)}(\\xi),O\\}=0$ and $(w_1,w_2)=0$, and the one-body term is a multiple of $\\mathcal B_1$.")
assert plane.endswith('\\end{proof}')
plane = plane[:-len('\\end{proof}')].rstrip('\n') + plane_add + '\n\\end{proof}'
twist = old(673, 690).replace('(Proposition \\ref{prop:c11-freedom} (iii))', '(Proposition \\ref{prop:c11-freedom} (iv))')
assert '(Proposition \\ref{prop:c11-freedom} (iv))' in twist
twist_add = (" Theorem \\ref{thm:c11-twistentries} gives the arity-one entries at $\\psi=0$, and Theorem \\ref{thm:c11-twistcheck} "
             "states what their comparison with [BK] tests.")
assert twist.endswith('\\end{remark}')
twist = twist[:-len('\\end{remark}')].rstrip('\n') + twist_add + '\n\\end{remark}'
twocomp = old(697, 715)
items = Path('twocomp_new_items.tex').read_text().rstrip('\n')
k = twocomp.index('\\item Weinberg.')
twocomp = twocomp[:k] + items + '\n' + twocomp[k:]
P = lambda f: Path(f).read_text().rstrip('\n')
tabs = (P('tab_defectpairs.tex') + '\n'
        '\\tabcaption{The leg sets of the defect cycles along $\\xi\\in S_+\\otimes W$ through the pairs of field and ghost letters '
        'without derivatives, with the fewest vertices of $I$ and the number of chains; along $\\xi\\in S_-\\otimes W^*$ the '
        'chiralities of $\\lambda$ and $\\lambda^*$ are exchanged. A leg set enters $\\Upsilon_\\xi$ when its output weight is at most '
        'the bound of the second column, and the empty leg set $\\emptyset$ enters the poles.\\label{tab:c11-defectpairs}}\n\n'
        + P('tab_defectcounts.tex') + '\n'
        '\\tabcaption{The number of leg sets of the defect cycles entering $\\Upsilon_\\xi$, $\\xi\\in S_+\\otimes W$, the '
        'translation-invariant part $\\Upsilon^0_\\eta$, $\\eta\\in S_+\\otimes W^*$, and the translation-invariant part '
        '$\\Upsilon^0_{K_b}$, for $n$ derivatives on the two letters and up to six vertices on the cycle, on the $45$ pairs of slots '
        'with defect cycles; the other ten pairs, of two letters '
        'among $A^*|_{W_{\\mathbb R}^\\perp}$, $\\lambda^*|_{S_+\\otimes W^*}$, $\\lambda^*|_{S_-\\otimes W}$, $c^*$, are those '
        'of Proposition \\ref{prop:c11-oddcensus} (ii).\\label{tab:c11-defectcounts}}')
parts = [s1to5, '', P('part6a.tex'), '', P('part6b.tex'), '', tabs, '', P('part6c.tex'), '', P('part6d.tex'), '',
         plane, '', twist, '', P('part6e.tex'), '', twocomp, '']
out = '\n'.join(parts)
while '\n\n\n' in out:
    out = out.replace('\n\n\n', '\n\n')
Path('ch11_new.tex').write_text(out)
print(len(out), out.count('\\section{'))
