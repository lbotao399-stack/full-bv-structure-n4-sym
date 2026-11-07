"""Build collab/smilezzm/A/report.pdf: a formula-only report (formulas, tables, a flowchart; no prose),
rendered with matplotlib's mathtext (no TeX installation needed).  Usage: python3 make_report.py"""
import os, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ra = json.load(open(os.path.join(ROOT, 'engine_a', 'results.json')))
ga = json.load(open(os.path.join(ROOT, 'engine_a', 'gates.json')))
rb_path = os.path.join(ROOT, 'engine_b', 'results.json')
rb = json.load(open(rb_path)) if os.path.exists(rb_path) else None
cmp_path = os.path.join(HERE, 'comparison.json')
cmp = json.load(open(cmp_path)) if os.path.exists(cmp_path) else None
ref_path = os.path.join(HERE, 'referee', 'findings.json')
ref = json.load(open(ref_path)) if os.path.exists(ref_path) else None

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.size'] = 11

def page(pdf, lines, title=None, fs=12, top=0.95, step=0.055):
    fig = plt.figure(figsize=(8.27, 11.69))
    y = top
    if title:
        fig.text(0.08, y, title, fontsize=14, weight='bold'); y -= 0.06
    for ln in lines:
        if isinstance(ln, tuple):
            txt, size = ln
        else:
            txt, size = ln, fs
        fig.text(0.08, y, txt, fontsize=size); y -= step
        if y < 0.05:
            pdf.savefig(fig); plt.close(fig)
            fig = plt.figure(figsize=(8.27, 11.69)); y = top
    pdf.savefig(fig); plt.close(fig)

def table(pdf, title, header, rows, colw=None, fs=9):
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.08, 0.95, title, fontsize=13, weight='bold')
    ax = fig.add_axes([0.05, 0.08, 0.9, 0.82]); ax.axis('off')
    t = ax.table(cellText=rows, colLabels=header, loc='upper center', cellLoc='left', colWidths=colw)
    t.auto_set_font_size(False); t.set_fontsize(fs); t.scale(1, 1.4)
    pdf.savefig(fig); plt.close(fig)

u = r'u := \frac{\hbar\,\mathrm{tr}_{\mathfrak{g}}(\mathrm{ad}_X \mathrm{ad}_Y)}{16\pi^2 (X,Y)} = -\frac{g^2 N}{32\pi^2}'
out = os.path.join(ROOT, 'report.pdf')
with PdfPages(out) as pdf:
    page(pdf, [
        (r'$O_{ab} = \int \omega\, (\iota_a A, \iota_b A),\quad K = \sum_I O_{e_I e_I},\quad ' + u + '$', 12),
        (r'$S_{12} = (a_1,a_2)(b_1,b_2) + (a_1,b_2)(b_1,a_2),\qquad T_{12} = (a_1,b_1)(a_2,b_2),\qquad \tau = T_{12}/S_{12}$', 12),
        (r'$D_\varepsilon(r) = \frac{1 - e^{-r^2/4\varepsilon}}{4\pi^2 r^2},\quad \hat D_\varepsilon(k) = \frac{e^{-\varepsilon k^2}}{k^2},\quad K_\varepsilon = -\Delta D_\varepsilon,\quad k^2 \hat D_\varepsilon = e^{-\varepsilon k^2}$', 12),
        (r'$F_0 \ni u\,\beta_2 \int (dA|_{W^\perp}, dA|_{W^\perp}),\qquad \psi = 0,\qquad T = \infty,\qquad \ell = 1$', 12),
        ('', 6),
        (r'$\langle O_1, O_2\rangle(x) = \frac{4\hbar\,\mathrm{dim}\mathfrak{g}\, S_{12}}{16\pi^4 |x|^4}\left[1 + u\left(\frac{27}{2} + 4\log\frac{3}{2} - 8\beta_2\right)\right] + \frac{32\hbar\,\mathrm{dim}\mathfrak{g}\, T_{12}}{16\pi^4 |x|^4}\, u\left(\log\frac{|x|^2}{\ell^2} + \gamma_E - 3\log 2\right)$', 11),
        ('', 6),
        (r'$\mathrm{(SE)} = \frac{43}{2} + 8\gamma_E - 24\log 2 + 8L,\qquad \mathrm{(B)} = -8\beta_2,\qquad L = \log|x|^2$', 11),
        (r'$\mathrm{(EX)} = -8 - 4\gamma_E + 8\log 2 + 4\log 3 - 4L + 4\log\varepsilon$', 11),
        (r'$\mathrm{(QU)} = (2\tau - 1)(4L - 4\log\varepsilon + 4\gamma_E - 12\log 2),\qquad \mathrm{(OC)} = -8(1-\tau)\log\varepsilon$', 11),
        (r'$\sum \log\varepsilon = 0,\qquad \sum L = 8\tau,\qquad \sum\gamma_E = 8\tau,\qquad \sum \log 2 = -4 - 24\tau,\qquad \sum \log 3 = 4,\qquad \mathrm{const} = \frac{27}{2}$', 11),
        ('', 6),
        (r'$\Sigma(p) = u\left[\frac{9}{2\varepsilon} - 2p^2\log\varepsilon - p^2\left(2\log p^2 + 2\gamma_E + 2\log 2 - \frac{43}{8}\right)\right]\quad (\mathrm{coefficient\ of\ }(A(p),A(-p)))$', 11),
        (r'$\Sigma^{su(N)}/\frac{g^2N}{16\pi^2}:\ \mathrm{seagull}_V: -4\varepsilon^{-1};\ \mathrm{seagull}_\perp: -5\varepsilon^{-1};\ \mathrm{bubble}: \frac{1}{2\varepsilon} - 2p^2(\log\varepsilon p^2 + \gamma_E + \log 2) + \frac{13}{8}p^2$', 10),
        (r'$\mathrm{gaugino}: 4\varepsilon^{-1} + 4p^2(\log\varepsilon p^2 + \gamma_E + \log 2) - 7p^2;\qquad \Sigma_{ren} = \frac{g^2N}{16\pi^2}p^2(2\log p^2 + 2\gamma_E + 2\log 2 - \frac{43}{8})$', 10),
        ('', 6),
        (r'$M[0,0,2] = 1 - \gamma_E - \log 2 - \log(\varepsilon p^2),\quad M[0,0,3] = \frac{1}{4\varepsilon} + p^2\left(\frac{\log(\varepsilon p^2) + \gamma_E + \log 2}{6} - \frac{25}{72}\right)$', 10),
        (r'$M[1,1,4] = -\frac{\log(\varepsilon p^2) + \gamma_E + \log 2}{6} + \frac{13}{72},\qquad c_{loc} = \int\frac{d^4k}{(2\pi)^4}\frac{e^{-2k^2}(e^{-k^2}-1)}{k^4} = \frac{\log(2/3)}{16\pi^2}$', 10),
        (r'$J = \frac{D_0^2}{8\pi^2}\left(\log\frac{x^2}{\varepsilon} + \gamma_E - 3\log 2\right),\qquad R = \frac{D_0}{16\pi^2}\left(\log\frac{x^2}{\varepsilon} + 1 + \gamma_E - 3\log 2\right)$', 10),
        (r'$(k+k^\prime)\cdot(2q-k-k^\prime) = 2q^2 - k^2 - (q-k)^2 - k^{\prime 2} - (q-k^\prime)^2 + (k-k^\prime)^2$', 10),
        (r'$\mathrm{(EX)} = 32\pi^2 u\,[\,0 - 4(R D_0 + D_0^2 c_{loc}) + J + 2 D_0^2 c_{loc}\,]/D_0^2,\qquad \mathrm{(QU)} = 32\pi^2 u (2\tau - 1) J / D_0^2$', 10),
    ], title='IP1: scalar bilinears at one loop (engine A, exact)')
    page(pdf, [
        (r'$(O_{ab}(x))^{\Theta_I} = (x,x)^{-2} O_{\bar a \bar b}(I(x)),\qquad \langle a, b\rangle_{rad} = \langle a^{\Theta_I,\hbar}, b\rangle$', 12),
        (r'$\langle O_{a_1 b_1}(0), O_{a_2 b_2}(0)\rangle_{rad} = \frac{4\hbar\mathrm{dim}\mathfrak{g}}{16\pi^4}\left[\bar S_{12}\left(1 + u(\frac{27}{2} + 4\log\frac{3}{2} - 8\beta_2)\right) + 8\bar T_{12}\, u(\gamma_E - 3\log 2)\right]$', 11),
        (r'$\mathrm{Gram}_{tree} = \mathrm{Gram}_{tree}^\dagger > 0:\quad \mathrm{diag}(2\ (\times 6),\ 1\ (\times 15))\cdot\frac{4\hbar\mathrm{dim}\mathfrak{g}}{16\pi^4}\ \mathrm{on}\ e_I e_J = \frac{1}{2}\mathrm{Frobenius};\qquad \langle K, \mathbf{20^\prime}\rangle = 0$', 11),
        (r'$\mathrm{Gram}_1 = \mathrm{Gram}_1^\dagger \Leftrightarrow \beta_2 \in \mathbb{R};\qquad \mathrm{Gram}_1 - \mathrm{Gram}_1^\dagger = -16 i\,\mathrm{Im}\beta_2\; u\,\mathrm{Gram}_{tree}/\hbar$', 11),
        ('', 6),
        (r'$J^{\Theta_n} = J,\ \ \psi(O^{\Theta_n}) = \psi(O)^{\Theta_n}\ \Rightarrow\ \langle b, a\rangle_n = \overline{\langle a, b\rangle_n},\qquad \langle a,b\rangle_n := \langle a^{\Theta_n}, b\rangle$', 11),
        (r'$W(P,J^\prime)^{\Theta_n} = W(P, J^{\prime\Theta_n}),\quad J^{\prime\Theta_n} = J + \bar t_1 a + \bar t_2 b^{\Theta_n}\ \Rightarrow\ \langle a, b^{\Theta_n}\rangle = (-1)^{|a||b|}\overline{\langle a,b\rangle_n}$', 10),
        (r'$\langle \Lambda^{(1)}(Y)F, G\rangle = -\langle F, \Lambda^{(1)}(Y) G\rangle,\ Y = P_n - K_n\ \Rightarrow\ \langle a,b\rangle_{rad} = \langle e^{\frac{\pi}{4}\delta^\hbar} a, e^{\frac{\pi}{4}\delta^\hbar} b\rangle_n\ \Rightarrow\ \langle b,a\rangle_{rad} = \overline{\langle a,b\rangle_{rad}}$', 10),
        (r'$\mathrm{Gram} = \overline{\mathrm{Gram}}^T,\qquad (X^{adj})^{adj} = X,\qquad \overline{X}^T\mathrm{Gram} = \mathrm{Gram}\,Y\ \Leftrightarrow\ \overline{Y}^T \mathrm{Gram} = \mathrm{Gram}\,X$', 11),
        ('', 8),
        (r'$\mathrm{flow:}\quad \mathrm{action} \to \mathrm{vertices}\ V_3, V_4,\ \mathrm{gaugino} \to \Sigma(p) \to \Sigma_{ren} \to \mathrm{FT} \to \mathrm{(SE)}$', 10),
        (r'$\qquad\quad\ \ \mathrm{EX} \to \mathrm{numerator\ identity} \to \{2q^2 I_5,\ K_\varepsilon\ \mathrm{contacts}\} \to \mathrm{localisation} \to \{J, R, c_{loc}\} \to \mathrm{(EX)}$', 10),
        (r'$\qquad\quad\ \ \mathrm{QU} \to \hat B(q)^2 \to J \to \mathrm{(QU)};\qquad \mathrm{OC} \to \partial_\delta I^{CT}(O) \to \mathrm{trees}$', 10),
        (r'$\qquad\quad\ \ \mathrm{gates\ G1-G8} \to \mathrm{engine\ A} = \mathrm{engine\ B}\ (\mathrm{per\ class}) \to \mathrm{referee} \to \mathrm{Thm}$', 10),
    ], title='IP0 and Hermiticity')
    # gate table
    rows = []
    for k, v in ga.items():
        if k == 'all_pass':
            continue
        rows.append([k, str(v.get('pass'))])
    table(pdf, 'Engine A gates', ['gate', 'pass'], rows, colw=[0.6, 0.2])
    g1 = ga['G1_selfenergy_singular_parts']['graphs']
    rows = [[k, v['eps_inv_units'], v['log_eps_units'], v['target'][0], v['target'][1]] for k, v in g1.items()]
    table(pdf, 'G1: self-energy singular parts (units)', ['graph', 'eps^-1', 'log eps', 'target eps^-1', 'target log eps'], rows)
    if cmp is not None and 'rows' in cmp:
        table(pdf, 'Engine A vs engine B (per class)', cmp['header'], cmp['rows'], fs=8)
    if ref is not None and 'findings' in ref:
        rows = [[f.get('id', ''), str(f.get('value_found', ''))[:40], str(f.get('value_claimed', ''))[:40], str(f.get('agreement', ''))[:30], f.get('severity', '')] for f in ref['findings']]
        table(pdf, 'Referee: ' + str(ref.get('verdict', {}).get('verdict', '')), ['id', 'found', 'claimed', 'agreement', 'severity'], rows, fs=8)
print('wrote', out)
