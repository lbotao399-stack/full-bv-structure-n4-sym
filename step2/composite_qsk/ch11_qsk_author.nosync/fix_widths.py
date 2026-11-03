from pathlib import Path
R = []
R.append(('part6a.tex',
r"""\[
[v.,Q^{GF}]=4(b,x)\,Q^{GF}+C_K,\qquad [v.,H]=4(b,x)\,H-4(w-1)(b,\partial)+4b_\nu\Sigma^{\rho\nu}\partial_\rho,\qquad C_K:=2\,\sigma^\rho(Q^{GF})\big(w\,b_\rho-b_\nu\Sigma^{\rho\nu}\big),
\]""",
r"""\[
\begin{aligned}
&[v.,Q^{GF}]=4(b,x)\,Q^{GF}+C_K,\qquad C_K:=2\,\sigma^\rho(Q^{GF})\big(w\,b_\rho-b_\nu\Sigma^{\rho\nu}\big),\\
&[v.,H]=4(b,x)\,H-4(w-1)(b,\partial)+4b_\nu\Sigma^{\rho\nu}\partial_\rho,
\end{aligned}
\]"""))
R.append(('part6a.tex',
r"""\[
[v.,D]=-\tfrac12(\operatorname{div}v)D-\tfrac w4\sigma^\rho(D)\,\partial_\rho(\operatorname{div}v)-\tfrac12\sigma^\rho(D)(\partial_\rho\omega_{\mu\nu})\Sigma^{\mu\nu},\qquad [v.,H]=-\tfrac12(\operatorname{div}v)H+\tfrac{w-1}2(\partial_\rho\operatorname{div}v)\partial_\rho+(\partial_\rho\omega_{\mu\nu})\Sigma^{\mu\nu}\partial_\rho,
\]""",
r"""\[
\begin{aligned}
&[v.,D]=-\tfrac12(\operatorname{div}v)D-\tfrac w4\sigma^\rho(D)\,\partial_\rho(\operatorname{div}v)-\tfrac12\sigma^\rho(D)(\partial_\rho\omega_{\mu\nu})\Sigma^{\mu\nu},\\
&[v.,H]=-\tfrac12(\operatorname{div}v)H+\tfrac{w-1}2(\partial_\rho\operatorname{div}v)\partial_\rho+(\partial_\rho\omega_{\mu\nu})\Sigma^{\mu\nu}\partial_\rho,
\end{aligned}
\]"""))
R.append(('part6a.tex',
r"""\hbar L_1(O)^{(1)}(x)=\sum_{1\le a<b\le n}(-1)^{|X_b|(|X_{a+1}|+\cdots+|X_{b-1}|)}\int\omega\,X_1\cdots X_{a-1}\,\Upsilon_x(X_a,X_b)\,X_{a+1}\cdots\widehat{X_b}\cdots X_n-\big\{\hbar F_1^{(1)}(x),O\big\},
\label{eq:c11-decompx}""",
r"""\begin{aligned}
\hbar L_1(O)^{(1)}(x)=\ &\sum_{1\le a<b\le n}(-1)^{|X_b|(|X_{a+1}|+\cdots+|X_{b-1}|)}\int\omega\,X_1\cdots X_{a-1}\,\Upsilon_x(X_a,X_b)\,X_{a+1}\cdots\widehat{X_b}\cdots X_n\\
&-\big\{\hbar F_1^{(1)}(x),O\big\},
\end{aligned}
\label{eq:c11-decompx}"""))
R.append(('part6b.tex',
r"""16\pi^2a_m=\frac1m\int_{\alpha_\sigma\text{ smallest}}\Big[\big(\min_{f\ne\sigma}\alpha_f\big)^m-\alpha_\sigma^m\Big]F_m\,d\alpha\quad(m\ge1),\qquad 16\pi^2a_0=\int_{\alpha_\sigma\text{ smallest}}F_0\,\log\frac{\min_{f\ne\sigma}\alpha_f}{\alpha_\sigma}\,d\alpha,
\label{eq:c11-oddsimplex}""",
r"""\begin{aligned}
16\pi^2a_m&=\frac1m\int_{\alpha_\sigma\text{ smallest}}\Big[\big(\min_{f\ne\sigma}\alpha_f\big)^m-\alpha_\sigma^m\Big]F_m\,d\alpha\qquad(m\ge1),\\
16\pi^2a_0&=\int_{\alpha_\sigma\text{ smallest}}F_0\,\log\frac{\min_{f\ne\sigma}\alpha_f}{\alpha_\sigma}\,d\alpha,
\end{aligned}
\label{eq:c11-oddsimplex}"""))
R.append(('part6c.tex',
r"""\[
\mathcal B_1:=\{S^{(1)}(\xi),O\},\quad \mathcal B_2:=(w_1,w_2)\!\int\!\omega\,\big(\rho(A|_{W_{\mathbb R}^\perp})\xi,\lambda\big),\quad \mathcal B_3:=(w_1,w_2)\!\int\!\omega\,\big(\rho(A|_V)\xi,\lambda\big),\quad \mathcal B_4:=(w_1,w_2)\!\int\!\omega\,\big(c,(\xi,\lambda^*)\big),
\]""",
r"""\[
\begin{aligned}
&\mathcal B_1:=\{S^{(1)}(\xi),O\},\qquad \mathcal B_2:=(w_1,w_2)\int\omega\,\big(\rho(A|_{W_{\mathbb R}^\perp})\xi,\lambda\big),\\
&\mathcal B_3:=(w_1,w_2)\int\omega\,\big(\rho(A|_V)\xi,\lambda\big),\qquad \mathcal B_4:=(w_1,w_2)\int\omega\,\big(c,(\xi,\lambda^*)\big),
\end{aligned}
\]"""))
R.append(('part6c.tex',
r"""\[
s\,Q_1(O)+L_1^{(0)}\big(\{S^{(1)}(\xi),O\}\big)+\big\{S^{(1)}(\xi),L_1^{(0)}(O)\big\}=\big(-4L-2+(2+4L)\big)\,\text{units}\cdot(w_1,w_2)\int\omega\,\big(\rho(dc)\xi,\lambda\big)=0 .
\]""",
r"""\[
\begin{aligned}
&s\,Q_1(O)+L_1^{(0)}\big(\{S^{(1)}(\xi),O\}\big)+\big\{S^{(1)}(\xi),L_1^{(0)}(O)\big\}\\
&\qquad=\big(-4L-2+(2+4L)\big)\,\text{units}\cdot(w_1,w_2)\int\omega\,\big(\rho(dc)\xi,\lambda\big)=0 .
\end{aligned}
\]"""))
R.append(('part6c.tex',
r"""\[
K_1(O)=-\big(\tfrac13+\theta\big)\,\text{units}\cdot\{S^{(1)}(K_b),O\}+\begin{cases}0&\text{(traceless part)},\\ 48\,\text{units}\cdot(b,x)\,O&\text{(trace)},\end{cases}\qquad \{S^{(1)}(K_b),O\}=v^\mu\partial_\mu O-4(b,x)\,O,
\]""",
r"""\[
\begin{aligned}
&K_1(O)=-\big(\tfrac13+\theta\big)\,\text{units}\cdot\{S^{(1)}(K_b),O\}+\begin{cases}0&\text{(traceless part)},\\ 48\,\text{units}\cdot(b,x)\,O&\text{(trace)},\end{cases}\\
&\{S^{(1)}(K_b),O\}=v^\mu\partial_\mu O-4(b,x)\,O,
\end{aligned}
\]"""))
R.append(('part6c.tex',
r"""\[
\sum_{a<b}\pm\Upsilon_\Delta(X_a,X_b)=\sum_{a<b}\pm\big(4\,\mathrm{lc}(X_aX_b)+\text{composite pair kernel along }\Delta\big),\qquad \sum_{a<b}\pm\Upsilon_{P_b}(X_a,X_b)=\sum_{a<b}\pm\big(\text{composite pair kernel along }P_b\big),
\]""",
r"""\[
\begin{aligned}
&\sum_{a<b}\pm\Upsilon_\Delta(X_a,X_b)=\sum_{a<b}\pm\big(4\,\mathrm{lc}(X_aX_b)+\text{composite pair kernel along }\Delta\big),\\
&\sum_{a<b}\pm\Upsilon_{P_b}(X_a,X_b)=\sum_{a<b}\pm\big(\text{composite pair kernel along }P_b\big),
\end{aligned}
\]"""))
R.append(('twocomp_new_items.tex',
r"""\[
Q_1(O)\big|_{\text{two-body}}=-2\log\tfrac32\,\mathcal B_1+4\log\tfrac32\,(\mathcal B_2+\mathcal B_3),\qquad L_1^{(0)}(O)=-\big(2+4\log\tfrac32\big)\tfrac14w_1^{ij}w_{2\,ij}\int\omega\,(\partial_mc,A_m),
\]""",
r"""\[
\begin{aligned}
&Q_1(O)\big|_{\text{two-body}}=-2\log\tfrac32\,\mathcal B_1+4\log\tfrac32\,(\mathcal B_2+\mathcal B_3),\\
&L_1^{(0)}(O)=-\big(2+4\log\tfrac32\big)\tfrac14w_1^{ij}w_{2\,ij}\int\omega\,(\partial_mc,A_m),
\end{aligned}
\]"""))
for f, a, b in R:
    p = Path(f); s = p.read_text()
    assert s.count(a) == 1, (f, a[:60])
    p.write_text(s.replace(a, b))
print('applied', len(R))
