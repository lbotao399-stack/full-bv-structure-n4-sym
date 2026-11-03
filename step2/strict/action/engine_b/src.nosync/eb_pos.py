"""Engine B, action side: position factors by momentum injection (Prop at-inject) and canonical densities.

A generator with Killing data of position degree one, S(x) = sum_mu x^mu Vdir(d_mu v) + V0, is evaluated by injecting a
momentum r at its vertex; legs 0..n-2 keep their momenta and the last leg carries -sum p_j - r (canonical density, Def
c11-canonical).  For each mu the cycle sum Gamma_mu(r = t e_mu) is a polynomial in t (Lemma at-special), interpolated
exactly on t in {0, +-1, +-2, ...}; then
   intrinsic part  alpha = Gamma_V0 + sum_mu (-i) d/dt Gamma_mu |_{t=0},
   position part   beta_mu = Gamma_mu(t=0)       (the functional int x^mu omega_mu with omega_mu of value beta_mu).
Both are dictionaries colour key -> {m: CLog}."""
from fractions import Fraction as Fr
import eb_core as E
import eb_cycle as Y


def lagrange_deriv0(ts):
    """weights w_t with f'(0) = sum w_t f(t) exactly for polynomials of degree < len(ts)."""
    ts = [Fr(t) for t in ts]
    ws = []
    for i, ti in enumerate(ts):
        # derivative at 0 of the Lagrange basis polynomial L_i
        others = [tj for j, tj in enumerate(ts) if j != i]
        den = Fr(1)
        for tj in others:
            den *= (ti - tj)
        # L_i(t) = prod (t - tj) / den ; L_i'(0) = sum_k prod_{j != k} (0 - tj) / den
        tot = Fr(0)
        for k in range(len(others)):
            pr = Fr(1)
            for jj, tj in enumerate(others):
                if jj != k:
                    pr *= (-tj)
            tot += pr
        ws.append(tot / den)
    return ws


def lagrange_top(ts):
    """weights giving the coefficient of t^{len(ts)-1} (degree check)."""
    ts = [Fr(t) for t in ts]
    ws = []
    for i, ti in enumerate(ts):
        den = Fr(1)
        for j, tj in enumerate(ts):
            if j != i:
                den *= (ti - tj)
        ws.append(1 / den)
    return ws


def add_into(acc, d, w):
    """acc[key][m] += w * d[key][m]; w Gaussian rational pair."""
    for ck, dm in d.items():
        for m, v in dm.items():
            a = acc.setdefault(ck, {})
            a[m] = a.get(m, E.CLog()) + v.scale(w)


def position_split(legs_at, gen_dir, gen_v0, extra_gens=(), tvals=(0, 1, -1, 2, -2), mlist=(0, 1), stats=None,
                   gen_first=True):
    """legs_at(r): legs with injected momentum r (list of 4 Fractions) at the generator vertex.
    gen_dir(mu, r): Gen object for Vdir(e_mu) with injection r; gen_v0: Gen for V0 (or None).
    extra_gens: further generator letters (arity two), placed after (gen_first) or before the position generator.
    Returns (alpha, beta, topcheck) with topcheck the coefficient of t^{len(tvals)-1} (must vanish)."""
    alpha, beta, top = {}, {}, {}
    zero = [Fr(0)] * 4

    def gl(g):
        return [g] + list(extra_gens) if gen_first else list(extra_gens) + [g]
    if gen_v0 is not None:
        r0 = Y.evaluate(legs_at(zero), gl(gen_v0), mlist=mlist, stats=stats)
        add_into(alpha, r0, (Fr(1), Fr(0)))
    wd = lagrange_deriv0(tvals)
    wt = lagrange_top(tvals)
    for mu in range(4):
        for t, w, w2 in zip(tvals, wd, wt):
            r = [Fr(t) if m == mu else Fr(0) for m in range(4)]
            val = Y.evaluate(legs_at(r), gl(gen_dir(mu, r)), mlist=mlist, stats=stats)
            # -i * w
            add_into(alpha, val, (Fr(0), -w))
            add_into(top, val, (w2, Fr(0)))
            if t == 0:
                add_into(beta.setdefault(mu, {}), val, (Fr(1), Fr(0)))
    return alpha, beta, top
