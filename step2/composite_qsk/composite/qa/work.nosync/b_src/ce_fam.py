"""Families: planar two-site reading of the composite kernels and the matrices of Conjecture P7 (b), (c).

Planar adjacent rule (Lemma c10-colour (iv), twoletter colour_rules.py): a kernel term c hbar/(16 pi^2) C_w with
C_w = (1 (x) ad_{Y_w1}...ad_{Y_wr})(-,-)^{-1} replaces the adjacent pair (X_a X_b) of a single trace by
(-1)^r (c/2) Y_wr ... Y_w1 in units u.  A labelled output form with two identical outputs is a functional after
division by 2.  Letter derivatives in the null mode: K((n.d)^a X, (n.d)^b Y) = a! b! [ss^a su^b]; an output
(n.d)^c Y_j corresponds to (i t_j)^c."""
import os, sys, math, json, time
from fractions import Fraction as Fr
import sympy as sp
from sympy.polys.domains import QQ_I
HERE = os.path.dirname(os.path.abspath(__file__))
import ce_job as CJ

I = sp.I


def qq(x):
    return sp.Rational(int(x.numerator), int(x.denominator))


def cplx(c):
    return qq(c.x) + I * qq(c.y)


def a_subst(job, P, avec):
    """substitute the line-factor direction a = avec (list of 4 sympy numbers) in a ring element."""
    R = job.R
    imgs = []
    for nm, gen in zip(job.names, R.gens):
        if nm.startswith('a') and nm[1:].isdigit():
            imgs.append((gen, R(QQ_I.from_sympy(sp.nsimplify(avec[int(nm[1:])])))))
    return P.compose(imgs) if imgs else P


def null_entries(job, P, nout, smax_in, extra=0):
    """P: ring element (null mode, a substituted).  Returns {(outs..., ins...): complex} with
    outs = derivative orders on the outputs (in label order), ins = (a, b) letter derivative orders,
    coefficient converted: x a! b! / i^{sum outs}; requires no w dependence (checked, returned flag)."""
    nm = job.names
    iss, isu = nm.index('ss'), nm.index('su')
    it = [nm.index(f't{j}') for j in range(nout)]
    iw = [nm.index(f'w{j}') for j in range(nout) if f'w{j}' in nm]
    out = {}
    wfree = True
    for mon, cf in P.items():
        if any(mon[i] for i in iw):
            wfree = False
            continue
        a, b = mon[iss], mon[isu]
        if a + b > smax_in:
            continue
        cs = tuple(mon[i] for i in it)
        v = cplx(cf) * math.factorial(a) * math.factorial(b)
        v = v / I ** sum(cs)
        out[(cs, (a, b))] = out.get((cs, (a, b)), 0) + v
    return out, wfree


def two_output_density(job, tab, which_m=0, avec=None, smax=4):
    """two identical outputs (Z): planar density {((c,d),(a,b)): units u} with sites (c,d) in trace order,
    from tab = res['Pb'] or res['D0'] (word -> m -> basis -> poly).  Factor 1/2 planar x 1/2 symmetrisation."""
    dens = {}
    logs_ok = True
    wfree_all = True
    for w, dm in tab.items():
        d = dm.get(which_m, {})
        for basis, P in d.items():
            if avec is not None:
                P = a_subst(job, P, avec)
            ent, wfree = null_entries(job, P, 2, smax)
            wfree_all &= wfree
            for (cs, ins), v in ent.items():
                v = sp.nsimplify(v)
                if v == 0:
                    continue
                if basis != '1':
                    logs_ok = False
                    continue
                c, dd = cs
                site = (dd, c) if w == (0, 1) else (c, dd)
                key = (site, ins)
                dens[key] = dens.get(key, 0) + v * sp.Rational(1, 4) * (-1) ** len(w)
    return {k: v for k, v in dens.items() if v != 0}, logs_ok, wfree_all


def one_output_density(job, tab, which_m=0, avec=None, smax=4):
    """one output: {((c,), (a,b)): units u}, word length 1 -> (-1)^1 c/2."""
    dens = {}
    logs_ok = True
    wfree_all = True
    for w, dm in tab.items():
        d = dm.get(which_m, {})
        for basis, P in d.items():
            if avec is not None:
                P = a_subst(job, P, avec)
            ent, wfree = null_entries(job, P, 1, smax)
            wfree_all &= wfree
            for (cs, ins), v in ent.items():
                v = sp.nsimplify(v)
                if v == 0:
                    continue
                if basis != '1':
                    logs_ok = False
                    continue
                dens[(cs, ins)] = dens.get((cs, ins), 0) + v * sp.Rational(1, 2) * (-1) ** len(w)
    return {k: v for k, v in dens.items() if v != 0}, logs_ok, wfree_all


def twist_two_map(h, AZ, ZA, a, b):
    """composite (or log) kernel summed over the pairs of O(a,b) = tr((n.D)^a Z (n.D)^b Z), two-field part:
    returns {(c,d): coefficient} in the ordered basis tr(d^c Z d^d Z).  h: ((c,d),(x,y)) two-output density;
    AZ[(cs, (r, s))]: one output for the pair (d^r A_n, d^s Z) [slot a = A_n]; ZA[(cs,(s, r))]: pair (d^s Z, d^r A_n)."""
    res = {}

    def add(c, d, v):
        if v != 0:
            res[(c, d)] = res.get((c, d), 0) + v

    for (site, ins), v in h.items():
        if ins == (a, b):
            add(site[0], site[1], v)
        if ins == (b, a):
            add(site[1], site[0], v)

    def kAZ(r, s_):     # pair (d^r A, d^s Z) -> {c: v}
        return {cs[0]: v for (cs, ins), v in AZ.items() if ins == (r, s_)}

    def kZA(r, s_):     # pair (d^s Z, d^r A): ZA keyed with ins = (s, r)
        return {cs[0]: v for (cs, ins), v in ZA.items() if ins == (s_, r)}

    # (n.D)^a Z = d^a Z + sum_r C(a, r+1) [d^r A_n, d^{a-1-r} Z] + ...
    for r in range(a):
        s_ = a - 1 - r; cf = math.comb(a, r + 1)
        # tr([A_r, Z_s] Z_b) = tr(A Zs Zb) - tr(Zs A Zb)
        for c, v in kAZ(r, s_).items():
            add(c, b, cf * v)                       # (A,Zs) adjacent -> tr(Y Zb)
        for c, v in kZA(r, b).items():
            add(c, s_, cf * v)                      # (Zb, A) adjacent (cyclic) -> tr(Y Zs)... order: tr(Zs Y)?
        for c, v in kZA(r, s_).items():
            add(c, b, -cf * v)                      # (Zs, A) -> tr(Y Zb)
        for c, v in kAZ(r, b).items():
            add(s_, c, -cf * v)                     # (A, Zb) -> tr(Zs Y)
    for r in range(b):
        s_ = b - 1 - r; cf = math.comb(b, r + 1)
        # tr(Za [A_r, Z_s]) = tr(Za A Zs) - tr(Za Zs A)
        for c, v in kZA(r, a).items():
            add(c, s_, cf * v)                      # (Za, A) -> tr(Y Zs)
        for c, v in kAZ(r, s_).items():
            add(a, c, cf * v)                       # (A, Zs) -> tr(Za Y)
        for c, v in kZA(r, s_).items():
            add(a, c, -cf * v)                      # (Zs, A) -> tr(Za Y)
        for c, v in kAZ(r, a).items():
            add(s_, c, -cf * v)                     # (A, Za) cyclic -> tr(Zs Y)
    return res
