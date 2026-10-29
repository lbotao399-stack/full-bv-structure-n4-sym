"""Audit #3 of the Part D theory (step2/theory/delta_theory.tex), item: the composite triangle Y1a of engine A.

Claim under test (engine A, delta_results.json task2_component; theory prop:d-T (iii), rem:d-Y1a):
  Q(Y1a) is a nonlocal Q-exact term, so [(d_CE+Q)O[infty]](Delta) is not of the form (local)[infty].

Engine A reports Y1a in its canonical form (field leg A|_Wperp at p1 fixed, dependent leg A^*|_Wperp at p2 = -p1-k):
  alpha(Y1a) = 0,  beta(Y1a) = unit T_O [ 8 p2/|p2|^2 + 8 |p1|^2 That(p1,p2) ]           (exact fit, 3 points)
  That_nu(a,b) = 16 pi^2 int d^4l/(2pi)^4 l_nu / (l^4 (l+b)^2 (l-a)^2) = - int_simplex x (y b - z a)_nu / Delta^2,
  Delta = y|b|^2 + z|a|^2 - |y b - z a|^2,
and Q(Y1a) as the kernel with the A^* leg replaced by its Q-image (-1/2 |p2|^2):
  alpha_Q = p2.beta,   beta_Q = -1/2 |p2|^2 beta.
After Q BOTH legs are A|_Wperp.  A functional of one even field sees only the part of the kernel that is symmetric
under exchange of the two legs.  For the canonical form D = alpha(p1) delta(S) + beta(p1).d delta(S), S = p1+p2+k,
the exchanged kernel is (f(y) delta'(y-a) = f(a) delta'(y-a) - f'(a) delta(y-a)):
  alpha^T(p1) = alpha(p2) - (div beta)(p2),   beta^T(p1) = beta(p2),
so the vertex value (second derivative) of the functional is VV = D + D^T.
With beta_Q(q) = 4(q+k) - 4|q+k|^2|q|^2 J(q),  J(q) := That(q,-q-k):
  div beta_Q(q) = 16 - 8 sum_nu [(q+k)_nu |q|^2 + |q+k|^2 q_nu] J_nu(q) - 4 |q+k|^2 |q|^2 div J(q).
Checks:
  (1) That(b,a) = -That(a,b)                                  (l -> -l)
  (2) div J(p2) = 2/(|p1|^2 |p2|^2)                           (d_nu(l_nu/l^4) = 2 pi^2 delta^4(l)), finite differences
  (3) VV[Q(Y1a)] = (8, 4k) unit T_O = -2 x VV[{S^(1)(Delta),O}],  VV[{S^(1)(Delta),O}] = (-4, -2k) T_O,
      while the antisymmetric part D - D^T varies between kinematic points (nonlocal): it is the term engine A saw.
  (4) consistency of the formula with engine A's raw alpha_Q / T_O at the logged point.
  (5) translation analogue (engine A: beta(Y1a) = -i x translation graph): symmetric part of Q(Y1a_v) = 4 i k.v.
"""
import json, os
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 20

PTS = [((1, 2, 0, 1), (0, -1, 1, 2)), ((2, -1, 1, 0), (1, 1, -2, 1)), ((1, 1, 1, -1), (-2, 0, 1, 1))]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def That(a, b, comps=range(4)):
    """engine A's Feynman-parameter form (d4_y1a_ward.py), arbitrary real momenta"""
    a = [mp.mpf(x) for x in a]; b = [mp.mpf(x) for x in b]
    aa, bb, ab = dot(a, a), dot(b, b), dot(a, b)
    out = []
    for nu in comps:
        def f(x, t, nu=nu):
            y = (1 - x) * t; z = (1 - x) * (1 - t)
            D = y * bb + z * aa - (y * y * bb - 2 * y * z * ab + z * z * aa)
            return -x * (y * b[nu] - z * a[nu]) / D ** 2 * (1 - x)
        out.append(mp.quad(f, [0, 1], [0, 1]))
    return out


def J(q, k, comps=range(4)):
    return That(q, [-q[m] - k[m] for m in range(4)], comps)


def divJ(q, k, h):
    s = mp.mpf(0)
    for nu in range(4):
        qp = list(q); qm = list(q); qp[nu] += h; qm[nu] -= h
        s += (J(qp, k, [nu])[0] - J(qm, k, [nu])[0]) / (2 * h)
    return s


def main():
    out = {'claim': 'Q(Y1a) nonlocal (engine A) vs symmetric part of the kernel', 'points': []}
    for (p1, k) in PTS:
        p1 = [mp.mpf(x) for x in p1]; k = [mp.mpf(x) for x in k]
        p2 = [-p1[m] - k[m] for m in range(4)]
        n1, n2 = dot(p1, p1), dot(p2, p2)
        r = {'p1': [int(x) for x in p1], 'k': [int(x) for x in k], 'p2': [int(x) for x in p2]}
        T12 = That(p1, p2); T21 = That(p2, p1)                   # J(p1) = T12, J(p2) = T21
        r['max|That(p2,p1)+That(p1,p2)|'] = mp.nstr(max(abs(T12[m] + T21[m]) for m in range(4)), 5)
        dJa = divJ(p2, k, mp.mpf('1e-4')); dJb = divJ(p2, k, mp.mpf('5e-5'))
        r['divJ(p2)*|p1|^2|p2|^2 (h=1e-4, 5e-5)'] = [mp.nstr(dJa * n1 * n2, 14), mp.nstr(dJb * n1 * n2, 14)]
        r['target'] = '2'
        dJ = (4 * dJb - dJa) / 3                                  # Richardson
        # unsymmetrized kernel of Q(Y1a), first leg p1, in units unit*T_O
        be1 = [8 * p2[m] / n2 + 8 * n1 * T12[m] for m in range(4)]
        aQ1 = dot(p2, be1); bQ1 = [-n2 / 2 * x for x in be1]
        # same kernel with first leg p2 (second leg p1)
        be2 = [8 * p1[m] / n1 + 8 * n2 * T21[m] for m in range(4)]
        aQ2 = dot(p1, be2); bQ2 = [-n1 / 2 * x for x in be2]
        # divergence of beta_Q at q = p2 (q + k = -p1)
        dvb = 16 - 8 * sum(((-p1[m]) * n2 + n1 * p2[m]) * T21[m] for m in range(4)) - 4 * n1 * n2 * dJ
        aT, bT = aQ2 - dvb, bQ2
        VVa = aQ1 + aT; VVb = [bQ1[m] + bT[m] for m in range(4)]
        ASa = aQ1 - aT; ASb = [bQ1[m] - bT[m] for m in range(4)]
        r['unsym alpha_Q/(unit T_O)'] = mp.nstr(aQ1, 12)
        r['unsym beta_Q/(unit T_O)'] = [mp.nstr(x, 12) for x in bQ1]
        r['VV alpha/(unit T_O)'] = mp.nstr(VVa, 12)
        r['VV beta/(unit T_O)'] = [mp.nstr(x, 12) for x in VVb]
        r['4k'] = [mp.nstr(4 * x, 5) for x in k]
        r['max |VV - (8,4k)|'] = mp.nstr(max([abs(VVa - 8)] + [abs(VVb[m] - 4 * k[m]) for m in range(4)]), 5)
        r['antisym alpha/(unit T_O)'] = mp.nstr(ASa, 12)
        r['antisym beta/(unit T_O)'] = [mp.nstr(x, 12) for x in ASb]
        # translation analogue: kernel -1/2|p2|^2 (i beta.v), no delta' term; symmetric part times 2
        sv = [-n2 / 2 * be1[m] - n1 / 2 * be2[m] for m in range(4)]
        r['translation VV/(i unit T_O v) (expect 4k)'] = [mp.nstr(x, 12) for x in sv]
        out['points'].append(r)
        print(json.dumps(r), flush=True)
    a_form = out['points'][0]['unsym alpha_Q/(unit T_O)']
    out['engine A raw alphaQ/T_O at point 1 (logs/d4_y1a_w44_su2.json)'] = '-2.10031347566'
    out['formula alphaQ/T_O at point 1'] = a_form
    out['VV[{S^(1)(Delta),O}]'] = '(-4, -2k) T_O  (2(E+2) on e^{ikx}: -(2 + k.d_k) x 2)'
    out['verdict'] = ('As a functional of the single field A|_Wperp, Q(Y1a) is local: vertex value (8, 4k) unit T_O = '
                      '-2 unit x vertex value of {S^(1)(Delta),O}.  The nonlocal part reported by engine A is the part of '
                      'the kernel antisymmetric under exchange of the two A|_Wperp legs (That(b,a) = -That(a,b)); it '
                      'contributes to no functional.  Q(Y1a) = -2 unit {S^(1)(Delta),O}.')
    json.dump(out, open(os.path.join(HERE, 'check_y1a_symmetry.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
