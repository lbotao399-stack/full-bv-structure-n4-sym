"""
checks.py -- machine checks for the gauge-fixing foundation (exact, sympy).  Writes checks.json.

Sections
  A  realisation of rho, (-,-), Gamma (two-component formulas of Definition c8-twocomp)
  B  formulations: general scale-free ansatz for Q^GF on the minimal second-order complex, the first-order
     complex, the non-minimal complex with b kept, and the reduced non-minimal complex (cbar, cbar^*)
  C  the recommended operator in frame Psi: axioms (1)-(5'), Hodge decomposition at p != 0, the R^4 form of
     [C07, Lemma 6.5.1], real slice, propagator components; frame C and the compression of the former operator
  D  engine calibration (Q from S_0 in every frame) and the su(2) master equation in frame Psi
  E  the classical equivariant master equation (abelian, every generator type, arities 1-2) in the frames
     min / C / Psi, the non-minimal sector (arities 1-4), the pair (cbar, b) and the elimination of b
  F  generators on the gauge fixing: [x., QGF], composite kernels, the s_x Psi vertices
Negative controls are marked NEG; their recorded value is True when the control fails as it must.
Usage: python3 checks.py [--fast]   (--fast skips section E's 26-component runs)
"""
import sys, json, time
import sympy as sp
import realization as R
import freecx as Fx
import gfsolve as GS
import firstorder as FO

p = Fx.p
mp = Fx.neg(p)
OUT = {}
T0 = time.time()


def rec(sec, key, val):
    OUT.setdefault(sec, {})[key] = val
    print('[%s] %s : %s' % (sec, key, val), flush=True)


# ------------------------------------------------------------------ A
for k, v in R.verify_realization().items():
    rec('A_realisation', k, v)

# two-component form of Gamma(e1,e2)|_V used in the displayed vertices:
# Gamma^m = i ( ebar2_{ad i} sigmabar^{m ad al} e1_al^i + ebar1_{ad i} sigmabar^{m ad al} e2_al^i )
import random as _rnd
_rnd.seed(11)
_ok = True
for _t in range(3):
    e1 = sp.Matrix([_rnd.randint(-3, 3) + sp.I * _rnd.randint(-2, 2) for _ in range(16)])
    e2 = sp.Matrix([_rnd.randint(-3, 3) + sp.I * _rnd.randint(-2, 2) for _ in range(16)])
    G = R.Gamma(e1, e2)
    for m in range(4):
        tc = 0
        for i in range(4):
            for ad in range(2):
                for a_ in range(2):
                    tc += sp.I * (e2[R.xib(ad, i)] * R.SIGB[m][ad, a_] * e1[R.xi(a_, i)] + e1[R.xib(ad, i)] * R.SIGB[m][ad, a_] * e2[R.xi(a_, i)])
        if sp.expand(tc - G[m]) != 0:
            _ok = False
rec('A_realisation', 'two-component Gamma|_V = i(ebar2 sigmabar e1 + ebar1 sigmabar e2)', _ok)

# ------------------------------------------------------------------ B
Q = Fx.Q_min('min')
Qm = Fx.Q_min('min', mp)
rec('B_formulations', 'min: Q^2=0', Fx.is_zero(Q * Q))
rec('B_formulations', 'min: Q preserves pairing', Fx.is_zero(Fx.adj_check_Q(Q, Qm, 'min')))
G, syms, sol, _ = GS.solve_ansatz(Q, Qm, 'min')
rec('B_formulations', 'min: scale-free H=p^2 & adjoint solutions', str(sol))
rec('B_formulations', 'min: unique solution equals former QGF', len(sol) == 1 and Fx.is_zero(G.subs(sol[0]) - Fx.QGF_old()))
rec('B_formulations', 'NEG min: (QGF)^2 != 0 for the unique solution', not Fx.is_zero((G.subs(sol[0])) ** 2))
fo = FO.check_first_order()
for k, v in fo.items():
    rec('B_formulations', 'first-order: ' + k, v)
# reduced non-minimal: symplectic condition, alpha-family, nilpotency forces alpha = -1/2
N = GS.Nmap_red(p)
Nm_ = GS.Nmap_red(mp)
X = GS.symplectic_condition(sp.eye(56) + N, sp.eye(56) + Nm_, 'red')
rec('B_formulations', 'red: Phi=1+N symplectic iff', str(set(sum([sp.Poly(x, *p).coeffs() for x in X if x != 0], []))))
al = GS.al
QP = GS.Q_Psi_red(p, al, -1, 1)
QPm = GS.Q_Psi_red(mp, al, -1, 1)
rec('B_formulations', 'red: Q_Psi^2=0 (all alpha)', Fx.is_zero(QP * QP))
rec('B_formulations', 'red: Q_Psi preserves pairing', Fx.is_zero(Fx.adj_check_Q(QP, QPm, 'red')))
G, syms, sol, _ = GS.solve_ansatz(QP, QPm, 'red')
rec('B_formulations', 'red: H=p^2 & adjoint, family in alpha', str(sol))
Gs = G.subs(sol[0])
nz = set(sp.factor(sp.expand(x)) for x in Gs * Gs if sp.expand(x) != 0)
rec('B_formulations', 'red: (QGF)^2 entries on the family', str(nz))
rec('B_formulations', 'red: (QGF)^2=0 iff alpha', str(sp.solve(list(nz), al)))
rec('B_formulations', 'red: at alpha=-1/2 QGF = antifield->field operator',
    Fx.is_zero(Gs.subs(al, -sp.Rational(1, 2)) - GS.QGF_Psi_red(p)))
# b kept
k_ = GS.kap
Qb = GS.Q_Psi_b(p, 1, -1, 1, k_, -k_)
Qbm = GS.Q_Psi_b(mp, 1, -1, 1, k_, -k_)
rec('B_formulations', 'b: Q^2=0 & pairing', Fx.is_zero(Qb * Qb) and Fx.is_zero(Fx.adj_check_Q(Qb, Qbm, 'b')))
Gb, sb = GS.ansatz('b', p)
Gbm, _ = GS.ansatz('b', mp)
H = (Qb * Gb + Gb * Qb).applyfunc(sp.expand) - Fx.P2 * sp.eye(58)
A = Fx.adj_check_GF(Gb, Gbm, 'b')
E = []
for Mx in (H, A):
    for x in Mx:
        if sp.expand(x) != 0:
            E += sp.Poly(sp.expand(x), *p).coeffs()
solb = sp.solve(E, sb + [k_], dict=True)
rec('B_formulations', 'b: solutions (frame parameter kappa free)', str(solb))
Gbs = Gb.subs(solb[0])
rec('B_formulations', 'b: (QGF)^2=0 for every kappa', Fx.is_zero(Gbs * Gbs))

# ------------------------------------------------------------------ C
half = -sp.Rational(1, 2)
QPsi = GS.Q_Psi_red(p, half, -1, 1)
QPsim = GS.Q_Psi_red(mp, half, -1, 1)
GPsi = GS.QGF_Psi_red(p)
GPsim = GS.QGF_Psi_red(mp)
n = 56
rec('C_axioms', '(0) Q^2=0', Fx.is_zero(QPsi * QPsi))
rec('C_axioms', '(0) Q preserves the pairing', Fx.is_zero(Fx.adj_check_Q(QPsi, QPsim, 'red')))
rec('C_axioms', '(1) QGF adjointness opposite to Q', Fx.is_zero(Fx.adj_check_GF(GPsi, GPsim, 'red')))
rec('C_axioms', '(2) QGF order <= 1', max(sp.Poly(x, *p).total_degree() for x in GPsi if x != 0) <= 1)
rec('C_axioms', '(3) (QGF)^2=0', Fx.is_zero(GPsi * GPsi))
rec('C_axioms', '(4) H=[Q,QGF]=|p|^2 on every slot', Fx.is_zero(QPsi * GPsi + GPsi * QPsi - Fx.P2 * sp.eye(n)))
Pinv = GPsi / Fx.P2
rec('C_axioms', 'R4 Lemma 6.5.1: Q P + P Q = 1 (P = QGF/|p|^2)', Fx.is_zero((QPsi * Pinv + Pinv * QPsi).applyfunc(sp.simplify) - sp.eye(n)))
rec('C_axioms', 'R4 Lemma 6.5.1: P^2 = 0, P Q P = P, Q P Q = Q',
    Fx.is_zero(Pinv * Pinv) and Fx.is_zero((Pinv * QPsi * Pinv - Pinv).applyfunc(sp.simplify)) and Fx.is_zero((QPsi * Pinv * QPsi - QPsi).applyfunc(sp.simplify)))
P1 = (QPsi * Pinv).applyfunc(sp.simplify)
P2_ = (Pinv * QPsi).applyfunc(sp.simplify)
rec('C_axioms', 'projectors QP, PQ: idempotent, orthogonal, sum 1',
    Fx.is_zero((P1 * P1 - P1).applyfunc(sp.simplify)) and Fx.is_zero((P2_ * P2_ - P2_).applyfunc(sp.simplify))
    and Fx.is_zero((P1 * P2_).applyfunc(sp.simplify)) and Fx.is_zero(P1 + P2_ - sp.eye(n)))
pnum = dict(zip(p, (sp.Rational(3, 2), -2, 1, sp.Rational(5, 3))))
Qn, Gn = QPsi.subs(pnum), GPsi.subs(pnum)
rec('C_axioms', '(5) at p != 0: rank Q, rank QGF, rank [Q|QGF], dim', (Qn.rank(), Gn.rank(), Qn.row_join(Gn).rank(), n))
Qo, Go = Fx.Q_min('min').subs(pnum), Fx.QGF_old().subs(pnum)
rec('C_axioms', 'NEG former operator at p != 0: rank Q + rank QGF > dim (Im Q cap Im QGF != 0)',
    (Qo.rank(), Go.rank(), Qo.row_join(Go).rank(), 54))
r = Fx.rng
blk = (Fx.Omega('red') * QPsi)
blkA = sp.Matrix([[blk[i, j] for j in list(r('AV')) + list(r('AW'))] for i in list(r('AV')) + list(r('AW'))])
rec('C_axioms', 'real slice: <A, Q A> block = -1/2 |p|^2 on all ten components', Fx.is_zero(blkA + Fx.P2 / 2 * sp.eye(10)))
blko = (Fx.Omega('min') * Fx.Q_min('min'))
blkAo = sp.Matrix([[blko[i, j] for j in r('AV')] for i in r('AV')])
rec('C_axioms', 'NEG former: <A|_V, Q A|_V> block has the longitudinal zero eigenvalue', sp.simplify((blkAo * sp.Matrix(p)).applyfunc(sp.expand)) == sp.zeros(4, 1))
Pk = (GPsi * Fx.K0hat('red')).applyfunc(sp.expand)


def slotname(i):
    for s in Fx.SLOTS:
        if s[1] <= i < s[1] + s[2]:
            return s[0]


pairs = sorted(set((slotname(i), slotname(j)) for i in range(n) for j in range(n) if Pk[i, j] != 0))
rec('C_propagator', 'nonzero slot pairs of P (times e^{-eps|p|^2}/|p|^2)', str(pairs))
rec('C_propagator', 'A(x)A, c(x)cbar, cbar(x)c fibre factors', str((Pk[1, 1], Pk[5, 5], Pk[0, 54], Pk[54, 0])))
rec('C_propagator', 'graded symmetry P_ij(p) = (-1)^{|i||j|} P_ji(-p)',
    Fx.is_zero(Pk - sp.Matrix(n, n, lambda i, j: (-1) ** (Fx.parity_vector('red')[i] * Fx.parity_vector('red')[j]) * (GPsim * Fx.K0hat('red'))[j, i])))
rec('C_propagator', 'K0 represents the identity (K0hat Omega Pi = 1)', Fx.is_zero(Fx.K0hat('red') * Fx.Omega('red') * Fx.Pi('red') - sp.eye(n)))
# frame C and compression
Nn = GS.Nmap_red(p, -1, 1)
Phi = sp.eye(n) + Nn
Phinv = sp.eye(n) - Nn
Phim = sp.eye(n) + GS.Nmap_red(mp, -1, 1)
GC = (Phi * GPsi * Phinv).applyfunc(sp.expand)
QC = GS.Qprime_red(p, half)
proj = sp.zeros(54, n)
proj[:, :54] = sp.eye(54)
iota = proj.T
j = proj * Phi
jdag = Fx.Omega('red').inv() * (proj * Phim).T * Fx.Omega('min')
rec('C_frames', 'Q_Psi = Phi^{-1} Q_C Phi', Fx.is_zero(Phinv * QC * Phi - QPsi))
rec('C_frames', 'frame C: (QGF_C)^2 = 0, [Q_C, QGF_C] = |p|^2', Fx.is_zero(GC * GC) and Fx.is_zero(QC * GC + GC * QC - Fx.P2 * sp.eye(n)))
rec('C_frames', 'frame C: QGF_C entry cbar -> cbar^* (order two, axiom (2) fails in frame C)', str(GC[55, 54]))
GCm = ((sp.eye(n) + GS.Nmap_red(mp, -1, 1)) * GPsim * (sp.eye(n) - GS.Nmap_red(mp, -1, 1))).applyfunc(sp.expand)
rec('C_frames', 'frame C: QGF_C adjointness (axiom (1))', Fx.is_zero(Fx.adj_check_GF(GC, GCm, 'red')))
rec('C_frames', 'frame C: ranks at p != 0 (axiom (5))', (QC.subs(pnum).rank(), GC.subs(pnum).rank(), QC.subs(pnum).row_join(GC.subs(pnum)).rank(), n))
rec('C_frames', 'chain map j = p Phi: j Q_Psi = Q j', Fx.is_zero(j * QPsi - Fx.Q_min('min') * j))
rec('C_frames', 'j j^dag = 1  ((j x j) K_t = K_t)', Fx.is_zero(j * jdag - sp.eye(54)))
rec('C_frames', 'j QGF_Psi j^dag = former QGF  ((j x j) P_Psi = P_former)', Fx.is_zero(j * GPsi * jdag - Fx.QGF_old(p)))
rec('C_frames', 'compression: p QGF_C iota = former QGF', Fx.is_zero(proj * GC * iota - Fx.QGF_old(p)))
pinm = sp.zeros(n, n)
pinm[54, 54] = 1
pinm[55, 55] = 1
rec('C_frames', '(former QGF)^2 = - p QGF_C pi_nm QGF_C iota', Fx.is_zero(Fx.QGF_old(p) ** 2 + proj * GC * pinm * GC * iota))

# ------------------------------------------------------------------ D
import models as Mo
import test_su2
Mmod = Mo.build_model(True)
Qmin_e = Fx.embed(Fx.Q_min('min'), 54, 56)
for frame, Qx in (('min', Qmin_e), ('C', GS.Qprime_red(p, half)), ('Psi', QPsi)):
    F = Mo.Fields(Mmod, frame)
    bad = Mo.calibrate_Q(Mmod, F, Qx, p)
    rec('D_engine', 'frame %s: {S_0, l} = l o Q on every coordinate (mismatches)' % frame, bad)
for k, v in test_su2.run().items():
    rec('D_engine', k, v)

# ------------------------------------------------------------------ E
import test_nm
for k, v in test_nm.run().items():
    rec('E_equivariant_nm_sector', k, v)
# gauge-parameter direction: G = 1/2 int (cbar^*+d^*A, cbar^*+d^*A) = - s int (cbar^*+d^*A, cbar)  (frame Psi)
_M = Mo.build_model(True)
_F = Mo.Fields(_M, 'Psi')
_S0 = Mo.S0(_F)
_N = _M.mul(_F.cbs, _F.cb)
_G = _M.scale(_M.mul(_F.cbs, _F.cbs), sp.Rational(1, 2))
rec('E_gauge_parameter', 'G = -s(int (cbar^*+d^*A, cbar)) (free part)', _M.is_zero_functional(_M.plus(_G, _M.bracket(_S0, _N))))
rec('E_gauge_parameter', 'NEG G = +s(...) fails', not _M.is_zero_functional(_M.plus(_G, _M.scale(_M.bracket(_S0, _N), -1))))
if '--fast' not in sys.argv:
    from test_equiv import run as erun
    for (frame, g2, ct, label) in [('min', None, True, 'min (former structure)'),
                                   ('C', -1, True, 'C, cbar weight 2, gamma2=-1'),
                                   ('Psi', -1, True, 'Psi, cbar weight 2, gamma2=-1'),
                                   ('Psi', 0, False, 'Psi, inert cbar, gamma2=0'),
                                   ('Psi', 0, True, 'NEG Psi, cbar weight 2, gamma2=0'),
                                   ('Psi', -1, False, 'NEG Psi, inert cbar, gamma2=-1')]:
        res = erun(frame, gamma2=g2, cbar_transforms=ct, verbose=False)
        fails = [k for k, v in res.items() if not v]
        rec('E_equivariant_abelian', label + ': components checked', len(res))
        rec('E_equivariant_abelian', label + ': failing components', fails)
        if label.startswith('NEG'):
            rec('E_equivariant_abelian', label + ': fails exactly on the four odd-odd arity-two components',
                sorted(fails) == sorted(['ar2_xi1_xi2', 'ar2_xi1_eta1', 'ar2_eta1_eta2', 'ar2_xi1_xi1']))
        else:
            rec('E_equivariant_abelian', label + ': all components vanish', fails == [])

# ------------------------------------------------------------------ F
import test_ops
res, tabs = test_ops.run(verbose=False)
for k, v in res.items():
    rec('F_generators', k, v)
json.dump(tabs, open('commutator_tables.json', 'w'), indent=1)
import runpy
vres = runpy.run_path('test_vertices.py')['res']
for k, v in vres.items():
    rec('F_generators', k, v)

OUT['_runtime_s'] = round(time.time() - T0, 1)
json.dump(OUT, open('checks.json', 'w'), indent=1, default=str)
print('written checks.json, %.0f s' % (time.time() - T0))
