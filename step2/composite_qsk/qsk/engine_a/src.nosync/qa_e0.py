"""Engine A, part Q: task E0 -- the operators on plane waves (exact Gaussian-integer arithmetic on the 54-component
single-colour field space of ea_core, column vectors e = (c, A_0..9, A*_0..9, lam_1..16, lam*_1..16, c*)).

Formula operators (index-free statements of the manuscript in the cliff.py basis):
  Q (Prop c1-Q):   A <- -dc,  A* <- -1/2 d*dA,  lam* <- dslash lam (rho on Sigma = sb),  c* <- d*(A*|_V);  d* = -iota_d;
  Q^GF (Lemma qs-arrows): A <- -2 A*,  c <- -d*(A|_V),  A*|_V <- d c*,  lam <- -dslash|_{Sigma*} lam* (rho = s);
  H = |p|^2;
  xi. (Lemma qs-arrows, sign corrected in the revision of 2026-10-30, finding Q3): A <- +Gamma(xi, lam),
       lam <- +1/2 rho(dA) xi,  plus the antifield arrows A* -> lam*,
       lam* -> A*, which are read off the S^(1)(xi) block (qa_s1) and checked against the formula field arrows;
  even v. (Def c7-vdot) for translations and the constant parts of rotations/dilatations.
On a plane wave e^{i(p,x)}, d_mu -> i p_mu.  Commutators are graded: odd xi. and odd Q, Q^GF anticommute in [,]."""
import numpy as np
import qa_s1 as S1
import ea_core as E
import cq_core as Q

NS = E.NS
O = E.OFS
SB, SS = S1.SB, S1.SS


def idx(t, a):
    return O[t] + a


def op_Q(p):
    M = np.zeros((NS, NS), complex)
    for mu in range(4):
        M[idx('A', mu), idx('c', 0)] = -1j * p[mu]                      # A = -dc
    pp = sum(p[m] * p[m] for m in range(4))
    for m in range(4):                                                  # A*_m = -1/2 (p^2 A_m - p_m p.A)
        for n in range(4):
            M[idx('As', m), idx('A', n)] = -0.5 * ((pp if m == n else 0) - p[m] * p[n])
    for I in range(4, 10):
        M[idx('As', I), idx('A', I)] = -0.5 * pp
    Dsl = sum(1j * p[mu] * SB[mu] for mu in range(4))                   # lam* = dslash lam
    M[O['ls']:O['ls'] + 16, O['l']:O['l'] + 16] = Dsl
    for mu in range(4):                                                 # c* = d*(A*_V) = -i p.A*
        M[idx('cs', 0), idx('As', mu)] = -1j * p[mu]
    return M


def op_QGF(p):
    M = np.zeros((NS, NS), complex)
    for Mi in range(10):
        M[idx('A', Mi), idx('As', Mi)] = -2
    for mu in range(4):
        M[idx('c', 0), idx('A', mu)] = 1j * p[mu]                       # -d*(A_V) = +i p.A
        M[idx('As', mu), idx('cs', 0)] = 1j * p[mu]                     # A*_V = d c*
    Dsl = sum(1j * p[mu] * SS[mu] for mu in range(4))
    M[O['l']:O['l'] + 16, O['ls']:O['ls'] + 16] = -Dsl
    return M


def op_H(p):
    return sum(p[m] * p[m] for m in range(4)) * np.eye(NS)


def PAR():
    return np.array([E.PAR[t] for t in E.TYPES for _ in range(E.SZ[t])])


def op_xi_fields(xi, p, sign=+1):
    """field arrows of Lemma qs-arrows (corrected sign +1: lam -> +Gamma(xi,lam), A -> +1/2 rho(dA) xi);
    sign=-1 gives the arrows as first printed in the theory note"""
    xi = np.asarray(xi, complex)
    M = np.zeros((NS, NS), complex)
    for Mi in range(10):
        M[idx('A', Mi), O['l']:O['l'] + 16] = (xi @ SB[Mi])            # A_M <- +xi^T sb[M] lam
    X = S1._X(xi)
    pv = np.zeros(10, complex); pv[:4] = p[:4]
    for Mi in range(10):
        for N in range(10):
            if Mi < N:
                # -1/2 F_MN X[M,N], F_MN = i(p_M A_N - p_N A_M)
                M[O['l']:O['l'] + 16, idx('A', N)] += 0.5 * 1j * pv[Mi] * X[Mi, N]
                M[O['l']:O['l'] + 16, idx('A', Mi)] += -0.5 * 1j * pv[N] * X[Mi, N]
    return sign * M


def block_operator(blk, p):
    """operator read off a two-leg block of a quadratic functional int (X e, e*): rows = target field slot,
    columns = source; (X)_{f <- g} = Kn[f, f*] (S h)[f*, g] at (kj, kl) = (-p, p) (the antifield leg absorbs -p)"""
    h = blk.h0 - sum(p[mu] * blk.Hj[mu] for mu in range(4)) + sum(p[mu] * blk.Hl[mu] for mu in range(4))
    Sd = Q.SD
    return Q.KN @ (Sd[:, None] * h)


def gcomm(A, B, parA, parB):
    return A @ B - (-1) ** (parA * parB) * B @ A


def run(log=print):
    rng = np.random.default_rng(3)
    out = {}
    p = np.array([1, -2, 3, 1], dtype=complex)
    Qm, Qg, H = op_Q(p), op_QGF(p), op_H(p)
    out['QQ=0'] = bool(np.allclose(Qm @ Qm, 0))
    out['[Q,QGF]=H'] = bool(np.allclose(Qm @ Qg + Qg @ Qm, H))
    # supercharge
    Bp = S1.chir_basis('lamP')
    xi = sum(int(c) * v for c, v in zip(rng.integers(-2, 3, 8), Bp)).astype(complex)
    gen = S1.quad_block(S1.odd_terms(xi, ('T1', 'T2', 'T3', 'T4')))
    Xb = block_operator(gen, p)
    Xf = op_xi_fields(xi, p)
    fld = [O['A'] + i for i in range(10)] + [O['l'] + a for a in range(16)] + [O['c']]
    # field-to-field arrows of the block versus the formula: common factor
    sub_b = Xb[np.ix_(fld, fld)]; sub_f = Xf[np.ix_(fld, fld)]
    nz = np.abs(sub_f) > 1e-9
    ratios = (sub_b[nz] / sub_f[nz]) if nz.any() else np.array([])
    rr = np.unique(np.round(ratios, 9))
    out['xi_block_field_arrows_over_formula'] = [complex(x) for x in rr]
    out['xi_block_zero_where_formula_zero'] = bool(np.allclose(sub_b[~nz], 0))
    c = complex(rr[0]) if len(rr) == 1 else None
    Xi = Xb / c if c is not None else Xb
    out['[xi.,Q]=0'] = bool(np.allclose(gcomm(Xi, Qm, 1, 1), 0))
    out['[xi.,H]=0'] = bool(np.allclose(gcomm(Xi, H, 1, 0), 0))
    Cxi = gcomm(Xi, Qg, 1, 1)
    comps = {}
    names = ['c', 'A', 'As', 'l', 'ls', 'cs']
    for t1 in names:
        for t2 in names:
            blk = Cxi[O[t1]:O[t1] + E.SZ[t1], O[t2]:O[t2] + E.SZ[t2]]
            if np.abs(blk).max() > 1e-9:
                comps['%s<-%s' % (t1, t2)] = True
    out['C_xi_components'] = sorted(comps)
    # lam -> c component = +i (p, Gamma(xi,lam)|_V)  (corrected Lemma qs-arrows)
    tgt = np.array([1j * sum(p[mu] * (xi @ SB[mu])[a] for mu in range(4)) for a in range(16)])
    out['C_xi(lam->c) = +i(p,Gamma(xi,lam)|_V)'] = bool(np.allclose(Cxi[O['c'], O['l']:O['l'] + 16], tgt))
    # finer components (finding Q9): vector part V (0..3) / perp (4..9), chirality of lam (G5S eigenvalue)
    G5 = np.asarray(S1.G5S); G5s = np.asarray(S1.G5SS)
    Pp, Pm = (np.eye(16) + G5) / 2, (np.eye(16) - G5) / 2
    Pps, Pms = (np.eye(16) + G5s) / 2, (np.eye(16) - G5s) / 2
    fine = {}
    def nzb(rows, cols):
        return bool(np.abs(Cxi[np.ix_(rows, cols)]).max() > 1e-9)
    Arow = {'V': [O['A'] + i for i in range(4)], 'perp': [O['A'] + i for i in range(4, 10)]}
    Asr = {'V': [O['As'] + i for i in range(4)], 'perp': [O['As'] + i for i in range(4, 10)]}
    lam = list(range(O['l'], O['l'] + 16)); lams = list(range(O['ls'], O['ls'] + 16))
    for k, rws in Arow.items():
        fine['A|_%s <- lam*' % k] = nzb(rws, lams)
    for k, cls in Asr.items():
        fine['lam <- A*|_%s' % k] = nzb(lam, cls)
    blk = Cxi[O['c'], O['l']:O['l'] + 16]
    fine['c <- lam_(+1 eigenspace of G5S)'] = bool(np.abs(blk @ Pp).max() > 1e-9)
    fine['c <- lam_(-1 eigenspace of G5S)'] = bool(np.abs(blk @ Pm).max() > 1e-9)
    fine['c <- A'] = nzb([O['c']], list(range(O['A'], O['A'] + 10)))
    out['C_xi_fine_components'] = fine
    # first order in p of every component
    def Cx(pp):
        return gcomm(block_operator(gen, pp) / c, op_QGF(pp), 1, 1)
    hom = np.allclose(Cx(2 * p), 2 * Cx(p))
    out['C_xi_homogeneous_degree_1'] = bool(hom)
    # composite line along xi: chain L S h L at line momentum k, zero block momentum (r = 0)
    k = p
    Lk = Q.LP0 + sum(k[mu] * Q.LP1[mu] for mu in range(4))
    h = gen.h0 - sum(-k[mu] * gen.Hj[mu] for mu in range(4)) * 0  # placeholder (overwritten below)
    hk = gen.h0 + sum(k[mu] * (gen.Hl[mu] - gen.Hj[mu]) for mu in range(4))   # kj = -k, kl = k
    Sd = Q.SD
    comp = Lk @ (Sd[:, None] * hk) @ (Sd[:, None] * Lk) / 2                   # vertex scale 2
    kk = complex(sum(k[m] * k[m] for m in range(4)))
    # components c x lam and lam x c only; multiple of xi * |k|^2 (numerator of xi e^{-2 eps k^2}/k^4 ... /k^2)
    nzc = np.argwhere(np.abs(comp) > 1e-9)
    slots = set()
    for r_, c_ in nzc:
        tr = [t for t in names if O[t] <= r_][-1]; tc = [t for t in names if O[t] <= c_][-1]
        slots.add((tr, tc))
    out['composite_line_components'] = sorted('%s x %s' % s for s in slots)
    row_c_lam = comp[O['c'], O['l']:O['l'] + 16]
    mult = row_c_lam / (kk * xi) if np.all(np.abs(xi) > 0) else None
    nzx = np.abs(xi) > 0
    m = np.unique(np.round(row_c_lam[nzx] / (kk * xi[nzx]), 9))
    out['composite_line_c_x_lam / (|k|^2 xi)'] = [complex(x) for x in m]
    m2 = np.unique(np.round(comp[O['l']:O['l'] + 16, O['c']][nzx] / (kk * xi[nzx]), 9))
    out['composite_line_lam_x_c / (|k|^2 xi)'] = [complex(x) for x in m2]
    # (C_xi (x) 1) K_0 graded symmetric: C Khat as a matrix M with M^T = (sign) M
    log(out)
    return out


if __name__ == '__main__':
    run()
