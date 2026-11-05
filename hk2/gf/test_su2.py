"""
Nonabelian check of the classical master equation in frame Psi for the gauge sector (g = su(2), basis with
[T_a,T_b] = eps_abc T_c, (T_a,T_b) = delta_ab):
  S_Psi = int ( -1/4 (F_A,F_A) + (d_A c, A^* - d cbar) + 1/2([c,c],c^*) - 1/4 (cbar^* + d^*A, cbar^* + d^*A) ),
with negative controls (wrong sign of the shift of A^*, wrong coefficient of the cbar^* term).
"""
import sympy as sp, time
import bvjet

X = bvjet.X
E4 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]


def model():
    M = bvjet.Model()
    for a in range(3):
        M.add('c%d' % a, 1)
    for a in range(3):
        for m in range(4):
            M.add('A%d_%d' % (a, m), 0)
    for a in range(3):
        M.add('cb%d' % a, 1)
    for a in range(3):
        for m in range(4):
            M.add('As%d_%d' % (a, m), 1)
    for a in range(3):
        M.add('cs%d' % a, 0)
    for a in range(3):
        M.add('cbs%d' % a, 0)
    for a in range(3):
        for m in range(4):
            M.pair('A%d_%d' % (a, m), 'As%d_%d' % (a, m), -1)
        M.pair('c%d' % a, 'cs%d' % a, -1)
        M.pair('cb%d' % a, 'cbs%d' % a, -1)
    return M


def eps3(a, b, c):
    return sp.LeviCivita(a, b, c)


def action(M, shift_sign=-1, cbs_coeff=-sp.Rational(1, 4), dstar_sign=1):
    j = M.jet
    A = [[j('A%d_%d' % (a, m)) for m in range(4)] for a in range(3)]
    c = [j('c%d' % a) for a in range(3)]
    cb = [j('cb%d' % a) for a in range(3)]
    cs = [j('cs%d' % a) for a in range(3)]
    As = [[M.plus(j('As%d_%d' % (a, m)), j('cb%d' % a, E4[m], shift_sign)) for m in range(4)] for a in range(3)]
    cbs = []
    for a in range(3):
        t = j('cbs%d' % a)
        for m in range(4):
            t = M.plus(t, j('A%d_%d' % (a, m), E4[m], -dstar_sign))
        cbs.append(t)
    S = {}
    # -1/4 (F,F) = -1/4 sum_{m<n} F_mn^a F_mn^a
    for m in range(4):
        for n in range(m + 1, 4):
            for a in range(3):
                Fmn = M.plus(M.Dtot(A[a][n], m), M.scale(M.Dtot(A[a][m], n), -1))
                for b in range(3):
                    for cc in range(3):
                        e = eps3(a, b, cc)
                        if e != 0:
                            Fmn = M.plus(Fmn, M.scale(M.mul(A[b][m], A[cc][n]), e))
                S = M.plus(S, M.scale(M.mul(Fmn, Fmn), -sp.Rational(1, 4)))
    # (d_A c, A^*_shifted)
    for m in range(4):
        for a in range(3):
            dc = M.Dtot(c[a], m)
            for b in range(3):
                for cc in range(3):
                    e = eps3(a, b, cc)
                    if e != 0:
                        dc = M.plus(dc, M.scale(M.mul(A[b][m], c[cc]), e))
            S = M.plus(S, M.mul(dc, As[a][m]))
    # 1/2 ([c,c], c^*)
    for a in range(3):
        for b in range(3):
            for cc in range(3):
                e = eps3(a, b, cc)
                if e != 0:
                    S = M.plus(S, M.scale(M.mul(M.mul(c[b], c[cc]), cs[a]), sp.Rational(1, 2) * e))
    for a in range(3):
        S = M.plus(S, M.scale(M.mul(cbs[a], cbs[a]), cbs_coeff))
    return M.simplify(S)


def run():
    out = {}
    M = model()
    t0 = time.time()
    S = action(M)
    out['su2_frame_Psi_{S,S}=0'] = M.is_zero_functional(M.bracket(S, S))
    S_bad = action(M, shift_sign=+1)
    out['NEG_su2_shift_A*+dcbar_{S,S}!=0'] = not M.is_zero_functional(M.bracket(S_bad, S_bad))
    S_bad2 = action(M, dstar_sign=-1)
    out['NEG_su2_cbar*-d^*A_{S,S}!=0'] = not M.is_zero_functional(M.bracket(S_bad2, S_bad2))
    out['_time_s'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    print(run())
