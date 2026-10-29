"""Family F2 su(1|1): tr(Z^{L-M} psi^M), L <= 6, letters without derivatives; and the su(2) check of F0.
Composite densities (route B) on collinear on-shell outputs p_j = t_j nbar, nbar = (1/2, -i/2, 0, 0) (engine comps),
P_b along b = n = (1, i, 0, 0); log densities from the twoletter route-B job (tl_engine.Job) in the same letters.
Z = A_4 + i A_5 (charges (1,0,0)); psi = lambda_+ letter of charges (1/2, 1/2, -1/2; 1/2, 1/2) (R; Lorentz),
on shell for p ~ nbar (checked)."""
import os, sys, json, time, pickle
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
import ce_core as C
import ce_letters as LT
import ce_job as CJ
import ce_planar as PL
import tl_engine as T

I = sp.I
n = [1, I, 0, 0]
nbar = [sp.Rational(1, 2), -I / 2, 0, 0]
ZP = {4: 1, 5: I}
ZO = {4: sp.Rational(1, 2), 5: -I / 2}
XP = {6: 1, 7: I}
XO = {6: sp.Rational(1, 2), 7: -I / 2}


def psi_data():
    lets = LT.lam_letters()
    Zm = sp.Matrix.hstack(*[v for v, ch, c in lets]).T          # rows = letters
    U = Zm.inv()                                                 # columns = dual basis (states)
    target = (sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(-1, 2), sp.Rational(1, 2), sp.Rational(1, 2))
    k = [i for i, (v, ch, c) in enumerate(lets) if ch == target][0]
    z = lets[k][0]
    u = U[:, k]
    # on shell for p ~ nbar
    rho = sum((nbar[m] * C.Rm[m] for m in range(4)), sp.zeros(16, 16))
    onshell = (rho * u).applyfunc(sp.expand).is_zero_matrix
    zp = {i: z[i] for i in range(16) if z[i] != 0}
    uo = {i: u[i] for i in range(16) if u[i] != 0}
    return zp, uo, onshell, (z.T * u)[0]


def run(Lmax=6, verbose=True):
    t0 = time.time()
    zpsi, upsi, onshell, norm = psi_data()
    info = dict(psi_letter=str(zpsi), psi_state=str(upsi), psi_onshell_for_p_nbar=onshell, z_dot_u=str(norm))
    LET = {'Z': ('AW', ZP, ZO, 0), 'psi': ('lamP', zpsi, upsi, 1)}
    odd = {'Z': 0, 'psi': 1}
    dens = {'D0': {}, 'Pb': {}, 'log': {}}
    nonrat = []
    for A in ('Z', 'psi'):
        for B in ('Z', 'psi'):
            sa, za, ua, pa = LET[A]; sb, zb, ub, pb = LET[B]
            outs_l = sorted([A, B], key=lambda x: (x != 'Z'))     # Z first
            outs = [LET[x][0] for x in outs_l]
            pols = [LET[x][2] for x in outs_l]
            oddo = [LET[x][3] for x in outs_l]
            job = CJ.CJob(sa, sb, outs, pol=pols, za=za, zb=zb, deriv=('null', 0, n, nbar), nullw=False)
            res = job.kernels(mlist=(0,), verbose=False)
            dP, nr1 = PL.density(job, res['Pb'], outs_l, oddo, 0, avec=n, smax=0)
            # other directions b: nbar, e_2, e_3 (b.nbar = 0): sum C2_Pb must vanish on collinear legs
            for nmb, bv in (('nbar', nbar), ('e2', [0, 0, 1, 0]), ('e3', [0, 0, 0, 1])):
                dX, _ = PL.density(job, res['Pb'], outs_l, oddo, 0, avec=bv, smax=0)
                info.setdefault('Pb_other_directions_nonzero_entries', {})[f'{A}|{B}|{nmb}'] = len(dX)
            dD, nr2 = PL.density(job, res['D0'], outs_l, oddo, 0, smax=0)
            nonrat += nr1 + nr2
            dens['Pb'][(A, B)] = dP; dens['D0'][(A, B)] = dD
            # log kernel (twoletter route-B job, same letters/outputs)
            tj = T.Job(sa, sb, outs, pol=pols, za=za, zb=zb, deriv=('null', 0, n, nbar))
            K = tj.kernel(keep_sequences=False)
            dl, nr3 = PL.density(tj, K['total'], outs_l, oddo, logkey='log', smax=0)
            dens['log'][(A, B)] = dl
            if verbose:
                print(A, B, 'done', round(time.time() - t0, 1), 's', flush=True)
    info['nonrational_entries'] = [str(x) for x in nonrat]
    return dens, odd, info


def matrices(dens, odd, Lmax=6, th=None):
    """cyclic matrices for tr(Z^{L-M} psi^M): composite D0 and Pb vs P7 targets; lc^natural and D_1."""
    thp, thl = sp.symbols('theta_perp theta_lambda') if th is None else th
    out = {}
    nF = {'Z': 2, 'psi': 4}
    cphi = {'Z': sp.Rational(2, 3), 'psi': sp.Rational(1, 3)}
    wX = {'Z': 1, 'psi': sp.Rational(3, 2)}
    thX = {'Z': thp, 'psi': thl}
    for L in range(2, Lmax + 1):
        for M in range(0, L + 1):
            content = {'Z': L - M, 'psi': M}
            B0 = PL.cyclic_basis(content, L, odd, 0)
            if not B0:
                out[f'L{L}_M{M}'] = dict(dim=0)
                continue
            B1 = PL.cyclic_basis(content, L, odd, 1)
            MD = PL.matrix(B0, B0, lambda w: PL.apply_pairs(w, dens['D0'], odd))
            MP = PL.matrix(B0, B1, lambda w: PL.apply_pairs(w, dens['Pb'], odd))
            TP = PL.matrix(B0, B1, lambda w: PL.deriv_map(w, odd)) * sp.Rational(-1, 3)
            TD = -sp.Rational(4, 3) * L * sp.eye(len(B0))
            # log kernel and lc^natural
            Mlog = PL.matrix(B0, B0, lambda w: PL.apply_pairs(w, dens['log'], odd))
            nFw = sum(nF[x[0]] for x in B0[0])
            lcn = Mlog - nFw * sp.eye(len(B0))
            # one-body -{F'_Delta, X}^0 = u (c_phi - 2 theta w_X) X ; Xi^0 = sum C^0 + one-body
            ob = sum(cphi[x[0]] - 2 * thX[x[0]] * wX[x[0]] for x in B0[0])
            Xi0 = MD + ob * sp.eye(len(B0))
            D1 = (4 * lcn + Xi0).applyfunc(sp.expand)
            D1c = D1.subs({thp: sp.Rational(-1, 3), thl: sp.Rational(-1, 3)})
            ev = D1c.eigenvals()
            out[f'L{L}_M{M}'] = dict(dim=len(B0), basis=[str(w) for w in B0],
                                    P7b_sumC_Pb=str(MP.tolist()), P7b_target=str(TP.tolist()), P7b_pass=(MP - TP).is_zero_matrix,
                                    P7c_sumC_D0=str(MD.tolist()), P7c_target=str(TD.tolist()), P7c_pass=(MD - TD).is_zero_matrix,
                                    lc_natural_units=str(lcn.tolist()),
                                    Xi_Delta_translation_invariant=str(Xi0.applyfunc(sp.expand).tolist()),
                                    D1_translation_invariant_symbolic=str(D1.tolist()),
                                    D1_at_theta_m13_eigenvalues={str(k): v for k, v in ev.items()})
    return out


if __name__ == '__main__':
    dens, odd, info = run()
    pickle.dump(dict(dens=dens, odd=odd, info=info), open(os.path.join(HERE, '..', 'out.nosync', 'su11_dens.pkl'), 'wb'))
    print(info)
    for k, v in dens.items():
        for pr, d in v.items():
            print(k, pr, {str(a): str(b) for a, b in d.items()})
    mats = matrices(dens, odd)
    for k, v in mats.items():
        print(k, v.get('dim'), v.get('P7b_pass'), v.get('P7c_pass'), v.get('D1_at_theta_m13_eigenvalues'))
    dstr = {k: {f'{pr[0]}|{pr[1]}': {str(a): str(b) for a, b in d.items()} for pr, d in v.items()} for k, v in dens.items()}
    json.dump(dict(info=info, densities_units_u=dstr, cyclic=mats), open(os.path.join(HERE, '..', 'out.nosync', 'su11.json'), 'w'), indent=1, default=str)
