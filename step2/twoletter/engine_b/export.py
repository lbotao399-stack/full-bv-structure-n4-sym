"""Assemble kernels.json from kernels_W0.json, the derivative tables (kernels_deriv.json, if present), the closed forms
verified by the gates, and the planar two-site densities (scalar SO(6), sl(2) null derivatives)."""
import os, sys, json, math
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main():
    W0 = json.load(open(os.path.join(HERE, 'kernels_W0.json')))
    try:
        DER = json.load(open(os.path.join(HERE, 'kernels_deriv.json')))
    except Exception:
        DER = {}
    try:
        TS = json.load(open(os.path.join(HERE, 'two_site.json')))
    except Exception:
        TS = {}
    gates = json.load(open(os.path.join(HERE, 'gates.json')))
    cov = {}
    g7_all = True
    for pair, e in W0.items():
        ls_tab = {}
        for ls, v in e['leg_sets'].items():
            g7_all &= bool(v['G7']['n_out_ok'] and v['G7']['momentum_degrees_ok'])
            terms = sorted(set(k for w in v['total'].values() for k in w))
            ls_tab[ls] = terms if terms else 'vanishes'
        cov[pair] = dict(W=e['W'], complete=e.get('complete'), failed=e.get('failed_leg_sets', []), leg_sets=ls_tab)
    for key, e in DER.items():
        cov[key] = dict(W=e['W'], complete=e.get('complete'),
                        leg_sets={ls: (sorted(set(k for w in v['total'].values() for k in w)) or 'vanishes')
                                  for ls, v in e['leg_sets'].items()})
    out = dict(
        coverage=cov,
        G7_output_bound_all_W0_kernels=g7_all,
        engine='B (route B: explicit momentum-space numerators, pinned-edge Feynman-parameter integrals, numerical eps residuals)',
        conventions=dict(
            kernel='K(X_a,X_b) = ct(int omega X_a X_b): coefficient of eps^-m and log eps of (hbar/2) sum_seq STr[P Pi h_1 ... P Pi H_M]; '
                   'equal to hbar Sing_eps (d^alpha x d^beta) P_e^{ab}(x,x) of Def def:tl-dressed (sigma_ab = 1)',
            units='every coefficient string is in units hbar/(16 pi^2); for su(N) a planar adjacent pair with coefficient c and word w '
                  'contributes (-1)^r c/2 units u = hbar tr(ad ad)/(16 pi^2 (,)) = -g^2 N/(32 pi^2), outputs in the order Y_wr ... Y_w1',
            colour='C_w = (1 (x) ad_{Y_w1} ... ad_{Y_wr}) (-,-)^{-1}; word "i,j" = (w1,w2) = output labels; X_a paired with the first factor',
            letters='Xa[comp], Xb[comp]: components of the two letters (dual coordinates); the kernel is linear in each',
            outputs='Y<j>[comp]: polarisation of output leg j, p<j>_m: its incoming momentum (i p_m <-> d_m acting on that output); '
                    'odd outputs: Grassmann parameters in front in increasing label order',
            components='two-component realisation of engine_b/realization.py: A_m (m=0..3) = A|_V, A_I1..A_I6 = A|_{W_R^perp} (SO(6) vector, '
                       'A_ij = A_I Sigma^I_ij with the t Hooft blocks of realization.py), lam_a<alpha>^<i> = lambda_+ in S_+ (x) W, '
                       'lambar_ad<alphadot>_<i> = lambda_- in S_- (x) W^*, lams^a<alpha>_<i> = lambda^*_+ in S_+ (x) W^*, '
                       'lambars^ad<alphadot><i> = lambda^*_- in S_- (x) W, As_m, As_I = A^*, c, cs = c^*; sigma = (1, i sigma_k), sigmabar = (1, -i sigma_k) '
                       '(Euclidean); Weinberg-Srednicki form via step2/notation/dictionary.json (Def c8-twocomp)',
            regulator='per-edge Schwinger cutoff t_e >= eps on every edge, T = infinity at generic momenta; no dimensional regularisation',
        ),
        closed_forms=dict(
            k0='K^(0) generating function (point splitting s,u): A|A: 2 (Xa.Xb) sum_m eps^-m (-1)^(m+1) (|u-s|^2/4)^(m-1)/m!; '
               'single powers eps^-(W/2), no log, no finite part (Lemma lem:tl-k0) -- verified to W<=6 (gate G6)',
            k1_bubble='N_V constant (quartic bubble): sum_{alpha,beta} s^alpha u^beta/(alpha! beta!) [log] K(d^alpha phi_a, d^beta phi_b) = '
                      '[log] K(phi_a,phi_b) x sum_{a,b,j} (a+j)!(b+j)!/(a! b! (j!)^2 (a+b+2j+1)!) (i(s,R))^a (i(u,R))^b |u-s|^{2j} (|R|^2/4)^j, '
                      'R the total output momentum (i R = d on the product of the outputs), eq:tl-bubbleclosed -- verified against the explicit '
                      'point-split tables to 4 derivatives (gate closed form)',
            scalar_two_site='planar density on letters A|_perp (outputs A|_perp A|_perp): h = u (2P - K) (SO(6)) -- gate G8',
            sl2_two_site='planar density on (n.d)^a Z (n null, Z = A_I1 + i A_I2): eigenvalues u(2 - 4 h(j)), j = 0..n; '
                         'h - 2u = -2u * 2h(J_12) (harmonic action) -- gate G8; matrices in two_site',
        ),
        pairs=W0,
        derivative_tables=DER,
        two_site=TS,
        gates_summary={k: v['pass_'] for k, v in gates.items()},
    )
    json.dump(out, open(os.path.join(HERE, 'kernels.json'), 'w'), indent=0)
    print('kernels.json written', os.path.getsize(os.path.join(HERE, 'kernels.json')) // 1024, 'kB')


if __name__ == '__main__':
    main()
