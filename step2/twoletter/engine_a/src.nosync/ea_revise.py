"""Engine A, revision 2026-10-29 (qa/disposition.json): gates G3-W2 (two-body Ward identity, ea_w2.py) and G8f (su(1|1)
fermion family, ea_su11.py); labels of the D_1 entries (one-body values without C^(1)_Delta; spectra conditional on
P7); D_1 and P_1 on the planar scalar sector (theory Prop prop:tl-scalarsector, input Thm c9-D1 / Thm c9-P1)."""
import json, os
import ea_w2 as W2, ea_su11 as SU

SRC = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(SRC)
OUT = os.path.join(HERE, 'out.nosync')


def gate(gid, name, value, target, ok, note=''):
    return {'gate': gid, 'name': name, 'value': value, 'target': target, 'pass': bool(ok), 'note': note}


def revise(kern, gates, rerun=True, log=print):
    w2path = os.path.join(OUT, 'w2_gate.json')
    if rerun or not os.path.exists(w2path):
        r = W2.run(w2path, deriv_pairs=[('AV', 'AV', 0, 1), ('AV', 'AV', 2, 2), ('AV', 'AV', 1, 1), ('lamP', 'lamM', 1, 3),
                                        ('lamP', 'lamM', 2, 2), ('AW', 'AW', 0, 0), ('AW', 'AW', 0, 3), ('AV', 'AW', 1, 2)],
                   log=lambda *a: None)
        r['negative_controls'] = W2.negative_controls()
        json.dump(r, open(w2path, 'w'), indent=1, default=str)
    r = json.load(open(w2path))
    sm, nc = r['summary'], r['negative_controls']
    gates = [g for g in gates if g['gate'] not in ('G3-W2', 'G3-W2n', 'G8f')]
    gates.append(gate('G3-W2', 'two-body Ward identity W2 as lc(sM) = s lc(M) + {L,M} on M = tr(zeta X_a X_b), field letters '
                      '(A_V, A_perp, lambda_+-, c), W_0, W_0+1, selected W_0+2, targets <= 3 field legs (ea_w2.py)',
                      {'tested': sm['tested'], 'pass': sm['pass'], 'nontrivial': sm['nontrivial'], 'skipped': sm['skipped']},
                      'all tested pass', sm['tested'] == sm['pass'] and sm['nontrivial'] > 0,
                      'skipped: targets reachable from an antifield output of lc(M) (s on antifield letters not implemented)'))
    gates.append(gate('G3-W2n', 'W2 negative controls (one convention flipped at a time)', nc, 'every flip detected',
                      all(v['detected'] for v in nc.values())))
    su = SU.run()
    json.dump(su, open(os.path.join(OUT, 'su11_gate.json'), 'w'), indent=1, default=str)
    gates.append(gate('G8f', 'su(1|1) sector {Z, psi}: planar lc^natural = -2u (1 - Pi) (graded) on (Z,psi), (Z,Z), (psi,psi)',
                      su['lc_natural_pair'], su['target_-2u(1-Pi)'], su['pass'], 'ea_su11.py'))
    # ------------------------------------------------------------------ labels and new entries
    D1 = kern['D1']
    for X in ('A_V', 'lambda', 'c'):
        D1['one_body_units'][X]['label'] = ('without C^(1)_Delta (nonlocal, Rem c9-Pi); parameters: -theta_%s v.X from '
                                            "F'_Delta" % ('V' if X in ('A_V', 'c') else 'lambda'))
        D1['one_body_units'][X]['v.X'] = '-theta_%s (+ part of C^(1)_Delta)' % ('V' if X in ('A_V', 'c') else 'lambda')
    D1['one_body_units']['A_perp']['v.X'] = '5/3 - theta_perp'
    D1['one_body_units']['A_perp']['label'] = 'one-output part, local (Pi_Delta|perp, Rem c9-Pi)'
    D1['two_body'] = ('4 K^log(X_a, X_b) (computed, tables) + C^(2)_Delta(X_a, X_b) (composite cycles). C^(2)_Delta is computed '
                      'by neither engine; on the planar scalar sector (letters A_perp without derivatives) it is '
                      '-u v.(X_a X_b) per adjacent pair modulo cyclically trivial terms (theory Prop prop:tl-scalarsector, '
                      'input Thm c9-D1 / Prop c9-triangle); elsewhere open (task T3, conjecture P7)')
    sp = D1.pop('spectrum_at_theta=-1/3_under_P7', None)
    if sp is not None:
        D1['spectrum_D1_conditional_on_P7_theta=-1/3'] = sp
        D1['spectrum_D1_conditional_on_P7_theta=-1/3']['label'] = (
            'D_1 eigenvalues in units u = 4 eig(lc^natural); proved on the planar scalar sector (theta_perp = -1/3, '
            'Thm c9-kappa), conditional on P7 elsewhere (twist two contains A_V letters); scheme-independent content: '
            'gamma = 2 eig(lc^natural)')
    D1['planar_scalar_sector'] = {
        'statement': 'letters A_perp without derivatives, planar single traces, part with A_perp outputs',
        'D1_pair_density_units_u': '-8 (1 - P + K/2) - (1/2)(1/3 + theta_perp) v.(X_a X_b)  [P: exchange, K: SO(6) trace]',
        'C2_Delta_pair_density': '-u v.(X_a X_b) (symmetric part; antisymmetric part cancels in cyclic sums)',
        'one_body': "u[-8 X + (5/3 - theta_perp) v.X] per letter",
        'Xi_Delta': '-(1/3 + theta_perp) u v.',
        'at_theta_perp=-1/3': 'D_1 = 4 lc^natural; Spec D_1 = 2 x (anomalous dimensions); su(2) spectra = P6 (gate G8)',
        'status': 'derived (theory), input Thm c9-D1 / Prop c9-triangle (Ch9); not an engine evaluation of C^(2)_Delta'}
    P1 = kern['P1']
    P1['planar_scalar_sector'] = {
        'C2_Pb_pair_density': '-u d_b(X_a X_b) (symmetric part)',
        'one_body': 'u (5/3 - theta_P) d_b X per letter',
        'P1': '-(1/3 + theta_P) u {S1(P_b), -}',
        'finite_part_parameter': ('psi\' = (a/2){S1(Delta),-} (Prop c9-scheme): P_1 -> P_1 - a {S1(P_b),-}; '
                                  'P_1 = 0 at a = -(1/3 + theta_P) u; then K_b component +a, Sigma -a/2, Sigma^* +a/2, D_1 unchanged'),
        'status': 'derived (theory), input Thm c9-P1 (Ch9); C^(2)_{P_b} elsewhere: open'}
    kern.setdefault('meta', {})['revision_2026-10-29'] = ('qa/disposition.json: gates G3-W2, G3-W2n, G8f added; D_1 labels; '
                                                          'planar scalar sector entries for D_1, P_1')
    return kern, gates


if __name__ == '__main__':
    kern = json.load(open(os.path.join(HERE, 'kernels.json')))
    gates = json.load(open(os.path.join(HERE, 'gates.json')))
    kern, gates = revise(kern, gates, rerun=False)
    for d in (HERE, OUT):
        json.dump(kern, open(os.path.join(d, 'kernels.json'), 'w'), indent=1, default=str)
        json.dump(gates, open(os.path.join(d, 'gates.json'), 'w'), indent=1, default=str)
    print('gates:', sum(g['pass'] for g in gates), '/', len(gates))
