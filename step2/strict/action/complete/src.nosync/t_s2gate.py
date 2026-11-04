"""gate GB-S2: Q S^(2)_2(xi,xi') + S_2(P_{Gamma|V}) - {S_2(xi), S_2(xi')} = 0 on two legs (ch02 line 173 / Lemma
c2-components), both chirality pairs."""
import numpy as np, random
import cb_common as B, cb_tree as T, cb_s2 as S
Q = B.Q
random.seed(11)
def rspin(rng):
    x = np.zeros(16, complex)
    for c in rng: x[c] = random.randint(-2, 2)
    return x
res = {}
for nm, r1, r2 in (('(+,-)', B.LAMP, B.LAMM), ('(+,+)', B.LAMP, B.LAMP), ('(-,-)', B.LAMM, B.LAMM), ('(-,+)', B.LAMM, B.LAMP)):
    x1, x2 = rspin(r1), rspin(r2)
    G12 = S.gam(x1, x2)
    lhs = Q.dsum((1, Q.apply_Q(S.s2_dens(x1, x2))), (1, S.translation_dens(G12[:4])),
                 (-1, T.bracket(S.S2xi_dens(x1), 0, S.S2xi_dens(x2))))
    alt = Q.dsum((1, Q.apply_Q(S.s2_dens(x1, x2))), (1, S.translation_dens(G12[:4])),
                 (1, T.bracket(S.S2xi_dens(x1), 0, S.S2xi_dens(x2))))
    w = wa = sz = 0
    for spcs in (['A', 'As'], ['As', 'A'], ['lam', 'lams'], ['lams', 'lam'], ['c', 'cs'], ['A', 'A'], ['lam', 'A'], ['c', 'As'], ['lam', 'lam'], ['A', 'cs']):
        for _ in range(3):
            legs = [Q.rand_leg(s) for s in spcs]; legs[1]['p'] = [-x for x in legs[0]['p']]
            w = max(w, abs(lhs(legs))); wa = max(wa, abs(alt(legs)))
            sz = max(sz, abs(S.translation_dens(G12[:4])(legs)), abs(Q.apply_Q(S.s2_dens(x1, x2))(legs)))
    res[nm] = (w, wa, sz, list(np.round(np.abs(G12), 3)))
    print(nm, 'residual', w, ' (opposite bracket sign:', wa, ')  scale', sz, ' |Gamma|', np.abs(G12[:4]).max(), np.abs(G12[4:]).max())
