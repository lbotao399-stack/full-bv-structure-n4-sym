"""QA (action side): replay engine B's T3 random stream (seed 23) without recomputing cycles, rebuild the leg
configurations, verify the replay on the stored V-sector candidate values, then fit B's stored raw lambda-lambda
values of Ob^(1)(M_omega)|_2 against Q of the spinor spin part of S(M_omega) in B's own realisation.
Read-only use of engine B modules."""
import sys, os, json, random
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
BSRC = os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b', 'src.nosync'))
sys.path.insert(0, BSRC)
import eb_qb as Q
import eb_t1

raw = json.load(open(os.path.join(BSRC, '..', 'out.nosync', 'T3.json')))['M']
VPAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

def omega_fn(om6):
    om = np.zeros((4, 4), complex)
    for v, (i, j) in zip(om6, VPAIRS):
        om[i, j] = v; om[j, i] = -v
    return om

random.seed(23)
cfg = {}
for sec in ['AVAV', 'ApAp', 'lamlam', 'cAsV']:
    cfg[sec] = []
    for ic in range(3):
        om6 = [random.randint(-2, 2) for _ in range(6)]
        legs_at = eb_t1.make_legs(sec)
        cfg[sec].append((om6, legs_at([Fr(0)] * 4)))

def cand_A(om, legs0):
    d = Q.dens(1, ('A', 'As'), lambda s: s[0]['u'][:4] @ om @ s[1]['u'][:4])
    return Q.to_fr(Q.apply_Q(d)(legs0))

ok = True
for sec in ['AVAV', 'cAsV']:
    for (om6, legs0), rec in zip(cfg[sec], raw[sec]):
        c = cand_A(omega_fn(om6), legs0)
        same = [str(c[0]), str(c[1])] == rec['cand']
        ok &= same
        print('replay check', sec, c, rec['cand'], same)
print('replay reproduces stored V-sector candidates:', ok)

R, RB, MP = Q.R, Q.RB, Q.MP
def spin_lam_density(om):
    # S(M_omega) spinor part in B's template: -( (1/4) sum r_MN RB[M] R[N] lam, lam* ), r = -omega
    r = -om
    M = sum(0.25 * r[a, b] * (RB[a] @ R[b]) for a in range(4) for b in range(4) if a != b)
    return Q.dens(1, ('lam', 'lams'), lambda s: -((M @ s[0]['u']) @ MP @ s[1]['u']))

rows = []
for (om6, legs0), rec in zip(cfg['lamlam'], raw['lamlam']):
    c = Q.to_fr(Q.apply_Q(spin_lam_density(omega_fn(om6)))(legs0))
    re = Fr(rec['raw_eps0']['re'].get('1', '0')); im = Fr(rec['raw_eps0']['im'].get('1', '0'))
    theory = (-re, -im)                     # theory component = -(raw) for even x (engine B normalisation)
    rows.append((c, theory))
    print('lamlam cand', c, 'theory', theory)
ratios = []
for (c, t) in rows:
    cz = complex(float(c[0]), float(c[1])); tz = complex(float(t[0]), float(t[1]))
    ratios.append(tz / cz if cz != 0 else None)
print('theory/cand ratios:', ratios)
json.dump({'replay_ok': ok, 'lamlam': [{'cand': [str(c[0]), str(c[1])], 'theory': [str(t[0]), str(t[1])]} for c, t in rows],
           'ratios': [str(r) for r in ratios]},
          open(os.path.join(HERE, '..', 'q_mspin_replay.json'), 'w'), indent=1)
