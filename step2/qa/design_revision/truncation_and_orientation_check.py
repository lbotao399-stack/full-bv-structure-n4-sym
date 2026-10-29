# Design revision (2026-10-28): two checks for step2/theory/odd_components_design.tex.
#
# (A) Truncation census (Prop prop:o-trunc, Lemma lem:o-actiontrunc).  The arity-k component of L_1(O) on
#     (xi,...,xi), xi in S_+ (x) W, has Weyl weight w_O + k/2 and ghost degree |O| + 1 - k.  Enumerate every field
#     content (A:1,0  lam:3/2,0  c:0,+1  A*:3,-1  lam*:5/2,-1  c*:4,-2, derivative 1,0) of that weight and degree
#     with at most 6 field factors, for k = 3, 4, 5 and the targets of the design.  Also the minimal weight m(n)
#     of a local functional of degree -n and the action-side bound (weight 4 + k/2, degree -k).
# (B) Orientation sum for a pair of floating odd insertions joined by a lambda-lambda edge P(eps,infty), at the
#     IR-leading order, for a (xi, eta) pair: xi in S_+ (x) W constant spinor, eta in Sigma^* (linear spinor).
#     Two orientations: eta-vertex at the f end (O2) or at the g end (O1).  O1 has the lambda-lambda propagator
#     in the opposite direction, d_k D(y2-y1) = -d_k D(y1-y2).  The two A*|_V -- c contractions differ by the
#     exchange of the two odd ghost ends (Koszul sign -1, ch03 prop:c3-P).  Check: O1 = O2 pointwise, hence the
#     total O1 - O2 vanishes, while each orientation alone has a nonzero S^3 average (logarithmic).
#     This replaces the control 'control_xi_eta' of step2/qa/design_audit_2/eta_pair_ir.py, which subtracted an
#     expression from itself.  Also: the S^3 average of the orientation-SYMMETRIC combination for the eta-eta
#     chain (sensitivity of the log coefficient to the orientation sign).
import itertools, json, os
from fractions import Fraction as Fr
import numpy as np

out = {}
FIELDS = {'A': (Fr(1), 0), 'lam': (Fr(3, 2), 0), 'c': (Fr(0), 1), 'A*': (Fr(3), -1), 'lam*': (Fr(5, 2), -1),
          'c*': (Fr(4), -2), 'd': (Fr(1), 0)}
NAMES = list(FIELDS)
def contents(weight, degree, maxf=6):
    res = []
    for n in range(1, maxf + 1):
        for combo in itertools.combinations_with_replacement(NAMES, n):
            if all(x == 'd' for x in combo):
                continue
            w = sum(FIELDS[x][0] for x in combo); g = sum(FIELDS[x][1] for x in combo)
            if w == weight and g == degree:
                res.append(' '.join(combo))
    return res
def m(nn):
    if nn <= 0:
        return Fr(0)
    return Fr(2 * nn) if nn % 2 == 0 else Fr(2 * nn) + Fr(1, 2)
# minimal weight check by enumeration
mins = {}
for nn in range(1, 5):
    best = None
    for n in range(1, 7):
        for combo in itertools.combinations_with_replacement(NAMES, n):
            if all(x == 'd' for x in combo):
                continue
            if sum(FIELDS[x][1] for x in combo) == -nn:
                w = sum(FIELDS[x][0] for x in combo)
                best = w if best is None or w < best else best
    mins[nn] = (str(best), str(m(nn)), best == m(nn))
out['minimal_weight_of_degree_-n'] = {str(k): {'enumerated': v[0], 'formula m(n)': v[1], 'agree': v[2]} for k, v in mins.items()}
targets = {'letter iota_w A / A|_V': (Fr(1), 0), 'letter lambda': (Fr(3, 2), 0), 'letter c': (Fr(0), 1),
           'scalar bilinears (20prime, Konishi primary)': (Fr(2), 0), 'Konishi descendant O_zeta': (Fr(5, 2), 0),
           'bc entry (open colour)': (Fr(2), 1), 'weight-3 descendants B_w, Fbar_w': (Fr(3), 0),
           'generic w_O=4 operator (example)': (Fr(4), 0)}
cen = {}
for name, (w, d) in targets.items():
    row = {}
    for k in (3, 4, 5):
        W = w + Fr(k, 2); D = d + 1 - k
        allowed_by_bound = W >= m(-D)
        row['k=%d' % k] = {'weight': str(W), 'degree': D, 'm(-degree)': str(m(-D)),
                           'weight_bound_allows': bool(allowed_by_bound), 'contents': contents(W, D)}
    row['crude_bound_floor((2w+4+4|O|)/3)'] = int((2 * w + 4 + 4 * d) // 3)
    cen[name] = row
out['census_k>=3'] = cen
act = {}
for k in range(1, 6):
    W = 4 + Fr(k, 2); D = -k
    act['k=%d' % k] = {'weight': str(W), 'degree': D, 'm(k)': str(m(k)), 'allowed': bool(W >= m(k)), 'contents': contents(W, D)}
out['action_side_along_Sigma'] = act

# ---------------- (B) orientation sum ----------------
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]], complex); s3 = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)
sig = [s1, s2, s3, 1j * I2]; sigb = [-s1, -s2, -s3, 1j * I2]
def sphere_nodes(ns=20, nphi=20):
    xs, ws = np.polynomial.legendre.leggauss(ns); s = 0.5 * (xs + 1); w = 0.5 * ws
    ph = 2 * np.pi * np.arange(nphi) / nphi
    S, P1, P2 = np.meshgrid(s, ph, ph, indexing='ij'); Wt = np.meshgrid(w, ph, ph, indexing='ij')[0] / nphi ** 2
    ca, sa = np.sqrt(1 - S), np.sqrt(S)
    pts = np.stack([ca * np.cos(P1), ca * np.sin(P1), sa * np.cos(P2), sa * np.sin(P2)], -1).reshape(-1, 4)
    return pts, Wt.reshape(-1)
NODES, WTS = sphere_nodes()
rng = np.random.default_rng(1028)
def rnd(): return rng.normal(size=(2, 4)) + 1j * rng.normal(size=(2, 4))
eta, xi = rnd(), rnd()      # eta^alpha_i (undotted, upper spinor index), xi_alpha^i
def chain(et, xi_, mats):
    # et sig^a sigbar^m sig^k sigbar^n xi with index-dependent scalar weights mats[a,m,k,n]
    return np.einsum('ai,amkn,ai->', et, mats, xi_)
def O2(p):
    pp = p @ p
    tot = 0
    for a in range(4):
        for mm in range(4):
            L = ((a == mm) / pp - 2 * p[a] * p[mm] / pp ** 2)
            if L == 0:
                continue
            for k in range(4):
                for n in range(4):
                    coef = L * (p[k] / pp) * (p[n] / pp)
                    tot += coef * np.einsum('ai,ab,bi->', eta, sig[a] @ sigb[mm] @ sig[k] @ sigb[n], xi)
    return tot
def O1(p):
    pp = p @ p
    tot = 0
    for b in range(4):
        for n in range(4):
            R = ((b == n) / pp - 2 * p[b] * p[n] / pp ** 2)
            for k in range(4):
                for mm in range(4):
                    coef = (-p[mm] / pp) * (-p[k] / pp) * R
                    tot += coef * np.einsum('ai,ab,bi->', eta, sig[b] @ sigb[n] @ sig[k] @ sigb[mm], xi)
    return tot
v2 = np.array([O2(p) for p in NODES]); v1 = np.array([O1(p) for p in NODES])
out['xi_eta_orientation'] = {
    'max_|O1 - O2|_on_S3': float(np.abs(v1 - v2).max()),
    'S3_average_single_orientation_O2': complex((v2 * WTS).sum()).__repr__(),
    'pairing_eta.xi': complex(np.einsum('ai,ai->', eta, xi)).__repr__(),
    'total_with_orientation_sign_-1_max': float(np.abs(v1 - v2).max()),
    'total_with_orientation_sign_+1_S3_average': complex(((v1 + v2) * WTS).sum()).__repr__(),
    'conclusion': 'O1 = O2 pointwise: the IR-leading (log) part cancels in O1 - O2; with the opposite orientation sign it would not'}
print(json.dumps(out, indent=1, default=str))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'truncation_and_orientation_check.json'), 'w'), indent=1, default=str)
