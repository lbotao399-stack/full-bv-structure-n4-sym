"""Part D4 (1b): the arity-one obstruction along Delta on the quadratic sector, from the UV boundary.

Flow of the arity-one defect at two legs, order hbar (right-ghost convention, J = I - I^CT[eps] - F_0 + (S_M - I^CT_Delta)xi):
   1PI one-loop graphs with S_M and one special edge Khat = K_eps - K_0 on the cycle      (class B1, exact via masters)
 + the K_0 vertices: -Q I^CT_Delta (= 0, Q X = 0) - {I^CT_2 + F_0, S_M} = +9 unit eps^-1 int(A,A)
 + S_M joined by P to the flowed arity-zero two-leg defect (-> 0 after F_0) + K̂-bridges (-> 0).
Canonical form on (leg 1 at p, leg 2 at -p-kappa, xi at kappa): alpha = G_w - i sum_mu d/dkappa_mu G_mu|_0,
beta^mu = -i G_mu|_0.  Exact coefficients of eps^-1 and eps^0, in Q + Q log(E/(E-1)) per graph."""
import json, sys, time, random
import sympy as sp
import delta_engine as D

I = sp.I
P = list(sp.symbols('P0:4')); KAP = list(sp.symbols('K0:4'))
SECT = {'AVAV': ('A', 'A', list(range(4)), list(range(4))),
        'ApAp': ('A', 'A', list(range(4, 10)), list(range(4, 10))),
        'lamlam': ('lam', 'lam', list(range(16)), list(range(16))),
        'cAsV': ('c', 'As', [0], list(range(4)))}

def rand_pol(sp_, comps, rnd):
    u = [sp.Integer(0)] * D.V.DIM[sp_]
    for c in comps:
        u[c] = sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3]))
    return u

def numP(rnd):
    return [sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3])) for _ in range(4)]

def B1(sec, u1, u2, pnum, verbose=False):
    """returns dict n -> (alpha, beta) with alpha, beta = {'1','L'} numbers (complex rationals) at the numeric momentum pnum."""
    s1, s2, c1, c2 = SECT[sec]
    out = {}
    def run(terms, kap):
        legs = [dict(sp=s1, u=u1, p=list(pnum)), dict(sp=s2, u=u2, p=[-pnum[m] - kap[m] for m in range(4)]),
                dict(sp='xi', u=[sp.Integer(1)], p=list(kap))]
        parts = [[('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC), ('SM', (2,), terms)]]
        if s1 == 'A' and s2 == 'A':
            parts.append([('quartic', (0, 1), D.V.QUARTIC), ('SM', (2,), terms)])
        recs = []
        for part in parts:
            recs += D.enumerate_terms(legs, part, 'K', [0, 1])
        return recs
    # weight part
    acc = {-1: {'1': 0, 'L2': 0, 'L32': 0}, 0: {'1': 0, 'L2': 0, 'L32': 0}}
    per = {}
    for r in run(D.SW_TERMS, [0] * 4):
        v = D.special_coeffs(r['Num'], r['curs'], r['s'])
        Lk = 'L2' if len(r['curs']) == 2 else 'L32'
        for n in (-1, 0):
            acc[n]['1'] += r['fac'] * v[n]['1']; acc[n][Lk] += r['fac'] * v[n]['L']
    alpha = {n: dict(acc[n]) for n in acc}
    beta = {n: [{'1': 0, 'L2': 0, 'L32': 0} for _ in range(4)] for n in (-1, 0)}
    for mu in range(4):
        recs = run(D.T_TERMS(mu), KAP)
        for r in recs:
            Lk = 'L2' if len(r['curs']) == 2 else 'L32'
            vd = D.special_coeffs(r['Num'], r['curs'], r['s'], D.kap_derivative_op(KAP, mu))
            vz = D.special_coeffs(r['Num'], r['curs'], r['s'], D.kap_zero_op(KAP))
            for n in (-1, 0):
                alpha[n]['1'] += r['fac'] * vd[n]['1']; alpha[n][Lk] += r['fac'] * vd[n]['L']
                beta[n][mu]['1'] += -I * r['fac'] * vz[n]['1']; beta[n][mu][Lk] += -I * r['fac'] * vz[n]['L']
    for n in (-1, 0):
        alpha[n] = {k: sp.expand(v) for k, v in alpha[n].items()}
        beta[n] = [{k: sp.expand(v) for k, v in b.items()} for b in beta[n]]
    return alpha, beta

if __name__ == '__main__':
    secs = sys.argv[1].split(',')
    nsamp = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    rnd = random.Random(23)
    out = {}
    for sec in secs:
        s1, s2, c1, c2 = SECT[sec]
        out[sec] = []
        for k in range(nsamp):
            t0 = time.time()
            u1 = rand_pol(s1, c1, rnd); u2 = rand_pol(s2, c2, rnd); pn = numP(rnd)
            al, be = B1(sec, u1, u2, pn)
            rec = dict(u1=[str(x) for x in u1], u2=[str(x) for x in u2], p=[str(x) for x in pn],
                       alpha={str(n): {kk: str(v) for kk, v in al[n].items()} for n in al},
                       beta={str(n): [{kk: str(v) for kk, v in b.items()} for b in be[n]] for n in be}, time=round(time.time() - t0, 1))
            out[sec].append(rec)
            print(sec, k, json.dumps(rec), flush=True)
            json.dump(out, open(f'delta_B1_{"_".join(secs)}.json', 'w'), indent=1)
