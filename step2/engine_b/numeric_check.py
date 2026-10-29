"""Engine B, Step 2: independent internal check -- the FULL regulated three-leg sum at finite eps, numerically, sector
(c, c, A*|_V), at two small eps and generic momenta, compared with the predicted limit Ob_3 + tree(Ob_2, I_3).
  I-a      : every special-edge triangle at finite eps: Schwinger integral (u_special in [0,1], others in [1,oo)) of the exact
             Gaussian-integrated integrand, NO Taylor expansion (tensor Gauss-Legendre in z = log u).
  I-b+III_Khat : -B_{Khat}(I_3, W_1^(2)(eps) - I^CT_2(eps)); W_1^(2) = regulated one-loop (c, A*) two-point function
             = i(k.a*) s_W(k^2, eps) (unique invariant), numerically (2-d Schwinger integral); I^CT_2 = -log eps int(dc,A*).
  I-d+II-b : B_{P(eps,oo)}(D_2(eps), I_3),  D_2(eps) = special bubbles - Q I^CT_2 = r(eps k^2) Ob_2 (unique invariant),
             r computed numerically from the special bubbles (2-d Schwinger integral).
  II-a, III_K0: exact (both zero in this sector: Q Y_3 = {I_3, Y_2} = 0 on (c, c, A*|_V)).
Units: unit * ([c_0,c_1],c_2)."""
import json, os, sys, time, math
import numpy as np
import sympy as sp
import ia_engine as IA
import loops as L
import fa
import sectors3 as S
import others as O
import vertices as V

HERE = os.path.dirname(os.path.abspath(__file__))
q = L.q
BET = IA.BET

def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def prep(Num, curs, special, with_k2):
    """Gaussian moments of Num (times -k_special^2 if with_k2) -> numpy callables C_k(beta), rho(beta)."""
    E = len(curs)
    Nt = Num
    if with_k2:
        ks = [q[m] - curs[special][m] for m in range(4)]
        Nt = sp.expand(-Num * sum(x * x for x in ks))
    shift = [sum(BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
    C = L.gaussian_moments(Nt, shift)
    rho = sp.expand(sum(BET[e] * sum(x * x for x in curs[e]) for e in range(E)) - sum(x * x for x in shift))
    Cf = {k: sp.lambdify(BET[:E], sp.expand(v), 'numpy') for k, v in C.items() if v != 0}
    rf = sp.lambdify(BET[:E], rho, 'numpy')
    return Cf, rf

def ev(f, *args):
    v = f(*args)
    return np.broadcast_to(np.asarray(v, dtype=complex), np.broadcast(*args).shape)

# ------------------------------------------------------------------ triangles with one special edge
def ia_numeric(terms, eps, n1=24, nz=160):
    """sum over collected I-a terms of fac * graph(eps), units (16 pi^2)^-1 (exact integrand, numerical quadrature)."""
    zmax = math.log(1.0 / eps) + 9.0
    u1, w1 = gl(n1, 0.0, 1.0)
    z, wz = gl(nz, 0.0, zmax)
    U1, Z2, Z3 = np.meshgrid(u1, z, z, indexing='ij')
    W = np.einsum('i,j,k->ijk', w1, wz, wz)
    U2 = np.exp(Z2); U3 = np.exp(Z3)
    W = W * U2 * U3
    T = U1 + U2 + U3
    tot = 0j
    per_class = {}
    # two-edge terms (quartic-vertex bubbles): u_special in [0,1], u_other in [1,oo)
    V1, Y2 = np.meshgrid(u1, z, indexing='ij')
    V2 = np.exp(Y2)
    W2 = np.einsum('i,j->ij', w1, wz) * V2
    T2 = V1 + V2
    for t in terms:
        if 'prep' not in t:
            t['prep'] = prep(t['Num'], t['curs'], t['s'], True)
        Cf, rf = t['prep']
        E = len(t['curs']); s = t['s']
        others = [e for e in range(E) if e != s]
        if E == 3:
            b = [None] * E
            b[s] = U1 / T; b[others[0]] = U2 / T; b[others[1]] = U3 / T
            TT, WW = T, W
        else:
            b = [None] * E
            b[s] = V1 / T2; b[others[0]] = V2 / T2
            TT, WW = T2, W2
        rho = ev(rf, *b).real
        acc = np.zeros_like(TT, dtype=complex)
        for k, f in Cf.items():
            acc += ev(f, *b) * (2 * eps * TT) ** (-k)
        val = complex(np.sum(WW * eps ** (E - 2) * TT ** (-2) * acc * np.exp(-eps * TT * rho)))
        val *= complex(t['fac'])
        tot += val
        per_class[t['class']] = per_class.get(t['class'], 0) + val
    return tot, per_class

# ------------------------------------------------------------------ two-point functions at finite eps
def bubble_terms(e1, e2):
    """the regulated one-loop two-point function W_1^(2) (all bubbles; no quartic tadpole for (c, A*)): list of
    (weight, Num, curs).  Same word formula as Step 1 (engine_b/sectors.py, bubble)."""
    out = []
    for (x, y, sgn) in [(e1, e2, 1), (e2, e1, (-1) ** (V.PAR[e1['sp']] * V.PAR[e2['sp']]))]:
        px = x['p']
        h1 = L.full_hessian(V.CUBIC, [x], [-v for v in q], [q[m] - px[m] for m in range(4)], {'f': 1})
        h2 = L.full_hessian(V.CUBIC, [y], [px[m] - q[m] for m in range(4)], list(q), {'f': 1})
        Num = L.str_product([L.prop_numerator(list(q)), L.PiH(h1), L.prop_numerator([q[m] - px[m] for m in range(4)]), L.PiH(h2)])
        if Num != 0:
            out.append((sp.Rational(1, 4) * sgn, Num, [[0] * 4, list(px)]))
    return out

def bubble_numeric(terms, eps, nz=200):
    zmax = math.log(1.0 / eps) + 9.0
    z, wz = gl(nz, 0.0, zmax)
    Z1, Z2 = np.meshgrid(z, z, indexing='ij')
    U1 = np.exp(Z1); U2 = np.exp(Z2)
    W = np.einsum('i,j->ij', wz, wz) * U1 * U2
    T = U1 + U2
    tot = 0j
    for (wt, Num, curs) in terms:
        Cf, rf = prep(Num, curs, 0, False)
        b = [U1 / T, U2 / T]
        rho = ev(rf, *b).real
        acc = np.zeros_like(T, dtype=complex)
        for k, f in Cf.items():
            acc += ev(f, *b) * (2 * eps * T) ** (-k)
        tot += complex(wt) * complex(np.sum(W * T ** (-2) * acc * np.exp(-eps * T * rho)))
    return tot

def special_bubble_numeric(terms, eps, n1=40, nz=200):
    zmax = math.log(1.0 / eps) + 9.0
    u1, w1 = gl(n1, 0.0, 1.0)
    z, wz = gl(nz, 0.0, zmax)
    U1, Z2 = np.meshgrid(u1, z, indexing='ij')
    U2 = np.exp(Z2)
    W = np.einsum('i,j->ij', w1, wz) * U2
    T = U1 + U2
    tot = 0j
    for t in terms:
        if 'prep' not in t:
            t['prep'] = prep(t['Num'], t['curs'], t['s'], True)
        Cf, rf = t['prep']
        s = t['s']; o = 1 - s
        b = [None, None]; b[s] = U1 / T; b[o] = U2 / T
        rho = ev(rf, *b).real
        acc = np.zeros_like(T, dtype=complex)
        for k, f in Cf.items():
            acc += ev(f, *b) * (2 * eps * T) ** (-k)
        tot += complex(t['fac']) * complex(np.sum(W * T ** (-2) * acc * np.exp(-eps * T * rho)))
    return tot

# ------------------------------------------------------------------ the check
def run(conf_ids=(0, 1), epss=(1e-4, 1e-5, 1e-6), out_json='numeric_check_ccAsV.json'):
    I = sp.I
    res = {}
    t0 = time.time()
    # reference two-point momentum for the scalar functions s_W and r (they depend on k^2 only)
    for cid in conf_ids:
        legs = S.config('ccAsV', cid)
        terms = []
        IA.evaluate(legs, collect=terms)
        exact_ia = sp.Integer(0)
        raw = json.load(open(os.path.join(HERE, 'ia_raw_ccAsV.json')))
        G4v = complex(S.G_value('G4', legs))
        tree_val = complex(sp.N(O.fval(fa.Bfun(fa.kernel_P(0), O.OB2, fa.I3), legs)))
        rc = {}
        for eps in epss:
            ia_val, ia_cls = ia_numeric(terms, eps)
            # --- the scalar functions s_W(k^2,eps) and r(eps k^2) at the momenta needed by the splittings
            cache_s = {}; cache_r = {}
            def s_W(k2v, kv):
                key = round(float(k2v), 9)
                if key not in cache_s:
                    a = [sp.Integer(x) for x in (1, 2, -1, 3)]
                    e1 = dict(sp='c', u=[sp.Integer(1)], p=list(kv)); e2 = dict(sp='As', u=a + [sp.Integer(0)] * 6, p=[-x for x in kv])
                    val = bubble_numeric(bubble_terms(e1, e2), eps)
                    norm = complex(sp.N(I * sum(kv[m] * a[m] for m in range(4))))      # int(dc, A*) on (c at k, A* at -k)
                    cache_s[key] = val / norm
                return cache_s[key]
            def r_D2(k2v, kv):
                key = round(float(k2v), 9)
                if key not in cache_r:
                    a = [sp.Integer(x) for x in (2, 1, 1, 3)]
                    l2 = [dict(sp='c', u=[sp.Integer(1)], p=list(kv)), dict(sp='A', u=a + [sp.Integer(0)] * 6, p=[-x for x in kv])]
                    tt = []
                    IA.evaluate(l2, collect=tt)
                    sb = special_bubble_numeric(tt, eps)
                    pa = complex(sp.N(sum(kv[m] * a[m] for m in range(4)))); k2 = float(k2v)
                    D2 = sb + 4.5 / eps * 1j * pa             # - Q I^CT_2 = + (9/2) eps^-1 int(dc, A|_V)
                    cache_r[key] = D2 / (-27.0 / 8 * 1j * pa * k2)
                return cache_r[key]
            # --- I-b + III_Khat:  -B_{(e^{-eps k^2}-1) F}(I_3, (s_W(k^2) + log eps) int(dc,A*))
            Kib = []
            for X, Y, f in fa.kernel_K0():
                def g(qv, f=f):
                    k2 = sum(x * x for x in qv)
                    return f(qv) * sp.Float(complex(math.expm1(-eps * float(k2)) * (s_W(k2, qv) + math.log(eps))).real) \
                        + f(qv) * sp.I * sp.Float(complex(math.expm1(-eps * float(k2)) * (s_W(k2, qv) + math.log(eps))).imag)
                Kib.append((X, Y, g))
            ib = fa.B(Kib, fa.I3, fa.local([(1, 'kingh')], 'int(dc,A*)', 1), legs).get('f', 0)
            ib = -complex(sp.N(ib))
            # --- I-d + II-b:  B_{P(eps) r(eps k^2)}(Ob_2, I_3)
            Pid = []
            for X, Y, f in fa.kernel_P(0):
                def g2(qv, f=f):
                    k2 = sum(x * x for x in qv)
                    rr = r_D2(k2, qv) * math.exp(-eps * float(k2))
                    return f(qv) * (sp.Float(rr.real) + sp.I * sp.Float(rr.imag))
                Pid.append((X, Y, g2))
            idv = complex(sp.N(fa.B(Pid, O.OB2, fa.I3, legs).get('f', 0)))
            total = ia_val + ib + idv
            rc[str(eps)] = dict(I_a=str(ia_val), I_a_per_class={k: str(v) for k, v in ia_cls.items()},
                                Ib_plus_III_Khat=str(ib), Id_plus_IIb=str(idv), IIa_plus_III_K0='0 (exact)',
                                total=str(total), total_over_G4=str(total / G4v),
                                s_W={str(k): str(v) for k, v in cache_s.items()}, r_D2={str(k): str(v) for k, v in cache_r.items()},
                                seconds=round(time.time() - t0, 1))
            print(cid, eps, rc[str(eps)], flush=True)
        res[str(cid)] = dict(G4_value=str(G4v), predicted_limit_total_over_G4='-27/4 (= Ob_3 [0] + tree of Ob_2 with I_3 [-27/4])',
                             tree_over_G4=str(tree_val / G4v), runs=rc)
    json.dump(res, open(os.path.join(HERE, out_json), 'w'), indent=1)
    return res

if __name__ == '__main__':
    run()
