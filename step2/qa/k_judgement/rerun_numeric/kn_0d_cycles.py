"""0-dim brute-force check of the one-loop cycle formulas used by kn_engine (independent of the partial files).
Fluctuation coordinates F (mixed parity) are contracted by an even propagator P (equal parities, graded symmetric)
and an odd kernel K (field-antifield, opposite parities); background coordinates B (mixed parity) and an odd constant
gamma (ghost of an even generator) are never contracted.  Every vertex J_i is quadratic in F, so every connected
graph without F-legs is a cycle.  With J = J1 + J2 + J3 (+ J4), H_i = Hessian of J_i in F (left derivatives):
  flow:          coefficient of J1 J2 J3 in hbar log(e^{hbar d_P} e^{J/hbar})  vs  (1/6) sum_{sigma in S3} str(P S H P S H P S H)
  special edge:  coefficient of J1 J2 J3 in e^{-W/hbar} e^{hbar d_P}( (1/2){J,J}_K e^{J/hbar} )  vs  c * sum_sigma str(K S H P S H P S H)
  same for n = 2 and n = 4.  {f,g}_K = Delta_K(fg)-Delta_K(f)g-(-1)^{|f|}f Delta_K(g), Delta_K = -d_K, d_phi = 1/2 sum phi_ab d_b d_a."""
import sys, random, itertools, json
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import sympy as sp
from grass import GA, PAR, newvar, dP
random.seed(2028)
fpar = [0, 0, 1, 1, 0, 1]
F = [newvar(700 + i, p) for i, p in enumerate(fpar)]
bpar = [1, 0, 1, 0, 1, 1]
B = [newvar(800 + i, p) for i, p in enumerate(bpar)]
GAM = newvar(899, 1)
def rnd(): return sp.Rational(random.randint(-4, 4) or 1, random.randint(1, 3))
P = {}; K = {}
for a in F:
    for b in F:
        if a <= b and PAR[a] == PAR[b] and not (PAR[a] == 1 and a == b):
            x = rnd(); P[(a, b)] = x; P[(b, a)] = x * (-1) ** (PAR[a] * PAR[b])
        if a < b and PAR[a] != PAR[b]:
            x = rnd(); K[(a, b)] = x; K[(b, a)] = x
def quadF(par_total, n=6):
    """random polynomial quadratic in F of given total parity (coefficients numbers)"""
    out = GA(); combos = [c for c in itertools.combinations_with_replacement(F, 2)
                          if not (PAR[c[0]] and c[0] == c[1]) and (PAR[c[0]] + PAR[c[1]]) % 2 == par_total]
    random.shuffle(combos)
    for c in combos[:n]:
        out = out + GA.var(c[0]) * GA.var(c[1]) * rnd()
    return out
def vertex(bmon, ghost=False):
    """J_i = (background monomial) * (quadratic in F) [* gamma], total even"""
    pb = sum(PAR[b] for b in bmon) + (1 if ghost else 0)
    q = quadF(pb % 2)
    t = GA.const(1)
    for b in bmon: t = t * GA.var(b)
    t = t * q
    if ghost: t = t * GA.var(GAM)
    return t
def dPg(Fx, M_): return dP(Fx, M_, F)
def brK(f, g, pf):
    DK = lambda X: dPg(X, K) * (-1)
    return DK(f * g) - DK(f) * g - (f * DK(g)) * ((-1) ** pf)
def hess(J): return {(b, c): J.dleft(c).dleft(b) for b in F for c in F}
def MSH(M_, H):
    out = {}
    for a in F:
        for c in F:
            s_ = GA()
            for b in F:
                if M_.get((a, b), 0) != 0:
                    s_ = s_ + H[(b, c)] * (M_[(a, b)] * (-1) ** PAR[b])
            out[(a, c)] = s_
    return out
def mm(X_, Y_): return {(a, c): sum((X_[(a, b)] * Y_[(b, c)] for b in F), GA()) for a in F for c in F}
def strace(X_): return sum((X_[(a, a)] * (-1) ** PAR[a] for a in F), GA())
def Fdeg0(X): return GA({k: sp.expand(v) for k, v in X.t.items() if not any(x in F for x in k) and sp.expand(v) != 0})
def same(A_, B_):
    keys = set(A_.t) | set(B_.t)
    return all(sp.simplify(A_.t.get(k, 0) - B_.t.get(k, 0)) == 0 for k in keys)
def fit(A_, B_):
    c = sp.symbols('c'); keys = set(A_.t) | set(B_.t)
    if not keys: return 'both sides zero (vacuous)'
    return sp.solve([sp.expand(A_.t.get(k, 0) - c * B_.t.get(k, 0)) for k in keys], [c], dict=True)
def multilinear(X, tags):
    """part of X containing each tag variable exactly once (tags = list of B variables / GAM marking the vertices)"""
    return GA({k: v for k, v in X.t.items() if all(k.count(t) == 1 for t in tags)})
res = {}
def run(n, ghost_at):
    # vertices tagged by distinct background variables; vertex i has background monomial (B[i],) or (B[i], B[i+3])
    bm = BM[:n]
    Js = [vertex(bm[i], ghost=(i == ghost_at)) for i in range(n)]
    tags = [m[0] for m in bm] + ([GAM] if ghost_at is not None else [])
    J = GA()
    for j in Js: J = J + j
    Hs = [hess(j) for j in Js]
    # brute flow: coefficient of J1..Jn in W = hbar log(...): connected part of d_P^n(J1...Jn)/n! * hbar^{n} hbar^{-n} -> order hbar^1 * (1/hbar)... 
    prod = GA.const(1)
    for j in Js: prod = prod * j
    Z = prod
    for k in range(n): Z = dPg(Z, P) * sp.Rational(1, k + 1)
    # e^{hbar dP} e^{J/hbar} multilinear part at F-degree 0: hbar^n/n! dP^n (J1..Jn / hbar^n) = dP^n(prod)/n!; W = hbar log -> hbar * (connected)
    # disconnected pieces multilinear in all J's at F-degree 0 are products of cycles over a partition: subtract products of lower cycles.
    # For the check we use n vertices on one cycle; subtract the products of sub-cycles (each vertex is quadratic: no tadpoles needed? d_P J_i can be nonzero)
    flow_brute = multilinear(Fdeg0(Z), tags)
    # formula for the full (possibly disconnected) Z: exp of cycles; compute the connected part via the formula and compare the
    # full Z with the exponential of the formula (multilinear part).
    def cyc_sum(seq_idx, first=None):
        tot = GA()
        for perm in itertools.permutations(seq_idx):
            Ms = [MSH(P if (first is None or i > 0) else first, Hs[perm[i]]) for i in range(len(perm))]
            X = Ms[0]
            for M_ in Ms[1:]: X = mm(X, M_)
            tot = tot + strace(X)
        return tot
    # connected cycle values C(S) = (1/(2|S|)) sum_perm str((PSH)...) for subsets S
    idx = list(range(n))
    def Cflow(S): return cyc_sum(S) * sp.Rational(1, 2 * len(S))
    # Z multilinear = sum over set partitions of idx of product of Cflow(block)
    def set_partitions(s):
        if not s: yield []; return
        first, rest = s[0], s[1:]
        for p in set_partitions(rest):
            yield [[first]] + p
            for i in range(len(p)):
                yield p[:i] + [[first] + p[i]] + p[i + 1:]
    def Zformula():
        tot = GA()
        for part in set_partitions(idx):
            t = GA.const(1)
            for blk in part: t = t * Cflow(blk)
            tot = tot + t
        return multilinear(Fdeg0(tot), tags)
    res['n=%d ghost=%s flow: brute Z == exp-of-cycles formula' % (n, ghost_at)] = same(flow_brute, Zformula())
    # special edge: e^{hbar dP}( D e^{J/hbar} ) with D = (1/2){J,J}_K, multilinear in J1..Jn, F-degree 0.
    # D multilinear part involves pairs (i,j): {J_i,J_j}_K ; then remaining vertices closed by P.
    D = brK(J, J, 0) * sp.Rational(1, 2)
    others = GA.const(1)
    Zs = D * prod   # contains J's twice; we need multilinear overall: D (uses 2 J's) times (n-2) remaining J's
    # build explicitly: sum over pairs i<j of {J_i,J_j}_K * prod_{k not in {i,j}} J_k, then dP^{n-1}/(n-1)! (one contraction is K)
    tot = GA()
    for i, j in itertools.combinations(idx, 2):
        t = brK(Js[i], Js[j], 0)
        for k in idx:
            if k not in (i, j): t = t * Js[k]
        tot = tot + t
    Zsp = tot
    for k in range(n - 1): Zsp = dPg(Zsp, P) * sp.Rational(1, k + 1)
    sp_brute_full = multilinear(Fdeg0(Zsp), tags)
    # connected special-edge cycle value (candidate) Cs(S) = c * sum_perm str(K S H P S H ...) over S; disconnected:
    # Zsp = sum over partitions with one block carrying the K edge: Cs(block) * prod Cflow(other blocks)
    def Cs(S, c): return cyc_sum(S, first=K) * c
    cc = sp.symbols('cc')
    def Zsp_formula(c):
        tot = GA()
        for part in set_partitions(idx):
            for bi, blk in enumerate(part):
                if len(blk) < 2: continue
                t = GA.const(1)
                for bj, b2 in enumerate(part):
                    t = t * (Cs(b2, c) if bj == bi else Cflow(b2))
                tot = tot + t
        return multilinear(Fdeg0(tot), tags)
    res['n=%d ghost=%s special: fit c in Zsp = Cs(c)*flows' % (n, ghost_at)] = str(fit(sp_brute_full, Zsp_formula(1)))
# background monomials: B[0] odd, B[1] even, B[2] odd, B[3] even, B[4] odd, B[5] odd
CHOICES = {
 (2, None): [(B[0],), (B[1],)],            # parities 1,0 -> total 1
 (2, 0):    [(B[1],), (B[3],)],            # 0,0 + ghost -> 1
 (2, 1):    [(B[0], B[2]), (B[1],)],       # 0,0 + ghost -> 1
 (3, None): [(B[0],), (B[2],), (B[4],)],   # 1,1,1 -> 1
 (3, 0):    [(B[0],), (B[1],), (B[2],)],   # 1,0,1 + ghost -> 1
 (3, 1):    [(B[1],), (B[0],), (B[3],)],   # 0,1,0 + ghost -> 0 ... adjust below
 (3, 2):    [(B[0],), (B[2],), (B[1],)],   # 1,1,0 + ghost -> 1
 (4, 3):    [(B[0],), (B[1],), (B[2],), (B[3],)],  # 1,0,1,0 + ghost -> 1
 (4, None): [(B[0],), (B[1],), (B[2],), (B[3],)],  # 1,0,1,0 -> 0 (vacuous control)
}
CHOICES[(3, 1)] = [(B[1],), (B[0],), (B[3], B[2])]   # 0,1,1 + ghost -> 1
for (n, g), bmc in CHOICES.items():
    BM = bmc
    run(n, g)
print(json.dumps(res, indent=1), flush=True)
json.dump(res, open('kn_0d_cycles.json', 'w'), indent=1)

