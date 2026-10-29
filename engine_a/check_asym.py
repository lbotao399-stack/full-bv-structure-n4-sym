"""numerical spot check of the exact eps-asymptotics (schwinger.py) for g(u)=1+2u+3u^2, c=5, s=p^2=1.3"""
import mpmath as mp
mp.mp.dps = 30
s = mp.mpf('1.3'); c = 5
g = lambda u: 1 + 2 * u + 3 * u ** 2
def J(eps):
    # t1 = T(1-u), t2 = T u, dt1 dt2 = T dT du ; region t1,t2 >= eps
    def inner(T):
        lo, hi = eps / T, 1 - eps / T
        if lo >= hi: return 0
        return mp.quad(lambda u: (g(u) / T ** 2 + c / T ** 3) * T * mp.e ** (-s * T * u * (1 - u)), [lo, 0.5, hi])
    return mp.quad(inner, [2 * eps, 10 * eps, 1, 10, mp.inf])
pred_inv = mp.mpf(c) / 4; pred_log = -(1 + 1 + 1) + s * c / 6
e1, e2 = mp.mpf('1e-3'), mp.mpf('1e-4')
J1, J2 = J(e1), J(e2)
# J = pred_inv/eps + pred_log*log(eps) + C + O(eps log eps)
C1 = J1 - pred_inv / e1 - pred_log * mp.log(e1); C2 = J2 - pred_inv / e2 - pred_log * mp.log(e2)
print('constant at eps=1e-3, 1e-4:', C1, C2, 'difference', C1 - C2)
