import mpmath as mp
mp.mp.dps=30
p2=1.0
def I(f,eps):
    # integrate over [eps,inf)^2 using substitution symmetric
    g=lambda t1,t2: f(t1,t2)
    return mp.quad(lambda t1: mp.quad(lambda t2: g(t1,t2),[eps,10*eps,1,mp.inf]),[eps,10*eps,1,mp.inf])
f3=lambda a,b: mp.e**(-a*b*p2/(a+b))/(a+b)**3
f2=lambda a,b: mp.e**(-a*b*p2/(a+b))/(a+b)**2
f4=lambda a,b: a*b*mp.e**(-a*b*p2/(a+b))/(a+b)**4
for name,f,a,b in [('T^-3',f3,0.25,1/6.),('T^-2',f2,0,-1),('t1t2T^-4',f4,0,-1/6.)]:
    vals=[]
    for eps in [1e-3,1e-4]:
        v=I(f,eps)-a/eps-b*mp.log(eps)
        vals.append(v)
    print(name,'remainder at eps=1e-3,1e-4:',[mp.nstr(x,8) for x in vals])
# Konishi conversion
from fractions import Fraction as F
print('unit su(N): hbar*2N/(16pi^2 * (-4hbar/g^2)) = -g^2N * ', F(2,16*4))
print('-24 * (-1/32) =', F(-24)*F(-1,32), ' i.e. 3/4 -> 3 g^2N/(4pi^2)')
