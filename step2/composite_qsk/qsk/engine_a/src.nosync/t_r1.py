import qa_s1 as S1, qa_r1 as R1, qa_util as U
from fractions import Fraction
B=S1.chir_basis('lamP')
xi=sum(int(c)*v for c,v in zip([1,2,-1,1,0,1,-2,1],B)).astype(complex)
a3=U.vzero(); a3[1]=(Fraction(-4),Fraction(0)); a3[2]=(Fraction(4),Fraction(0))
R1.run(xi,a3,npts=2)
