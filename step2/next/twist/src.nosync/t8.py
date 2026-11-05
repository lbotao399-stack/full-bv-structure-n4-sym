import env, twdriver as D, sector as S, sympy as sp
for d, lab in ((D.D1, 'd1'), (D.D2, 'd2')):
    B0, pi, R = D.projection(*D.sec_db(d))
    print('db', lab, 'dimZ', len(R['Z']), 'dimB', len(R['Bim']), 'basis/pi:', [(S.name(b), p) for b, p in zip(B0, pi) if p != 0])
# the (xi, lam*) sector: charge q_xi (= 0 mod q_xi), t = 4
B0, Bm = S.sector_basis(S.wt5(), 4, 3), S.sector_basis(S.wt5(), 3, 3)
R = S.cohomology(S.wt5(), 4, 2)
print('charge-0 t=4 sector:', [S.name(b) for b in R['B0']], 'dimZ', len(R['Z']), 'dimB', len(R['Bim']))
for z in R['Z']: print('  Z', {S.name(R['B0'][i]): z[i] for i in range(len(z)) if z[i] != 0})
for z in R['Bim']: print('  B', {S.name(R['B0'][i]): z[i] for i in range(len(z)) if z[i] != 0})
