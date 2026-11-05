import env, sector as S, sympy as sp, time
t0 = time.time()
for lab, chi, t in (('b', S.wt5((-1, 1)), 4),):
    R = S.cohomology(chi, t, maxdeg=2)
    print(lab, 'basis t:', [S.name(b) for b in R['B0']]); print('  t-1:', [S.name(b) for b in R['Bm']]); print('  t+1:', [S.name(b) for b in R['Bp']])
    print('  dim Z', len(R['Z']), 'dim B', len(R['Bim']))
    for z in R['Z']: print('  Z:', {S.name(R['B0'][i]): z[i] for i in range(len(z)) if z[i] != 0})
    for z in R['Bim']: print('  B:', {S.name(R['B0'][i]): z[i] for i in range(len(z)) if z[i] != 0})
print(time.time() - t0)
