import time, sys, json, pickle
import sympy as sp
import ce_job as CJ, ce_fam as CF
I = sp.I
Smax = int(sys.argv[1]) if len(sys.argv) > 1 else 2
OPT = int(sys.argv[2]) if len(sys.argv) > 2 else 3
NW = (sys.argv[3] == 'w') if len(sys.argv) > 3 else False
ZP = {4: 1, 5: I}; ZO = {4: sp.Rational(1, 2), 5: -I / 2}
n = [1, I, 0, 0]; nbar = [sp.Rational(1, 2), -I / 2, 0, 0]
NV = {0: 1, 1: I}
t0 = time.time()
job = CJ.CJob('AW', 'AW', ['AW', 'AW'], pol=[ZO, ZO], za=ZP, zb=ZP, deriv=('null', Smax, n, nbar), nullw=NW)
res = job.kernels(mlist=(0, 1), option=OPT, verbose=True)
print('ZZ done', round(time.time() - t0, 1), flush=True)
hP, okP, wP = CF.two_output_density(job, res['Pb'], 0, avec=n, smax=Smax)
hD, okD, wD = CF.two_output_density(job, res['D0'], 0, smax=Smax)
print('logs absent', okP, okD, 'w-free', wP, wD)
poles = {m: sum(len(d) for dm in res['Pb'].values() for mm, d in dm.items() if mm == m) for m in (1, 2)}
print('pole entries (Pb) m=1,2:', poles)
one = {}
for tag, sa, sb, za, zb in (('AZ', 'AV', 'AW', NV, ZP), ('ZA', 'AW', 'AV', ZP, NV)):
    j1 = CJ.CJob(sa, sb, ['AW'], pol=[ZO], za=za, zb=zb, deriv=('null', Smax, n, nbar), nullw=NW)
    r1 = j1.kernels(mlist=(0, 1), option=OPT, verbose=False)
    one[tag + 'P'], ok1, w1 = CF.one_output_density(j1, r1['Pb'], 0, avec=n, smax=Smax - 1)
    one[tag + 'D'], ok2, w2 = CF.one_output_density(j1, r1['D0'], 0, smax=Smax - 1)
    print(tag, 'logs absent', ok1, ok2, 'w-free', w1, w2, round(time.time() - t0, 1), flush=True)
pickle.dump(dict(hP=hP, hD=hD, one=one), open(f'../out.nosync/tw2_S{Smax}.pkl', 'wb'))
import subprocess
print(subprocess.run(['python3', 't_ibp.py', f'../out.nosync/tw2_S{Smax}.pkl', str(Smax)], capture_output=True, text=True).stdout)
for S in range(0, Smax + 1):
    sym = [(a, S - a) for a in range(S // 2 + 1)]
    symo = [(a, S + 1 - a) for a in range((S + 1) // 2 + 1)]
    MP = sp.zeros(len(symo), len(sym)); MD = sp.zeros(len(sym), len(sym))
    TP = sp.zeros(len(symo), len(sym))
    for j, (a, b) in enumerate(sym):
        rP = CF.twist_two_map(hP, one['AZP'], one['ZAP'], a, b)
        rD = CF.twist_two_map(hD, one['AZD'], one['ZAD'], a, b)
        for i, (c, d) in enumerate(symo):
            MP[i, j] = sp.nsimplify(rP.get((c, d), 0) + (rP.get((d, c), 0) if c != d else 0))
        for i, (c, d) in enumerate(sym):
            MD[i, j] = sp.nsimplify(rD.get((c, d), 0) + (rD.get((d, c), 0) if c != d else 0))
        # target -1/3 d_n : tr(d^a Z d^b Z) -> tr(d^{a+1} Z d^b Z) + tr(d^a Z d^{b+1} Z)
        tg = {}
        for (c, d) in ((a + 1, b), (a, b + 1)):
            key = (min(c, d), max(c, d))
            tg[key] = tg.get(key, 0) + sp.Rational(-1, 3)
        for i, (c, d) in enumerate(symo):
            TP[i, j] = tg.get((c, d), 0)
    print('S\'=', S, '\n  sum C_Pn =', MP.tolist(), '\n  target  =', TP.tolist(), '\n  sum C_D0 =', MD.tolist(), ' target', -sp.Rational(2, 3) * (4 + S))
print('total time', round(time.time() - t0, 1))
