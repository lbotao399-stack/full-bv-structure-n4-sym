import cProfile, pstats, time
import sympy as sp
import ce_job as CJ
I = sp.I
ZP = {4: 1, 5: I}; ZO = {4: sp.Rational(1, 2), 5: -I / 2}
n = [1, I, 0, 0]; nbar = [sp.Rational(1, 2), -I / 2, 0, 0]
job = CJ.CJob('AW', 'AW', ['AW', 'AW'], pol=[ZO, ZO], za=ZP, zb=ZP, deriv=('null', 2, n, nbar))
seq = [(0,), (1,)]
def run():
    t0 = time.time()
    types, curs, words = job.config(seq, 1, 0, 1, 'D2')
    t1 = time.time()
    for w, Num in words.items():
        print('Num terms', len(Num))
        co = job.integrate(Num, curs, types, (0,))
    print('config', t1 - t0, 'integrate', time.time() - t1)
cProfile.run('run()', '../out.nosync/prof.out')
p = pstats.Stats('../out.nosync/prof.out'); p.sort_stats('cumulative').print_stats(18)
