import env, sys, json, time, resource
import twdriver as D, sympy as sp
OUT = 'entries.json'
try: done = json.load(open(OUT))
except Exception: done = {}
def save(): json.dump(done, open(OUT, 'w'), indent=1, default=str)
def log(s): print(s, flush=True)
R = sp.Rational
specs = []
for (dU, dV) in ((D.D1, D.D2), (D.D2, D.D1)):
    tag = 'd1d2' if dU == D.D1 else 'd2d1'
    specs.append(('betagamma_00_' + tag, D.L_beta(0), D.L_gamma(0), R(5, 2), 0, D.sec_dc(dU), D.sec_dc(dV)))
    specs.append(('betabeta_01_c_g2_' + tag, D.L_beta(0), D.L_beta(1), R(3), 0, D.sec_dc(dU), D.sec_dgamma(2, dV)))
    specs.append(('bgamma_0_c_g0_' + tag, D.L_b(), D.L_gamma(0), R(3), 0, D.sec_dc(dU), D.sec_dgamma(0, dV)))
    specs.append(('betab_0_b0_c_' + tag, D.L_beta(0), D.L_b(), R(7, 2), 0, D.sec_dbeta(0, dU), D.sec_dc(dV)))
    specs.append(('betab_0_g1_g2_' + tag, D.L_beta(0), D.L_b(), R(7, 2), 0, D.sec_dgamma(1, dU), D.sec_dgamma(2, dV)))
which = sys.argv[1:] or [s[0] for s in specs]
for nm, La, Lb, w, d, sU, sV in specs:
    if nm not in which or nm in done: continue
    t0 = time.time(); log('== ' + nm)
    r = D.assemble2(La, Lb, sU, sV, log=log)
    done[nm] = {'UX': str(r['UX']), 'VX': str(r['VX']), 'parts': {k: [str(x) for x in v] for k, v in r['parts'].items()},
                'time_s': round(time.time() - t0, 1)}
    save()
    log('   UX = %s   VX = %s   (%.1fs, rss %.0f MB)' % (r['UX'], r['VX'], time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6))
# ---- arity-two entries
specs2 = []
for (dU, dV) in ((D.D1, D.D2), (D.D2, D.D1)):
    tag = 'd1d2' if dU == D.D1 else 'd2d1'
    specs2.append(('bc_c_c_' + tag, D.L_b(), D.L_c(), D.sec_dc(dU), D.sec_dc(dV)))
for nm, La, Lb, sU, sV in specs2:
    if (sys.argv[1:] and nm not in sys.argv[1:]) or nm in done: continue
    t0 = time.time(); log('== ' + nm)
    r = D.assemble2(La, Lb, sU, sV, log=log)
    done[nm] = {'UX': str(r['UX']), 'VX': str(r['VX']), 'parts': {k: [str(x) for x in v] for k, v in r['parts'].items()},
                'time_s': round(time.time() - t0, 1)}
    save()
    log('   UX = %s   VX = %s   (%.1fs, rss %.0f MB)' % (r['UX'], r['VX'], time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6))
