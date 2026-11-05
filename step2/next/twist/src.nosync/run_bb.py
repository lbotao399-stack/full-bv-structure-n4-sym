import env, sys, json, time, resource
import twdriver as D, sympy as sp
OUT = 'entries_bb.json'
try: done = json.load(open(OUT))
except Exception: done = {}
def log(s): print(s, flush=True)
specs = []
for (dU, dV) in ((D.D1, D.D2), (D.D2, D.D1)):
    tag = 'd1d2' if dU == D.D1 else 'd2d1'
    specs.append(('bb_c_b_' + tag, D.sec_dc(dU), D.sec_db(dV), (False, False)))
    specs.append(('bb_b0_g0_' + tag, D.sec_dbeta(0, dU), D.sec_dgamma(0, dV), (False, False)))
specs.append(('bc_c_c_d1d2_check', D.sec_dc(D.D1), D.sec_dc(D.D2), (False, False)))   # targeted-interpolation recheck
which = sys.argv[1:] or [s[0] for s in specs]
for nm, sU, sV, rev in specs:
    if nm not in which or nm in done: continue
    t0 = time.time(); log('== ' + nm)
    Lb = D.L_c() if nm.startswith('bc') else D.L_b()
    r = D.assemble3(D.L_b(), Lb, sU, sV, log=log, rev=rev)
    done[nm] = {'UX': str(r['UX']), 'VX': str(r['VX']), 'parts': {k: [str(x) for x in v] for k, v in r['parts'].items()}, 'time_s': round(time.time() - t0, 1)}
    json.dump(done, open(OUT, 'w'), indent=1)
    log('   UX = %s   VX = %s   (%.1fs, rss %.0f MB)' % (r['UX'], r['VX'], time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6))
# reversed projection (alternative homotopy data) on the d-b output
for nm, sU, sV, rev in (('bb_c_b_d1d2_revp', D.sec_dc(D.D1), D.sec_db(D.D2), (False, True)),
                        ('bb_c_b_d2d1_revp', D.sec_dc(D.D2), D.sec_db(D.D1), (False, True))):
    if nm in done or (sys.argv[1:] and nm not in sys.argv[1:]): continue
    t0 = time.time(); log('== ' + nm)
    r = D.assemble3(D.L_b(), D.L_b(), sU, sV, log=log, rev=rev)
    done[nm] = {'UX': str(r['UX']), 'VX': str(r['VX']), 'parts': {k: [str(x) for x in v] for k, v in r['parts'].items()}, 'time_s': round(time.time() - t0, 1)}
    json.dump(done, open(OUT, 'w'), indent=1)
    log('   UX = %s   VX = %s   (%.1fs)' % (r['UX'], r['VX'], time.time() - t0))
