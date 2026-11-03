#!/usr/bin/env python3
"""Engine A (route A), strict Costello route, action side: single entry point.
Usage: python3 run_all.py [--only T1,T4] [--force]
Resumable: every task result is cached in out.nosync/task_<id>.json and skipped when present (unless --force).
Writes results.json, gates.json (and copies in out.nosync/), appends to PROGRESS.txt, prints wall time and max RSS.
A watchdog thread aborts the run when the resident set exceeds 1.4 GB (hard rule < 1.5 GB)."""
import os, sys, json, time, threading, resource, subprocess, datetime
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'src.nosync'); OUT = os.path.join(HERE, 'out.nosync')
sys.path.insert(0, SRC)
os.makedirs(OUT, exist_ok=True)
LIMIT = 1.4e9
T0 = time.time()


def log(msg):
    line = '%s %s' % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'), msg)
    print(line, flush=True)
    with open(os.path.join(HERE, 'PROGRESS.txt'), 'a') as f:
        f.write(line + '\n')


def rss():
    try:
        return int(subprocess.check_output(['ps', '-o', 'rss=', '-p', str(os.getpid())]).strip()) * 1024
    except Exception:
        return 0


def watchdog():
    while True:
        r = rss()
        if r > LIMIT:
            log('WATCHDOG: RSS %.2f GB > 1.4 GB, aborting' % (r / 1e9))
            os._exit(3)
        time.sleep(5)


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items() if k not in ('raw', 'pole_raw')}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, (bool, int, float, str)) or x is None:
        return x
    return str(x)


def main():
    args = sys.argv[1:]
    only = None; force = '--force' in args
    if '--only' in args:
        only = args[args.index('--only') + 1].split(',')
    threading.Thread(target=watchdog, daemon=True).start()
    import sa_tasks as T
    import sa_tasks2 as T2
    tasks = [('T0', T.task_T0), ('T1', T.task_T1), ('T4', T.task_T4), ('B1', T.task_B1),
             ('T3', T2.task_T3), ('T6', T2.task_T6), ('B2', T2.task_B2), ('T2', T2.task_T2), ('T5', T2.task_T5),
             ('T8', T2.task_T8)]
    results, gates = {}, []
    for tid, fn in tasks:
        cache = os.path.join(OUT, 'task_%s.json' % tid)
        if only is not None and tid not in only and os.path.exists(cache):
            d = json.load(open(cache)); results[tid] = d['results']; gates += d['gates']; continue
        if only is not None and tid not in only:
            continue
        if os.path.exists(cache) and not force:
            d = json.load(open(cache)); results[tid] = d['results']; gates += d['gates']
            log('task %s: cached' % tid); continue
        t = time.time(); log('task %s: start' % tid)
        try:
            r, g = fn()
        except Exception as e:
            import traceback; traceback.print_exc()
            log('task %s: FAILED %r' % (tid, e)); continue
        r, g = jsonable(r), jsonable(g)
        json.dump({'results': r, 'gates': g}, open(cache, 'w'), indent=1)
        results[tid] = r; gates += g
        log('task %s: done in %.1f s, RSS %.0f MB, gates %s' % (tid, time.time() - t, rss() / 1e6,
                                                                 ['%s:%s' % (x['id'], 'pass' if x['pass'] else 'FAIL') for x in g]))
        write(results, gates)
    write(results, gates)
    mx = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    log('run_all: wall %.1f s, max RSS %.0f MB, gates %d/%d pass' % (time.time() - T0, mx / 1e6 if mx > 1e7 else mx / 1e3,
                                                                     sum(1 for x in gates if x['pass']), len(gates)))


def write(results, gates):
    meta = {'engine': 'A (route A: closed special-edge cycles in scaled-Schwinger/simplex form)', 'unit': 'u',
            'date': datetime.date.today().isoformat(), 'spec': 'step2/strict/action/theory/engine_spec.json'}
    for nm, obj in (('results.json', {'meta': meta, 'results': results}), ('gates.json', gates)):
        for d in (HERE, OUT):
            json.dump(obj, open(os.path.join(d, nm), 'w'), indent=1)


if __name__ == '__main__':
    main()
