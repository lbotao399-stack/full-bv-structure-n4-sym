"""Rerun of the qa/compare entry-by-entry comparison (cmp_W0.py, read-only) on the regenerated engine A tables of the
2026-10-29 revision run; the comparator source is executed unchanged except for its output folder (this folder)."""
import os, sys, runpy
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', '..', 'qa', 'compare', 'src.nosync')
sys.path.insert(0, SRC)
name = sys.argv[1]
code = open(os.path.join(SRC, name)).read()
code = code.replace("OUT = os.path.join(L.QA, 'out.nosync')", "OUT = %r" % HERE)
assert "OUT = %r" % HERE in code
sys.argv = [name] + sys.argv[2:]
g = {'__name__': '__main__', '__file__': os.path.join(SRC, name)}
exec(compile(code, os.path.join(SRC, name), 'exec'), g)
