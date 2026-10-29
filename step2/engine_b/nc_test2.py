import json, time, sympy as sp, numpy as np, mpmath
import numeric_check as NC, ia_engine as IA, sectors3 as S
legs = S.config('ccAsV', 0)
terms = []
IA.evaluate(legs, collect=terms)
raw = json.load(open('ia_raw_ccAsV.json'))['0']['classes']
def exval(d):
    return complex(sum(complex(sp.N(sp.sympify(v))) * float({'1':1,'log2':mpmath.log(2),'log3':mpmath.log(3)}[b]) for b, v in d.items()))
for eps in (1e-3, 1e-4, 1e-5, 1e-6):
    tot, pc = NC.ia_numeric(terms, eps, 20, 200)
    print('eps', eps, 'total', tot, ' diffs per class:', [round((v - exval(raw[cl].get('0', {}))).imag, 6) for cl, v in pc.items()], flush=True)
