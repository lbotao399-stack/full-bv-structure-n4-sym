import os, sys
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QA = os.path.join(ROOT, 'step2/composite_qsk/qsk/engine_a/src.nosync')
for p in (QA, os.path.join(ROOT, 'step2/twoletter/engine_a/src.nosync'), os.path.join(ROOT, 'step2/composite_qsk/composite/engine_a/src.nosync')):
    if p not in sys.path:
        sys.path.insert(0, p)
