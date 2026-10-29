"""Run reader_build/prep_supplements.py with its output (supplement_blocks.json) redirected to step2/qa/ch11_review."""
import os, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
REV = HERE.parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / 'reader_build'))
import prep_supplements
prep_supplements.HERE = REV
ALT = '/Users/libotao/.npm/_npx/4b4c857f6efdfb61/node_modules'
try:
    open(os.path.join(prep_supplements.NODE_PATH, 'uc.micro', 'package.json')).read(10); print('node_modules: live')
except OSError:
    prep_supplements.NODE_PATH = ALT; print('node_modules: npx cache')
prep_supplements.main()
