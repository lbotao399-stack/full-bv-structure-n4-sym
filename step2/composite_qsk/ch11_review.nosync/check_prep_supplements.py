"""Run reader_build/prep_supplements.py with its output (supplement_blocks.json) redirected to this folder."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(ROOT / 'reader_build'))
import ch00_guard
print('ch00:', ch00_guard.install())
import prep_supplements
prep_supplements.HERE = HERE
import os
ALT = '/Users/libotao/.npm/_npx/4b4c857f6efdfb61/node_modules'
if os.path.exists(os.path.join(prep_supplements.NODE_PATH, 'uc.micro', 'package.json')):
    try:
        open(os.path.join(prep_supplements.NODE_PATH, 'uc.micro', 'package.json')).read(10); print('node_modules: live')
    except OSError:
        prep_supplements.NODE_PATH = ALT; print('node_modules: npx cache (live copy evicted)')
prep_supplements.main()
