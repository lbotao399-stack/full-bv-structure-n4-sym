"""Use the live manuscript/ch00_overview.tex when it can be read; when iCloud has evicted it and the read fails,
substitute manuscript/build/backup_pre_fix3/ch00_overview.tex (read-only use, outputs go to this folder only)."""
import subprocess
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
CH00 = ROOT / 'manuscript' / 'ch00_overview.tex'
BACKUP = ROOT / 'manuscript' / 'build' / 'backup_pre_fix3' / 'ch00_overview.tex'

def ch00_readable():
    try:
        r = subprocess.run(['head', '-c', '64', str(CH00)], capture_output=True, timeout=30)
        return r.returncode == 0 and len(r.stdout) > 0
    except subprocess.TimeoutExpired:
        return False

def install():
    if ch00_readable():
        return 'live ch00'
    orig = Path.read_text
    def rt(self, *a, **k):
        if self.resolve() == CH00.resolve():
            return orig(BACKUP, *a, **k)
        return orig(self, *a, **k)
    Path.read_text = rt
    return 'backup ch00 (live file evicted)'
