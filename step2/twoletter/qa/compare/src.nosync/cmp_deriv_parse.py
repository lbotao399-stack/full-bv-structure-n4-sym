"""parser of engine B derivative-table strings (copy of the one in cmp_deriv.py)"""
from fractions import Fraction
import cmp_lib as L


def parse_Bd(s, n):
    """like cmp_lib.parse_B, plus point-splitting variables s0..s3, u0..u3 (appended after the 4n momenta)"""
    out = {}
    if not s or s.strip() == '0':
        return out
    for term in s.split(' + '):
        m = L._TERM.match(term.strip())
        coef = Fraction(m.group(1)); rest = m.group(2)
        im = False; ia = ib = None; opos = [None] * n; mon = [0] * (4 * n + 8)
        facs = rest[1:].split('*') if rest.startswith('*') else rest.split('*')
        k = 0
        while k < len(facs):
            f = facs[k]; e = 1
            if not f:
                k += 1; continue
            if k + 2 < len(facs) and facs[k + 1] == '':
                e = int(facs[k + 2]); k += 2
            if f == 'I':
                im = not im
            elif f.startswith('Xa['):
                ia = f[3:-1]
            elif f.startswith('Xb['):
                ib = f[3:-1]
            elif f.startswith('Y'):
                lab = int(f[1:f.index('[')]); opos[lab] = f[f.index('[') + 1:-1]
            elif f.startswith('p'):
                lab, mm = f[1:].split('_'); mon[4 * int(lab) + int(mm)] += e
            elif f[0] in 'su' and f[1:].isdigit():
                mon[4 * n + (0 if f[0] == 's' else 4) + int(f[1:])] += e
            else:
                raise ValueError(f)
            k += 1
        val = complex(coef) * (1j if im else 1)
        d = out.setdefault((ia, ib) + tuple(opos), {})
        d[tuple(mon)] = d.get(tuple(mon), 0) + val
    return out


