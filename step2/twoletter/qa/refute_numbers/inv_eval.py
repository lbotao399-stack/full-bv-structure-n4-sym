"""evaluate engine A index-free strings (invariant_forms_bosonic_W0) on explicit 10-vectors (data, not code)"""
import re
import sympy as sp

def ev(form, slots, rs):
    """form: string like '(3*I)*d(b,o1) a.r1 + (-3*I)*d(a,o1) b.r1'; slots: name -> 10-vector; rs: list of 4-vectors"""
    def d(x, y):
        return sum(slots[x][m] * slots[y][m] for m in range(4))
    def dW(x, y):
        return sum(slots[x][m] * slots[y][m] for m in range(4, 10))
    def sdot(tok):
        u, v = tok.split('.')
        def vec(n):
            if re.fullmatch(r'r\d+', n):
                return list(rs[int(n[1:]) - 1]) + [0]*6
            return slots[n]
        U, V = vec(u), vec(v)
        return sum(U[m] * V[m] for m in range(4))
    expr = re.sub(r'\bI\b', 'IMAG', form)
    expr = re.sub(r'(\d+)/(\d+)', r'RAT(\1,\2)', expr)
    expr = re.sub(r'\((-?)(\d+)\)', r'(\1INT(\2))', expr)
    # tokenise terms: replace d(..), dW(..), x.y with python calls
    expr = re.sub(r'dW\((\w+),(\w+)\)', r'dW("\1","\2")', expr)
    expr = re.sub(r'(?<![\w"])d\((\w+),(\w+)\)', r'd("\1","\2")', expr)
    expr = re.sub(r'\b([a-z]+\d*)\.([a-z]+\d*)\b', r'sdot("\1.\2")', expr)
    expr = re.sub(r'\)\s+(d|dW|sdot)\(', r')*\1(', expr)
    return sp.expand(eval(expr, {'d': d, 'dW': dW, 'sdot': sdot, 'IMAG': sp.I, 'RAT': sp.Rational, 'INT': sp.Integer}))
