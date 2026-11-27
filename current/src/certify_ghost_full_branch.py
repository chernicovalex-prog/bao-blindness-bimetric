#!/usr/bin/env python3
"""Exact interval certificate for the reduced benchmark branch (epsilon=0).

The reduced polynomial has no radiation and no Lambda correction. All sign
decisions use exact rationals, not floats. This is not an exact proof for the
full benchmark model with finite epsilon and radiation.
"""
from fractions import Fraction
import json
import mpmath as mp
import sympy as sp

mp.mp.dps = 100
t_mp = mp.tan(mp.mpf('6.1') * mp.pi / 180) ** 2
scale = 10 ** 35
lo = Fraction(int(mp.floor(t_mp * scale)), scale)
hi = Fraction(int(mp.ceil(t_mp * scale)), scale)
assert mp.mpf(lo.numerator) / lo.denominator < t_mp < mp.mpf(hi.numerator) / hi.denominator and lo > 0

y, t = sp.symbols('y t')
F2 = 3*y**2 - 6*y + 4
g = 8*y**3 - 9*y**2 + 4
F3 = 3*t*y**2*F2 + g
A = 9*t*y**4 - 18*t*y**3 + 12*t*y**2 + 8*y**3 - 9*y**2 + 4
assert sp.expand(F3-A) == 0
# F2 = 3(y-1)^2 + 1 > 0; g >= g(3/4) = 37/16 > 0.
assert sp.simplify(F2 - (3*(y-1)**2 + 1)) == 0
assert sp.simplify(g.subs(y, sp.Rational(3, 4)) - sp.Rational(37, 16)) == 0

# N4 is derived directly from the reduced quartic, not from mapper shorthand.
u = sp.Symbol('u')
P = (3*t*y**4 + (-9*t+4)*y**3 + (12*t-9)*y**2
     + (t*(-6+3*u)+9)*y - 4)
N4 = sp.expand(3*t*u*y - P)
F4 = 4 + (6*t-5)*y + (4-6*t)*y**2 + 3*t*y**3
assert sp.expand(N4 - (1-y)*F4) == 0

# F4 = Q(y) + 3 t y^3 >= Q(y), Q = 4 + (6t-5)y + (4-6t)y^2.
y0 = (5-6*t)/(2*(4-6*t))
Qmin = sp.simplify(4 - (5-6*t)**2/(4*(4-6*t)))
for endpoint in (lo, hi):
    assert 0 < endpoint * y0.subs(t, endpoint) < 1
assert hi < sp.Rational(1, 2)
# Qmin increases here: its derivative is positive for 0<t<1/2.
dQmin = sp.diff(Qmin, t)
assert sp.simplify(dQmin - 3*(5-6*t)*(3-6*t)/(2*(4-6*t)**2)) == 0
qmin_lo = sp.simplify(Qmin.subs(t, lo))
assert qmin_lo > 0

result = {
    'status': 'PASS',
    'claim': 'For the benchmark branch, F2>0, F3>0, F4>0 on 0<y<1; hence N4=(1-y)F4>=0 and dy/dN=3yN4/F3>=0 on that open branch.',
    't_tan2_6.1deg_decimal': mp.nstr(t_mp, 45),
    't_interval_exact': [str(lo), str(hi)],
    'g_min': '37/16',
    'Qmin_at_t_lo': str(qmin_lo),
    'Qmin_at_t_lo_decimal': mp.nstr(mp.mpf(qmin_lo.numerator) / qmin_lo.denominator, 30),
    'proof_scope': 'reduced algebraic benchmark branch (epsilon=0, radiation-free) 0<y<1; not an exact full-model proof, not full cosmological viability and not X>=0 everywhere',
}
print(json.dumps(result, indent=2, ensure_ascii=False))
