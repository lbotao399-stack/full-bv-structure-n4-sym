"""Radial derivatives of D_t(r) = (1 - e^{-r^2/4t})/(4 pi^2 r^2) (checked against mp.diff in localisation.py)."""
import mpmath as mp


def dD(r, t):
    a = 1 / (4 * t)
    E = mp.e**(-a * r * r)
    return (2 * a * E / r - 2 * (-mp.expm1(-a * r * r)) / r**3) / (4 * mp.pi**2)


def ddD(r, t):
    a = 1 / (4 * t)
    E = mp.e**(-a * r * r)
    return (-4 * a * a * E - 6 * a * E / r**2 + 6 * (-mp.expm1(-a * r * r)) / r**4) / (4 * mp.pi**2)
