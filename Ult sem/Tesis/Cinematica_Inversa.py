import numpy as np
import sympy as sp


def rrrIK(L1, L2, L3, p):
    px = p[0]
    py = p[1]
    phi = p[2]

    pwx = px - (L3 * np.cos(phi))
    pwy = py - (L3 * np.sin(phi))

    c2 = ((pwx ** 2) + (pwy ** 2) - (L1 ** 2) - (L2 ** 2)) / (2 * L1 * L2)

    if c2 >= 1 or c2 <= -1:
        print('not reachable')
        return None

    s2 = -np.sqrt(1 - c2 ** 2)
    theta2 = np.arctan2(s2, c2)

    s1 = ((L1 + L2 * c2) * pwy - L2 * s2 * pwx) / (pwx ** 2 + pwy ** 2)
    c1 = ((L1 + L2 * c2) * pwx + L2 * s2 * pwy) / (pwx ** 2 + pwy ** 2)
    theta1 = np.arctan2(s1, c1)

    theta3 = phi - theta1 - theta2

    return np.array([theta1, theta2, theta3])

