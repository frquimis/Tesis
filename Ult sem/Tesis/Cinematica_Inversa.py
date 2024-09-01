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


def forward_kinematics_3R(theta1, theta2, theta3, l1, l2, l3):
    # Posición del primer eslabón
    x1 = l1 * np.cos(theta1)
    y1 = l1 * np.sin(theta1)

    # Posición del segundo eslabón
    x2 = x1 + l2 * np.cos(theta1 + theta2)
    y2 = y1 + l2 * np.sin(theta1 + theta2)

    # Posición del efector final
    x3 = x2 + l3 * np.cos(theta1 + theta2 + theta3)
    y3 = y2 + l3 * np.sin(theta1 + theta2 + theta3)

    # Ángulo de orientación del efector final
    phi = theta1 + theta2 + theta3

    return x3, y3, phi
