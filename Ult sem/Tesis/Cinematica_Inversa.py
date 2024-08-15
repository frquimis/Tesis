import numpy as np
import sympy as sp
import scipy.optimize as opt


def Cinematica_Inversa_3R(xd, yd):
    # Longitudes de los eslabones
    L1 = 0.5
    L2 = 0.35
    L3 = 0.3

    # Definir la función de error basada en la cinemática directa
    def error_func(thetas):
        theta1, theta2, theta3 = thetas
        x = L1 * np.cos(theta1) + L2 * np.cos(theta1 + theta2) + L3 * np.cos(theta1 + theta2 + theta3)
        y = L1 * np.sin(theta1) + L2 * np.sin(theta1 + theta2) + L3 * np.sin(theta1 + theta2 + theta3)
        return np.array([x - xd, y - yd])

    # Solución inicial para los ángulos
    initial_guess = np.array([1.659, -1.979, 1.105])

    # Resolvemos el problema de optimización
    solution = opt.least_squares(error_func, initial_guess)

    theta1, theta2, theta3 = solution.x
    arreglo = sp.Matrix([[theta1], [theta2], [theta3]])

    return arreglo
