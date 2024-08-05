import sympy as sp

from Tesis.Matriz_Inercia import forward_kinematics_dh

Kc = sp.Matrix([
    [sp.symbols('Kx'), 0],
    [0, sp.symbols('Ky')]
])
Cc = sp.Matrix([
    [sp.symbols('Cx'), 0],
    [0, sp.symbols('Cy')]
])


def matriz_jaco_planar(DH, a):
    n = len(a)  # Número de coordenadas generalizadas
    J = sp.zeros(3, n)  # Inicializar la matriz Jacobiana con ceros (3 filas, n columnas)
    H1 = forward_kinematics_dh(DH, 4)  # Matriz identidad de 4x4 para transformaciones homogéneas
    prev = sp.simplify(H1[0:3, 3])
    for i in range(n):
        # Vector z de la articulación i (eje z de la articulación anterior)
        z = sp.Matrix([0, 0, 1])  # z es fijo (perpendicular al plano)

        # Calcular las derivadas
        if i == 0:
            p_prev = sp.Matrix([0, 0, 0])
        else:
            H_prev = forward_kinematics_dh(DH, i + 1)
            p_prev = H_prev[0:3, 3]

        J[0:3, i] = z.cross(sp.simplify(prev - p_prev))  # Parte traslacional (solo x e y)
        J[2, i] = z[2]  # Parte rotacional (alrededor de z)

    return J


def espacio_N_Normali(vector, matrizIner):
    vector_normalizado = vector / ((vector.T * matrizIner * vector)[0])
    return vector_normalizado


def vector_N(Kg, angulos, valores):
    espacioNulo = Kg.nullspace()
    if espacioNulo:
        return espacioNulo[0].subs(angulos).subs(valores)

    else:
        print("La matriz no tiene un espacio nulo no trivial.")


def rigidez(J):
    nuva_matriz = (J.T * Kc * J)
    return nuva_matriz


def Amortiguamiento(J):
    nuva_matriz = J.transpose() * Cc * J
    return nuva_matriz


def simplify_expression(expr, tol=1e-10):
    """Simplifica una expresión simbólica eliminando términos insignificantes."""
    simplified_expr = expr.evalf()
    if abs(simplified_expr) < tol:
        return 0
    return simplified_expr


def simplify_vector(vector, tol=1e-10):
    """Simplifica cada componente de un vector simbólico."""
    return sp.Matrix([simplify_expression(comp, tol) for comp in vector])
