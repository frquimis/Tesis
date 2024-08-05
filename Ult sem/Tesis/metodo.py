import sympy as sp

from Tesis.wf import simplify_expression


def Matriz_S(Vector, Matriz):
    q1, q2, q3 = sp.symbols('q1 q2 q3')
    q = sp.Matrix([q1, q2, q3])
    u0T_M = Vector.T * Matriz
    simplificado = sp.sympify(u0T_M * q)

    equation = sp.Eq(simplificado[0], 0)

    # Resolver para q3 en términos de q1 y q2
    solution = sp.solve(equation, q3)
    Matriz_S1 = sp.zeros(3, 2)
    Matriz_S1[0, 0] = 1
    Matriz_S1[1, 1] = 1
    valorq1 = solution[0].subs({sp.symbols('q2'): 0}) / sp.symbols('q1')
    valorq2 = solution[0].subs({sp.symbols('q1'): 0}) / sp.symbols('q2')
    Matriz_S1[2, 0] = valorq1
    Matriz_S1[2, 1] = valorq2
    return Matriz_S1


def qzb(impulso, Uo, C):
    t, s = sp.symbols('t s')
    qimpulse = Uo.T*impulso

    x=(Uo.T * C * Uo)[0]
    x=simplify_expression(x)

    denominator = s ** 2 + x * s
    sp.pprint(qimpulse[0])

    laplace_function = qimpulse[0] / denominator

    # Aplicar la transformada inversa de Laplace
    beta_t = sp.inverse_laplace_transform(laplace_function, s, t)
    sp.pprint(beta_t)
    qzbr = sp.simplify(beta_t) * Uo
    return qzbr


# NZP motions
def NZP(Mm, R, A, S1):
    Mp = (S1.T * Mm * S1)
    Cp = (S1.T * R * S1)
    Kp = (S1.T * A * S1)

    O = sp.zeros(2, 2)
    I = sp.eye(2)
    nueva = O.row_join(I)
    MMinv = Mp.inv()

    mm3 = (-MMinv) * Kp
    mm4 = (-MMinv) * Cp

    partInfe = mm3.row_join(mm4)
    A = nueva.col_join(partInfe)
    B = O.col_join(MMinv)
    X1 = A.eigenvects()
    eigen_matrix = sp.Matrix.hstack(*[vects[0] for val, mult, vects in X1])
    Y1 = (eigen_matrix.inv()).T

    return eigen_matrix, Y1, X1, B


def matrize(valores):
    valores_propios = [val for val, mult, vects in valores]
    exp_eigenvalues = [sp.exp(val * sp.symbols('t')) for val in valores_propios]
    matriz = sp.diag(*exp_eigenvalues)
    return matriz


def xt(X1, Y1, valores, B, S, Q):
    matrizval = matrize(valores)
    x = X1 * matrizval * Y1.T * B * S.T * Q
    return x
