import sympy as sp
from tensor import calcular_tensor_inercia
from sympy import nsimplify
from rigidez import rigidez, Amortiguamiento, eliminar_filas_redundantes


def denavit_hartenberg_com(a, Lcom):
    """
    Create the Denavit-Hartenberg parameter matrix for the center of mass.

    Parameters:
    a (list): List of link lengths.
    Lcom (list): List of center of mass distances.

    Returns:
    comDH (sympy Matrix): Denavit-Hartenberg parameter matrix for the center of mass.
    """
    comDH = sp.Matrix([
        [0, 0, 0, 0],
        [a[0], 0, Lcom[0], 0],
        [a[1], 0, Lcom[1], 0],
        [a[2], 0, Lcom[2], 0]
    ])

    return comDH


def dh_transform(theta, d, a, alpha):
    """
    Create the Denavit-Hartenberg transformation matrix.
    """
    return sp.Matrix([
        [sp.cos(theta), -sp.sin(theta) * sp.cos(alpha), sp.sin(theta) * sp.sin(alpha), a * sp.cos(theta)],
        [sp.sin(theta), sp.cos(theta) * sp.cos(alpha), -sp.cos(theta) * sp.sin(alpha), a * sp.sin(theta)],
        [0, sp.sin(alpha), sp.cos(alpha), d],
        [0, 0, 0, 1]
    ])


def forward_kinematics_dh(DH, up_to_joint):
    """
    Calculate the forward kinematics for a chain up to a given joint.
    """
    T = sp.eye(4)
    for i in range(up_to_joint):
        theta, d, a, alpha = DH.row(i)
        T *= dh_transform(theta, d, a, alpha)
    return T


def forward_kinematics_dh_com(DH, comDH, rb):
    """
    Calculate the forward kinematics to the center of mass of a requested rigid body.
    """
    Lcom = sp.symbols(f'Lcom{rb}')

    frame = None
    for i in range(comDH.shape[0]):
        if comDH[i, 2] == Lcom:
            frame = i
            break

    if frame is None:
        raise ValueError(f"Lcom{rb} not found in comDH matrix.")

    # Forward kinematics to the previous reference frame
    H = forward_kinematics_dh(DH, frame - 1)

    # Forward kinematics to the requested center of mass
    Rz = dh_transform(comDH[frame, 0], 0, 0, 0)
    Tz = dh_transform(0, comDH[frame, 1], 0, 0)
    Tx = dh_transform(0, 0, comDH[frame, 2], 0)

    Hcom = sp.simplify(H * Rz * Tz * Tx)

    return Hcom


def matriz_jaco(DH, comDH, a, rb):
    Lcom1 = sp.symbols(f'Lcom{rb}')
    comFrame = None
    for i in range(comDH.shape[0]):
        if comDH[i, 2] == Lcom1:
            comFrame = i
            break
    # Calcular la cinemática directa al centro de masa solicitado
    Hcom = forward_kinematics_dh_com(DH, comDH, rb)
    # Obtener el número de coordenadas generalizadas
    n = len(a)

    # Inicializar la matriz Jacobiana con ceros
    Jcom = sp.zeros(6, n)
    com1 = comDH[:, 0]

    for i in range(1, n + 1):  # i toma valores 1, 2 y 3
        symbol_to_find = sp.Symbol(f'a{i}')

        jointFrame = next((j for j, val in enumerate(com1) if val == symbol_to_find), None)
        # print(f'Para i = {i}, jointFrame = {jointFrame + 1}')

        # Si la articulación actual está en un marco de referencia que no afecta el COM
        if jointFrame > comFrame:
            # Detener la construcción
            break
        Rz = dh_transform(comDH[jointFrame, 0], 0, 0, 0)
        H = forward_kinematics_dh(DH, jointFrame - 1)
        H = H * Rz
        z = sp.Matrix(H[0:3, 2])
        r = sp.Matrix(Hcom[0:3, 3] - H[0:3, 3])
        # Jcom[0:3, i - 1] = np.cross(z, r)  # Ajuste de índices para Python (índice base 0)
        Jcom[0:3, i - 1] = z.cross(r)
        Jcom[3:6, i - 1] = z
    return Jcom


def MatrizInercia(m, dh, comDh, a, steps):
    rb = len(m)
    n = len(a)
    D = sp.zeros(n, n)

    for j in range(1, rb + 1):
        hcom = forward_kinematics_dh_com(dh, comDh, j)
        tensor1 = calcular_tensor_inercia(m[j - 1], dh[j - 1, 2])
        hcom = sp.Matrix(hcom)

        I = hcom[0:3, 0:3].transpose() * tensor1 * hcom[0:3, 0:3]
        jcom = matriz_jaco(DH, comDH, a, j)

        jcom = sp.Matrix(jcom)
        ele1 = sp.simplify(m[j - 1] * jcom[0:3, :].transpose() * jcom[0:3, :])
        ele2 = sp.simplify(jcom[3:, :].transpose() * I * jcom[3:, :])
        D = D + sp.simplify(ele1 + ele2)
    return nsimplify(D)


def espacio_N_Normali(matriz_rigidez):
    espacioNulo = matriz_rigidez.nullspace()
    if espacioNulo:
        primer_vector = espacioNulo[0]
        primer_vector_normalizado = primer_vector / primer_vector.norm()
        return primer_vector_normalizado
    else:
        print("La matriz no tiene un espacio nulo no trivial.")


def Matriz_S(Vector, Matriz):
    q1, q2, q3 = sp.symbols('q1 q2 q3')
    q = sp.Matrix([q1, q2, q3])
    u0T_M = Vector.T * Matriz
    simplificado = sp.sympify(u0T_M * q).evalf(5)

    equation = sp.Eq(simplificado[0], 0)

    # Resolver para q3 en términos de q1 y q2
    solution = sp.solve(equation, q3)
    Matriz_S1 = sp.zeros(3, 2)
    Matriz_S1[0, 0] = 1
    Matriz_S1[1, 1] = 1
    valorq1 = solution[0].subs({sp.symbols('q2'): 0}) / sp.symbols('q1')
    valorq2 = solution[0].subs({sp.symbols('q1'): 0}) / sp.symbols('q2')
    Matriz_S1[2, 0] = valorq1.evalf(4)
    Matriz_S1[2, 1] = valorq2.evalf(4)
    return Matriz_S1


def qzb(impulso, Uo, C):
    t, s = sp.symbols('t s')
    qimpulse = impulso.dot(Uo.T)

    denominator = s ** 2 + (Uo.T * C * Uo)[0] * s

    laplace_function = qimpulse / denominator

    # Aplicar la transformada inversa de Laplace
    beta_t = sp.inverse_laplace_transform(laplace_function, s, t)

    qzbr = beta_t.simplify().evalf(4) * Uo
    return qzbr


# NZP motions
def NZP(q0, Mm, u0, qdot0, R, A, S):
    Mp = (sp.transpose(S) * Mm * S).evalf(4)
    Cp = (sp.transpose(S) * R * S).evalf(4)
    Kp = (sp.transpose(S) * A * S).evalf(4)
    q0NRB = (q0 - (u0.T * Mm * q0)[0] * u0).evalf(4)
    qdot0NRB = qdot0 - (u0.T * Mm * qdot0)[0] * u0

    O = sp.zeros(2, 2)
    I = sp.eye(2)
    nueva = O.row_join(I)
    MMinv = Mp.inv()
    mm3 = -MMinv * Kp
    mm4 = -MMinv*Cp
    partInfe= mm3.row_join(mm4)
    A=nueva.col_join(partInfe)
    sp.pprint(A)



# Example usage
if __name__ == "__main__":
    # Example link lengths and center of mass distances
    a = [sp.symbols('a1'), sp.symbols('a2'), sp.symbols('a3')]

    # centro de masa de los eslabones
    Lcom = [sp.symbols('Lcom1'), sp.symbols('Lcom2'), sp.symbols('Lcom3')]

    # parametros DH de centros de masa
    comDH = denavit_hartenberg_com(a, Lcom)

    # Example DH parameters for a 3-joint manipulator
    DH = sp.Matrix([
        [sp.symbols('a1'), 0, 1, 0],
        [sp.symbols('a2'), 0, 1, 0],
        [sp.symbols('a3'), 0, 1, 0],
    ])
    # condiciones iniciales de los eslabones masas  etc
    m = [1, 1, 1]
    valores = {sp.symbols('Lcom1'): 1 / 2, sp.symbols('Lcom2'): 1 / 2, sp.symbols('Lcom3'): 1 / 2}
    angulos = {sp.symbols('a1'): sp.pi / 4, sp.symbols('a2'): sp.pi / 4, sp.symbols('a3'): sp.pi / 4}
    k = {sp.symbols('Kx'): 200, sp.symbols('Ky'): 300}
    c = {sp.symbols('Cx'): 300, sp.symbols('Cy'): 500}

    d1 = MatrizInercia(m, DH, comDH, a, 1)
    matriz_evaluada = d1.subs(valores).subs(angulos)
    matriz = matriz_jaco(DH, comDH, a, 3)
    jac_eval = sp.Matrix(matriz).subs(valores).subs(angulos)
    # descomentar para ver las matrices de inercia
    # print("Inertia Matrix:")
    # sp.pprint(matriz_evaluada.evalf(5))
    Kq = rigidez(jac_eval[0:2, :]).evalf(5).subs(k)
    # print("AMORTIGUAMIENTO Matrix:")
    Cq = Amortiguamiento(jac_eval[0:2, :]).evalf(5).subs(c)
    normalizado = espacio_N_Normali(Kq)
    S = Matriz_S(normalizado, matriz_evaluada)
    M12 = eliminar_filas_redundantes(jac_eval.evalf(3))
    sp.pprint(Amortiguamiento(M12[0:2, :]).evalf(5).subs(c))

    # matrices de rigidez amortiguamieno,e inercia

    sp.pprint(normalizado.T)
    # sp.pprint(qzb(sp.Matrix([[1, 0, 0]]),normalizado,Cq))
    NZP(sp.Matrix([1, 1, 1]), matriz_evaluada, normalizado, sp.Matrix([0, 0, 0]), Kq, Cq, S)
