import sympy as sp
from sympy import nsimplify

from Tesis.rigidez import rigidez, Amortiguamiento
from Tesis.tensor import calcular_tensor_inercia


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
    H = forward_kinematics_dh(DH, frame)

    # Forward kinematics to the requested center of mass
    Rz = dh_transform(comDH[frame, 0], 0, 0, 0)
    Tz = dh_transform(0, comDH[frame, 1], 0, 0)
    Tx = dh_transform(0, 0, comDH[frame, 2], 0)

    Hcom = sp.simplify(H * Rz * Tz * Tx)
    #sp.pprint(Hcom)
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

        # Si la articulación actual está en un marco de referencia que no afecta el COM
        if jointFrame > comFrame:
            # Detener la construcción
            break
        H = forward_kinematics_dh(DH, jointFrame)
        Rz = dh_transform(comDH[jointFrame, 0], 0, 0, 0)

        H = H * Rz
        z = sp.Matrix(H[0:3, 2])

        r = sp.Matrix(sp.simplify(Hcom[0:3, 3] - H[0:3, 3]))

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

        I = hcom[0:3, 0:3].transpose() * tensor1 * hcom[0:3, 0:3]
        jcom = matriz_jaco(DH, comDH, a, j)
        ele1 = sp.simplify(m[j - 1] * jcom[0:3, :].transpose() * jcom[0:3, :])
        ele2 = sp.simplify(jcom[3:, :].transpose() * I * jcom[3:, :])
        D = D + sp.simplify(ele1 + ele2)
    return nsimplify(D)


import sympy as sp


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
            H_prev = forward_kinematics_dh(DH,i+1)
            p_prev = H_prev[0:3, 3]

        J[0:3, i] = z.cross(sp.simplify(prev - p_prev))  # Parte traslacional (solo x e y)
        J[2, i] = z[2]  # Parte rotacional (alrededor de z)

    return J


#nuevas funciones
def espacio_N_Normali(matriz_rigidez):
    espacioNulo = matriz_rigidez.nullspace()
    if espacioNulo:
        primer_vector = espacioNulo[0]
        primer_vector_normalizado = primer_vector / primer_vector.norm()
        return primer_vector_normalizado
    else:
        print("La matriz no tiene un espacio nulo no trivial.")

if __name__ == "__main__":
    a = [sp.symbols('a1'), sp.symbols('a2'), sp.symbols('a3')]

    # centro de masa de los eslabones
    Lcom = [sp.symbols('Lcom1'), sp.symbols('Lcom2'), sp.symbols('Lcom3')]

    # parametros DH de centros de masa
    comDH = denavit_hartenberg_com(a, Lcom)

    # Example DH parameters for a 3-joint manipulator
    DH = sp.Matrix([
        [0, 0, 0, 0],
        [sp.symbols('a1'), 0, sp.symbols('l1'), 0],
        [sp.symbols('a2'), 0, sp.symbols('l2'), 0],
        [sp.symbols('a3'), 0, sp.symbols('l3'), 0],
    ])
    # condiciones iniciales de los eslabones masas  etc
    #valores iniciales
    m = [20, 10, 5]
    valores = {sp.symbols('Lcom1'): 0.5 / 2, sp.symbols('Lcom2'):  0.35 / 2, sp.symbols('Lcom3'): 0.3 / 2}
    angulos = {sp.symbols('a1'):  1.659, sp.symbols('a2'): -1.979, sp.symbols('a3'):1.105}
    k = {sp.symbols('Kx'): 100, sp.symbols('Ky'): 100}
    c = {sp.symbols('Cx'): 15, sp.symbols('Cy'): 15}
    valores1 = {sp.symbols('l1'): 0.5, sp.symbols('l2'): 0.35, sp.symbols('l3'): 0.3}
    jac=matriz_jaco_planar(DH,a).subs(valores1).subs(angulos).evalf(5)
    #sp.pprint(matriz_jaco(DH,comDH, a, 3))
    d1 = MatrizInercia(m, DH, comDH, a, 1).subs(valores).subs(angulos).subs(valores1).evalf(5)
    Kq = rigidez(jac[0:2, :]).evalf(4).subs(k)
    Cq = Amortiguamiento(jac[0:2, :]).evalf(3).subs(c)
    sp.pprint(d1)
    sp.pprint(Kq)
    sp.pprint(Cq)


    espacioNulo = Kq.nullspace()
    sp.pprint(espacioNulo)
