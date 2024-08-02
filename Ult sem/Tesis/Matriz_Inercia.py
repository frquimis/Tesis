import sympy as sp
from sympy import nsimplify
import pickle


# Funciones importadas (deben estar definidas en rigidez.py y metodo.py)
# from Tesis.rigidez import rigidez, Amortiguamiento
# from Tesis.tensor import calcular_tensor_inercia

def denavit_hartenberg_com(a, Lcom):
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
        jcom = matriz_jaco(dh, comDh, a, j)
        ele1 = sp.simplify(m[j - 1] * jcom[0:3, :].transpose() * jcom[0:3, :])
        ele2 = sp.simplify(jcom[3:, :].transpose() * I * jcom[3:, :])
        D = D + sp.simplify(ele1 + ele2)
    return nsimplify(D)


# Función para calcular el tensor de inercia simbólicamente
def calcular_tensor_inercia(m, l):
    Ixx = Iyy = (1 / 12) * m * l ** 2
    Izz = (1 / 6) * m * l ** 2
    tensor = sp.Matrix([
        [Ixx, 0, 0],
        [0, Iyy, 0],
        [0, 0, Izz]
    ])
    return tensor.evalf(2)


def save_variable(variable, filename):
    with open(filename, 'wb') as file:
        pickle.dump(variable, file)


def load_variable(filename):
    with open(filename, 'rb') as file:
        return pickle.load(file)
