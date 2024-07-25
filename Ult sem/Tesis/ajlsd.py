import sympy as sp


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
    H = forward_kinematics_dh(DH, frame )


    # Forward kinematics to the requested center of mass
    Rz = dh_transform(comDH[frame, 0], 0, 0, 0)
    Tz = dh_transform(0, comDH[frame, 1], 0, 0)
    Tx = dh_transform(0, 0, comDH[frame, 2], 0)

    Hcom = sp.simplify(H * Rz * Tz * Tx)
    sp.pprint(Hcom)
    return Hcom


def matriz_jaco(DH, comDH, a, rb):
    Lcom1 = sp.symbols(f'Lcom{rb}')
    comFrame = None
    for i in range(comDH.shape[0]):
        if comDH[i, 2] == Lcom1:
            comFrame = i
            break
    # Calcular la cinemática directa al centro de masa solicitado
    Hcom = forward_kinematics_dh_com(DH, comDH, 2)
    sp.pprint(Hcom)
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
        H = forward_kinematics_dh(DH, jointFrame - 1)
        Rz = dh_transform(comDH[jointFrame, 0], 0, 0, 0)

        H = H * Rz
        sp.pprint(H)
        z = sp.Matrix(H[0:3, 2])

        r = sp.Matrix(Hcom[0:3, 3] - H[0:3, 3])
        sp.pprint(r)
        sp.pprint(z.cross(r))
        # Jcom[0:3, i - 1] = np.cross(z, r)  # Ajuste de índices para Python (índice base 0)
        Jcom[0:3, i - 1] = z.cross(r)
        Jcom[3:6, i - 1] = z
    return Jcom


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
    m = [20, 10, 5]

forward_kinematics_dh_com(DH, comDH, 2)