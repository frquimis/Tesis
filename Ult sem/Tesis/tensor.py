import sympy as sp

import sympy as sp

# Función para calcular el tensor de inercia simbólicamente
def calcular_tensor_inercia(m, l):
    Ixx = Iyy = (1/12) * m * l**2
    Izz = (1/6) * m * l**2
    tensor = sp.Matrix([
        [Ixx, 0, 0],
        [0, Iyy, 0],
        [0, 0, Izz]
    ])
    return tensor.evalf(2)


