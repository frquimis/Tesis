import numpy as np
import sympy as sp
from matplotlib import pyplot as plt

from Tesis.Matriz_Inercia import load_variable
from Tesis.metodo import NZP, xt, matrize

Mamort = load_variable('amortiguamiento.pkl')
Mrestric = load_variable('RESTRICCIONES.pkl')
Mrigi = load_variable('rigidez.pkl')
Miner = load_variable('inercia.pkl')
qo = sp.Matrix([0, 0, 1])
derechos, izquierdos, valores, B = NZP(Miner, Mrigi, Mamort, Mrestric)
x1t = xt(derechos, izquierdos, valores, B, Mrestric, qo)
QNBZ = sp.simplify(Mrestric * x1t[0:2, :])
QNBZ = sp.simplify(QNBZ.evalf(6))

qNzb1_funcs = [sp.lambdify(sp.symbols('t'), comp, 'numpy') for comp in QNBZ]

# Evaluar las funciones en un rango de valores para t
t_values = np.linspace(0, 15, 400)
qzb1_values = [func(t_values) for func in qNzb1_funcs]

# Graficar las componentes de qzb1
plt.figure(figsize=(10, 6))
for i, qzb1_val in enumerate(qzb1_values):
    plt.plot(t_values, qzb1_val, label=f'qzb1_{i + 1}')
plt.xlabel('t')
plt.ylabel('Valor')
plt.title('Componentes de qzb1')
plt.legend()
plt.grid(True)
plt.show()
