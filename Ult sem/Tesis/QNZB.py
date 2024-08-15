import numpy as np
import sympy as sp
from matplotlib import pyplot as plt

from Tesis.Matriz_Inercia import load_variable, save_variable
from Tesis.metodo import NZP, xt, matrize

Mamort = load_variable('amortiguamiento.pkl')
Mrestric = load_variable('RESTRICCIONES.pkl')
Mrigi = load_variable('rigidez.pkl')
Miner = load_variable('inercia.pkl')
qo = load_variable('xcero.pkl')
derechos, izquierdos, valores, B = NZP(Miner, Mrigi, Mamort, Mrestric)
x1t = xt(derechos, izquierdos, valores, B, Mrestric, qo)
sp.pprint(x1t)
'''QNBZ = sp.simplify(Mrestric * x1t[0:2, :])
QNBZ = sp.simplify(QNBZ.evalf(6))
save_variable(QNBZ, 'cuerpoNoRigido.pkl')

qNzb1_funcs = [sp.lambdify(sp.symbols('t'), comp, 'numpy') for comp in QNBZ]

# Evaluar las funciones en un rango de valores para t
t_values = np.linspace(0, 10, 400)
qzb1_values = [func(t_values) for func in qNzb1_funcs]

# Graficar las componentes de qzb1
plt.figure(figsize=(10, 6))
for i, qzb1_val in enumerate(qzb1_values):
    plt.plot(t_values, qzb1_val, label=f'QNBZ_{i + 1}')
plt.xlabel('t')
plt.ylabel('Valor')
plt.title('Componentes de QNBZ')
plt.legend()
plt.grid(True)
plt.show()'''
