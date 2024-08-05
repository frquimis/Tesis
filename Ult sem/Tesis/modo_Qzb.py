from Tesis.Matriz_Inercia import load_variable, save_variable
import sympy as sp
import numpy as np

import matplotlib.pyplot as plt

from Tesis.metodo import Matriz_S, qzb
from Tesis.wf import espacio_N_Normali

impulso = sp.Matrix([1, 0, 0])

matriz_nula_cargada = load_variable('matriz_nula.pkl')
m1 = load_variable('amortiguamiento.pkl')
m2 = load_variable('rigidez.pkl')
m3 = load_variable('inercia.pkl')
#print("matrices")
#output = f"{sp.pretty(m1.evalf(5))}\n{sp.pretty(m2.evalf(5))}\n{sp.pretty(m3.evalf(5))}"
#print(output)
#sp.pprint(matriz_nula_cargada)
#vector normalizado
normalizado = espacio_N_Normali(matriz_nula_cargada, m3)
#matriz de restriccion

S = Matriz_S(normalizado, m3)

qzb1 = qzb(impulso, normalizado, m1)
save_variable(qzb1, 'cuerpo_rigido.pkl')
save_variable(S, 'RESTRICCIONES.pkl')
sp.pprint(qzb1)

qzb1_funcs = [sp.lambdify(sp.symbols('t'), comp, 'numpy') for comp in qzb1]

# Evaluar las funciones en un rango de valores para t
t_values = np.linspace(0, 10, 1000)
qzb1_values = [func(t_values) for func in qzb1_funcs]

# Graficar las componentes de qzb1
plt.figure(figsize=(10, 6))
for i, qzb1_val in enumerate(qzb1_values):
    plt.plot(t_values, qzb1_val, label=f'qzb1_{i+1}')
plt.xlabel('t')
plt.ylabel('Valor')
plt.title('Componentes de qzb1')
plt.legend()
plt.grid(True)
plt.show()
