import numpy as np
from matplotlib import pyplot as plt
import sympy as sp

from Tesis.Matriz_Inercia import load_variable

qzb = load_variable('cuerpo_rigido.pkl')
qznb = load_variable('cuerpoNoRigido.pkl')

qz=qzb+qznb
#sp.pprint(qz)
qNzb1_funcs = [sp.lambdify(sp.symbols('t'), comp, 'numpy') for comp in qz]

# Evaluar las funciones en un rango de valores para t
t_values = np.linspace(0, 0.5, 400)
qzb1_values = [func(t_values) for func in qNzb1_funcs]

# Graficar las componentes de qzb1
plt.figure(figsize=(10, 6))
for i, qzb1_val in enumerate(qzb1_values):
    plt.plot(t_values, qzb1_val, label=f'Qimp_{i + 1}')
plt.xlabel('t')
plt.ylabel('Valor')
plt.title('Qimp')
plt.legend()
plt.grid(True)
plt.show()
