import numpy as np
from matplotlib import pyplot as plt
import sympy as sp

from Tesis.Matriz_Inercia import load_variable, save_variable

qzb = load_variable('cuerpo_rigido.pkl')
qznb = load_variable('cuerpoNoRigido.pkl')

qz = qzb + qznb

save_variable(qz, 'qt.pkl')
#sp.pprint(qz)

qNzb1_funcs = [sp.lambdify(sp.symbols('t'), comp, 'numpy') for comp in qz]

# Evaluar las funciones en un rango de valores para t
t_values = np.linspace(0, 15, 400)
qzb1_values = [func(t_values) for func in qNzb1_funcs]

# Graficar las componentes de qzb1

plt.figure(figsize=(10, 6))
for i, qzb1_val in enumerate(qzb1_values):
    plt.plot(t_values, qzb1_val, label=f'q(t){i + 1}')
    np.savetxt("Qimp{}.csv".format(i + 1), np.column_stack((qzb1_val.real, t_values)), delimiter=',')
    print(qzb1_val)
plt.xlabel('t(s)')
plt.ylabel('Amplitud')
plt.title('freeVibration')
plt.legend()
plt.grid(True)
plt.show()
