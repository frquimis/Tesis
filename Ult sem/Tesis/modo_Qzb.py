from Tesis.Cinematica_Inversa import rrrIK
from Tesis.Matriz_Inercia import load_variable, save_variable
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from Tesis.metodo import Matriz_S, qzb
from Tesis.wf import espacio_N_Normali, funcion

#elementos cargados
matriz_nula_cargada = load_variable('matriz_nula.pkl')
m1 = load_variable('amortiguamiento.pkl')
jaco = load_variable('Jacobiano.pkl')
m2 = load_variable('rigidez.pkl')
m3 = load_variable('inercia.pkl')

#creacion de distintos escenarios
qx = 0.20
qy = 0.00
vx = 0
vy = 0
q0 = rrrIK(0.5, 0.35, 0.3, [(0.685 - qx), (0.84 + qy), 1.23])
condicio = sp.Matrix(q0)
qdot0 = funcion(vx, vy, jaco[0:2, :])
normalizado = espacio_N_Normali(matriz_nula_cargada, m3)

qzb1, qrb, qrbdot = qzb(condicio, qdot0, normalizado, m3)
sp.pprint(qzb1)
QNRB = condicio - qrb
QNRB_dot = qdot0 - qrbdot

x0 = QNRB[0:2, :].col_join(QNRB_dot[0:2, :])

S = Matriz_S(normalizado, m3)
save_variable(qzb1, 'cuerpo_rigido.pkl')
save_variable(x0, 'xcero.pkl')
save_variable(S, 'RESTRICCIONES.pkl')

'''
qzb1_funcs = [sp.lambdify(sp.symbols('t'), comp, 'numpy') for comp in qzb1]

# Evaluar las funciones en un rango de valores para t
t_values = np.linspace(0, 10, 400)
qzb1_values = [func(t_values) for func in qzb1_funcs]

# Graficar las componentes de qzb1
plt.figure(figsize=(10, 6))
for i, qzb1_val in enumerate(qzb1_values):
    plt.plot(t_values, qzb1_val, label=f'QZB_{i + 1}')
plt.xlabel('t')
plt.title('Componentes de qzb1')
plt.legend()
plt.grid(True)
plt.show()'''
'''
q = {sp.symbols('Q'): 7, sp.symbols('Fy'): 5}#

Qo = sp.Matrix([
    [sp.symbols('Fx')],
    [sp.symbols('Fy')]
])

impulso = funcion_Impulso(Q, jaco[0:2, :]).subs(F)
sp.pprint(impulso)

#print("matrices")
#output = f"{sp.pretty(m1.evalf(5))}\n{sp.pretty(m2.evalf(5))}\n{sp.pretty(m3.evalf(5))}"
#sp.pprint(matriz_nula_cargada)
#vector normalizado
normalizado = espacio_N_Normali(matriz_nula_cargada, m3)
#matriz de restriccion

S = Matriz_S(normalizado, m3)
sp.pprint(S)

qzb1 = qzb(impulso, normalizado, m1)
save_variable(qzb1, 'cuerpo_rigido.pkl')
save_variable(S, 'RESTRICCIONES.pkl')
sp.pprint(qzb1)
save_variable(impulso, 'funcionImpulso.pkl')

qzb1_funcs = [sp.lambdify(sp.symbols('t'), comp, 'numpy') for comp in qzb1]

# Evaluar las funciones en un rango de valores para t
t_values = np.linspace(0, 10, 400)
qzb1_values = [func(t_values) for func in qzb1_funcs]

# Graficar las componentes de qzb1
plt.figure(figsize=(10, 6))
for i, qzb1_val in enumerate(qzb1_values):
    plt.plot(t_values, qzb1_val, label=f'QZB_{i + 1}')
plt.xlabel('t')
plt.title('Componentes de qzb1')
plt.legend()
plt.grid(True)
plt.show()'''
