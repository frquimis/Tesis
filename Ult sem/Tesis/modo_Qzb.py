from Tesis.Matriz_Inercia import load_variable, save_variable
import sympy as sp

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
sp.pprint(matriz_nula_cargada)
#vector normalizado
normalizado = espacio_N_Normali(matriz_nula_cargada, m3)
#matriz de restriccion

S = Matriz_S(normalizado, m3)

qzb1 = qzb(impulso, normalizado, m1)
save_variable(qzb1, 'cuerpo_rigido.pkl')
save_variable(S, 'RESTRICCIONES.pkl')
