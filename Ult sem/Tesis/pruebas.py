import sympy as sp

from Tesis.Matriz_Inercia import denavit_hartenberg_com, MatrizInercia, save_variable
from Tesis.wf import matriz_jaco_planar, rigidez, Amortiguamiento, vector_N

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
m = [10, 5, 2]
valores = {sp.symbols('Lcom1'): 0.4 / 2, sp.symbols('Lcom2'): 0.25 / 2, sp.symbols('Lcom3'): 0.2 / 2}
angulos = {sp.symbols('a1'): 0.529, sp.symbols('a2'): 0.6, sp.symbols('a3'): 0.105}
k = {sp.symbols('Kx'): 300, sp.symbols('Ky'): 300}
c = {sp.symbols('Cx'): 20, sp.symbols('Cy'): 15}
longitudes = {sp.symbols('l1'): 0.4, sp.symbols('l2'): 0.25, sp.symbols('l3'): 0.2}

jac = matriz_jaco_planar(DH, a)
jacEvalf = matriz_jaco_planar(DH, a).subs(angulos).subs(longitudes)

Kq = rigidez(jac[0:2, :]).subs(k)

Cq = Amortiguamiento(jac[0:2, :]).subs(c).subs(angulos).subs(longitudes)

d1 = MatrizInercia(m, DH, comDH, a, 1)
'''matrizInercia=sp.Matrix([
        [2.013, 0.4634, 0.2542],
        [0.4634, 0.6215, 0.1339],
        [0.2542, 0.1339, 0.0750]
    ])
'''
matriz_evaluada = d1.subs(valores).subs(angulos).subs(longitudes)
sp.pprint(matriz_evaluada)

m1 = vector_N(Kq, angulos, longitudes)

save_variable(m1, 'matriz_nula.pkl')

save_variable(jacEvalf, 'Jacobiano.pkl')

Kqeva = Kq.subs(angulos).subs(longitudes)

save_variable(Cq, 'amortiguamiento.pkl')

save_variable(Kqeva, 'rigidez.pkl')
save_variable(matriz_evaluada, 'inercia.pkl')
