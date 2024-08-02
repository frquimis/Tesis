import sympy as sp

from Tesis.Matriz_Inercia import load_variable
from Tesis.metodo import NZP, xt, matrize

Mamort = load_variable('amortiguamiento.pkl')
Mrestric = load_variable('RESTRICCIONES.pkl')
Mrigi = load_variable('rigidez.pkl')
Miner = load_variable('inercia.pkl')
qo = sp.Matrix([0, 0, 1])
derechos, izquierdos, valores, B = NZP(Miner, Mrigi, Mamort, Mrestric)
x1t = xt(derechos, izquierdos, valores, B, Mrestric, qo)

