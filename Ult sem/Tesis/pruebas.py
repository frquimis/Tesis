import sympy as sp

Kc = sp.Matrix([
    [61.00, 33.33, 23.33],
    [33.33, 30.67, 13.71],
    [23.33, 13.71, 9.000],
])
valores_propios = Kc.eigenvals()
# Verificar si la matriz es semidefinida positiva
es_semidefinida_positiva = all(val >= 0 for val in valores_propios)

# Encontrar el espacio nulo de la matriz K
espacio_nulo = Kc.nullspace()

print(f"Valores propios: {valores_propios}")
print(f"La matriz es semidefinida positiva: {es_semidefinida_positiva}")
print(f"Espacio nulo: {espacio_nulo}")

# Si hay un espacio nulo no trivial, normalizar un vector en el espacio nulo
if espacio_nulo:
    vector_nulo = espacio_nulo[0]
    vector_nulo_normalizado = vector_nulo / vector_nulo.norm()
    print(f"Vector en el espacio nulo (normalizado): {vector_nulo_normalizado}")
else:
    print("El espacio nulo es trivial (solo contiene el vector cero).")
