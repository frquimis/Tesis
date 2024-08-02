import sympy as sp

def eliminar_filas_redundantes(A):
    filas_no_redundantes = []

    rango_original = A.rank()

    for i in range(A.rows):
        # Crear una nueva matriz excluyendo la fila i
        A_sin_fila_i = sp.Matrix([A.row(j) for j in range(A.rows) if j != i])
        # Calcular el rango de la nueva matriz
        rango_sin_fila_i = A_sin_fila_i.rank()

        # Si el rango de la nueva matriz es menor que el rango original,
        # la fila i es redundante y se puede eliminar.
        if rango_sin_fila_i < rango_original:
            filas_no_redundantes.append(A.row(i))

        # Convertir las filas no redundantes en una nueva matriz de SymPy
    A_simplificad1 = sp.Matrix(filas_no_redundantes)

    # Eliminar filas que consisten únicamente de ceros
    filas_no_nulas = [fila for fila in A_simplificad1.tolist() if any(fila)]

    # Convertir de nuevo a una mat
    # riz de SymPy
    M = sp.Matrix(filas_no_nulas)

    return M
