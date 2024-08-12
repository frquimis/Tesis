import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

from Tesis.Matriz_Inercia import load_variable

qt = load_variable('qt.pkl')
sampling_rate = 25  # Frecuencia de muestreo (Hz)
time_array = np.linspace(0, 1, sampling_rate)  # Vector de tiempo
vector_funcs = [sp.lambdify(sp.symbols('t'), qt[i], 'numpy') for i in range(3)]
vector_numeric = np.array([func(time_array) for func in vector_funcs])


# Definir la función para la FFT de un vector 3x1 en función del tiempo
def fft_3x1(vector_3x1, sampling_rate):
    fft_result = []
    freqs = np.fft.fftfreq(len(vector_3x1[0]), 1 / sampling_rate)

    for i in range(3):
        fft_vector = np.fft.fft(vector_3x1[i])
        fft_result.append(fft_vector)

    return np.array(fft_result), freqs


fft_result1, freqs1 = fft_3x1(vector_numeric, sampling_rate)
plt.figure(figsize=(12, 8))

for i in range(3):
    plt.subplot(3, 1, i + 1)
    plt.plot(freqs1[:len(freqs1) // 2], np.abs(fft_result1[i])[:len(freqs1) // 2])
    plt.title(f'FFT del vector componente {i + 1}')
    plt.xlabel('Frecuencia (Hz)')
    plt.ylabel('Amplitud')

plt.tight_layout()
plt.show()
