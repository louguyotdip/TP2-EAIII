import numpy as np

# --- FUNCIONES DE DISEÑO DE MICROTIRA ---

def calcularMicrotiraLambda4(Z, er, h, t, fp):
    # Cálculo de A y B
    A = (Z / 60) * np.sqrt((er + 1) / 2) + ((er - 1) / (er + 1)) * (0.23 + 0.11 / er)
    B = (377 * np.pi) / (2 * Z * np.sqrt(er))

    # Cálculo de W/H (Usando la fórmula de tu función para B > 0.66)
    W_H = (2 / np.pi) * ((B - 1) - np.log(2 * B - 1) + ((er - 1) / (2 * er)) * (np.log(B - 1) + 0.39 - 0.61 / er))

    # Ancho W
    W = h * W_H

    # Corrección de ancho efectivo
    We = W + (t / np.pi) * (1 + np.log(2 * h / t))

    # Cálculo de epsilon efectivo
    er_eff = (er + 1) / 2 + (er - 1) / 2 * (1 / np.sqrt(1 + 12 * h / We))

    # Impedancia característica corregida
    Zo_e = (120 * np.pi) / (np.sqrt(er_eff) * (W / h + 1.393 + 0.667 * np.log(W / h + 1.444)))

    # Longitud de onda en el sustrato
    lambda_val = 300 / (fp * 1e-6)  # fp en Hz, convertimos para obtener lambda en m
    lambda_p = lambda_val / np.sqrt(er_eff)

    # Largo de la microtira λ/4
    l_microtira = lambda_p / 4

    # Ángulo eléctrico
    beta = (2 * np.pi) / lambda_p
    angulo_electrico = beta * l_microtira * (180 / np.pi)

    print('\nResultados del diseño de microtira \u03bb/4:')
    print('------------------------------------------')
    print(f'A(Zin) = {A:.4f}')
    print(f'B(Zin) = {B:.4f}')
    print(f'Impedancia adaptador (Zo \u03bb/4) = {Z:.2f} Ohms')
    print(f'Ancho W = {W*1e3:.4f} mm')
    print(f'Ancho efectivo We = {We*1e3:.4f} mm')
    print(f'Epsilon efectivo = {er_eff:.4f}')
    print(f'Impedancia corregida Zo_e = {Zo_e:.2f} Ohms')
    print(f'Longitud de microtira = {l_microtira*1e3:.4f} mm')
    print(f'\u00c1ngulo el\u00e9ctrico = {angulo_electrico:.2f} grados')


def calcularMicrotiraCap(Z, er, h, t, fp, C):
    # Cálculo de A y B
    A = (Z / 60) * np.sqrt((er + 1) / 2) + ((er - 1) / (er + 1)) * (0.23 + 0.11 / er)
    B = (377 * np.pi) / (2 * Z * np.sqrt(er))

    # Cálculo de W/H
    W_H = (2 / np.pi) * ((B - 1) - np.log(2 * B - 1) + ((er - 1) / (2 * er)) * (np.log(B - 1) + 0.39 - 0.61 / er))

    # Ancho W
    W = h * W_H

    # Corrección de ancho efectivo
    We = W + (t / np.pi) * (1 + np.log(2 * h / t))

    # Cálculo de epsilon efectivo
    er_eff = (er + 1) / 2 + (er - 1) / 2 * (1 / np.sqrt(1 + 12 * h / We))

    # Impedancia característica corregida
    Zo_e = (120 * np.pi) / (np.sqrt(er_eff) * (W / h + 1.393 + 0.667 * np.log(W / h + 1.444)))

    # Longitud de onda en el sustrato
    lambda_val = 300 / (fp * 1e-6)
    lambda_p = lambda_val / np.sqrt(er_eff)

    # Calculo de beta y reactancia
    beta = (2 * np.pi) / lambda_p
    Xc = 1 / (2 * np.pi * fp * C)

    # Cálculo del largo para el capacitor (stub)
    d = np.arctan(Zo_e / Xc) / beta

    # Ángulo eléctrico
    angulo_electrico = beta * d * (180 / np.pi)

    print('\nResultados del diseño de microtira capacitiva:')
    print('----------------------------------------------')
    print(f'A(Zin) = {A:.4f}')
    print(f'B(Zin) = {B:.4f}')
    print(f'Ancho W = {W*1e3:.4f} mm')
    print(f'Ancho efectivo We = {We*1e3:.4f} mm')
    print(f'Epsilon efectivo = {er_eff:.4f}')
    print(f'Impedancia corregida Zo_e = {Zo_e:.2f} Ohms')
    print(f'Longitud de microtira (capacitor) = {d*1e3:.4f} mm')
    print(f'\u00c1ngulo el\u00e9ctrico = {angulo_electrico:.2f} grados')


# --- SCRIPT PRINCIPAL (MAIN) ---

# Datos

#s_11 = 0.8734 * np.exp(1j * 142.5 * np.pi / 180)
#s_21 = 1.097 * np.exp(1j * 45.4 * np.pi / 180)
#s_12 = 0.0961 * np.exp(1j * 47 * np.pi / 180)
#s_22 = 0.701 * np.exp(1j * 145.3 * np.pi / 180)

s_11 =-0.776386 + 0.234405j
s_21 = 1.627867 + 3.969150j
s_12 = 0.033230 + 0.021170j
s_22 = -0.256968 + 0.004037j

Zo = 50
fp = 2.2e9
er = 4
h = 1.546e-3
t = 29.5e-6

# Cálculos de Estabilidad y Gamma
delta = s_11 * s_22 - s_12 * s_21
k = (1 + np.abs(delta)**2 - np.abs(s_22)**2 - np.abs(s_11)**2) / (2 * np.abs(s_12 * s_21))

C1 = s_11 - delta * np.conj(s_22)
B1 = 1 + np.abs(s_11)**2 - np.abs(s_22)**2 - np.abs(delta)**2
C2 = s_22 - delta * np.conj(s_11)
B2 = 1 + np.abs(s_22)**2 - np.abs(s_11)**2 - np.abs(delta)**2

gamma_in_mag = (B1 - np.sqrt(B1**2 - 4 * np.abs(C1)**2)) / (2 * np.abs(C1))
gamma_in = gamma_in_mag * np.exp(1j * np.angle(C1))

gamma_out_mag = (B2 - np.sqrt(B2**2 - 4 * np.abs(C2)**2)) / (2 * np.abs(C2))
gamma_out = gamma_out_mag * np.exp(1j * np.angle(C2))

# Impedancias
Zin = Zo * ((1 + gamma_in) / (1 - gamma_in))
Zout = Zo * ((1 + gamma_out) / (1 - gamma_out))

# Ganancia MAG
MAG_dB = 10 * np.log10((np.abs(s_21) / np.abs(s_12)) * (k - np.sqrt(k**2 - 1)))

# Conversión a paralelo
Rin_s, Xin_s = np.real(Zin), np.imag(Zin)
Rin_p = Rin_s * (1 + (Xin_s / Rin_s)**2)
Xin_p = Xin_s * (1 + (Rin_s / Xin_s)**2)

# Componentes entrada
Cin = 1 / (2 * np.pi * fp * np.abs(Xin_p))
Zo_in = np.sqrt(Rin_p * Zo)

# Componentes salida
Rout_s, Xout_s = np.real(Zout), np.imag(Zout)
Zo_out = np.sqrt(Zo * Rout_s)
Xout_p = np.abs(Zo_out**2 / Xout_s)
Cout = 1 / (2 * np.pi * fp * Xout_p)

# --- IMPRESIÓN DE RESULTADOS GENERALES ---
print(f"Zin = {Zin:.4f}")
print(f"Zout = {Zout:.4f}")
print(f"MAG = {MAG_dB:.2f} dB")
print(f"Cin = {Cin:.2e} F")
print(f"Cout = {Cout:.2e} F")

print("\n" + "="*42)
print("ADAPTADOR DE ENTRADA")
calcularMicrotiraLambda4(Zo_in, er, h, t, fp)

print("\n" + "="*42)
print("ADAPTADOR DE SALIDA")
calcularMicrotiraLambda4(Zo_out, er, h, t, fp)

print("\n" + "="*46)
print("CAPACITOR DE ENTRADA")
calcularMicrotiraCap(Zo, er, h, t, fp, Cin)

print("\n" + "="*46)
print("CAPACITOR DE SALIDAkk")
calcularMicrotiraCap(60, er, h, t, fp, Cout)
