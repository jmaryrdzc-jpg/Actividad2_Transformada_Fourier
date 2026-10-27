"""Actividad formativa 2: análisis de señales con FFT. Ejecutar: python simulacion_fourier.py"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent / 'graficas'
OUT.mkdir(exist_ok=True)
FS = 1000  # muestras por segundo
N = 2000  # ventana de 2 s; resolución espectral 0.5 Hz
t = np.arange(N) / FS - 1.0

def espectro(x):
    """FFT bilateral centrada y normalizada por N (amplitud por bin)."""
    return np.fft.fftshift(np.fft.fft(x)) / len(x)

f = np.fft.fftshift(np.fft.fftfreq(N, d=1/FS))

def graficar(nombre, x, rango=60):
    X = espectro(x)
    fig, ax = plt.subplots(3, 1, figsize=(10, 8), constrained_layout=True)
    ax[0].plot(t, x, lw=1.5); ax[0].set(xlabel='Tiempo (s)', ylabel='Amplitud', title=f'{nombre}: dominio del tiempo'); ax[0].grid(alpha=.3)
    sel = np.abs(f) <= rango
    ax[1].plot(f[sel], np.abs(X[sel]), lw=1.3); ax[1].set(xlabel='Frecuencia (Hz)', ylabel='|FFT|/N', title='Magnitud bilateral'); ax[1].grid(alpha=.3)
    phase = np.angle(X[sel]); phase[np.abs(X[sel]) < 1e-5] = np.nan
    ax[2].scatter(f[sel], phase, s=7); ax[2].set(xlabel='Frecuencia (Hz)', ylabel='Fase (rad)', title='Fase (solo bins con magnitud significativa)', ylim=(-3.4, 3.4)); ax[2].grid(alpha=.3)
    fig.savefig(OUT / (nombre.lower().replace(' ', '_') + '.png'), dpi=160); plt.close(fig)
    return X

rect = (np.abs(t) < .1).astype(float)
escalon = (t >= 0).astype(float)
seno = np.sin(2*np.pi*10*t)
for name, x in [('Pulso rectangular', rect), ('Escalon unitario', escalon), ('Senoidal 10 Hz', seno)]: graficar(name, x)

# Linealidad: FFT(a*x + b*y) = a*FFT(x) + b*FFT(y)
a, b = 1.7, -0.4
err_lin = np.max(np.abs(espectro(a*rect+b*seno) - (a*espectro(rect)+b*espectro(seno))))
# Desplazamiento circular discreto: roll por k muestras -> exp(-j2*pi*m*k/N)
k = 120
m = np.fft.fftshift(np.fft.fftfreq(N)*N)
err_shift = np.max(np.abs(espectro(np.roll(rect,k))-espectro(rect)*np.exp(-2j*np.pi*m*k/N)))
# Escalamiento temporal x(a*t): seno a 10 Hz -> 20 Hz si a=2
seno_doble = np.sin(2*np.pi*20*t)
graficar('Senoidal 20 Hz', seno_doble)

fig, ax = plt.subplots(2, 1, figsize=(10, 6), constrained_layout=True)
ax[0].plot(t, seno, label='10 Hz'); ax[0].plot(t, seno_doble, label='20 Hz', alpha=.75); ax[0].set_xlim(-.2,.2); ax[0].set(xlabel='Tiempo (s)',ylabel='Amplitud',title='Escalamiento temporal x(2t)'); ax[0].legend(); ax[0].grid(alpha=.3)
for x, label in [(seno,'10 Hz'),(seno_doble,'20 Hz')]:
    ax[1].plot(f, np.abs(espectro(x)),label=label)
ax[1].set_xlim(-30,30); ax[1].set(xlabel='Frecuencia (Hz)',ylabel='|FFT|/N',title='Picos espectrales ±10 Hz y ±20 Hz'); ax[1].legend(); ax[1].grid(alpha=.3)
fig.savefig(OUT/'comparacion_escalamiento.png',dpi=160); plt.close(fig)

fig, ax = plt.subplots(2,1,figsize=(10,6),constrained_layout=True)
ax[0].plot(t,rect,label='Original'); ax[0].plot(t,np.roll(rect,k),label='Desplazado 0.12 s',alpha=.8); ax[0].legend(); ax[0].set(xlabel='Tiempo (s)',ylabel='Amplitud',title='Desplazamiento circular de 120 muestras'); ax[0].grid(alpha=.3)
sel = np.abs(f)<=40
ax[1].plot(f[sel],np.abs(espectro(rect))[sel],label='Original'); ax[1].plot(f[sel],np.abs(espectro(np.roll(rect,k)))[sel],ls='--',label='Desplazado'); ax[1].legend(); ax[1].set(xlabel='Frecuencia (Hz)',ylabel='|FFT|/N',title='Magnitud invariante al desplazamiento circular'); ax[1].grid(alpha=.3)
fig.savefig(OUT/'comparacion_desplazamiento.png',dpi=160); plt.close(fig)
print(f'Error máximo linealidad: {err_lin:.3e}')
print(f'Error máximo desplazamiento circular: {err_shift:.3e}')
print('Frecuencia dominante seno original: 10 Hz; escalado: 20 Hz')
