# Simulación y análisis de señales con la transformada de Fourier

**Autora:** Juana María Rodríguez Camarillo · **Materia:** Señales y Sistemas · **UCNL**

Se estudian pulso rectangular, escalón unitario y señales senoidales en los dominios temporal y frecuencial mediante la FFT de NumPy. Se muestran magnitud y fase, y se comprueban linealidad, desplazamiento circular y escalamiento temporal/frecuencial.

## Requisitos
Python 3.10+; `numpy`, `matplotlib` (SciPy es opcional para ampliar el trabajo).

```bash
python -m pip install -r requirements.txt
python simulacion_fourier.py
```

Los archivos PNG se guardan en `graficas/`. La ventana es de 2 s, con Fs=1000 Hz y resolución Δf=0.5 Hz. La FFT está centrada con `fftshift` y normalizada por N. La fase en bins de magnitud casi nula se oculta por no ser interpretable. El escalón y el pulso se observan sobre una ventana finita: su FFT no es exactamente la transformada continua ideal.

## Publicación en GitHub
Crea un repositorio vacío en GitHub, abre una terminal en esta carpeta y ejecuta (reemplaza USUARIO y REPOSITORIO):

```bash
git init
git add .
git commit -m "Simulacion y analisis de Fourier"
git branch -M main
git remote add origin https://github.com/USUARIO/REPOSITORIO.git
git push -u origin main
```

No se ha publicado automáticamente ni se afirma que el repositorio exista.
