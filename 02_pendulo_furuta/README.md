# Péndulo de Furuta (QUBE-Servo 2) — Koopman Estimation

> **Planta:** Péndulo de Furuta (Quanser QUBE-Servo 2) — sistema MIMO no lineal subactuado.  
> **Tipo:** 4 estados / 1 entrada — $[\alpha, \theta, \dot{\alpha}, \dot{\theta}]$, voltaje $V$.  
> **Objetivo:** Identificación del operador de Koopman mediante EDMD y SINDy, con diseño de observador de estado para MPC en lazo cerrado.

---

## Modelo de la planta

El péndulo de Furuta consiste en un brazo horizontal (ángulo $\theta$) que impulsa un péndulo libre (ángulo $\alpha$). Las ecuaciones de movimiento son no lineales:

$$\frac{d}{dt}\begin{bmatrix}\alpha \\ \theta \\ \dot{\alpha} \\ \dot{\theta}\end{bmatrix} = f\!\left(\alpha, \theta, \dot{\alpha}, \dot{\theta}, V\right)$$

con inerciales del brazo ($J_{\text{arm}}$), inercia del péndulo ($J_p$), masas, longitudes y fricciones como parámetros (ver `parameters.m`).

**Señales:**

| Señal | Descripción | Unidades |
|-------|-------------|---------|
| $\alpha$ | Ángulo del péndulo (posición vertical = 0) | rad |
| $\theta$ | Ángulo del brazo | rad |
| $\dot{\alpha}$, $\dot{\theta}$ | Velocidades angulares | rad/s |
| $V$ | Tensión aplicada al motor | V |

**Datos de entrenamiento:** 10 000 muestras a $\Delta t = 0.001$ s, excitación PRBS en voltaje.

---

## Archivos incluidos

```
src/
├── parameters.m          — Parámetros físicos del sistema (masas, inercias, longitudes)
├── koopman_edmd.m        — Identificación Koopman con EDMD (algoritmo principal)
├── koopman_sindy.m       — Identificación Koopman con SINDy (alternativa dispersa)
└── observer_design.m     — Diseño del observador de estado (filtro de Kalman / Luenberger)
```

**Datos** (descargar de Google Drive):
- `datos/pendulum_time_series.mat` — Serie temporal de entrenamiento (PRBS)
- `datos/pendulum_time_series_test.mat` — Serie temporal de validación

---

## Descripción de los algoritmos

### `koopman_edmd.m` — EDMD (Extended Dynamic Mode Decomposition)

Estima el operador de Koopman finito $\mathbf{A} \in \mathbb{R}^{p \times p}$ y la matriz de salida $\mathbf{C}$ resolviendo:

$$\min_{\mathbf{A},\mathbf{C}} \sum_{i=0}^{N-1} \bigl\|\boldsymbol{\Phi}_{i+1} - \mathbf{A}\boldsymbol{\Phi}_i\bigr\|^2 + \bigl\|y_i - \mathbf{C}\boldsymbol{\Phi}_i\bigr\|^2$$

donde el diccionario de funciones de elevación $\boldsymbol{\Phi}_i = \Psi(x_i)$ incluye términos polinomiales, sinusoidales y productos cruzados de los estados.

**Flujo del script:**
1. Carga y pre-procesa las series temporales (conversión deg → rad, selección de ventana)
2. Construye los snapshots $X = [\Phi_0 \cdots \Phi_{N-1}]$ y $X' = [\Phi_1 \cdots \Phi_N]$
3. Resuelve EDMD vía pseudoinversa: $\hat{A} = X' X^+$
4. Calcula el espectro del operador (eigenvalores de $\hat{A}$)
5. Valida la predicción $k$-pasos sobre el conjunto de prueba
6. Grafica errores de predicción y modos de Koopman

**Parámetros clave:**

| Variable | Descripción |
|----------|-------------|
| `Npoints` | Longitud del conjunto de entrenamiento (default: 10 001) |
| `N_test` | Longitud del conjunto de prueba (default: 10 000) |
| Diccionario | Ajustable en la sección `Lifting function` del script |

### `koopman_sindy.m` — SINDy (Sparse Identification of Nonlinear Dynamics)

Alternativa dispersa que identifica cuáles términos de una biblioteca de candidatos son activos en la dinámica:

$$\dot{x} \approx \boldsymbol{\Theta}(x)\,\boldsymbol{\Xi}$$

usando la minimización LASSO/sequentialThresholding para obtener $\boldsymbol{\Xi}$ con pocos términos no nulos.

**Ventaja sobre EDMD:** El modelo resultante es interpretable (muestra los términos físicamente activos en la dinámica).

### `observer_design.m` — Diseño del observador

Diseña un observador de estado (Luenberger o Kalman) sobre el modelo linealizado del péndulo para:
- Estimar $\dot{\alpha}$ y $\dot{\theta}$ cuando solo se miden $\alpha$ y $\theta$
- Proveer el vector de estado completo al controlador MPC en lazo cerrado

---

## Cómo ejecutar

### Requisitos

- MATLAB R2020b o superior
- Toolboxes: **Control System**, **Signal Processing**
- Hardware Quanser QUBE-Servo 2 (opcional — solo para experimentos físicos; los scripts de identificación funcionan con los datos `.mat` incluidos)

### Pasos

```matlab
% 1. Navegar a la carpeta:
cd('ruta/a/codigo/02_pendulo_furuta/src')

% 2. Descargar datos y colocarlos en 02_pendulo_furuta/datos/
%    (ver enlace Google Drive en el README principal)

% 3. Ejecutar los parámetros primero:
parameters

% 4. Ejecutar la identificación con EDMD:
koopman_edmd
% Genera figuras de entrenamiento, espectro de eigenvalores y validación.

% 5. (Opcional) Comparar con SINDy:
koopman_sindy

% 6. Diseñar el observador sobre el modelo identificado:
observer_design
```

### Salidas esperadas

- **Figura 1:** Dataset de entrenamiento — $\alpha$, $\theta$ vs. tiempo y voltaje aplicado
- **Figura 2:** Error de predicción $k$-pasos del modelo Koopman
- **Figura 3:** Eigenvalores del operador estimado en el plano complejo
- **Figura 4 (SINDy):** Coeficientes dispersos identificados — qué términos son activos

---

## Referencia

Publicación relacionada:
- L. I. Minchala-Avila *et al.*, "Data-Driven Linear Representations of Forced Nonlinear MIMO Systems via Hankel Dynamic Mode Decomposition with Lifting," *Mathematics* (MDPI), vol. 14, no. 4, 2026. DOI: [10.3390/math14040625](https://doi.org/10.3390/math14040625)
