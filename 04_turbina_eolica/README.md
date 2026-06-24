# Turbina Eólica — Predicción de Desviación de Frecuencia con Koopman/SINDy

> **Planta:** Sistema de potencia con generación eólica — modelo de desviación de frecuencia.  
> **Tipo:** SISO (predicción de desviación de frecuencia $\Delta f$ bajo perturbaciones de potencia eólica).  
> **Objetivo:** Modelamiento data-driven de la dinámica de frecuencia mediante el operador de Koopman y SINDy, para soporte en diseño de controladores de frecuencia robustos.

---

## Descripción del problema

En sistemas de potencia con alta penetración de energía eólica, la variabilidad del recurso genera desviaciones de frecuencia que pueden comprometer la estabilidad de la red. La dinámica del sistema puede modelarse como:

$$\Delta\dot{f} = g(\Delta f, P_{\text{wind}})$$

donde $g(\cdot)$ es no lineal e incierta. Los algoritmos incluidos identifican esta dinámica directamente desde datos medidos, sin necesidad de un modelo físico detallado.

---

## Archivos incluidos

```
src/
├── sindy_prediction.m    — Identificación SINDy de la dinámica de frecuencia (principal)
└── koopman_estimation.m  — Estimación Koopman alternativa
```

**Datos** (descargar de Google Drive):
- `datos/Datta.mat` — Dataset de desviación de frecuencia (múltiples escenarios de viento)

---

## Descripción de los algoritmos

### `sindy_prediction.m` — SINDy para predicción de frecuencia

**SINDy (Sparse Identification of Nonlinear Dynamics)** identifica la estructura del modelo seleccionando el mínimo conjunto de términos de una biblioteca de candidatos:

$$\dot{x} = \boldsymbol{\Theta}(x)\,\boldsymbol{\Xi}, \quad \text{con } \boldsymbol{\Xi} \text{ disperso}$$

La biblioteca $\boldsymbol{\Theta}$ incluye: polinomios $[1,\, x,\, x^2,\, x^3,\, \ldots]$, funciones trigonométricas y productos cruzados con la entrada eólica $P_{\text{wind}}$.

**Flujo del script:**

1. **Carga de datos:** Lee `Datta.mat`, extrae la desviación de frecuencia y la potencia eólica normalizada.
2. **Construcción de la biblioteca $\Theta(x)$:** Evalúa cada función candidata en los datos.
3. **Regresión dispersa:** Resuelve iterativamente el LASSO con umbral (Sequential Thresholding Least Squares, STLS):
   $$\boldsymbol{\Xi}^* = \arg\min_{\boldsymbol{\Xi}} \|\dot{X} - \boldsymbol{\Theta}(X)\boldsymbol{\Xi}\|_F^2 + \lambda\|\boldsymbol{\Xi}\|_1$$
4. **Identificación del modelo:** Muestra los términos activos — el modelo es explícito e interpretable.
5. **Validación:** Simula el modelo identificado $k$-pasos hacia adelante y compara con datos reales.

**Ventajas de SINDy sobre EDMD aquí:**
- La dinámica de frecuencia tiene estructura física conocida → SINDy recupera términos interpretables
- Dataset relativamente pequeño → SINDy requiere menos datos que EDMD con diccionario grande

### `koopman_estimation.m` — Koopman para desviación de frecuencia

Implementación alternativa basada en la representación lineal global del operador de Koopman, útil cuando la estructura del modelo no es conocida a priori.

---

## Cómo ejecutar

### Requisitos

- MATLAB R2020b o superior
- No requiere toolboxes adicionales

### Pasos

```matlab
% 1. Navegar a la carpeta:
cd('ruta/a/codigo/04_turbina_eolica/src')

% 2. Descargar datos y colocarlos en 04_turbina_eolica/datos/

% 3. Ejecutar la identificación SINDy:
sindy_prediction
% El script carga Datta.mat, identifica el modelo y valida la predicción.

% 4. (Opcional) Comparar con Koopman puro:
koopman_estimation
```

### Salidas esperadas

- **Figura 1:** Dataset de entrenamiento — desviación de frecuencia y potencia eólica
- **Figura 2:** Términos del modelo identificado (coeficientes $\boldsymbol{\Xi}^*$)
- **Figura 3:** Predicción $k$-pasos del modelo vs. datos reales de prueba
- **Consola:** Ecuación del modelo identificado impresa en forma legible

---

## Nota sobre los datos

El dataset `Datta.mat` contiene mediciones reales/simuladas de:
- Desviación de frecuencia $\Delta f$ (%) para múltiples turbinas
- Potencia eólica normalizada $P_{\text{wind}}$ como entrada exógena
- Muestreo: extraído según el paso de tiempo del dataset (variable `deltaT`)

Los datos de múltiples turbinas están organizados en columnas; el script usa por defecto la primera columna para entrenamiento y la segunda para validación.
