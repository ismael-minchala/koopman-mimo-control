# Sistema Multitanque — Koopman MIMO Control

> **Planta:** Sistema de tres tanques interconectados (laboratorio de automatización industrial, Universidad de Cuenca, Campus Balzay).  
> **Tipo:** MIMO 3 entradas / 2 salidas (niveles $h_1$, $h_2$; bombas $u_1$, $u_2$, $u_3$).  
> **Objetivo:** Control de nivel mediante MPC basado en el operador de Koopman (EDMD) con detección y diagnóstico de fallas (FDD).

---

## Modelo de la planta

El sistema está descrito por las ecuaciones diferenciales no lineales:

$$\dot{h}_1 = \frac{q_0\,u_1}{A_1} - \frac{a_2(u_2+\alpha_1)\sqrt{2g\,h_1}}{A_1}$$

$$\dot{h}_2 = \frac{a_2(u_2+\alpha_1)\sqrt{2g\,h_1}}{A_2} - \frac{a_3(u_3+\alpha_2)\sqrt{2g\,h_2}}{A_2}$$

**Parámetros:**

| Símbolo | Valor | Descripción |
|---------|-------|-------------|
| $A_1, A_2$ | 0.04 m² | Sección transversal de los tanques |
| $a_2$ | 4.17 × 10⁻⁵ m² | Área del orificio entre tanques |
| $a_3$ | 1.81 × 10⁻⁵ m² | Área del orificio de salida |
| $q_0$ | 3.37 × 10⁻⁵ m³/s/u | Ganancia de la bomba |
| $g$ | 9.81 m/s² | Gravedad |
| $\Delta t$ | 1 s | Período de muestreo |

**Punto de operación nominal:** $h_{e1} = 0.0879$ m, $h_{e2} = 0.1153$ m.

---

## Archivos incluidos

```
src/
├── koopman_edmd_mpc.m       — Identificación Koopman-EDMD + control MPC (script principal)
├── mpc_koopman_lineal.m     — MPC sobre modelo linealizado con capacidades Koopman
├── fdd_espacio_paridad.m    — Detección de fallas por espacios de paridad (MPC baseline)
├── fdd_filtro_kalman.m      — Estimación y aislamiento de fallas con filtro de Kalman
└── utils/
    ├── hildreths.m          — Solver QP de Hildreth para MPC sin librerías externas
    ├── mat_f_phi.m          — Construcción del vector de elevación φ(x)
    ├── mod_aumentado_mimo.m — Modelo aumentado MIMO para diseño del controlador
    └── calc_accion_control.m — Cálculo de la acción de control óptima
```

**Datos** (descargar de Google Drive — ver enlace principal):
- `datos/Dataset_13.mat` — Dataset simulado de 300 000 pasos para entrenamiento
- `datos/linear_control_Results.mat` — Resultados de referencia del control lineal
- `datos/data_*.csv` — Datos experimentales reales del sistema multitanque

---

## Descripción de los algoritmos

### `koopman_edmd_mpc.m` — Algoritmo principal

Implementa el flujo completo **EDMD → Identificación Koopman → MPC**:

1. **Recolección de datos:** Carga el `Dataset_13.mat` (simulado con PRBS).
2. **Construcción del diccionario:** Selecciona funciones de elevación $\Psi(x)$ de tipo polinomial.
3. **Estimación EDMD:** Resuelve el problema de mínimos cuadrados:
   $$\min_{A,C} \sum_{i} \|\Phi_{i+1} - A\Phi_i\|^2 + \|y_i - C\Phi_i\|^2$$
   donde $\Phi_i = \Psi(x_i)$ es el vector de variables elevadas.
4. **Diseño del MPC:** Usa el modelo bilineal de Koopman $\Phi_{k+1} = A\Phi_k + Bu_k$ para formular el problema de optimización con horizonte $N_p$.
5. **Validación:** Simula el lazo cerrado y grafica el seguimiento de referencia, el error y las entradas de control.

**Parámetros ajustables (parte superior del script):**

| Variable | Descripción |
|----------|-------------|
| `Np` | Horizonte de predicción del MPC |
| `Nc` | Horizonte de control |
| `Q`, `R` | Matrices de ponderación del MPC |
| `N_efun` | Número de funciones de elevación (diccionario) |

### `fdd_espacio_paridad.m` — FDD por espacios de paridad

Detecta fallas en actuadores y sensores construyendo un residuo de paridad:

$$r_k = V_s \begin{bmatrix} y_k \\ y_{k-1} \\ \vdots \\ y_{k-s} \end{bmatrix} - V_s H_s \begin{bmatrix} u_k \\ \vdots \\ u_{k-s} \end{bmatrix}$$

El umbral de detección $J_{th}$ se calibra estadísticamente a partir de datos sin falla.

### `fdd_filtro_kalman.m` — FDD con filtro de Kalman

Estima el vector de estado aumentado $[x; f]$ (estado + falla) mediante un filtro de Kalman discreto sobre el modelo linealizado, permitiendo detección y estimación simultánea de la magnitud de la falla.

---

## Cómo ejecutar

### Requisitos

- MATLAB R2020b o superior (probado hasta R2023b)
- Toolboxes: **Control System**, **Optimization** (o solver QP propio `hildreths.m` incluido)
- No se requiere Simulink para los scripts `.m`

### Pasos

```matlab
% 1. Desde MATLAB, navegar a la carpeta de la planta:
cd('ruta/a/codigo/01_sistema_multitanque/src')

% 2. Agregar las utilidades al path:
addpath('utils')

% 3. Descargar los datos (ver sección "Datos") y colocarlos en:
%    01_sistema_multitanque/datos/

% 4. Ajustar la ruta de datos en el script (buscar 'AJUSTAR_RUTA'):
%    load(fullfile('..','datos','Dataset_13.mat'))
%    (ya está ajustada por defecto)

% 5. Ejecutar el algoritmo principal:
koopman_edmd_mpc

% 6. Para probar detección de fallas:
fdd_espacio_paridad
% o:
fdd_filtro_kalman
```

### Salidas esperadas

Al ejecutar `koopman_edmd_mpc.m` se generan:

- **Figura 1:** Dataset de entrenamiento ($h_1$, $h_2$, entradas $u$)
- **Figura 2:** Valor singular de la matriz de Gramm (selección del número de modos)
- **Figura 3:** Validación del modelo Koopman — predicción vs. datos de prueba
- **Figura 4:** Simulación en lazo cerrado — seguimiento de referencias escalón
- **Figura 5:** Señales de control

---

## Referencia

Script desarrollado por Juan Francisco Durán Siguenza bajo la dirección del Dr. Luis Ismael Minchala Ávila como parte del proyecto *"Control activo tolerante a fallas basado en el operador de Koopman"* (UC-DEET, 2024–2026).

Publicaciones relacionadas:
- L. I. Minchala-Avila *et al.*, "Control based on the Koopman operator: A comprehensive review," *J. Franklin Institute*, 2025. DOI: [10.1016/j.jfranklin.2025.108256](https://doi.org/10.1016/j.jfranklin.2025.108256)
- L. I. Minchala-Avila *et al.*, "An experimental comparison of model-based fault detection and isolation techniques in a multi-tank system," *IEEE ETCM*, 2024. DOI: [10.1109/ETCM63562.2024.10746149](https://doi.org/10.1109/ETCM63562.2024.10746149)
