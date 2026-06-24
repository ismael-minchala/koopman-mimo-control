# TCL Plant (Temperature Control Lab) — Koopman Eigenfunction Discovery

> **Planta:** Temperature Control Lab (TCLab) — sistema térmico MIMO 2×2.  
> **Tipo:** 2 entradas ($Q_1$, $Q_2$ — potencia de calentadores) / 2 salidas ($T_1$, $T_2$ — temperaturas).  
> **Objetivo:** Descubrimiento de eigenfunciones del operador de Koopman mediante conjuntos no recurrentes, con control MPC en coordenadas intrínsecas.

---

## Modelo de la planta

El TCLab es un sistema de laboratorio con dos calentadores y dos sensores de temperatura que interactúan térmicamente. Su dinámica es no lineal (pérdidas radiativas, acoplamiento cruzado):

$$\tau_1 \dot{T}_1 = U_a(T_{amb} - T_1) + U_{12}(T_2 - T_1) + \alpha_1 Q_1$$
$$\tau_2 \dot{T}_2 = U_a(T_{amb} - T_2) + U_{21}(T_1 - T_2) + \alpha_2 Q_2$$

con parámetros obtenidos experimentalmente (ver `Plant_parameters.xlsx`).

**Señales:**

| Señal | Rango | Descripción |
|-------|-------|-------------|
| $Q_1$, $Q_2$ | 0–100 % | Potencia de los calentadores |
| $T_1$, $T_2$ | ≈ 25–75 °C | Temperatura de los sensores |
| $\Delta t$ | 1 s | Período de muestreo |

**Dataset:** Series de tiempo de 10–12 horas, filtradas y normalizadas para entrenamiento.

---

## Archivos incluidos

```
src/
├── dataset_generator.m          — Generación del dataset de entrenamiento (PRBS)
├── koopman_nonrecurrent_sets.m  — Descubrimiento de eigenfunciones (algoritmo principal)
├── mpc_koopman.m                — Control MPC basado en el modelo Koopman (Hildreth QP)
└── utils/
    ├── hildreths.m              — Solver QP de Hildreth (MPC sin toolboxes)
    ├── mat_f_phi.m              — Construcción del vector de elevación φ(x)
    ├── mod_aumentado_mimo.m     — Modelo aumentado para el MPC
    ├── calc_accion_control.m    — Cálculo de la acción de control
    ├── energy_bal.m             — Balance de energía del calentador
    └── heat.m                   — Ecuaciones de calor del TCLab
```

**Datos** (descargar de Google Drive):
- `datos/Data4_12h_filter.mat` — 12 horas de datos experimentales filtrados (entrenamiento)
- `datos/Data1_10h_filter.mat` — Dataset de 10 horas (validación)
- `datos/parameters_controller.mat` — Parámetros del controlador identificado

---

## Descripción de los algoritmos

### `koopman_nonrecurrent_sets.m` — Eigenfunciones via conjuntos no recurrentes

Este es el algoritmo más novedoso del proyecto. A diferencia de EDMD, que fija un diccionario de funciones *a priori*, este método **descubre** las eigenfunciones $\varphi_j$ del operador de Koopman aprovechando la propiedad de que las eigenfunciones son constantes a lo largo de las órbitas del sistema:

$$\mathcal{K}\varphi_j = \lambda_j\,\varphi_j \quad \Leftrightarrow \quad \varphi_j(x_{k+1}) = \lambda_j\,\varphi_j(x_k)$$

**Algoritmo:**

1. **Preparación de datos:** Normaliza las series temporales y organiza $N_{\text{traj}}$ trayectorias de longitud fija.
2. **Identificación de conjuntos no recurrentes:** Detecta regiones del espacio de estados que no son visitadas repetidamente (donde las órbitas se separan). Estas regiones son informativas para identificar eigenfunciones.
3. **Estimación de eigenfunciones:** Ajusta una red de funciones radiales (RBF) o polinomiales que satisfacen la ecuación de eigenvalor en los datos.
4. **Validación:** Calcula el error de predicción en datos reservados ($y_{\text{hom}}$ vs. $y_{\text{for}}$) y el espectro del operador estimado.

**Parámetros clave:**

| Variable | Descripción |
|----------|-------------|
| `N_efun` | Número de eigenfunciones a identificar (default: 20) |
| `n_traj` | Número de trayectorias de entrenamiento (default: 40) |
| `n`, `m` | Estados (2) y entradas (2) |

### `mpc_koopman.m` — Control MPC con Hildreth

Implementa un MPC sobre el modelo Koopman linealizado usando el **algoritmo de Hildreth** como solver QP (sin necesidad de Optimization Toolbox):

$$\min_{\Delta U} \; \frac{1}{2}\Delta U^\top H\,\Delta U + f^\top \Delta U$$
$$\text{s.a.} \quad M\,\Delta U \leq \gamma$$

El problema se resuelve iterativamente mediante el método de Hildreth, eficiente para horizontes de predicción pequeños a medianos.

**Restricciones manejadas:**
- Límites de las entradas: $0 \le Q_i \le 100$
- Tasa de cambio máxima: $|\Delta Q_i| \le \Delta Q_{\max}$
- Límites de temperatura (opcionales)

---

## Cómo ejecutar

### Requisitos

- MATLAB R2021a o superior
- Hardware TCLab (opcional — con los datos `.mat` se puede simular sin hardware)
- **No requiere** Optimization Toolbox (solver Hildreth incluido)

### Pasos

```matlab
% 1. Navegar a la carpeta:
cd('ruta/a/codigo/03_tcl_plant/src')

% 2. Agregar utilidades al path:
addpath('utils')

% 3. Descargar datos y colocarlos en 03_tcl_plant/datos/

% 4. Generar un nuevo dataset (opcional):
dataset_generator
% Esto requiere el hardware TCLab conectado.
% Sin hardware: usar los datos pre-grabados en /datos/

% 5. Ejecutar la identificación de eigenfunciones:
koopman_nonrecurrent_sets
% Genera: eigenfunciones identificadas, espectro, error de validación.

% 6. Ejecutar el control MPC:
mpc_koopman
% Requiere los parámetros identificados (parameters_controller.mat).
```

### Salidas esperadas

Al ejecutar `koopman_nonrecurrent_sets.m`:
- **Figura 1:** Datos de entrenamiento (temperaturas y potencias)
- **Figura 2:** Trayectorias libres en el espacio de estados (sin control)
- **Figura 3:** Eigenfunciones identificadas proyectadas sobre el espacio de estados
- **Figura 4:** Error de predicción — validación cruzada

Al ejecutar `mpc_koopman.m`:
- **Figura 1:** Seguimiento de referencias de temperatura ($T_1^{ref}$, $T_2^{ref}$)
- **Figura 2:** Señales de control ($Q_1$, $Q_2$)
- **Figura 3:** Índices de desempeño ITAE/IAE

---

## Referencia

Este método corresponde al trabajo de identificación mediante conjuntos no recurrentes descrito en el artículo sometido a revisión:

- L. I. Minchala-Avila *et al.*, "Data-Driven Model Predictive Control based on Koopman Identification Via Constrained Non-Recurrent Sets," *IEEE Transactions on Control Systems Technology* (bajo revisión, 2026).
