# Koopman Operator-Based Control for Nonlinear MIMO Systems
## Repositorio de Código — UC-DEET Research Project (2024–2026)

[![Universidad de Cuenca](https://img.shields.io/badge/Universidad%20de%20Cuenca-DEET-blue)](https://www.ucuenca.edu.ec)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2020b%2B-orange)](https://www.mathworks.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)

> **Proyecto:** Control activo tolerante a fallas basado en el operador de Koopman para sistemas MIMO no lineales  
> **Institución:** Universidad de Cuenca, Departamento de Ingeniería Eléctrica, Electrónica y Telecomunicaciones (DEET)  
> **Período:** Marzo 2024 – Febrero 2026  
> **Director:** Dr. Luis Ismael Minchala Ávila — ismael.minchala@ucuenca.edu.ec  
> **Código de proyecto:** XX Concurso Universitario de Proyectos de Investigación

---

## Descripción

Este repositorio contiene las implementaciones MATLAB de los algoritmos de identificación y control basados en el **operador de Koopman** desarrollados durante el proyecto de investigación. El operador de Koopman transforma sistemas dinámicos no lineales en representaciones lineales de dimensión elevada, habilitando técnicas de control lineal (MPC, LQR) para sistemas intrínsecamente no lineales.

La representación central es:

$$\dot{x} = f(x) + Bu \quad \xrightarrow{\quad\mathcal{K}\quad} \quad \dot{\varphi} = \mathcal{K}\varphi$$

donde $\varphi = \Psi(x)$ son las **eigenfunciones** del operador, y $\mathcal{K}$ actúa linealmente sobre el espacio de observables.

---

## Estructura del repositorio

```
codigo/
├── 01_sistema_multitanque/   — Sistema de 3 tanques interconectados (MIMO 3u/2y)
│   ├── README.md             — Documentación detallada + modelo + instrucciones
│   ├── src/
│   │   ├── koopman_edmd_mpc.m      ← ALGORITMO PRINCIPAL: EDMD + MPC
│   │   ├── mpc_koopman_lineal.m    ← MPC sobre modelo linealizado Koopman
│   │   ├── fdd_espacio_paridad.m   ← Detección de fallas: espacios de paridad
│   │   ├── fdd_filtro_kalman.m     ← Detección y estimación de fallas: Kalman
│   │   └── utils/                  ← Funciones auxiliares (Hildreth QP, etc.)
│   └── datos/                      ← Colocar aquí los .mat/.csv (ver abajo)
│
├── 02_pendulo_furuta/        — Péndulo de Furuta / QUBE-Servo 2 (MIMO 4x/1u)
│   ├── README.md
│   ├── src/
│   │   ├── parameters.m            ← Parámetros físicos del sistema
│   │   ├── koopman_edmd.m          ← EDMD — estimación del operador Koopman
│   │   ├── koopman_sindy.m         ← SINDy — identificación no lineal dispersa
│   │   └── observer_design.m       ← Diseño del observador de estado
│   └── datos/
│
├── 03_tcl_plant/             — Temperature Control Lab 2×2 (MIMO 2u/2y)
│   ├── README.md
│   ├── src/
│   │   ├── dataset_generator.m         ← Generación de datos de entrenamiento
│   │   ├── koopman_nonrecurrent_sets.m ← Eigenfunciones vía conjuntos no recurrentes
│   │   ├── mpc_koopman.m               ← MPC con solver Hildreth
│   │   └── utils/                      ← Solver QP, modelos de calor, etc.
│   └── datos/
│
├── 04_turbina_eolica/        — Predicción de desviación de frecuencia (SISO)
│   ├── README.md
│   ├── src/
│   │   ├── sindy_prediction.m      ← SINDy para dinámica de frecuencia
│   │   └── koopman_estimation.m    ← Koopman alternativo
│   └── datos/
│
└── 05_exoesqueleto/          — Péndulo doble (referencia / trabajo futuro)
    ├── README.md
    └── src/
        ├── pendulo_doble_params.m
        └── interpolacion.m
```

---

## Inicio rápido

### Prerrequisitos

- **MATLAB R2020b** o superior (probado hasta R2023b)
- Toolboxes recomendados: **Control System Toolbox**, **Signal Processing Toolbox**
- **No se requiere** Optimization Toolbox — el solver QP de Hildreth está incluido

### Instalación

```bash
# Clonar el repositorio
git clone https://github.com/ismael-minchala/koopman-mimo-control
cd koopman-mimo-control
```

En MATLAB:

```matlab
% Agregar todas las carpetas al path
addpath(genpath('codigo'))
```

### Datos de entrada

Los datos de experimentación (archivos `.mat` y `.csv`) **no están incluidos** en el repositorio debido a su tamaño. Descárgarlos desde:

📁 **[Google Drive — Repositorio del Proyecto](https://drive.google.com/drive/folders/1PJfk9FfFzaYXfAU78d7J0IUnML7zHtIU?usp=sharing)**

Colocar cada archivo en la subcarpeta `datos/` correspondiente a cada planta.

| Planta | Archivos necesarios |
|--------|-------------------|
| `01_sistema_multitanque/datos/` | `Dataset_13.mat`, `linear_control_Results.mat`, archivos `data_*.csv` |
| `02_pendulo_furuta/datos/` | `pendulum_time_series.mat`, `pendulum_time_series_test.mat` |
| `03_tcl_plant/datos/` | `Data4_12h_filter.mat`, `Data1_10h_filter.mat` |
| `04_turbina_eolica/datos/` | `Datta.mat` |

### Ejecutar el ejemplo principal

```matlab
% Sistema multitanque — flujo completo EDMD → Koopman → MPC
cd('codigo/01_sistema_multitanque/src')
addpath('utils')
koopman_edmd_mpc   % Ejecuta identificación + simulación en lazo cerrado
```

---

## Algoritmos implementados

| Algoritmo | Planta | Descripción |
|-----------|--------|-------------|
| **EDMD** | Multitanque, Furuta | Descomposición de Modos Dinámicos Extendida para estimación de Koopman |
| **SINDy** | Furuta, Turbina | Identificación No Lineal Dispersa — modelos interpretables |
| **Eigenfunciones (conjuntos no recurrentes)** | TCL | Descubrimiento automático de eigenfunciones del operador |
| **MPC Koopman (Hildreth QP)** | Multitanque, TCL | Control Predictivo sobre modelo Koopman sin Optimization Toolbox |
| **FDD — Espacios de Paridad** | Multitanque | Detección y aislamiento de fallas basada en residuos de paridad |
| **FDD — Filtro de Kalman** | Multitanque | Estimación simultánea de estado y magnitud de falla |
| **Observador de Estado** | Furuta | Estimación de velocidades angulares no medidas |

---

## Publicaciones relacionadas

1. L. I. Minchala-Avila *et al.*, **"Control based on the Koopman operator: A comprehensive review"**, *Journal of the Franklin Institute*, vol. 362, no. 8, 2025. [DOI: 10.1016/j.jfranklin.2025.108256](https://doi.org/10.1016/j.jfranklin.2025.108256)

2. L. I. Minchala-Avila *et al.*, **"Data-Driven Linear Representations of Forced Nonlinear MIMO Systems via Hankel Dynamic Mode Decomposition with Lifting"**, *Mathematics* (MDPI), vol. 14, no. 4, art. 625, 2026. [DOI: 10.3390/math14040625](https://doi.org/10.3390/math14040625)

3. L. I. Minchala-Avila *et al.*, **"Digital Twin Framework for Multi-Tank Systems with OPC UA Interoperability and Real-Time Control and Visualization"**, *IEEE*, 2025. [IEEE Xplore](https://ieeexplore.ieee.org/abstract/document/11304457)

4. L. I. Minchala-Avila *et al.*, **"An Experimental Comparison of Model-Based Fault Detection and Isolation Techniques in a Multi-Tank System"**, *IEEE ETCM*, 2024. [DOI: 10.1109/ETCM63562.2024.10746149](https://doi.org/10.1109/ETCM63562.2024.10746149)

5. L. I. Minchala-Avila *et al.*, **"Data-Driven Model Predictive Control based on Koopman Identification Via Constrained Non-Recurrent Sets"**, *IEEE Transactions on Control Systems Technology* (bajo revisión, 2026).

---

## Equipo

| Rol | Nombre |
|-----|--------|
| Director | Dr. Luis Ismael Minchala Ávila |
| Investigador | Alcides Fabián Araujo Pacheco |
| Investigador asociado | Dr. Luis Eduardo Garza Castañón |
| Técnico de investigación / Desarrollador principal | Juan Francisco Durán Siguenza |
| Tesista de Maestría | Ing. Marcos Lenin Villarreal Esquivel |
| Tesista de Pregrado | Luis Antonio Andrade Matute |

---

## Licencia

MIT License — ver [LICENSE](LICENSE). Al citar el código, por favor referenciar las publicaciones listadas arriba.

---

## Contacto

Dr. Luis Ismael Minchala Ávila  
Universidad de Cuenca — DEET  
✉ ismael.minchala@ucuenca.edu.ec  
🌐 [Página del Proyecto](https://ismael-minchala.github.io/koopman-mimo-control/)
