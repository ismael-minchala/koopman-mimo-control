# Exoesqueleto Robótico — Modelo de Péndulo Doble

> **Planta:** Exoesqueleto de extremidades inferiores modelado como péndulo doble articulado.  
> **Tipo:** Sistema mecánico no lineal MIMO con 4 estados y 2 entradas de torque.  
> **Estado:** Material de referencia y punto de partida para trabajo futuro. El exoesqueleto físico no estuvo disponible durante el período de ejecución del proyecto (2024–2026); la actividad fue reorientada hacia el gemelo digital del sistema multitanque.  
> **Nota:** Los scripts de esta carpeta son punto de partida para investigación futura sobre control de exoesqueletos con el operador de Koopman.

---

## Descripción del sistema

El exoesqueleto se modela como un péndulo doble planar con dos eslabones (muslo y pierna) conectados por articulaciones rotacionales (cadera y rodilla):

$$M(q)\ddot{q} + C(q,\dot{q})\dot{q} + G(q) = \tau$$

donde:
- $q = [\theta_1, \theta_2]^\top$ son los ángulos de cadera y rodilla
- $M(q)$ es la matriz de inercia (dependiente de configuración)
- $C(q,\dot{q})$ es la matriz de Coriolis/centrífuga
- $G(q)$ es el vector de gravedad
- $\tau = [\tau_1, \tau_2]^\top$ son los torques de los actuadores

**Parámetros del modelo (péndulo doble):**

| Parámetro | Descripción |
|-----------|-------------|
| $m_1, m_2$ | Masas de muslo y pierna |
| $l_1, l_2$ | Longitudes de los eslabones |
| $I_1, I_2$ | Inercias de cada segmento |
| $g$ | Gravedad (9.81 m/s²) |

---

## Archivos incluidos

```
src/
├── pendulo_doble_params.m  — Parámetros físicos del péndulo doble (punto de partida)
└── interpolacion.m         — Rutinas de interpolación de trayectorias de marcha
```

---

## Descripción de los archivos

### `pendulo_doble_params.m`

Define los parámetros físicos del péndulo doble (masas, inercias, longitudes) y configura el modelo en Simulink/MATLAB. Este archivo es el punto de partida para cargar el modelo dinámico.

### `interpolacion.m`

Implementa la interpolación de trayectorias de referencia de marcha (ángulos de cadera y rodilla) a partir de datos de captura de movimiento. Útil para generar referencias suaves para el controlador.

---

## Hoja de ruta para trabajo futuro

Este módulo está pensado como punto de partida. Los pasos sugeridos para una implementación completa con Koopman:

1. **Modelado:** Completar la parametrización física del exoesqueleto específico (Exoesqueleto UC).
2. **Generación de datos:** Recolectar datos de operación (simulados o experimentales) con excitación PRBS en los torques.
3. **Identificación Koopman:** Adaptar `koopman_edmd.m` del péndulo de Furuta para 4 estados y 2 entradas.
4. **Control MPC:** Aplicar `mpc_koopman.m` con las restricciones físicas de las articulaciones ($\theta_{\min} \le \theta_i \le \theta_{\max}$, $\tau_{\max}$).
5. **Validación:** Comparar contra controladores PID y espacio de estados en trayectorias de marcha.

---

## Cómo ejecutar los scripts disponibles

```matlab
% 1. Navegar a la carpeta:
cd('ruta/a/codigo/05_exoesqueleto/src')

% 2. Cargar parámetros del péndulo doble:
pendulo_doble_params
% Esto define las variables del espacio de trabajo para el modelo.

% 3. (Opcional) Generar y visualizar trayectorias interpoladas:
interpolacion
```

---

## Referencia

El modelo de péndulo doble para exoesqueleto fue desarrollado como parte de la actividad inicial del proyecto antes de la reorientación hacia el gemelo digital del sistema multitanque.

Para el desarrollo futuro, se recomienda consultar:
- L. I. Minchala-Avila *et al.*, "Control based on the Koopman operator: A comprehensive review," *J. Franklin Institute*, 2025. DOI: [10.1016/j.jfranklin.2025.108256](https://doi.org/10.1016/j.jfranklin.2025.108256) — Sección sobre robótica y exoesqueletos.
