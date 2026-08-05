# TP2 Macroeconometría — Exchange Rate Pass-Through (Colombia)

Trabajo Práctico 2 — Maestría en Economía, Universidad de San Andrés (2026).
Análisis del traspaso de tipo de cambio a precios (ERPT) para Colombia.

## Estructura
- `code/01_datos.R` — procesamiento de datos: lectura, limpieza, transformaciones y construcción del panel.
- `data_raw/` — bases crudas descargadas de las fuentes oficiales.
- `data/` — panel procesado (`panel_final.csv`, `datos_ts.rds`).
- `Gráficos/` — figuras del trabajo.

## Cómo correr
1. Abrir `TP2_macrometrics.Rproj` (o cualquier IDE; el script usa `here` para las rutas).
2. Instalar paquetes: `install.packages(c("readxl","dplyr","readr","lubridate","tidyr","seasonal","here"))`
3. Correr `code/01_datos.R`.

## Datos
Ventana muestral: octubre 2003 – marzo 2023 (234 obs. mensuales).
Fuentes detalladas en `Fichas_Variables_TP2.docx`.

## Integrantes
- Martín [apellido]
- Felicitas Blum
- Mateo [apellido]