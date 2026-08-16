# TP2 Macroeconometría — Exchange Rate Pass-Through (Colombia)

Trabajo Práctico 2 — Maestría en Economía, Universidad de San Andrés (2026).
Estimación del traspaso de tipo de cambio a precios (ERPT) en Colombia, con
datos mensuales de octubre 2003 a marzo 2023.

El análisis combina VAR estructurales (identificación de Cholesky con bloque
exógeno) y Local Projections (LP y LP-IV), incluyendo extensiones no lineales
por signo y por magnitud del shock, y una comparación entre ERPT esperado
(expectativas) y realizado.

## Estructura del repositorio

- `code/` — scripts del trabajo, numerados en orden de ejecución (`01`–`10`).
- `tools/` — funciones provistas por la cátedra (LP, SVAR, bootstrap, gráficos).
  Material de referencia; no es necesario correrlos por separado.
- `data_raw/` — bases crudas descargadas de las fuentes oficiales.
- `data/` — objetos intermedios generados por los scripts (`.rds`).
- `output/` — figuras (`.pdf`) y tablas (`.csv`) que produce el análisis.
- `TP2_macrometrics.Rproj` — proyecto de RStudio (fija la raíz para `here()`).
- `MacroMetrics_TP2_2026.pdf` — consigna del trabajo.
- `Ficha_Variables.pdf` — detalle de variables y fuentes.
- `TP2_Macroeconometria_Arechavala_Blum_Fradkin.pdf` — informe final.

## Requisitos

Todas las rutas se resuelven con el paquete `here` a partir de la raíz del
proyecto, de modo que el repositorio corre en cualquier ordenador sin editar
rutas: basta con abrir `TP2_macrometrics.Rproj`.

Paquetes de R utilizados:

```r
install.packages(c(
  "here", "readxl", "readr", "dplyr", "tidyr", "lubridate", "seasonal",
  "ggplot2", "patchwork", "cowplot", "sandwich", "vars"
))
```

`seasonal` requiere el motor X-13ARIMA-SEATS; en instalaciones recientes viene
incluido con el paquete. `sandwich` se usa para los errores estándar HAC de las
Local Projections.

## Cómo correr

1. Abrir `TP2_macrometrics.Rproj` en RStudio.
2. Instalar los paquetes de la lista anterior (una sola vez).
3. Correr los scripts de `code/` **en orden**, de `01` a `10`.

**Correr cada script en una sesión de R limpia** (en RStudio: *Session → Restart
R*, o `Ctrl+Shift+F10`, antes de cada uno). Esto evita colisiones de nombres
entre paquetes (por ejemplo, `seasonal`/`MASS` enmascaran `dplyr::select`). Cada
script es autónomo: hace su propia limpieza inicial, carga sus librerías, lee sus
insumos desde `data/` y guarda sus resultados en `data/` y `output/`.

Los scripts posteriores reutilizan objetos guardados por los anteriores (no se
reestima nada dos veces). Por eso el orden importa: para correr un script del
medio, primero deben existir los `.rds` que consume, generados por los scripts
previos. Si se conserva la carpeta `data/` provista, cualquier script puede
correrse de forma aislada.

## Scripts y qué produce cada uno

| Script | Contenido | Insumos (`data/`) | Productos principales |
|---|---|---|---|
| `01_Datos.R` | Lee `data_raw/`, limpia y transforma las series, arma el panel mensual. | — (lee `data_raw/`) | `panel_final.rds`, `datos_ts.rds` |
| `02_Ejercicios1.2.R` | Incisos 1–2: gráficos de dos ejes, tests de raíz unitaria, VAR bivariado (TCN–IPC) e IRF/ERPT base. | `datos_ts.rds`, `panel_final.rds` | `objetos_02.rds`; `graficos_dos_ejes.pdf`, `tests_raiz_unitaria.csv`, `irf_punto2_tcn_ipc.pdf`, `erpt_punto2.pdf` |
| `03_Ejercicio3.R` | Inciso 3: IRF y ERPT del sistema bivariado por VAR y por LP; identificación del shock estructural al TCN. | `datos_ts.rds`, `panel_final.rds`, `objetos_02.rds` | `objetos_03.rds`; `irf_punto3_tcn_var_lp.pdf`, `irf_punto3_ipc_var_lp.pdf`, `erpt_punto3_var_lp.pdf` |
| `04_Ejercicio4.R` | Inciso 4: IRF y ERPT asimétricos por **signo** del shock (LP piecewise, dos especificaciones de controles). | `objetos_03.rds` | `irf_punto4_tcn_pw.pdf`, `irf_punto4_ipc_pw.pdf`, `erpt_punto4_pw.pdf`, 4 anexos de robustez (`anexo_*_spec.pdf`) |
| `05_Ejercicio5.R` | Inciso 5: histograma del shock estructural y ERPT condicional por **magnitud** (tres regímenes a ±1 desvío). | `objetos_03.rds` | `hist_punto5_shock.pdf`, `erpt_punto5_size.pdf` |
| `06_Ejercicio6.R` | Inciso 6: VAR de 7 variables (Fed, petróleo, EBP, TOT, EMBI, TCN, IPC), Cholesky con bloque exógeno; los 6 shocks estructurales alimentan una LP-IV. | `datos_ts.rds` | `objetos_06.rds`; `erpt_punto6_LPIV.pdf`, `tabla_erpt_punto6.csv` |
| `07_Ejercicio7.R` | Inciso 7: ERPT asimétrico por **signo**, para los 6 shocks del VAR7 (extensión piecewise del LP-IV). | `objetos_06.rds` | `erpt_punto7_asimetrico.pdf`, `tabla_erpt_punto7.csv` |
| `08_Ejercicio8.R` | Inciso 8: IRF de expectativas por LP (mismo shock y controles del inciso 3, cambiando la variable dependiente). | `datos_ts.rds`, `objetos_02.rds`, `objetos_03.rds` | `objetos_04.rds`; `08_irf_inflacion.pdf`, `08_irf_tcn.pdf`, `08_irf_diferencias.pdf` |
| `09_Ejercicio9.R` | Inciso 9: ERPT esperado vs. realizado (forward) para el shock al TCN. | `objetos_04.rds` | `09_erpt_esp_vs_real.pdf` |
| `10_Ejercicio10.R` | Inciso 10 (extra): ERPT esperado vs. realizado para los 6 shocks del VAR7. | `objetos_06.rds`, `datos_ts.rds` | `erpt_punto10_esp_vs_real.pdf`, `tabla_erpt_punto10.csv` |

> Nota: los nombres de los `.rds` intermedios siguen la numeración del inciso que
> los genera (`objetos_02.rds` proviene del script `02`, etc.). El `objetos_04.rds`
> lo produce el inciso 8 por la numeración histórica de la sección de expectativas.

## Fuentes de datos (`data_raw/`)

- `ipc.xlsx` — IPC de Colombia.
- `trm_historico.csv` — tipo de cambio nominal (TRM).
- `imf_pctot_xm_RW_IX.csv` — términos de intercambio (IMF/PCTOT).
- `Daily_Embi_Panel_wide_selected.xlsx` — spread EMBI.
- `eme_expectativas.xlsx` — expectativas de inflación y de tipo de cambio.
- `fed_jk_shocks.csv` — shocks de política monetaria de la Fed (Jarociński–Karadi).
- `oil_supply_bh.xlsx` — shocks de oferta de petróleo (Baumeister–Hamilton).
- `ebp.csv` — Excess Bond Premium (Gilchrist–Zakrajšek).

## Integrantes

- Mateo Arechavala
- Felicitas Blum
- Martín Fradkin