#------------------------------------------------------------------------------#
# Universidad de San Andrés - Macroeconometría 2026
# TP2: Exchange Rate Pass-Through - Colombia
# 01_Procesamiento_Datos.R
#
# Este script lee las bases crudas desde 'data_raw/' y construye las series
# transformadas para el análisis. Todas las series locales quedan en frecuencia
# mensual.
#
# Convención de tipo de cambio: pesos colombianos por dólar (COP/USD).
#   -> un aumento del índice = DEPRECIACIÓN del peso.
#
# Transformaciones (según enunciado TP2):
#   - TCN, IPC, términos de intercambio, exp. de TCN  -> 100 * log(x)
#   - Expectativa de inflación a 12m                  -> 100 * log(1 + x/100)
#     (en el Excel de la EME x viene en decimal, i.e. x/100 ya está aplicado)
#   - EMBI                                            -> niveles (sin log)
#------------------------------------------------------------------------------#

#Una vez:
#install.packages(c("readxl","dplyr","readr","lubridate","tidyr","seasonal"))
#install.packages("here")

#Si "seasonal" genera problemas:
#install.packages("seasonal", repos="https://cloud.r-project.org")



remove(list = ls(all.names = TRUE))
gc()

# install.packages(c("readxl", "dplyr", "readr", "lubridate", "seasonal"))
library(readxl)
library(dplyr)
library(readr)
library(lubridate)
library(here)

#------------------------------------------------------------------------------#
# 0. Directorio de trabajo y rutas
#    Preferimos rutas RELATIVAS. Abrir el proyecto de RStudio en la raíz del TP,
#    o setear aquí manualmente la raíz UNA sola vez. No usar rutas absolutas de
#    una máquina en particular.
#------------------------------------------------------------------------------#

# Si usás RStudio Project, comentá la línea de setwd y descomentá here::here().
# setwd("<raíz del proyecto TP2>")   # <- cada integrante ajusta SOLO esto si hace falta

ruta_raw <- here("data_raw")   # carpeta con las bases crudas
ruta_out <- here("data")      # carpeta con las bases crudas
if (!dir.exists(ruta_out)) dir.create(ruta_out, showWarnings = FALSE)






#------------------------------------------------------------------------------#
# 1. Expectativas de inflación a 12 meses (EME - Banco de la República)
#    Fuente: eme_expectativas.xlsx, hoja INFLACION_TOTAL
#    Columna 10 = MEDIANA del bloque "EXPECTATIVAS DE INFLACIÓN A DOCE MESES"
#    (12 meses MÓVILES, no fin de año).  x viene en decimal (0.061 = 6.1%).
#    Nota fechas: la EME registra en años recientes el día exacto de la encuesta
#    (ej. 2026-06-12); normalizamos al primer día del mes con floor_date para
#    tener una clave temporal homogénea de cara al merge.
#------------------------------------------------------------------------------#

archivo_eme <- file.path(ruta_raw, "eme_expectativas.xlsx")

exp_inflacion <- read_excel(
  path = archivo_eme,
  sheet = "INFLACION_TOTAL",
  skip = 4,
  col_names = FALSE
) |>
  dplyr::select(
    fecha = 1,
    exp_inf_12m_orig = 10
  ) |>
  filter(!is.na(fecha), !is.na(exp_inf_12m_orig)) |>
  mutate(
    fecha = floor_date(as.Date(fecha), unit = "month"),
    # x/100 ya está aplicado (dato en decimal), por eso 100*log(1 + x):
    exp_inf_12m = 100 * log(1 + exp_inf_12m_orig)
  ) |>
  arrange(fecha)


#chequeos:deberías ver 273 observaciones, desde octubre 2003 hasta mediados de 2026, con las fechas ya normalizadas a día 01
str(exp_inflacion)
head(exp_inflacion)
tail(exp_inflacion)
nrow(exp_inflacion)
range(exp_inflacion$fecha)
summary(exp_inflacion$exp_inf_12m)


#------------------------------------------------------------------------------#
# 2. Expectativa de tipo de cambio (TRM) a 12 meses (EME - Banco de la República)
#    Fuente: eme_expectativas.xlsx, hoja TRM
#    Columna 10 = MEDIANA del bloque "+ 12 meses". Viene en NIVEL (pesos/dólar).
#    Nota metodológica: hasta dic-2014 se preguntó TRM promedio mensual;
#    desde ene-2015, TRM del último día del mes.
#    Misma normalización de fechas (floor_date) que la sección 1.
#------------------------------------------------------------------------------#

exp_tipo_cambio <- read_excel(
  path = archivo_eme,
  sheet = "TRM",
  skip = 5,
  col_names = FALSE
) |>
  select(
    fecha = 1,
    exp_tcn_12m_orig = 10
  ) |>
  filter(!is.na(fecha), !is.na(exp_tcn_12m_orig)) |>
  mutate(
    fecha = floor_date(as.Date(fecha), unit = "month"),
    exp_tcn_12m = 100 * log(exp_tcn_12m_orig)
  ) |>
  arrange(fecha)

#chequeos: deberías ver 274 obs, desde 2003-09-01 hasta 2026-06-01, con las fechas ya normalizadas a día 01
str(exp_tipo_cambio)
head(exp_tipo_cambio)
tail(exp_tipo_cambio)
range(exp_tipo_cambio$fecha)
summary(exp_tipo_cambio$exp_tcn_12m)


#------------------------------------------------------------------------------#
# 3. Tipo de cambio nominal (TRM diaria -> promedio mensual)
#    Fuente: trm_historico.csv (Banco de la República / Datos Abiertos)
#    Formato crudo: "VALOR" viene como texto "$3,132.42"; fechas dd/mm/yyyy;
#    orden descendente. Colapsamos a PROMEDIO MENSUAL.
#    Convención: pesos por dólar (un aumento = depreciación del peso).
#el TCN se hizo con "promedio de días hábiles"
#------------------------------------------------------------------------------#

trm_raw <- read_csv(
  file.path(ruta_raw, "trm_historico.csv"),
  show_col_types = FALSE
)

tcn <- trm_raw |>
  mutate(
    # limpiar "$" y separador de miles "," -> numérico
    valor = as.numeric(gsub("[$,]", "", VALOR)),
    fecha_dia = dmy(VIGENCIADESDE)
  ) |>
  filter(!is.na(valor), !is.na(fecha_dia)) |>
  mutate(fecha = floor_date(fecha_dia, unit = "month")) |>
  group_by(fecha) |>
  summarise(tcn_nivel = mean(valor, na.rm = TRUE), .groups = "drop") |>
  arrange(fecha) |>
  mutate(tcn = 100 * log(tcn_nivel))


#Chequeos: Deberías ver ~416 obs, desde dic-1991 hasta jul-2026.
# El último mes (jul-2026) es parcial (solo hasta el 31/7), y arranca mucho antes que las demás (1991).
str(tcn)
head(tcn)
tail(tcn)
range(tcn$fecha)
summary(tcn$tcn)


#------------------------------------------------------------------------------#
# 4. IPC - índice desestacionalizado
#    Fuente: ipc.xlsx (DANE vía Banrep), hoja "Datos", índice base 2018=100,
#    "Dato fin de mes".
#    Detalles del archivo crudo:
#      - 2 filas de encabezado (nombres + "dd/mm/aaaa"/"índice") -> skip = 2
#      - decimal con COMA ("159,53") -> reemplazar por punto
#      - filas de basura al final (vacías + "Descargado de sistema...") -> filtrar
#    El índice viene SIN desestacionalizar -> aplicamos X-13ARIMA-SEATS.
#    Orden: (1) leer nivel -> (2) desestacionalizar -> (3) 100*log(.)
#------------------------------------------------------------------------------#

library(seasonal)   # interfaz a X-13ARIMA-SEATS

ipc_raw <- read_excel(
  path = file.path(ruta_raw, "ipc.xlsx"),
  sheet = "Datos",
  skip = 2,
  col_names = FALSE
) |>
  select(fecha = 1, ipc_txt = 2) |>
  # conservar solo filas con fecha dd/mm/aaaa válida (descarta basura final)
  filter(grepl("^\\d{2}/\\d{2}/\\d{4}$", fecha)) |>
  mutate(
    fecha     = dmy(fecha),
    ipc_nivel = as.numeric(gsub(",", ".", ipc_txt))   # coma decimal -> punto
  ) |>
  filter(!is.na(ipc_nivel)) |>
  arrange(fecha) |>
  select(fecha, ipc_nivel)

# Recorte a la ventana de trabajo (colchón desde ene-2003; la EME arranca oct-2003)
ipc_raw <- ipc_raw |> filter(fecha >= as.Date("2003-01-01"))

# Serie ts mensual
anio0  <- year(min(ipc_raw$fecha)); mes0 <- month(min(ipc_raw$fecha))
ipc_ts <- ts(ipc_raw$ipc_nivel, start = c(anio0, mes0), frequency = 12)

# Desestacionalización con X-13ARIMA-SEATS (maneja Semana Santa, días hábiles,
# outliers y elige modelo ARIMA automáticamente).
# Si diera error "x13binary not found": install.packages("x13binary")
ipc_seas <- seas(ipc_ts)
ipc_sa   <- final(ipc_seas)              # serie ajustada estacionalmente

# Transformación del TP sobre la serie desestacionalizada.
# Normalizamos la fecha al primer día del mes (floor_date) para que la clave
# temporal sea homogénea con el resto de las series en el merge.
ipc <- data.frame(
  fecha  = floor_date(ipc_raw$fecha, unit = "month"),
  ipc_sa = as.numeric(ipc_sa),
  ipc    = 100 * log(as.numeric(ipc_sa))
)

# (Opcional) inspección del ajuste:
summary(ipc_seas); plot(ipc_seas)



#------------------------------------------------------------------------------#
# 5. Términos de intercambio de commodities (CTOT - FMI, Gruss & Kebhaj)
#    Fuente: imf_pctot_xm_RW_IX.csv (serie M.CO.xm.R_RW_IX, mensual, ROLLING
#    weights - la apropiada para Colombia por el peso creciente del petróleo).
#    Formato: col 'period' = "YYYY-MM"; col 2 = índice. Transformación 100*log(x).
#------------------------------------------------------------------------------#

tot_raw <- read_csv(
  file.path(ruta_raw, "imf_pctot_xm_RW_IX.csv"),
  show_col_types = FALSE
)

tot <- tot_raw |>
  rename(period = 1, tot_nivel = 2) |>
  mutate(
    fecha = as.Date(paste0(period, "-01")),   # "1980-01" -> 1980-01-01
    tot = 100 * log(tot_nivel)
  ) |>
  filter(!is.na(tot_nivel)) |>
  arrange(fecha) |>
  select(fecha, tot_nivel, tot)

#Chequeos
str(tot)
head(tot)
tail(tot)
range(tot$fecha)
summary(tot$tot)



#------------------------------------------------------------------------------#
# 6. EMBI Colombia (spread soberano)
#    Fuente: Daily_Embi_Panel_wide_selected.xlsx (panel diario, varios países).
#    Usamos la columna 'colombia'. Es DIARIO -> promedio mensual.
#    Queda en NIVELES (puntos básicos), sin log.
#    Nota: la serie de Colombia cubre feb-1997 a mar-2023 -> acota la ventana
#    del panel por el extremo derecho.
#    (Se fuerza lectura numérica: las ~1155 primeras filas son NA y confunden
#     la detección automática de tipo de read_excel.)
#------------------------------------------------------------------------------#

embi_raw <- read_excel(
  path = file.path(ruta_raw, "Daily_Embi_Panel_wide_selected.xlsx"),
  sheet = "Sheet1",
  guess_max = 20000        # mira más filas para tipar bien las columnas
)

embi <- embi_raw |>
  select(fecha_dia = date, embi_col = colombia) |>
  mutate(
    embi_col  = as.numeric(embi_col),        # asegura numérico
    fecha_dia = as.Date(fecha_dia)
  ) |>
  filter(!is.na(embi_col)) |>
  mutate(fecha = floor_date(fecha_dia, unit = "month")) |>
  group_by(fecha) |>
  summarise(embi = mean(embi_col, na.rm = TRUE), .groups = "drop") |>
  arrange(fecha)


#Chequeos: 
summary(embi$embi)
head(embi)
tail(embi)


#------------------------------------------------------------------------------#
# 7. Merge de las series locales en un panel mensual
#    Se unen por 'fecha' (todas ya normalizadas a día 01).
#    inner_join sucesivo -> el panel queda en la VENTANA COMÚN, donde todas las
#    series tienen dato. Binding: EME (oct-2003) por izquierda, EMBI (mar-2023)
#    por derecha. Resultado esperado: oct-2003 a mar-2023 (234 meses, sin huecos).
#------------------------------------------------------------------------------#

datos_locales <- exp_inflacion |>
  select(fecha, exp_inf_12m) |>
  inner_join(exp_tipo_cambio |> select(fecha, exp_tcn_12m), by = "fecha") |>
  inner_join(tcn  |> select(fecha, tcn),  by = "fecha") |>
  inner_join(ipc  |> select(fecha, ipc),  by = "fecha") |>
  inner_join(tot  |> select(fecha, tot),  by = "fecha") |>
  inner_join(embi |> select(fecha, embi), by = "fecha") |>
  arrange(fecha)

# Chequeos del panel local
cat("Panel local:", nrow(datos_locales), "meses, de",
    format(min(datos_locales$fecha)), "a", format(max(datos_locales$fecha)), "\n")

# ¿hay huecos? (meses esperados vs. observados)
meses_esperados <- length(seq(min(datos_locales$fecha),
                              max(datos_locales$fecha), by = "month"))
cat("Meses esperados:", meses_esperados,
    "| observados:", nrow(datos_locales),
    "| huecos:", meses_esperados - nrow(datos_locales), "\n")

# ¿algún NA?
cat("NAs por columna:\n"); print(colSums(is.na(datos_locales)))

head(datos_locales)
tail(datos_locales)






#==============================================================================#
#==============================================================================#
# VARIABLES GLOBALES: Fed shocks (Jarocinski-Karadi), Oil supply (Baumeister-
# Hamilton), Excess Bond Premium (EBP).
# Todas en NIVELES (son shocks/premios ya estandarizados): NO llevan log*100.
# Se leen desde data_raw/ (bases descargadas, no desde URL).
#==============================================================================#
#==============================================================================#

library(readxl)

# --- 1) Fed shock: Jarocinski-Karadi, columna MP_median -----------------------
fed_raw <- read_csv(file.path(ruta_raw, "fed_jk_shocks.csv"),
                    show_col_types = FALSE)
str(fed_raw)   # confirmar: year, month, MP_median

fed_shock <- fed_raw |>
  transmute(
    fecha = make_date(year, month, 1),
    fed_mp = MP_median
  ) |>
  arrange(fecha)

# --- 2) Excess Bond Premium ---------------------------------------------------
ebp_raw <- read_csv(file.path(ruta_raw, "ebp.csv"),
                    show_col_types = FALSE)
str(ebp_raw)   # confirmar nombres: date, ebp (y ojo si date viene "YYYY-MM-DD")

ebp <- ebp_raw |>
  transmute(
    fecha = floor_date(as.Date(date), unit = "month"),
    ebp = ebp
  ) |>
  arrange(fecha)

# --- 3) Oil supply shock: Baumeister-Hamilton --------------------------------
oil_raw <- read_excel(file.path(ruta_raw, "oil_supply_bh.xlsx"))
str(oil_raw)   # estructura irregular: fila 1 basura, cols 1-2 = fecha/valor

oil <- oil_raw[-1, 1:2]
names(oil) <- c("date", "oil_shock")
oil <- oil |>
  mutate(
    oil_shock = as.numeric(oil_shock),
    fecha = floor_date(as.Date(date), unit = "month")
  ) |>
  filter(!is.na(oil_shock)) |>
  select(fecha, oil_shock) |>
  arrange(fecha)


# Chequeos de las 3 globales
cat("Fed (JK):  ", nrow(fed_shock), "obs |",
    format(min(fed_shock$fecha)), "a", format(max(fed_shock$fecha)), "\n")
cat("EBP:       ", nrow(ebp), "obs |",
    format(min(ebp$fecha)), "a", format(max(ebp$fecha)), "\n")
cat("Oil (BH):  ", nrow(oil), "obs |",
    format(min(oil$fecha)), "a", format(max(oil$fecha)), "\n")

summary(fed_shock$fed_mp)
summary(ebp$ebp)
summary(oil$oil_shock)

tail(fed_shock); tail(ebp); tail(oil)



#------------------------------------------------------------------------------#
# 8. Merge final: panel local + variables globales
#    Las globales cubren de sobra la ventana local, así que el binding sigue
#    siendo EME (oct-2003) y EMBI (mar-2023). inner_join recorta a la común.
#    Globales en NIVELES (no llevan log*100).
#------------------------------------------------------------------------------#

datos <- datos_locales |>
  inner_join(fed_shock |> select(fecha, fed_mp),  by = "fecha") |>
  inner_join(ebp       |> select(fecha, ebp),     by = "fecha") |>
  inner_join(oil       |> select(fecha, oil_shock), by = "fecha") |>
  arrange(fecha)

# Chequeos del panel final
cat("Panel final:", nrow(datos), "meses, de",
    format(min(datos$fecha)), "a", format(max(datos$fecha)), "\n")

meses_esp <- length(seq(min(datos$fecha), max(datos$fecha), by = "month"))
cat("Meses esperados:", meses_esp, "| observados:", nrow(datos),
    "| huecos:", meses_esp - nrow(datos), "\n")

cat("NAs por columna:\n"); print(colSums(is.na(datos)))
cat("\nVariables del panel:\n"); print(names(datos))

head(datos)
tail(datos)




#------------------------------------------------------------------------------#
# 9. Objetos de serie de tiempo y guardado
#    Guardamos el data.frame 'datos' (con fecha) y el mts 'datos_ts'.
#    Las series ts individuales NO se guardan: se derivan al vuelo en el
#    análisis con as.list(datos_ts), evitando duplicar la fuente de verdad.
#------------------------------------------------------------------------------#

datos_ts <- ts(
  datos |> dplyr::select(-fecha),
  start     = c(2003, 10),   # panel arranca oct-2003
  frequency = 12
)

# Chequeos
stopifnot(identical(dim(datos_ts)[1], nrow(datos)))
cat("datos_ts:", paste(start(datos_ts), collapse = "-"), "a",
    paste(end(datos_ts), collapse = "-"),
    "| vars:", paste(colnames(datos_ts), collapse = ", "), "\n")

# Guardado (una sola vez; el análisis lee de acá)
saveRDS(datos,    file.path(ruta_out, "panel_final.rds"))
saveRDS(datos_ts, file.path(ruta_out, "datos_ts.rds"))


