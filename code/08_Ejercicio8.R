#------------------------------------------------------------------------------#
# TP2 Macroeconometría 2026 — Traspaso cambiario (ERPT), Colombia
# INCISO 8: IRF de expectativas por Local Projections
#
# Qué hace este script:
#   - Reconstruye el shock estructural al TCN (û_t) del VAR del inciso 2/3.
#   - Estima 3 objetos por LP (mismo shock, mismos controles que inciso 3),
#     cambiando SOLO la variable dependiente (LHS):
#       Paso 1: IRF de la expectativa de inflación a 12m en nivel  (π^{e,12}_{t+h})
#       Paso 2: IRF de la inflación realizada 12m forward          (IPC_{t+h+12}-IPC_{t+h})
#       Paso 3: IRF de la DIFERENCIA (test)                        (paso1 - paso2)
#     y sus análogos para el TCN (depreciación esperada vs realizada forward).
#
# Qué objetos consume:
#   - datos_ts.rds  (mts, start=2003.10, freq=12): series en niveles/log×100
#   - objetos_02.rds: SVAR, VAR, Y (diferencias), p, m, H, gamma
#
# Muestra efectiva (forma B): se recortan las últimas 12 obs. para que el
# forward realizado exista en todas las fechas y los 3 pasos sean comparables.
#------------------------------------------------------------------------------#

remove(list = ls()); gc()

library(here)
library(vars)
library(sandwich)
library(ggplot2)

#---------------------------------------------------------------------------#
# Carga de datos y objetos previos
#---------------------------------------------------------------------------#
datos_ts <- readRDS(here("data", "datos_ts.rds"))
for (v in colnames(datos_ts)) assign(v, datos_ts[, v])

obj <- readRDS(here("data", "objetos_02.rds"))
for (nm in names(obj)) assign(nm, obj[[nm]])

obj3 <- readRDS(here("data", "objetos_03.rds"))
u.tcn <- obj3$u.tcn
Y     <- obj3$Y

source(here("tools", "PS3_LP_Tools.R"))

#---------------------------------------------------------------------------#
# 1. Paso 1: IRF de la expectativa de inflación a 12m (nivel)
#---------------------------------------------------------------------------#

# --- 1.a Insumos heredados del inciso 3 (NO se recalculan) -------------------
# u.tcn : shock estructural al TCN (vector), generado en el script 03.
# Y     : matriz de diferencias (dtcn, dipc) del VAR del inciso 2.
# Ambos se cargan vía objetos_03.rds en la Sección 0.


# --- 1.b Alineación de las expectativas con las filas de Y -------------------
# Y (diferencias) tiene 233 filas y arranca UNA obs. después que datos_ts
# (perdió la primera fila al diferenciar). Las expectativas viven en datos_ts
# (234 filas, en nivel). Para que "fila t de Y" y "fila t del LHS" sean la
# MISMA fecha, descartamos la primera fila de la serie de nivel.
einf <- as.numeric(exp_inf_12m)[-1]   # π^{e,12} alineada a Y  -> largo 233
etcn <- as.numeric(exp_tcn_12m)[-1]   # Δ^{e,12} alineada a Y  -> largo 233

Tfull <- nrow(Y)                      # 233 (largo antes de recortar)


# --- 1.c Muestra efectiva (forma B): recorte común de 12 obs. ----------------
# POR QUÉ acá: el paso 2 usará realizado "12 meses hacia adelante", que no
# existe para las últimas 12 fechas de la muestra. Para que los 3 pasos se
# estimen sobre las MISMAS observaciones (única forma de que el test del paso 3
# sea limpio), definimos AHORA una ventana común y la usamos en todo el inciso.
h12  <- 12                            # horizonte del forward (fijo por la encuesta)
keep <- 1:(Tfull - h12)              # filas 1:221 de Y  (descartamos las 12 últimas)

# Aplicamos el recorte a los controles (Y) y al shock:
Y.e <- Y[keep, ]                      # controles: mismas filas para todos los pasos

# El shock u.tcn (largo 232) arranca en la fila 2 de Y. Al recortar Y a las
# filas 1:221, el shock correspondiente son las filas 2:221 de Y, que en el
# índice de u.tcn (que empieza en la fila 2) son las posiciones 1:220.
u.e <- u.tcn[1:(length(keep) - 1)]    # shock alineado a Y.e -> largo 220


# --- 1.d LHS del paso 1: expectativa en nivel, recortada ---------------------
einf.e <- einf[keep]                  # π^{e,12} sobre la muestra efectiva (221)
etcn.e <- etcn[keep]                  # Δ^{e,12} sobre la muestra efectiva (221)


# --- 1.e Estimación LP del paso 1 --------------------------------------------
# Mismo shock, mismos controles (rezagos de Y.e vía get.controls dentro de la
# función). Solo cambia el LHS: la expectativa en nivel.
lp.einf <- lp.shock.y(Y.e, u.e, einf.e, p = p, H = H, gamma = gamma)  # inflación esp.
lp.etcn <- lp.shock.y(Y.e, u.e, etcn.e, p = p, H = H, gamma = gamma)  # deprec. esp.


# --- VERIFICACIÓN (correr y pegar output) ---
cat("nrow(Y.e):", nrow(Y.e), "\n")
cat("length(u.e):", length(u.e), "\n")
cat("length(einf.e):", length(einf.e), "\n")
cat("length(etcn.e):", length(etcn.e), "\n")
cat("length(lp.einf$pe):", length(lp.einf$pe), "\n")


#---------------------------------------------------------------------------#
# 2. Paso 2: IRF del realizado forward 12m (comparable con el paso 1)
#---------------------------------------------------------------------------#

# --- 2.a Niveles alineados a Y -----------------------------------------------
# Igual que las expectativas: IPC y TCN están en nivel en datos_ts (234 filas)
# y hay que descartar la primera para alinearlos con las filas de Y.
ipc.lvl <- as.numeric(ipc)[-1]        # nivel IPC alineado a Y -> largo 233
tcn.lvl <- as.numeric(tcn)[-1]        # nivel TCN alineado a Y -> largo 233

# --- 2.b LHS forward 12m: realizado de los próximos 12 meses -----------------
# π^{real,12}_t = IPC_{t+12} - IPC_t. Para cada fila t tomamos el IPC 12 filas
# más adelante y le restamos el actual. Las últimas 12 filas no tienen "t+12",
# por eso las completamos con NA... pero como la muestra efectiva (keep=1:221)
# YA excluye esas 12, al recortar no queda ningún NA.
ipc.fwd <- c(ipc.lvl[(h12 + 1):Tfull], rep(NA, h12)) - ipc.lvl   # largo 233
tcn.fwd <- c(tcn.lvl[(h12 + 1):Tfull], rep(NA, h12)) - tcn.lvl   # largo 233

# Recorte a la muestra efectiva (las mismas filas que el paso 1):
ipc.fwd.e <- ipc.fwd[keep]            # realizado inflación forward -> 221, sin NA
tcn.fwd.e <- tcn.fwd[keep]            # realizado depreciación forward -> 221, sin NA

# --- 2.c Estimación LP del paso 2 --------------------------------------------
# Mismo shock (u.e), mismos controles (Y.e). Solo cambia el LHS.
lp.ipc.fwd <- lp.shock.y(Y.e, u.e, ipc.fwd.e, p = p, H = H, gamma = gamma)  # inflación real.
lp.tcn.fwd <- lp.shock.y(Y.e, u.e, tcn.fwd.e, p = p, H = H, gamma = gamma)  # deprec. real.


# --- VERIFICACIÓN (correr y pegar output) ---
cat("length(ipc.fwd.e):", length(ipc.fwd.e), " NAs:", sum(is.na(ipc.fwd.e)), "\n")
cat("length(tcn.fwd.e):", length(tcn.fwd.e), " NAs:", sum(is.na(tcn.fwd.e)), "\n")
cat("length(lp.ipc.fwd$pe):", length(lp.ipc.fwd$pe), "\n")
# Chequeo manual del forward en la fila 1: IPC[13]-IPC[1] debe igualar ipc.fwd.e[1]
cat("check fila 1:", ipc.lvl[13] - ipc.lvl[1], "vs", ipc.fwd.e[1], "\n")




#---------------------------------------------------------------------------#
# 3. Paso 3: IRF de la DIFERENCIA (expectativa − realizado) + test
#---------------------------------------------------------------------------#

# --- 3.a Variable diferencia (LHS del test) ----------------------------------
# D_t = expectativa − realizado forward, sobre la MISMA muestra efectiva.
# Como einf.e y ipc.fwd.e ya están recortados a keep (221) y alineados, la
# resta es directa, posición por posición.
dif.inf.e <- einf.e - ipc.fwd.e       # diferencia inflación  -> 221
dif.tcn.e <- etcn.e - tcn.fwd.e       # diferencia depreciación -> 221

# --- 3.b Estimación LP de la diferencia --------------------------------------
# El coeficiente de cada horizonte ES β^e_h − β^real_h (la diferencia de IRFs).
# La banda (HAC, vía vcovHAC dentro de la función) es el test: si excluye 0,
# la diferencia es significativa en ese horizonte.
lp.dif.inf <- lp.shock.y(Y.e, u.e, dif.inf.e, p = p, H = H, gamma = gamma)
lp.dif.tcn <- lp.shock.y(Y.e, u.e, dif.tcn.e, p = p, H = H, gamma = gamma)



# --- VERIFICACIÓN (correr y pegar output) ---
cat("length(dif.inf.e):", length(dif.inf.e), " NAs:", sum(is.na(dif.inf.e)), "\n")
cat("length(lp.dif.inf$pe):", length(lp.dif.inf$pe), "\n")
# Identidad de linealidad: la IRF de la diferencia debe igualar la diferencia
# de IRFs (punto a punto). Si esto da ~0, el método está bien planteado.
cat("max |dif − (paso1 − paso2)|:",
    max(abs(lp.dif.inf$pe - (lp.einf$pe - lp.ipc.fwd$pe))), "\n")
# ¿En qué horizontes la diferencia es significativa? (banda excluye 0)
sig <- which(lp.dif.inf$lb > 0 | lp.dif.inf$ub < 0) - 1  # horizontes 0..24
cat("horizontes con diferencia significativa (inflación):", sig, "\n")





#---------------------------------------------------------------------------#
# 4. Gráficos (ggplot2, mismo formato que el inciso 3)
#---------------------------------------------------------------------------#

# Paletas nombradas (líneas y bandas). Serie "Esperada" en azul pizarra,
# "Realizada" en terracota, "Diferencia" en verde musgo.
colores_serie <- c("Esperada" = "#2C5F8A", "Realizada" = "#C2410C",
                   "Diferencia" = "#4D7C5A")
colores_banda <- colores_serie   # mismas tonalidades para el ribbon



# --- 4.a Figura A: inflación (paneles separados: esperada | realizada) --------
df_inf <- bind_rows(
  data.frame(horizonte = 0:H, panel = "Inflación esperada (12m)",
             pe = lp.einf$pe,    lb = lp.einf$lb,    ub = lp.einf$ub,
             serie = "Esperada"),
  data.frame(horizonte = 0:H, panel = "Inflación realizada (fwd 12m)",
             pe = lp.ipc.fwd$pe, lb = lp.ipc.fwd$lb, ub = lp.ipc.fwd$ub,
             serie = "Realizada")
)
df_inf$panel <- factor(df_inf$panel,
                       levels = c("Inflación esperada (12m)",
                                  "Inflación realizada (fwd 12m)"))

p_inf <- ggplot(df_inf, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = serie), alpha = 0.35) +
  geom_line(aes(y = pe, color = serie), linewidth = 0.9) +
  scale_color_manual(values = colores_serie, name = NULL) +
  scale_fill_manual(values = colores_banda, name = NULL) +
  labs(x = "Meses", y = "Respuesta (log × 100)") +
  facet_wrap(~ panel) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
p_inf

# --- 4.b Figura B: TCN (paneles separados: esperada | realizada) --------------
df_tcn <- bind_rows(
  data.frame(horizonte = 0:H, panel = "Depreciación esperada (12m)",
             pe = lp.etcn$pe,    lb = lp.etcn$lb,    ub = lp.etcn$ub,
             serie = "Esperada"),
  data.frame(horizonte = 0:H, panel = "Depreciación realizada (fwd 12m)",
             pe = lp.tcn.fwd$pe, lb = lp.tcn.fwd$lb, ub = lp.tcn.fwd$ub,
             serie = "Realizada")
)
df_tcn$panel <- factor(df_tcn$panel,
                       levels = c("Depreciación esperada (12m)",
                                  "Depreciación realizada (fwd 12m)"))

p_tcn <- ggplot(df_tcn, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = serie), alpha = 0.35) +
  geom_line(aes(y = pe, color = serie), linewidth = 0.9) +
  scale_color_manual(values = colores_serie, name = NULL) +
  scale_fill_manual(values = colores_banda, name = NULL) +
  labs(x = "Meses", y = "Respuesta (log × 100)") +
  facet_wrap(~ panel) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
p_tcn

# --- 4.c Figura C: diferencias (test) ----------------------------------------
df_dif <- bind_rows(
  data.frame(horizonte = 0:H, panel = "Diferencia — Inflación",
             pe = lp.dif.inf$pe, lb = lp.dif.inf$lb, ub = lp.dif.inf$ub),
  data.frame(horizonte = 0:H, panel = "Diferencia — TCN",
             pe = lp.dif.tcn$pe, lb = lp.dif.tcn$lb, ub = lp.dif.tcn$ub)
)
df_dif$serie <- "Diferencia"

p_dif <- ggplot(df_dif, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = serie), alpha = 0.35) +
  geom_line(aes(y = pe, color = serie), linewidth = 0.9) +
  scale_color_manual(values = colores_serie, name = NULL) +
  scale_fill_manual(values = colores_banda, name = NULL) +
  labs(x = "Meses", y = "Diferencia de respuestas (log × 100)") +
  facet_wrap(~ panel) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
p_dif

# --- 4.d Guardado -------------------------------------------------------------
ggsave(file.path("output", "08_irf_inflacion.pdf"),   p_inf, width = 6, height = 4, units = "in")
ggsave(file.path("output", "08_irf_tcn.pdf"),         p_tcn, width = 6, height = 4, units = "in")
ggsave(file.path("output", "08_irf_diferencias.pdf"), p_dif, width = 6, height = 4, units = "in")




#---------------------------------------------------------------------------#
# Guardado de objetos para el inciso 9
#---------------------------------------------------------------------------#
saveRDS(list(
  lp.einf = lp.einf, lp.ipc.fwd = lp.ipc.fwd,
  lp.etcn = lp.etcn, lp.tcn.fwd = lp.tcn.fwd,
  lp.dif.inf = lp.dif.inf, lp.dif.tcn = lp.dif.tcn,
  einf.e = einf.e, etcn.e = etcn.e,
  ipc.fwd.e = ipc.fwd.e, tcn.fwd.e = tcn.fwd.e,
  Y.e = Y.e, u.e = u.e,
  keep = keep, h12 = h12
), here("data", "objetos_04.rds"))



