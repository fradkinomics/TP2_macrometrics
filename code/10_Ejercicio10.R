#------------------------------------------------------------------------------#
# Universidad de San Andrés - Macroeconometría 2026
# TP2: Exchange Rate Pass-Through - Colombia
# INCISO 10 [Extra]: ERPT esperado vs. estimado (realizado forward 12m) para
# cada uno de los 6 shocks estructurales del inciso 6. Generalización del
# inciso 9 (par esperado/realizado, sin cumsum) a los 6 shocks del VAR7.
# Controles: sistema de 7 variables (Opción 1). Shocks e insumos se cargan
# de objetos guardados; no se reestima nada.
#------------------------------------------------------------------------------#

remove(list = ls(all.names = TRUE))
gc()

library(sandwich)
library(here)

#---( Insumos guardados )---#
# Del inciso 6: los 6 shocks estructurales + el sistema Y (7 var) + p
objetos_06 <- readRDS(here("data", "objetos_06.rds"))
struct.shocks <- objetos_06$struct.shocks   # matriz (T_var x 7), shocks del VAR7
Y             <- objetos_06$Y                # mts de 7 columnas (controles)
p             <- objetos_06$p                # orden de rezagos del VAR7 (=2)

# Datos en niveles para reconstruir las LHS forward y esperadas sobre las
# fechas del VAR7 (ancladas por ts, no por posición)
datos_ts <- readRDS(here("data", "datos_ts.rds"))

#---( Tools de la cátedra )---#
source(here("tools", "PS3_LP_Tools.R"))  # lp.multiplier.shock.y, lp.shock.y, get.controls, first.stage.F

#---( Verificación de carga: fechas y largos de los insumos )---#
cat("--- struct.shocks ---\n")
cat("dim:", dim(struct.shocks)[1], "x", dim(struct.shocks)[2], "\n")
cat("colnames:", paste(colnames(struct.shocks), collapse = ", "), "\n\n")

cat("--- Y (controles VAR7) ---\n")
cat("dim:", nrow(Y), "x", ncol(Y), "| p =", p, "\n")
cat("start:", paste(start(Y), collapse = ":"),
    "| end:", paste(end(Y), collapse = ":"),
    "| freq:", frequency(Y), "\n\n")

cat("--- datos_ts ---\n")
cat("colnames:", paste(colnames(datos_ts), collapse = ", "), "\n")
cat("start:", paste(start(datos_ts), collapse = ":"),
    "| end:", paste(end(datos_ts), collapse = ":"), "\n")




#---------------------------------------------------------------------------#
# BLOQUE 2: Reconstrucción de las 4 LHS + alineación del shock + ventana común
#---------------------------------------------------------------------------#
# Reloj de referencia: las fechas de Y (VAR7): 2003:11 .. 2023:03, 233 filas.
# - struct.shocks empieza 2 filas después (2004:01) por los p=2 rezagos del VAR.
# - las LHS forward pierden las últimas 12 filas por el forward a 12m.
# Ventana común contigua resultante: 2004:01 .. 2022:03.

#---( 2.a Niveles/expectativas alineados a Y, como en el inciso 8 )---#
# datos_ts tiene 234 filas (2003:10); Y tiene 233 (2003:11). Descartamos la
# primera fila del nivel para que "fila t de Y" y "fila t del LHS" sean la
# misma fecha (idéntico criterio que el script 08).
einf    <- as.numeric(datos_ts[, "exp_inf_12m"])[-1]  # π^{e,12}, largo 233
etcn    <- as.numeric(datos_ts[, "exp_tcn_12m"])[-1]  # Δ^{e,12}, largo 233
ipc.lvl <- as.numeric(datos_ts[, "ipc"])[-1]          # nivel IPC, largo 233
tcn.lvl <- as.numeric(datos_ts[, "tcn"])[-1]          # nivel TCN, largo 233

Tfull <- nrow(Y)   # 233
h12   <- 12

#---( 2.b LHS forward 12m: realizado de los próximos 12 meses )---#
# π^{real,12}_t = IPC_{t+12} - IPC_t ; Δ^{real,12}_t = TCN_{t+12} - TCN_t.
# Últimas 12 filas sin "t+12" -> NA (se descartan luego en la ventana común).
ipc.fwd <- c(ipc.lvl[(h12 + 1):Tfull], rep(NA, h12)) - ipc.lvl   # largo 233
tcn.fwd <- c(tcn.lvl[(h12 + 1):Tfull], rep(NA, h12)) - tcn.lvl   # largo 233

#---( 2.c Convertir todo a ts sobre el reloj de Y, para alinear por FECHA )---#
sY  <- start(Y); fY <- frequency(Y)
einf.ts    <- ts(einf,    start = sY, frequency = fY)
etcn.ts    <- ts(etcn,    start = sY, frequency = fY)
ipc.fwd.ts <- ts(ipc.fwd, start = sY, frequency = fY)
tcn.fwd.ts <- ts(tcn.fwd, start = sY, frequency = fY)

# struct.shocks arranca en 2004:01 (2 filas después de Y). Lo anclamos a su
# fecha real: fila 1 de struct.shocks = start(Y) + p meses.
s.shock <- c(sY[1] + (sY[2] - 1 + p) %/% fY, (sY[2] - 1 + p) %% fY + 1)
shocks.ts <- ts(struct.shocks, start = s.shock, frequency = fY)  # 231 x 7, 2004:01..

#---( 2.d Ventana común contigua: intersección temporal )---#
# Inicio: max(inicio shock, inicio Y) = 2004:01 (manda el shock).
# Fin: última fecha con forward definido = end(Y) - 12 meses = 2022:03.
win.start <- start(shocks.ts)                                  # 2004:01
win.end   <- c(end(Y)[1] - (h12 - end(Y)[2]) %/% fY - (end(Y)[2] <= h12),
               (end(Y)[2] - h12 - 1) %% fY + 1)                # 2022:03

# Recortamos TODO a esa ventana (mismas fechas exactas en cada objeto):
Y.c        <- window(Y,          start = win.start, end = win.end)
shocks.c   <- window(shocks.ts,  start = win.start, end = win.end)
einf.c     <- window(einf.ts,    start = win.start, end = win.end)
etcn.c     <- window(etcn.ts,    start = win.start, end = win.end)
ipc.fwd.c  <- window(ipc.fwd.ts, start = win.start, end = win.end)
tcn.fwd.c  <- window(tcn.fwd.ts, start = win.start, end = win.end)

#---( 2.e VERIFICACIÓN: fechas, largos y ausencia de NA (correr y pegar) )---#
cat("--- Ventana común ---\n")
cat("win.start:", paste(win.start, collapse = ":"),
    "| win.end:", paste(win.end, collapse = ":"), "\n\n")

chk <- function(x, nm) cat(sprintf("%-12s start=%s end=%s  n=%d  NAs=%d\n",
                                   nm, paste(start(x), collapse=":"), paste(end(x), collapse=":"),
                                   NROW(x), sum(is.na(as.matrix(x)))))
chk(Y.c, "Y.c"); chk(shocks.c, "shocks.c")
chk(einf.c, "einf.c"); chk(etcn.c, "etcn.c")
chk(ipc.fwd.c, "ipc.fwd.c"); chk(tcn.fwd.c, "tcn.fwd.c")



#---------------------------------------------------------------------------#
# BLOQUE 3: Loop sobre los 6 shocks -> par (ERPT esperado, ERPT realizado)
#---------------------------------------------------------------------------#
# Para cada shock s del VAR7 se estiman DOS cocientes por LP-IV con
# lp.multiplier.shock.y (misma función del inciso 9, SIN cumsum):
#   ERPT esperado_s  = IRF(einf | s) / IRF(etcn | s)      [encuesta, nivel 12m]
#   ERPT realizado_s = IRF(ipc.fwd | s) / IRF(tcn.fwd | s) [forward 12m realizado]
# El shock instrumenta el denominador; el numerador se proyecta sobre la
# 1ra etapa. Controles: Y.c (sistema VAR7, Opción 1).

H     <- 36
gamma <- 0.95

shocks.nombres <- c("Shock Fed (JK)", "Shock Petróleo (BH)", "Excess Bond Premium",
                    "Términos de Intercambio", "EMBI", "TCN")

erpt.esp10  <- vector("list", 6)   # ERPT esperado  por shock
erpt.real10 <- vector("list", 6)   # ERPT realizado por shock
names(erpt.esp10)  <- shocks.nombres
names(erpt.real10) <- shocks.nombres

for (s in 1:6) {
  sh <- shocks.c[, s]   # shock estructural s, alineado a Y.c
  
  # ERPT esperado: numerador inflación esperada, denominador deprec. esperada
  erpt.esp10[[s]] <- lp.multiplier.shock.y(
    Y = Y.c, shock = sh,
    y.num = einf.c, y.den = etcn.c,
    p = p, H = H, gamma = gamma)
  
  # ERPT realizado: numerador inflación fwd, denominador deprec. fwd
  erpt.real10[[s]] <- lp.multiplier.shock.y(
    Y = Y.c, shock = sh,
    y.num = ipc.fwd.c, y.den = tcn.fwd.c,
    p = p, H = H, gamma = gamma)
}

#---( VERIFICACIÓN: que los 12 objetos existan y tengan largo H+1 )---#
cat("--- Chequeo de estructura de resultados ---\n")
cat("Largo esperado por serie (H+1):", H + 1, "\n")
for (s in 1:6) {
  cat(sprintf("%-24s esp: pe[0]=%7.2f pe[36]=%7.2f | real: pe[0]=%8.2f pe[36]=%8.2f\n",
              shocks.nombres[s],
              erpt.esp10[[s]]$pe[1],  erpt.esp10[[s]]$pe[H + 1],
              erpt.real10[[s]]$pe[1], erpt.real10[[s]]$pe[H + 1]))
}



#---------------------------------------------------------------------------#
# BLOQUE 4: Diagnóstico de instrumento débil (F de primera etapa)
#---------------------------------------------------------------------------#
# Para cada shock, el F de regresar el DENOMINADOR (limpiado de controles)
# sobre el shock, por horizonte. Dos denominadores por shock:
#   esperado  -> etcn.c   (deprec. esperada 12m)
#   realizado -> tcn.fwd.c (deprec. realizada forward 12m)
# Regla de lectura (Staiger-Stock): F < 10 = instrumento débil.
# first.stage.F NO está en las tools: se define en el script 09. La replicamos.

first.stage.F <- function(Y, shock, y.den, p, H) {
  m <- ncol(Y); T <- nrow(Y)
  Wt.full <- get.controls(Y, p, m)
  rsp.r <- y.den[(p + 1):T]
  Fh <- rep(NA, H + 1)
  for (h in 0:H) {
    yh.r <- rsp.r[(h + 1):(T - p)]
    Wt   <- Wt.full[1:(T - p - h), ]
    st   <- shock[1:(T - p - h)]
    yh.r <- resid(lm(yh.r ~ -1 + Wt))
    fit  <- lm(yh.r ~ -1 + st)
    Fh[h + 1] <- summary(fit)$fstatistic[1]
  }
  Fh
}

#---( F por shock, para los dos denominadores )---#
F.esp10  <- vector("list", 6)   # F del denominador esperado  por shock
F.real10 <- vector("list", 6)   # F del denominador realizado por shock
names(F.esp10)  <- shocks.nombres
names(F.real10) <- shocks.nombres

for (s in 1:6) {
  sh <- shocks.c[, s]
  F.esp10[[s]]  <- first.stage.F(Y.c, sh, etcn.c,   p, H)
  F.real10[[s]] <- first.stage.F(Y.c, sh, tcn.fwd.c, p, H)
}

#---( VERIFICACIÓN: mediana del F por shock y denominador )---#
cat("--- F de primera etapa (mediana sobre h=0..36) ---\n")
cat(sprintf("%-24s %8s %10s\n", "Shock", "F esp", "F real"))
for (s in 1:6) {
  cat(sprintf("%-24s %8.2f %10.2f\n",
              shocks.nombres[s],
              median(F.esp10[[s]],  na.rm = TRUE),
              median(F.real10[[s]], na.rm = TRUE)))
}
cat("\nUmbral Staiger-Stock: F > 10 para instrumento fuerte.\n")





#---------------------------------------------------------------------------#
# BLOQUE 5: Figura principal - 6 paneles, ERPT esperado vs. realizado por shock
#---------------------------------------------------------------------------#
# Un panel por shock; en cada uno se superponen las dos definiciones
# (esperado / realizado forward) con sus bandas. Zoom por panel (mediana +/-
# 4*IQR combinando ambas series), sin descartar datos. Estilo del inciso 7.

library(ggplot2)
library(dplyr)
library(patchwork)

h <- 0:H

# Data frame largo: dos series (Esperado/Realizado) por shock
df_erpt10 <- lapply(1:6, function(s) {
  bind_rows(
    data.frame(horizonte = h, shock = shocks.nombres[s], tipo = "Esperado",
               pe = erpt.esp10[[s]]$pe,
               lb = erpt.esp10[[s]]$lb,
               ub = erpt.esp10[[s]]$ub),
    data.frame(horizonte = h, shock = shocks.nombres[s], tipo = "Realizado",
               pe = erpt.real10[[s]]$pe,
               lb = erpt.real10[[s]]$lb,
               ub = erpt.real10[[s]]$ub)
  )
}) %>% bind_rows()

df_erpt10$shock <- factor(df_erpt10$shock, levels = shocks.nombres)
df_erpt10$tipo  <- factor(df_erpt10$tipo,  levels = c("Esperado", "Realizado"))

colores_tipo       <- c("Esperado" = "steelblue4", "Realizado" = "darkorange")
colores_banda_tipo <- c("Esperado" = "lightsteelblue", "Realizado" = "navajowhite2")

# Zoom por panel: mediana +/- 4*IQR combinando ambas series (no elimina datos)
limites10 <- df_erpt10 %>%
  group_by(shock) %>%
  summarise(
    centro = median(pe, na.rm = TRUE),
    escala = IQR(c(pe, lb, ub), na.rm = TRUE),
    ymin = centro - 4 * escala,
    ymax = centro + 4 * escala,
    .groups = "drop"
  )

hacer_panel10 <- function(shock.nombre) {
  lims <- limites10 %>% dplyr::filter(shock == shock.nombre)
  
  ggplot(df_erpt10 %>% dplyr::filter(shock == shock.nombre), aes(x = horizonte)) +
    geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = tipo), alpha = 0.25) +
    geom_line(aes(y = pe, color = tipo), linewidth = 0.9) +
    scale_color_manual(values = colores_tipo, name = NULL) +
    scale_fill_manual(values = colores_banda_tipo, name = NULL) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax)) +
    labs(x = NULL, y = NULL, title = shock.nombre) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          panel.grid.minor = element_blank(),
          legend.position = "none")
}

paneles10 <- lapply(shocks.nombres, hacer_panel10)

# En lugar de extraer la leyenda con cowplot, dejamos que patchwork la
# recolecte y la ubique en un area propia (guide_area) al pie del layout.
hacer_panel10_leg <- function(shock.nombre) {
  hacer_panel10(shock.nombre) + theme(legend.position = "bottom")
}
paneles10 <- lapply(shocks.nombres, hacer_panel10_leg)

p_erpt10 <- wrap_plots(paneles10, nrow = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = "Meses (eje x) | ERPT condicional (%) (eje y)") &
  theme(legend.position = "bottom")

p_erpt10

p_erpt10

if (!dir.exists(here("output"))) dir.create(here("output"))
ggsave(here("output", "erpt_punto10_esp_vs_real.pdf"), p_erpt10,
       width = 10, height = 6.5, units = "in")



#---------------------------------------------------------------------------#
# BLOQUE 6: Tabla resumen - ERPT esperado vs. realizado (h=36) + F de 1ra etapa
#---------------------------------------------------------------------------#
# Resultado central: los F (mediana sobre h) muestran que NINGUN shock
# identifica ninguno de los dos cocientes. Los point estimates a h=36 se
# reportan solo para ilustrar la inestabilidad, no como magnitudes leibles.

tabla.erpt10 <- data.frame(
  Shock       = shocks.nombres,
  ESP_pe36    = sapply(erpt.esp10,  function(x) round(x$pe[H + 1], 2)),
  ESP_lb36    = sapply(erpt.esp10,  function(x) round(x$lb[H + 1], 2)),
  ESP_ub36    = sapply(erpt.esp10,  function(x) round(x$ub[H + 1], 2)),
  ESP_Fmed    = sapply(F.esp10,     function(x) round(median(x, na.rm = TRUE), 2)),
  REAL_pe36   = sapply(erpt.real10, function(x) round(x$pe[H + 1], 2)),
  REAL_lb36   = sapply(erpt.real10, function(x) round(x$lb[H + 1], 2)),
  REAL_ub36   = sapply(erpt.real10, function(x) round(x$ub[H + 1], 2)),
  REAL_Fmed   = sapply(F.real10,    function(x) round(median(x, na.rm = TRUE), 2))
)
print(tabla.erpt10, row.names = FALSE)

write.csv(tabla.erpt10, here("output", "tabla_erpt_punto10.csv"), row.names = FALSE)