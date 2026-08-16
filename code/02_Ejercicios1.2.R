#------------------------------------------------------------------------------#
# Universidad de San Andrés - Macroeconometría 2026
# TP2: Exchange Rate Pass-Through - Colombia
#------------------------------------------------------------------------------#

#===============================================================================
#===============================INCISO 1========================================
#===============================================================================
#   (a) Gráficos de dos ejes verticales: TCN (eje izq., en 100*log) vs. cada
#       una de las demás series ya transformadas (eje der.).
#   (b) Tests de raíz unitaria (ADF y KPSS) en niveles y en diferencias para
#       cada serie, como insumo para decidir niveles vs. diferencias en el VAR.
#-------------------------------------------------------------------------------

remove(list = ls(all.names = TRUE))
gc()


library(here); library(ggplot2); library(urca)

datos_ts <- readRDS(here("data", "datos_ts.rds"))   # el mts, lectura instantánea
datos     <- readRDS(here("data", "panel_final.rds")) # df con 'fecha', para ggplot

# Series individuales al vuelo (para tests de raíz unitaria univariados):
for (v in colnames(datos_ts)) assign(v, datos_ts[, v], envir = environment())
# ahora existen tcn, ipc, tot, embi, exp_inf_12m, exp_tcn_12m, fed_mp, ebp, oil_shock



#Parte de Felicitas
# Carpeta de salida para gráficos y tablas (se crea sola si no existe)
if (!dir.exists(here("output"))) dir.create(here("output"))

#------------------------------------------------------------------------------#
# 1.a Gráficos de dos ejes verticales
#------------------------------------------------------------------------------#
# Todas las series ya comparten el mismo índice temporal (vienen del mismo
# panel mergeado en 01_Procesamiento_Datos.R), así que no hace falta alinear
# ventanas: se grafican directamente.

plot.two.axis <- function(tcn, y, name.y, main.title,
                          col.tcn = "black", col.y = "steelblue4") {
  
  # Margen inferior más grande para dejarle lugar a la leyenda por fuera
  # del cuadro del gráfico (debajo de las etiquetas del eje X).
  par(mar = c(6.5, 4.5, 3, 4.5), xpd = FALSE)
  
  plot(tcn, type = "l", col = col.tcn, lwd = 2,
       xlab = "", ylab = "TCN (100*log)", main = main.title,
       cex.main = 0.95, cex.lab = 0.85)
  
  grid(NULL, NULL, lty = 3, col = "gray80")
  lines(tcn, col = col.tcn, lwd = 2) # redibujar arriba de la grilla
  
  par(new = TRUE)
  
  plot(y, type = "l", col = col.y, lwd = 2, lty = 2,
       axes = FALSE, xlab = "", ylab = "")
  axis(4, col = col.y, col.axis = col.y)
  mtext(name.y, side = 4, line = 3, col = col.y, cex = 0.85)
  
  # Leyenda AFUERA del área de datos: xpd = TRUE habilita dibujar fuera de
  # la región del gráfico, e inset con y negativo la empuja por debajo del
  # eje X (no tapa ninguna serie, sea cual sea su forma).
  par(xpd = TRUE)
  legend("bottom", inset = c(0, -0.32), legend = c("TCN", name.y),
         col = c(col.tcn, col.y), lty = c(1, 2), lwd = 2,
         horiz = TRUE, bty = "n", cex = 0.8, seg.len = 2.5)
  par(xpd = FALSE)
  
  par(mar = c(5, 4, 4, 2) + 0.1) # reset margins por las dudas
}

# Lista de series a graficar contra el TCN (nombre para el eje derecho)
series.plot <- list(
  IPC     = list(serie = ipc,         label = "IPC (100*log, SA)"),
  EXP_INF = list(serie = exp_inf_12m, label = "Exp. Inflación 12m (100*log(1+x))"),
  EXP_TCN = list(serie = exp_tcn_12m, label = "Exp. TCN 12m (100*log)"),
  TOT     = list(serie = tot,         label = "Términos de Intercambio (100*log)"),
  EMBI    = list(serie = embi,        label = "EMBI Colombia (pb, nivel)"),
  FED     = list(serie = fed_mp,      label = "Shock Fed (JK, MP median)"),
  EBP     = list(serie = ebp,         label = "Excess Bond Premium"),
  OIL     = list(serie = oil_shock,   label = "Shock Oferta Petróleo (BH)")
)

# En pantalla (una por vez)
for (nm in names(series.plot)) {
  plot.two.axis(tcn, series.plot[[nm]]$serie, series.plot[[nm]]$label,
                paste("TCN vs.", series.plot[[nm]]$label))
}

# Guardado en PDF (una página por gráfico) para el informe
pdf(here("output", "graficos_dos_ejes.pdf"), width = 8, height = 5.5)
for (nm in names(series.plot)) {
  plot.two.axis(tcn, series.plot[[nm]]$serie, series.plot[[nm]]$label,
                paste("TCN vs.", series.plot[[nm]]$label))
}
dev.off()



#------------------------------------------------------------------------------#
# 1.b Tests de raíz unitaria (ADF), en niveles y diferencias
#------------------------------------------------------------------------------#
# Procedimiento secuencial (top-down, Dickey-Fuller/Enders) para elegir qué
# términos determinísticos incluir en la regresión de test:
#
#   Paso 1: modelo con TENDENCIA y constante (type = "trend")
#     - phi3 testea la hipótesis conjunta H0: gamma = 0 Y beta.tendencia = 0.
#       Si SE RECHAZA phi3 -> la tendencia es relevante, nos quedamos acá y
#       usamos tau3 para la conclusión de raíz unitaria.
#       Si NO se rechaza -> la tendencia sobra, bajamos al Paso 2.
#   Paso 2: modelo solo con CONSTANTE (type = "drift")
#     - PISO del procedimiento: nunca bajamos al modelo sin constante
#       (type = "none"). Para series económicas (inflación, tasas de
#       crecimiento, etc.) casi nunca tiene sentido testear con media
#       forzada a cero.
#     - Usamos tau2 para la conclusión de raíz unitaria, y de paso miramos
#       el t-stat de la constante (informativo, no cambia la especificación).
#
# En todos los casos, la conclusión final de raíz unitaria es con el
# estadístico tau correspondiente: H0 = raíz unitaria, se rechaza (=> serie
# estacionaria) si tau ESTÁ POR DEBAJO del valor crítico (test de cola
# izquierda).
#
# NOTA TÉCNICA: en urca, el slot @testreg de un objeto ur.df YA es el
# resultado de summary(lm(...)), no el lm crudo. Por eso acá se llama
# coef(@testreg) directamente, SIN envolver de nuevo en summary().

adf.sequential <- function(x, alpha = 0.05) {
  
  x <- na.omit(x)
  lag.max <- trunc(12 * (length(x) / 100)^(1 / 4)) # regla de Schwert
  cv.col <- paste0(alpha * 100, "pct")             # "5pct" para alpha = 0.05
  
  # --- Paso 1: modelo con tendencia y constante ---
  test.trend <- ur.df(x, type = "trend", lags = lag.max, selectlags = "AIC")
  phi3.stat  <- unname(test.trend@teststat[1, "phi3"])
  phi3.cval  <- unname(test.trend@cval["phi3", cv.col])
  trend.util <- phi3.stat > phi3.cval # TRUE = se rechaza H0 -> tendencia relevante
  
  if (isTRUE(trend.util)) {
    
    modelo      <- "trend"
    tau.stat    <- unname(test.trend@teststat[1, "tau3"])
    tau.cval    <- unname(test.trend@cval["tau3", cv.col])
    const.tstat <- unname(coef(test.trend@testreg)["(Intercept)", "t value"])
    
  } else {
    
    # --- Paso 2: modelo solo con constante (piso del procedimiento) ---
    test.drift  <- ur.df(x, type = "drift", lags = lag.max, selectlags = "AIC")
    modelo      <- "drift"
    tau.stat    <- unname(test.drift@teststat[1, "tau2"])
    tau.cval    <- unname(test.drift@cval["tau2", cv.col])
    const.tstat <- unname(coef(test.drift@testreg)["(Intercept)", "t value"])
    
  }
  
  estacionaria     <- tau.stat < tau.cval        # H0 raíz unitaria; se rechaza si tau < cval
  const.significat <- abs(const.tstat) > 1.96    # informativo, no cambia la especificación
  
  list(modelo = modelo, tau.stat = tau.stat, tau.cval = tau.cval,
       estacionaria = estacionaria, const.tstat = const.tstat,
       const.significat = const.significat)
}

# Lista de series a testear
vars.list <- list(
  TCN = tcn, IPC = ipc, EXP_INF = exp_inf_12m, EXP_TCN = exp_tcn_12m,
  TOT = tot, EMBI = embi, FED = fed_mp, EBP = ebp, OIL = oil_shock
)

resultados <- data.frame(
  Variable     = character(),
  Modelo_Nivel = character(),
  ADF_Nivel    = character(),
  Const_Nivel  = character(),
  Modelo_Dif   = character(),
  ADF_Dif      = character(),
  Const_Dif    = character(),
  Sugerencia   = character(),
  stringsAsFactors = FALSE
)

for (nm in names(vars.list)) {
  
  x  <- vars.list[[nm]]
  dx <- diff(x)
  
  adf.n <- adf.sequential(x)
  adf.d <- adf.sequential(dx)
  
  sugerencia <- if (!adf.n$estacionaria && adf.d$estacionaria) {
    "Diferencias (I(1))"
  } else if (adf.n$estacionaria) {
    "Niveles (I(0))"
  } else {
    "No estacionaria ni en niveles ni en dif.: revisar (¿I(2)? ¿mala especificación?)"
  }
  
  resultados <- rbind(resultados, data.frame(
    Variable     = nm,
    Modelo_Nivel = adf.n$modelo,
    ADF_Nivel    = sprintf("%.2f (cv5%%=%.2f) -> %s", adf.n$tau.stat, adf.n$tau.cval,
                           ifelse(adf.n$estacionaria, "estac.", "no estac.")),
    Const_Nivel  = sprintf("t=%.2f -> %s", adf.n$const.tstat,
                           ifelse(adf.n$const.significat, "signif.", "no signif.")),
    Modelo_Dif   = adf.d$modelo,
    ADF_Dif      = sprintf("%.2f (cv5%%=%.2f) -> %s", adf.d$tau.stat, adf.d$tau.cval,
                           ifelse(adf.d$estacionaria, "estac.", "no estac.")),
    Const_Dif    = sprintf("t=%.2f -> %s", adf.d$const.tstat,
                           ifelse(adf.d$const.significat, "signif.", "no signif.")),
    Sugerencia   = sugerencia,
    stringsAsFactors = FALSE
  ))
}

print(resultados, row.names = FALSE)

write.csv(resultados, here("output", "tests_raiz_unitaria.csv"), row.names = FALSE)

# Para inspeccionar el detalle completo de un test puntual (ej. TCN, niveles):
# summary(ur.df(na.omit(tcn), type = "trend", lags = 12, selectlags = "AIC"))

#------------------------------------------------------------------------------#
# Test de raíz unitaria con quiebre estructural (Zivot-Andrews) - IPC
#------------------------------------------------------------------------------#
# El ADF secuencial dio un resultado borderline para el IPC en diferencias
# (inflación mensual): estadístico -2.88 vs. cv5% -2.88, prácticamente un
# empate. Un motivo típico es un quiebre estructural en la muestra (ej. el
# shock inflacionario global de 2021-2022), que reduce la potencia del ADF
# estándar. Zivot-Andrews (urca::ur.za) permite que la alternativa incluya un
# quiebre endógeno -en el intercepto, la tendencia, o ambos- y es el remedio
# habitual en estos casos.
#
# H0: raíz unitaria (sin quiebre).
# H1: estacionaria, con un quiebre estructural en la fecha que el propio test
#     elige de forma endógena (la que da MÁS evidencia en contra de H0).
# Se rechaza H0 si el estadístico es MENOR al valor crítico (igual que ADF,
# pero con su propia tabla de valores críticos, específica de este test).

# --- IPC en diferencias (inflación mensual) ---
dipc <- diff(ipc)

za.dipc <- ur.za(dipc, model = "intercept", lag = NULL)
summary(za.dipc)

# Traducir el punto de quiebre (índice de la serie) a fecha calendario
fecha.quiebre.dipc <- time(dipc)[za.dipc@bpoint]
cat("Quiebre estimado en IPC (diferencias):", fecha.quiebre.dipc, "\n")

# --- IPC en niveles (para comparar, model = "both" porque el nivel de
#     precios sí puede tener quiebre tanto en intercepto como en tendencia) ---
za.ipc <- ur.za(ipc, model = "both", lag = NULL)
summary(za.ipc)

fecha.quiebre.ipc <- time(ipc)[za.ipc@bpoint]
cat("Quiebre estimado en IPC (niveles):", fecha.quiebre.ipc, "\n")


#------------------------------------------------------------------------------#
# Chequeo de robustez de TOT: batería completa (ERS, KPSS, PP, ZA)
# El ADF secuencial dio TOT I(0), pero de forma marginal (-3.14 vs -2.88).
# Confirmamos con tests de nula opuesta (KPSS) y mayor potencia (ERS).
#------------------------------------------------------------------------------#

lag_tot <- trunc(12 * (length(na.omit(tot)) / 100)^(1/4))  # regla de Schwert

cat("\n---- TOT: ERS (DF-GLS, const) ----\n")
print(summary(ur.ers(tot, type = "DF-GLS", model = "const", lag.max = lag_tot)))

cat("\n---- TOT: KPSS (mu, long) ----\n")
print(summary(ur.kpss(tot, type = "mu", lags = "long")))

cat("\n---- TOT: Phillips-Perron (constant, long) ----\n")
print(summary(ur.pp(tot, type = "Z-tau", model = "constant", lags = "long")))

cat("\n---- TOT: Zivot-Andrews (intercepto, lag AIC) ----\n")
print(summary(ur.za(tot, model = "intercept", lag = NULL)))










#===============================================================================
#===============================INCISO 2========================================
#===============================================================================

#------------------------------------------------------------------------------#
# Datos: cambio en TCN y cambio en IPC
#------------------------------------------------------------------------------#
# Según el inciso 1, ambas series son I(1) en niveles -> se trabaja con la
# primera diferencia (ya en 100*log, así que diff() da una var. porcentual
# aprox. mensual).

dtcn <- diff(tcn)
dipc <- diff(ipc)

# ORDEN DE CHOLESKY: TCN primero, IPC segundo.
# Un shock al TCN puede afectar al IPC dentro del mismo período (traspaso
# contemporáneo); un shock al IPC no mueve al TCN en el mismo período.
Y <- cbind(dtcn, dipc)
colnames(Y) <- c("dtcn", "dipc")

#------------------------------------------------------------------------------#
# VAR Estimation (Reduced Form)
#------------------------------------------------------------------------------#

library(vars)

# Selección de rezagos
pmax <- 12 # Rezago máximo

popt <- VARselect(Y, lag.max = pmax, type = "const")
popt
p <- popt$selection["HQ(n)"]


VAR <- VAR(Y, p = p, type = "const")

m <- VAR$K   # Número de variables (2: dtcn, dipc)
T <- VAR$obs # Observaciones efectivas, netas de los p rezagos iniciales

# Chequeos de especificación
roots(VAR, modulus = TRUE) # Estabilidad: todas las raíces deben estar < 1

h.BG <- 6
serial.test(VAR, lags.bg = h.BG, type = "ES") # Autocorrelación residual

#------------------------------------------------------------------------------#
# SVAR Estimation (Structural Form) - Identificación Cholesky (recursiva)
#------------------------------------------------------------------------------#
# Se implementa el Cholesky vía el modelo AB de vars::SVAR(): Amat
# triangular inferior con 1 en la diagonal y NA por debajo (a estimar),
# Bmat diagonal con NA (a estimar). Con esta estructura, la matriz de
# impacto implícita P = A^{-1} B coincide exactamente con la factorización
# de Cholesky de la matriz de var-cov de los residuos reducidos.

# A Matrix (triangular inferior, orden: fila/columna 1 = dtcn, 2 = dipc)
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}
Amat

# B Matrix (diagonal)
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}
Bmat

# SVAR Estimation (AB model configuration)
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR

# Verificación: la matriz de impacto implícita por el modelo AB debe
# coincidir con la Cholesky directa de la matriz de var-cov de residuos.
P.SVAR <- solve(SVAR$A, SVAR$B) # inv(A) %*% B

S <- t(resid(VAR)) %*% resid(VAR) / (T - m * p - 1)
P.chol <- t(chol(S))

P.SVAR
P.chol # Deberían ser (casi) idénticas

#------------------------------------------------------------------------------#
# SVAR Analysis: IRF ante shock al TCN
#------------------------------------------------------------------------------#

source(here("tools", "PS2_Miscellaneous.R"))
source(here("tools", "PS2_Bootstrap.R"))
source(here("tools", "PS2_SVAR_Tools.R"))
source(here("tools", "PS2_SVAR_Plots.R"))

# FIX: boot.estimate() de PS2_Bootstrap.R no devuelve el VAR reestimado
# cuando resmat = NULL (nuestro caso, ya que no restringimos el VAR
# original). Se corrige agregando la línea final explícita.
boot.estimate <- function(var.names, Y, m, p, lag.max, ic, resmat = NULL) {
  colnames(Y) <- var.names
  VAR <- VAR(Y, p = p, type = "const")
  if (!is.null(resmat)) { VAR <- restrict(VAR, method = "man", resmat = resmat) }
  VAR
}



H <- 24        # Horizonte para IRF/FEVD
H.ERPT <- 36   # Horizonte para ERPT

# IRF (no acumulada): respuesta de Delta-TCN y Delta-IPC a shocks
# estructurales S.1 (TCN) y S.2 (IPC)
IRF <- SVAR.sirf(SVAR, H)
plot.sirf(IRF, m, H)

# IRF acumulada (nivel del TCN y del IPC ante los mismos shocks)
IRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
plot.sirf(IRF.c, m, H)

# Extracción puntual: IRF ante shock al TCN (S.1) únicamente, en un solo
# gráfico con las dos respuestas (nivel, ya que es lo relevante para ERPT)
par(mfrow = c(2, 1))

plot(0:H, IRF.c["DTCN", "S.1", ], type = "o", lwd = 2,
     main = "Respuesta del TCN (nivel) ante shock al TCN",
     xlab = "Horizonte", ylab = "")
grid(NULL, NULL, lty = 1)

plot(0:H, IRF.c["DIPC", "S.1", ], type = "o", lwd = 2,
     main = "Respuesta del IPC (nivel) ante shock al TCN",
     xlab = "Horizonte", ylab = "")
grid(NULL, NULL, lty = 1)

# FEVD
FEVD <- SVAR.fevd(SVAR, H)
plot.fevd(FEVD, m, H)

# HD (Historical Decomposition)
HD <- SVAR.hd(SVAR)
plot.hd(Y, HD, m)

#------------------------------------------------------------------------------#
# ERPT incondicional (promedio)
#------------------------------------------------------------------------------#
# ERPT_h(TCN) = CIRF_h(IPC | shock TCN) / CIRF_h(TCN | shock TCN)
# En SVAR.erpt(SVAR, H, vx, vy): vx = índice de la variable numerador (IPC),
# vy = índice de la variable denominador (TCN). Con el orden Y = (dtcn,
# dipc): dtcn = fila/shock 1, dipc = fila/shock 2 -> vx = 2, vy = 1.

ERPT <- SVAR.erpt(SVAR, H.ERPT, 2, 1)
plot.erpt(ERPT, H.ERPT)

# Tabla de ERPT para algunos horizontes de referencia
horizontes.tabla <- c(0, 1, 3, 6, 12, 24, 36)
data.frame(Horizonte = horizontes.tabla, ERPT_pct = round(ERPT[horizontes.tabla + 1], 2))

#------------------------------------------------------------------------------#
# Bootstrap Inference
#------------------------------------------------------------------------------#

R <- 500 # Número de replicaciones bootstrap
type <- "nonparametric"
gamma <- 0.95 # Nivel de confianza

# lag.max / ic: no se usan efectivamente dentro de boot.estimate() (solo
# reestima un VAR(p) con el p ya elegido), pero deben existir en el entorno
# porque SVAR.sirf.boot() y SVAR.erpt.boot() las referencian.
lag.max <- pmax
ic <- "HQ"
set.seed(123)

# Replicaciones bootstrap de los datos
Y.boot <- boot.replicate(VAR, R, type)

# IRF (Bootstrap, no acumulada)
IRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, H, gamma, Y.boot)
plot.sirf.boot(IRF.boot, m, H)

# IRF (Bootstrap, acumulada / nivel)
IRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, H, gamma, Y.boot, cumulative = TRUE)
plot.sirf.boot(IRF.c.boot, m, H)

# FEVD (Bootstrap)
FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, H, gamma, Y.boot)
plot.fevd.boot(FEVD.boot, m, H)

# ERPT (Bootstrap)
ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, H.ERPT, 2, 1, gamma, Y.boot)
plot.erpt.boot(ERPT.boot, H.ERPT)

# Tabla de ERPT con bandas de confianza para los horizontes de referencia
data.frame(
  Horizonte = horizontes.tabla,
  ERPT_pe   = round(ERPT.boot$pe[horizontes.tabla + 1], 2),
  ERPT_lb   = round(ERPT.boot$lb[horizontes.tabla + 1], 2),
  ERPT_ub   = round(ERPT.boot$ub[horizontes.tabla + 1], 2)
)


#------------------------------------------------------------------------------#
# Punto 2: Gráficos para el informe (formato TP1)
#------------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)

# --- Panel de IRF (nivel) ante shock al TCN, con bandas bootstrap ---

h <- 0:H
vars_orden <- c("DTCN", "DIPC")
vars_label <- c(DTCN = "TCN (nivel, 100*log)", DIPC = "IPC (nivel, 100*log)")

df_panel <- lapply(vars_orden, function(v) {
  data.frame(
    horizonte = h,
    variable  = vars_label[v],
    pe = IRF.c.boot$pe[v, "S.1", ],
    lb = IRF.c.boot$lb[v, "S.1", ],
    ub = IRF.c.boot$ub[v, "S.1", ]
  )
}) %>% bind_rows()

df_panel$variable <- factor(df_panel$variable, levels = vars_label[vars_orden])

p_panel <- ggplot(df_panel, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "steelblue", alpha = 0.18) +
  geom_line(aes(y = pe), color = "black", linewidth = 0.9) +
  facet_wrap(~ variable, nrow = 1, scales = "free_y") +
  labs(x = "Meses", y = "Respuesta acumulada (niveles)") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank())

p_panel

ggsave(here("output", "irf_punto2_tcn_ipc.pdf"), p_panel, width = 7, height = 3.8, units = "in")

# --- ERPT incondicional con banda bootstrap ---

df_erpt <- data.frame(
  horizonte = 0:H.ERPT,
  pe = ERPT.boot$pe,
  lb = ERPT.boot$lb,
  ub = ERPT.boot$ub
)

p_erpt <- ggplot(df_erpt, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub), fill = "steelblue", alpha = 0.18) +
  geom_line(aes(y = pe), color = "black", linewidth = 0.9) +
  labs(x = "Meses", y = "ERPT acumulado (%)") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())

p_erpt

ggsave(here("output", "erpt_punto2.pdf"), p_erpt, width = 5, height = 3.6, units = "in")






# ---- Guardado de objetos ----
saveRDS(list(VAR = VAR, SVAR = SVAR, Y = Y, p = p, m = m,
             Amat = Amat, Bmat = Bmat, H = H, H.ERPT = H.ERPT, gamma = gamma,
             IRF.c.boot = IRF.c.boot, ERPT.boot = ERPT.boot, Y.boot = Y.boot),
        here("data", "objetos_02.rds"))



