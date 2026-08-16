#------------------------------------------------------------------------------#
# Universidad de San Andrés - Macroeconometría 2026
# TP2: Exchange Rate Pass-Through - Colombia
# INCISO 6: VAR de 7 variables (shock Fed, shock petróleo, EBP, TOT, EMBI,
# TCN, IPC), identificación Cholesky con ese orden + restricción de bloque
# exógeno (Cushman-Zha, 1997). Los 6 shocks estructurales resultantes
# (se excluye el propio del IPC) se usan como instrumento en una LP-IV
# para estimar el ERPT condicional a cada uno.
#------------------------------------------------------------------------------#

remove(list = ls(all.names = TRUE))
gc()

library(vars)
library(sandwich)
library(here)

datos_ts <- readRDS(here("data", "datos_ts.rds"))

tcn       <- datos_ts[, "tcn"]
ipc       <- datos_ts[, "ipc"]
tot       <- datos_ts[, "tot"]
embi      <- datos_ts[, "embi"]
fed_mp    <- datos_ts[, "fed_mp"]
ebp       <- datos_ts[, "ebp"]
oil_shock <- datos_ts[, "oil_shock"]

#------------------------------------------------------------------------------#
# Datos: 3 shocks globales en niveles + TOT, EMBI, TCN, IPC en diferencias
#------------------------------------------------------------------------------#

dtot  <- diff(tot)
dembi <- diff(embi)
dtcn  <- diff(tcn)
dipc  <- diff(ipc)

fed_mp.a    <- window(fed_mp,    start = start(dtcn))
oil_shock.a <- window(oil_shock, start = start(dtcn))
ebp.a       <- window(ebp,       start = start(dtcn))

Y <- cbind(fed_mp.a, oil_shock.a, ebp.a, dtot, dembi, dtcn, dipc)
colnames(Y) <- c("fed_mp", "oil_shock", "ebp", "dtot", "dembi", "dtcn", "dipc")

#------------------------------------------------------------------------------#
# Selección de rezagos
#------------------------------------------------------------------------------#

pmax <- 12
popt <- VARselect(Y, lag.max = pmax, type = "const")
popt

h.BG <- 6

Y.trim <- Y[2:nrow(Y), ]
VAR.p1 <- VAR(Y.trim, p = 1, type = "const")
serial.test(VAR.p1, lags.bg = h.BG, type = "ES")

VAR.p2 <- VAR(Y, p = 2, type = "const")
serial.test(VAR.p2, lags.bg = h.BG, type = "ES")

p <- 2

#------------------------------------------------------------------------------#
# VAR Estimation (Reduced Form)
#------------------------------------------------------------------------------#

VAR <- VAR(Y, p = p, type = "const")
m <- VAR$K
T <- VAR$obs

roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = h.BG, type = "ES")

#------------------------------------------------------------------------------#
# SVAR Estimation - Cholesky recursiva + bloque exógeno (Cushman-Zha)
#------------------------------------------------------------------------------#

Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

resmat <- matrix(1, m, m * p + 1)
for (lag in 0:(p - 1)) {
  cols.domesticas <- (lag * m + 4):(lag * m + 7)
  resmat[1:3, cols.domesticas] <- 0
}

VAR.restricted <- restrict(VAR, method = "man", resmat = resmat)

SVAR <- SVAR(VAR.restricted, Amat = Amat, Bmat = Bmat, lrtest = FALSE, max.iter = 1000)
SVAR

# Verificación de identificación: P %*% t(P) debe reproducir S
P.SVAR <- solve(SVAR$A, SVAR$B)
S <- t(resid(VAR.restricted)) %*% resid(VAR.restricted) / (T - m * p - 1)
S.SVAR <- P.SVAR %*% t(P.SVAR)
max(abs(S.SVAR - S) / abs(S))

# Diagnósticos sobre el VAR restringido
roots(VAR.restricted, modulus = TRUE)
serial.test(VAR.restricted, lags.bg = h.BG, type = "ES")

res.restricted <- residuals(VAR.restricted)
cat("--- Ljung-Box por ecuación (VAR restringido) ---\n")
for (v in colnames(res.restricted)) {
  bt <- Box.test(res.restricted[, v], lag = h.BG, type = "Ljung-Box")
  cat(v, ": p-value =", round(bt$p.value, 4), "\n")
}

#------------------------------------------------------------------------------#
# Extracción de los 6 shocks estructurales (para usarlos como instrumento)
#------------------------------------------------------------------------------#

source(here("tools", "PS3_LP_Tools.R"))  # get.controls, lp.shock, lp.multiplier.shock, etc.

E.hat <- resid(VAR.restricted)
struct.shocks <- t(solve(SVAR$B) %*% SVAR$A %*% t(E.hat))
colnames(struct.shocks) <- colnames(Y)

# Persistir los 6 shocks estructurales (+ Y y p) para el inciso 10
saveRDS(list(struct.shocks = struct.shocks, Y = Y, p = p),
        here("data", "objetos_06.rds"))

#------------------------------------------------------------------------------#
# ERPT condicional (LP-IV) para los primeros 6 shocks estructurales
#------------------------------------------------------------------------------#

H.ERPT <- 36
gamma <- 0.95

shocks.nombres <- c("Shock Fed (JK)", "Shock Petróleo (BH)", "Excess Bond Premium",
                    "Términos de Intercambio", "EMBI", "TCN")

ERPT.iv <- vector("list", 6)
names(ERPT.iv) <- shocks.nombres

for (s in 1:6) {
  ERPT.iv[[s]] <- lp.multiplier.shock(Y, struct.shocks[, s], p,
                                      idx.rl = 7, idx.rr = 6,   # 7=dipc, 6=dtcn
                                      H.ERPT, gamma)
}

#------------------------------------------------------------------------------#
# Diagnóstico: F-stat de primera etapa por shock y horizonte
#------------------------------------------------------------------------------#

Wt.full <- get.controls(Y, p, m)
rsp.r <- Y[(p + 1):nrow(Y), "dtcn"]

cat("--- F-stat de primera etapa por shock ---\n")
for (s in 1:6) {
  st <- struct.shocks[, s]
  Fstats <- rep(NA, H.ERPT + 1)
  for (h in 0:H.ERPT) {
    n <- length(st) - h
    Wt.h <- Wt.full[1:n, ]
    st.h <- st[1:n]
    yh.r <- cumsum.h(rsp.r, h)
    yh.r <- resid(lm(yh.r ~ -1 + Wt.h))
    first.stage <- lm(yh.r ~ -1 + st.h)
    Fstats[h + 1] <- summary(first.stage)$fstatistic[1]
  }
  cat(shocks.nombres[s], ": F-stat mínimo =", round(min(Fstats), 1),
      "| mediana =", round(median(Fstats), 1), "\n")
}

#------------------------------------------------------------------------------#
# Gráfico: panel de ERPT condicional por shock (LP-IV), con eje y acotado
# por panel (zoom, sin descartar datos) via patchwork
#------------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)
library(patchwork)

h <- 0:H.ERPT

df_erpt6 <- lapply(1:6, function(s) {
  data.frame(
    horizonte = h,
    shock     = shocks.nombres[s],
    pe = ERPT.iv[[s]]$pe,
    lb = ERPT.iv[[s]]$lb,
    ub = ERPT.iv[[s]]$ub
  )
}) %>% bind_rows()

df_erpt6$shock <- factor(df_erpt6$shock, levels = shocks.nombres)

# Límites de zoom por panel: mediana +/- 4*IQR de esa serie. Solo define
# el rango visible del eje y - ningún dato se elimina del data frame.
limites6 <- df_erpt6 %>%
  group_by(shock) %>%
  summarise(
    centro = median(pe, na.rm = TRUE),
    escala = IQR(c(pe, lb, ub), na.rm = TRUE),
    ymin = centro - 4 * escala,
    ymax = centro + 4 * escala,
    .groups = "drop"
  )

hacer_panel6 <- function(shock.nombre) {
  lims <- limites6 %>% dplyr::filter(shock == shock.nombre)
  
  ggplot(df_erpt6 %>% dplyr::filter(shock == shock.nombre), aes(x = horizonte)) +
    geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
    geom_ribbon(aes(ymin = lb, ymax = ub), fill = "steelblue", alpha = 0.18) +
    geom_line(aes(y = pe), color = "black", linewidth = 0.9) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax)) +
    labs(x = NULL, y = NULL, title = shock.nombre) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          panel.grid.minor = element_blank())
}

paneles6 <- lapply(shocks.nombres, hacer_panel6)

p_erpt6 <- wrap_plots(paneles6, nrow = 2) +
  plot_annotation(caption = "Meses (eje x) | ERPT condicional (%) (eje y)")

p_erpt6

if (!dir.exists(here("output"))) dir.create(here("output"))
ggsave(here("output", "erpt_punto6_LPIV.pdf"), p_erpt6,
       width = 10, height = 6.5, units = "in")

#------------------------------------------------------------------------------#
# Tabla resumen: ERPT condicional de largo plazo (h=36) por shock, con F-stat
#------------------------------------------------------------------------------#

Fstat.resumen <- sapply(1:6, function(s) {
  st <- struct.shocks[, s]
  n <- length(st) - H.ERPT
  Wt.h <- Wt.full[1:n, ]
  st.h <- st[1:n]
  yh.r <- cumsum.h(rsp.r, H.ERPT)
  yh.r <- resid(lm(yh.r ~ -1 + Wt.h))
  summary(lm(yh.r ~ -1 + st.h))$fstatistic[1]
})

tabla.erpt6 <- data.frame(
  Shock   = shocks.nombres,
  ERPT_pe = sapply(ERPT.iv, function(x) round(x$pe[H.ERPT + 1], 2)),
  ERPT_lb = sapply(ERPT.iv, function(x) round(x$lb[H.ERPT + 1], 2)),
  ERPT_ub = sapply(ERPT.iv, function(x) round(x$ub[H.ERPT + 1], 2)),
  Fstat_h36 = round(Fstat.resumen, 2)
)
print(tabla.erpt6, row.names = FALSE)

write.csv(tabla.erpt6, here("output", "tabla_erpt_punto6.csv"), row.names = FALSE)