#------------------------------------------------------------------------------#
# Universidad de San Andrés - Macroeconometría 2026
# TP2: Exchange Rate Pass-Through - Colombia
#
# Inciso 5: (a) distribución del shock estructural al TCN (histograma) y
#           (b) ERPT condicional por MAGNITUD del shock: tres regímenes según
#           el shock supere +1 desvío, caiga bajo -1 desvío, o quede en el medio.
#
# Este script NO recalcula los incisos previos: LEE los objetos ya guardados
# desde 'data/objetos_03.rds' (u.tcn, Y, p, H, gamma).
#
# Herramientas de la cátedra (viven en 'tools/'):
#   - PS3_LP_Tools.R : get.controls(), cumsum.h(), HAC
#------------------------------------------------------------------------------#

remove(list = ls(all.names = TRUE))
gc()

library(here)
library(sandwich)   # vcovHAC, usado internamente por las funciones LP
library(ggplot2)
library(dplyr)

#------------------------------------------------------------------------------#
# 0. Carga de objetos del inciso 3 (NO se recalcula el shock ni el VAR)
#------------------------------------------------------------------------------#

obj03 <- readRDS(here("data", "objetos_03.rds"))
for (nm in names(obj03)) assign(nm, obj03[[nm]])
remove(obj03)

p <- unname(p)   # p venía como named int; lo dejamos como escalar limpio

#------------------------------------------------------------------------------#
# 1. Herramientas LP del curso (tools)
#------------------------------------------------------------------------------#

source(here("tools", "PS3_LP_Tools.R"))


df_hist <- data.frame(shock = u.tcn)

media <- mean(u.tcn)
sd.shock <- sd(u.tcn)

p_hist <- ggplot(df_hist, aes(x = shock)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", color = "white", alpha = 0.75) +
  stat_function(fun = dnorm, args = list(mean = media, sd = sd.shock), color = "black", linewidth = 0.8, linetype = "dashed") +
  geom_vline(xintercept = media, color = "gray30", linewidth = 0.5) +
  geom_vline(xintercept = media + sd.shock, color = "gray30", linewidth = 0.5, linetype = "dotted") +
  geom_vline(xintercept = media - sd.shock, color = "gray30", linewidth = 0.5, linetype = "dotted") +
  labs(x = "Shock estructural al TCN", y = "Densidad", title = "Distribucion del shock estructural al TCN") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 11))

p_hist

ggsave(here("output", "hist_punto5_shock.pdf"), p_hist, 
       width = 7.5, height = 3, units = "in")


#------------------------------------------------------------------------------#
# Inciso 5 (continuación): ERPT por magnitud del shock (3 regímenes)
#------------------------------------------------------------------------------#

# Dummies de magnitud, sobre el shock estructural (umbral: 1 sd, ya calculado)
D <- ifelse(u.tcn > sd.shock, 1, 0)   # muy positivo ("grande")
B <- ifelse(u.tcn < -sd.shock, 1, 0)  # muy negativo
# El régimen "normal" es (1 - D - B), no necesita vector propio

lp.multiplier.size.shock <- function(Y, shock, D, B, p, idx.rl, idx.rr, H, gamma) {
  
  m <- ncol(Y); T <- nrow(Y)
  Wt.full <- get.controls(Y, p, m)
  rsp.l <- Y[(p + 1):T, idx.rl]  # numerador (IPC)
  rsp.r <- Y[(p + 1):T, idx.rr]  # denominador (TCN)
  
  irf.pe <- array(NA, c(H + 1, 3))
  irf.lb <- array(NA, dim(irf.pe))
  irf.ub <- array(NA, dim(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    yh.l <- cumsum.h(rsp.l, h)
    yh.r <- cumsum.h(rsp.r, h)
    
    Wt <- Wt.full[1:(T - p - h), ]
    st <- shock[1:(T - p - h)]
    Dt <- D[1:(T - p - h)]
    Bt <- B[1:(T - p - h)]
    Nt <- 1 - Dt - Bt  # régimen "normal"
    
    # Controles totalmente interactuados por régimen
    Xc <- cbind(Dt * cbind(1, Wt), Bt * cbind(1, Wt), Nt * cbind(1, Wt))
    
    yh.l <- resid(lm(yh.l ~ -1 + Xc))
    yh.r <- resid(lm(yh.r ~ -1 + Xc))
    
    # Proyección del TCN (residual) sobre el shock, interactuado por régimen
    st.D <- Dt * st
    st.B <- Bt * st
    st.N <- Nt * st
    yh.r.fit <- fitted(lm(yh.r ~ -1 + st.D + st.B + st.N))
    
    # Regresión final: coeficientes separados por régimen
    proj <- lm(yh.l ~ -1 + I(Dt * yh.r.fit) + I(Bt * yh.r.fit) + I(Nt * yh.r.fit))
    
    b <- 100 * unname(coef(proj))
    se <- 100 * unname(sqrt(diag(vcovHAC(proj))))
    
    irf.pe[h + 1, ] <- b
    irf.lb[h + 1, ] <- b - z * se
    irf.ub[h + 1, ] <- b + z * se
  }
  
  dimnames(irf.pe) <- list(0:H, c("Grande (+)", "Grande (-)", "Normal"))
  dimnames(irf.lb) <- dimnames(irf.pe)
  dimnames(irf.ub) <- dimnames(irf.pe)
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
}

LP.ERPT.size <- lp.multiplier.size.shock(Y, u.tcn, D, B, p, idx.rl = 2, idx.rr = 1, H, gamma)

#------------------------------------------------------------------------------#
# Gráfico: ERPT por magnitud del shock (Normal, Grande +, Grande -)
#------------------------------------------------------------------------------#

colores_size <- c("Normal" = "black", "Grande (+)" = "firebrick3", "Grande (-)" = "darkgoldenrod2")
colores_banda_size <- c("Normal" = "gray70", "Grande (+)" = "lightpink1", "Grande (-)" = "lightgoldenrod1")

df_erpt5 <- bind_rows(
  data.frame(horizonte = 0:H, regimen = "Normal",
             pe = LP.ERPT.size$pe[, "Normal"], lb = LP.ERPT.size$lb[, "Normal"], ub = LP.ERPT.size$ub[, "Normal"]),
  data.frame(horizonte = 0:H, regimen = "Grande (+)",
             pe = LP.ERPT.size$pe[, "Grande (+)"], lb = LP.ERPT.size$lb[, "Grande (+)"], ub = LP.ERPT.size$ub[, "Grande (+)"]),
  data.frame(horizonte = 0:H, regimen = "Grande (-)",
             pe = LP.ERPT.size$pe[, "Grande (-)"], lb = LP.ERPT.size$lb[, "Grande (-)"], ub = LP.ERPT.size$ub[, "Grande (-)"])
)
df_erpt5$regimen <- factor(df_erpt5$regimen, levels = c("Normal", "Grande (+)", "Grande (-)"))

p_erpt5 <- ggplot(df_erpt5, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = regimen), alpha = 0.45) +
  geom_line(aes(y = pe, color = regimen), linewidth = 0.9) +
  scale_color_manual(values = colores_size, name = NULL) +
  scale_fill_manual(values = colores_banda_size, name = NULL) +
  coord_cartesian(ylim = c(-120, 120)) +
  labs(x = "Meses", y = "ERPT acumulado (%)",
       title = "ERPT por magnitud del shock") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 11))
p_erpt5

ggsave(here("output", "erpt_punto5_size.pdf"), p_erpt5, width = 6, height = 4, units = "in")
