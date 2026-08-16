#------------------------------------------------------------------------------#
# Universidad de San Andrés - Macroeconometría 2026
# TP2: Exchange Rate Pass-Through - Colombia
#------------------------------------------------------------------------------#

#===============================================================================
#===============================INCISO 3========================================
#===============================================================================

#   Estimar las IRF del TCN y del IPC ante un shock estructural al TCN mediante
#   PROYECCIONES LOCALES (Jordà, 2005), y comparar contra las IRF y el ERPT del
#   VAR del inciso 2.
#
#   Esquema en dos pasos (enunciado):
#     Paso 1 (inferencia del shock): se reutiliza el VAR bivariado del inciso 2
#       (dtcn, dipc) con identificación Cholesky TCN -> IPC. El shock estructural
#       al TCN, u^b_t, se recupera ortogonalizando los residuos reducidos.
#     Paso 2 (estimación de IRF): para cada horizonte h se corre la proyección
#       local   x_{t+h} - x_{t-1} = a_h + b_h * u^b_t + g_h(L) Y_{t-1} + w_{t+h},
#       con x = TCN o IPC. El coeficiente b_h ES la IRF acumulada en h.
#
# Este script NO recalcula el inciso 2: LEE los objetos ya estimados desde
# 'data/objetos_02.rds' (guardados al final de 02_Ejercicios1_2.R). Los .rds de
# datos se leen igual que en 02, por consistencia.
#
# Herramientas de la cátedra (viven en 'tools/'):
#   - PS3_LP_Tools.R : lp.shock(), lp.multiplier.shock(), get.shock(), HAC
#   - PS3_LP_Plots.R : plot.lp()
#------------------------------------------------------------------------------#

remove(list = ls(all.names = TRUE))
gc()

library(here)
library(sandwich)   # vcovHAC, usado internamente por las funciones LP
library(ggplot2)
library(dplyr)



#------------------------------------------------------------------------------#
# 0. Carga de objetos (NO se recalcula el inciso 2)
#------------------------------------------------------------------------------#
# Datos base (mismos .rds que usa el 02). No son estrictamente necesarios para
# las LP, pero se cargan por consistencia y para eventuales chequeos.
datos_ts <- readRDS(here("data", "datos_ts.rds"))
datos    <- readRDS(here("data", "panel_final.rds"))

# Objetos del inciso 2: VAR, SVAR, Y, p, m, Amat, Bmat, H, H.ERPT, gamma,
# IRF.c.boot, ERPT.boot, Y.boot.
obj02 <- readRDS(here("data", "objetos_02.rds"))
for (nm in names(obj02)) assign(nm, obj02[[nm]])
remove(obj02)

p <- unname(p)   # p venía como named int ("HQ(n)"); lo dejamos como escalar limpio



#------------------------------------------------------------------------------#
# 1. Herramientas LP del curso (tools)
#------------------------------------------------------------------------------#

source(here("tools", "PS3_LP_Tools.R"))
source(here("tools", "PS3_LP_Plots.R"))  # ya trae plot.lp()
       # vcovHAC
#--------------------------------------------------------------------------



#------------------------------------------------------------------------------#
# 2. Shock estructural al TCN (u^b_t)
#------------------------------------------------------------------------------#
# Shocks estructurales (T x m): columna 1 = shock a dtcn, columna 2 = shock a dipc
U <- t(solve(SVAR$B) %*% SVAR$A %*% t(resid(VAR)))
colnames(U) <- colnames(Y)

u.tcn <- U[, "dtcn"]  # u^b_t del enunciado

saveRDS(list(u.tcn = u.tcn, Y = Y, p = p, H = H, gamma = gamma),
        here("data", "objetos_03.rds"))

#------------------------------------------------------------------------------#
# 3. Paso 2 - IRF acumuladas vía LP ante el shock al TCN
#------------------------------------------------------------------------------#
# IRF vía LP: respuesta del TCN (nivel) y del IPC (nivel) ante el shock al TCN
LP.tcn <- lp.shock(Y, u.tcn, p, idx.r = 1, H, gamma, cumulative = TRUE)
LP.ipc <- lp.shock(Y, u.tcn, p, idx.r = 2, H, gamma, cumulative = TRUE)

# plot.lp(LP.tcn)
# plot.lp(LP.ipc)

# ERPT vía LP: ratio directo, sin dividir IRFs
LP.ERPT <- lp.multiplier.shock(Y, u.tcn, p, idx.rl = 2, idx.rr = 1, H.ERPT, gamma)
#plot.lp(LP.ERPT)


#------------------------------------------------------------------------------#
# Gráficos ggplot para el informe: VAR (inciso 2) vs. LP (inciso 3)
#------------------------------------------------------------------------------#
library(ggplot2)
library(dplyr)

if (!dir.exists("output")) dir.create("output")

colores_metodo <- c("VAR" = "steelblue", "LP" = "firebrick")
colores_banda  <- c("VAR" = "lightblue", "LP" = "lightpink")

# --- Gráfico 1: IRF del TCN, VAR vs. LP ------------------------------------- #
df_tcn <- bind_rows(
  data.frame(horizonte = 0:H, metodo = "VAR",
             pe = IRF.c.boot$pe["DTCN", "S.1", ],
             lb = IRF.c.boot$lb["DTCN", "S.1", ],
             ub = IRF.c.boot$ub["DTCN", "S.1", ]),
  data.frame(horizonte = 0:H, metodo = "LP",
             pe = LP.tcn$pe, lb = LP.tcn$lb, ub = LP.tcn$ub)
)
df_tcn$metodo <- factor(df_tcn$metodo, levels = c("VAR", "LP"))

df_tcn$panel <- "Respuesta del TCN ante shock al TCN: VAR vs. LP"

p_tcn <- ggplot(df_tcn, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = metodo), alpha = 0.35) +
  geom_line(aes(y = pe, color = metodo), linewidth = 0.9) +
  scale_color_manual(values = colores_metodo, name = NULL) +
  scale_fill_manual(values = colores_banda, name = NULL) +
  labs(x = "Meses", y = "Respuesta acumulada (nivel)") + facet_wrap(~ panel) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11))
p_tcn
ggsave(file.path("output", "irf_punto3_tcn_var_lp.pdf"), p_tcn, width = 6, height = 4, units = "in")

# --- Gráfico 2: IRF del IPC, VAR vs. LP ------------------------------------- #
df_ipc <- bind_rows(
  data.frame(horizonte = 0:H, metodo = "VAR",
             pe = IRF.c.boot$pe["DIPC", "S.1", ],
             lb = IRF.c.boot$lb["DIPC", "S.1", ],
             ub = IRF.c.boot$ub["DIPC", "S.1", ]),
  data.frame(horizonte = 0:H, metodo = "LP",
             pe = LP.ipc$pe, lb = LP.ipc$lb, ub = LP.ipc$ub)
)
df_ipc$metodo <- factor(df_ipc$metodo, levels = c("VAR", "LP"))

df_ipc$panel <- "Respuesta del IPC ante shock al TCN: VAR vs. LP"

p_ipc <- ggplot(df_ipc, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = metodo), alpha = 0.35) +
  geom_line(aes(y = pe, color = metodo), linewidth = 0.9) +
  scale_color_manual(values = colores_metodo, name = NULL) +
  scale_fill_manual(values = colores_banda, name = NULL) +
  labs(x = "Meses", y = "Respuesta acumulada (nivel)") + facet_wrap(~ panel) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11))
p_ipc
ggsave(file.path("output", "irf_punto3_ipc_var_lp.pdf"), p_ipc, width = 6, height = 4, units = "in")

# --- Gráfico 3: ERPT, VAR vs. LP --------------------------------------------- #
df_erpt3 <- bind_rows(
  data.frame(horizonte = 0:H.ERPT, metodo = "VAR",
             pe = ERPT.boot$pe, lb = ERPT.boot$lb, ub = ERPT.boot$ub),
  data.frame(horizonte = 0:H.ERPT, metodo = "LP",
             pe = LP.ERPT$pe, lb = LP.ERPT$lb, ub = LP.ERPT$ub)
)
df_erpt3$metodo <- factor(df_erpt3$metodo, levels = c("VAR", "LP"))

df_erpt3$panel <- "ERPT incondicional: VAR vs. LP"

p_erpt3 <- ggplot(df_erpt3, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = metodo), alpha = 0.35) +
  geom_line(aes(y = pe, color = metodo), linewidth = 0.9) +
  scale_color_manual(values = colores_metodo, name = NULL) +
  scale_fill_manual(values = colores_banda, name = NULL) +
  labs(x = "Meses", y = "ERPT acumulado (%)") +
  facet_wrap(~ panel) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 11))
p_erpt3
ggsave(file.path("output", "erpt_punto3_var_lp.pdf"), p_erpt3, width = 6, height = 4, units = "in")
