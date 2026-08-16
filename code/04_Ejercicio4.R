#------------------------------------------------------------------------------#
# Universidad de San Andrés - Macroeconometría 2026
# TP2: Exchange Rate Pass-Through - Colombia
#
# Inciso 4: IRF y ERPT asimétricos por signo del shock estructural al TCN.
#   Extensión piecewise del LP: se separa la respuesta según el signo del shock
#   (positivo vs. negativo), en dos especificaciones:
#     Spec 1 = controles totalmente interactuados por régimen;
#     Spec 2 = controles compartidos, solo el shock interactuado.
#
# Este script NO recalcula los incisos previos: LEE los objetos ya guardados
# desde 'data/objetos_03.rds' (u.tcn, Y, p, H, gamma).
#
# Herramientas de la cátedra (viven en 'tools/'):
#   - PS3_LP_Tools.R : lp.pw.shock(), lp.multiplier.pw.shock(), get.controls(), HAC
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

#------------------------------------------------------------------------------#
# 2. Dummy de signo (sobre el shock estructural ya identificado en el inciso 3) y estimación piecewise
#------------------------------------------------------------------------------#

D <- ifelse(u.tcn > 0, 1, 0)

LP.ERPT.pw <- lp.multiplier.pw.shock(Y, u.tcn, D, p, idx.rl = 2, idx.rr = 1, H, gamma)

# --- Estimación: ambas specs, para TCN e IPC --------------------------------- #

LP.tcn.s1 <- lp.pw.shock(Y, u.tcn, D, p, idx.r = 1, H, gamma, 
                         cumulative = TRUE, interact.controls = TRUE)
LP.tcn.s2 <- lp.pw.shock(Y, u.tcn, D, p, idx.r = 1, H, gamma, 
                         cumulative = TRUE, interact.controls = FALSE)

LP.ipc.s1 <- lp.pw.shock(Y, u.tcn, D, p, idx.r = 2, H, gamma, 
                         cumulative = TRUE, interact.controls = TRUE)
LP.ipc.s2 <- lp.pw.shock(Y, u.tcn, D, p, idx.r = 2, H, gamma, 
                         cumulative = TRUE, interact.controls = FALSE)

# --- Comparación: ¿qué tan distintos son los point estimates? --------------- #
# Diferencia absoluta punto a punto entre Spec 1 y Spec 2, por horizonte
diff.tcn <- abs(LP.tcn.s1$pe - LP.tcn.s2$pe)
diff.ipc <- abs(LP.ipc.s1$pe - LP.ipc.s2$pe)

round(diff.tcn, 3)
round(diff.ipc, 3)

# Referencia útil: comparar esa diferencia contra el ancho típico del IC de
# Spec 1 (si la diferencia es chica en relación al IC, son "similares")
ancho.ic.tcn <- LP.tcn.s1$ub - LP.tcn.s1$lb
ancho.ic.ipc <- LP.ipc.s1$ub - LP.ipc.s1$lb

round(diff.tcn / ancho.ic.tcn, 3)  # si esto es << 1 en la mayoría de horizontes, similares
round(diff.ipc / ancho.ic.ipc, 3)


#------------------------------------------------------------------------------#
# Gráficos: IRF y ERPT asimétricos por signo del shock (Spec 1)
#------------------------------------------------------------------------------#

colores_regimen <- c("Positivo" = "darkgreen", "Negativo" = "darkorange")
colores_banda_regimen <- c("Positivo" = "palegreen3", "Negativo" = "navajowhite2")

ancho_panel <- 3.6
alto_panel  <- 3.9

# --- Gráfico 1: IRF del TCN, por régimen ------------------------------------ #
df_tcn4 <- bind_rows(
  data.frame(horizonte = 0:H, regimen = "Positivo",
             pe = LP.tcn.s1$pe[, "Positivo"], lb = LP.tcn.s1$lb[, "Positivo"], ub = LP.tcn.s1$ub[, "Positivo"]),
  data.frame(horizonte = 0:H, regimen = "Negativo",
             pe = LP.tcn.s1$pe[, "Negativo"], lb = LP.tcn.s1$lb[, "Negativo"], ub = LP.tcn.s1$ub[, "Negativo"])
)
df_tcn4$regimen <- factor(df_tcn4$regimen, levels = c("Positivo", "Negativo"))

p_tcn4 <- ggplot(df_tcn4, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = regimen), alpha = 0.35) +
  geom_line(aes(y = pe, color = regimen), linewidth = 0.9) +
  scale_color_manual(values = colores_regimen, name = NULL) +
  scale_fill_manual(values = colores_banda_regimen, name = NULL) +
  labs(x = "Meses", y = "Respuesta acumulada (nivel)",
       title = "TCN: shock positivo vs. negativo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10))
p_tcn4
ggsave(here("output", "irf_punto4_tcn_pw.pdf"), p_tcn4, width = ancho_panel, height = alto_panel, units = "in")

# --- Gráfico 2: IRF del IPC, por régimen ------------------------------------ #
df_ipc4 <- bind_rows(
  data.frame(horizonte = 0:H, regimen = "Positivo",
             pe = LP.ipc.s1$pe[, "Positivo"], lb = LP.ipc.s1$lb[, "Positivo"], ub = LP.ipc.s1$ub[, "Positivo"]),
  data.frame(horizonte = 0:H, regimen = "Negativo",
             pe = LP.ipc.s1$pe[, "Negativo"], lb = LP.ipc.s1$lb[, "Negativo"], ub = LP.ipc.s1$ub[, "Negativo"])
)
df_ipc4$regimen <- factor(df_ipc4$regimen, levels = c("Positivo", "Negativo"))

p_ipc4 <- ggplot(df_ipc4, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = regimen), alpha = 0.35) +
  geom_line(aes(y = pe, color = regimen), linewidth = 0.9) +
  scale_color_manual(values = colores_regimen, name = NULL) +
  scale_fill_manual(values = colores_banda_regimen, name = NULL) +
  labs(x = "Meses", y = "Respuesta acumulada (nivel)",
       title = "IPC: shock positivo vs. negativo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10))
p_ipc4
ggsave(here("output", "irf_punto4_ipc_pw.pdf"), p_ipc4, width = ancho_panel, height = alto_panel, units = "in")

# --- Gráfico 3: ERPT, por régimen ------------------------------------------- #
df_erpt4 <- bind_rows(
  data.frame(horizonte = 0:H, regimen = "Positivo",
             pe = LP.ERPT.pw$pe[, "Positivo"], lb = LP.ERPT.pw$lb[, "Positivo"], ub = LP.ERPT.pw$ub[, "Positivo"]),
  data.frame(horizonte = 0:H, regimen = "Negativo",
             pe = LP.ERPT.pw$pe[, "Negativo"], lb = LP.ERPT.pw$lb[, "Negativo"], ub = LP.ERPT.pw$ub[, "Negativo"])
)
df_erpt4$regimen <- factor(df_erpt4$regimen, levels = c("Positivo", "Negativo"))

p_erpt4 <- ggplot(df_erpt4, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = regimen), alpha = 0.35) +
  geom_line(aes(y = pe, color = regimen), linewidth = 0.9) +
  scale_color_manual(values = colores_regimen, name = NULL) +
  scale_fill_manual(values = colores_banda_regimen, name = NULL) +
  coord_cartesian(ylim = c(-130, 130)) +   # ajustá el rango según lo que veas al probar
  labs(x = "Meses", y = "ERPT acumulado (%)",
       title = "ERPT: shock positivo vs. negativo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10))
p_erpt4
ggsave(here("output", "erpt_punto4_pw.pdf"), p_erpt4, width = ancho_panel, height = alto_panel, units = "in")


#------------------------------------------------------------------------------#
# Anexo: Comparación Spec 1 (interactuada) vs. Spec 2 (controles compartidos)
#------------------------------------------------------------------------------#

colores_spec <- c("Spec. 1 (interactuada)" = "darkorchid4", "Spec. 2 (compartida)" = "darkcyan")
colores_banda_spec <- c("Spec. 1 (interactuada)" = "plum2", "Spec. 2 (compartida)" = "paleturquoise2")

ancho_panel <- 3.6
alto_panel  <- 3.9

# --- Gráfico A1: TCN, régimen Positivo, Spec 1 vs. Spec 2 ------------------- #
df_tcn_pos <- bind_rows(
  data.frame(horizonte = 0:H, spec = "Spec. 1 (interactuada)",
             pe = LP.tcn.s1$pe[, "Positivo"], lb = LP.tcn.s1$lb[, "Positivo"], ub = LP.tcn.s1$ub[, "Positivo"]),
  data.frame(horizonte = 0:H, spec = "Spec. 2 (compartida)",
             pe = LP.tcn.s2$pe[, "Positivo"], lb = LP.tcn.s2$lb[, "Positivo"], ub = LP.tcn.s2$ub[, "Positivo"])
)
df_tcn_pos$spec <- factor(df_tcn_pos$spec, levels = c("Spec. 1 (interactuada)", "Spec. 2 (compartida)"))

p_tcn_pos <- ggplot(df_tcn_pos, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = spec), alpha = 0.35) +
  geom_line(aes(y = pe, color = spec), linewidth = 0.9) +
  scale_color_manual(values = colores_spec, name = NULL) +
  scale_fill_manual(values = colores_banda_spec, name = NULL) +
  labs(x = "Meses", y = "Respuesta acumulada (nivel)",
       title = "TCN, régimen positivo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10))
ggsave(here("output", "anexo_tcn_pos_spec.pdf"), p_tcn_pos, width = ancho_panel, height = alto_panel, units = "in")

# --- Gráfico A2: TCN, régimen Negativo, Spec 1 vs. Spec 2 ------------------- #
df_tcn_neg <- bind_rows(
  data.frame(horizonte = 0:H, spec = "Spec. 1 (interactuada)",
             pe = LP.tcn.s1$pe[, "Negativo"], lb = LP.tcn.s1$lb[, "Negativo"], ub = LP.tcn.s1$ub[, "Negativo"]),
  data.frame(horizonte = 0:H, spec = "Spec. 2 (compartida)",
             pe = LP.tcn.s2$pe[, "Negativo"], lb = LP.tcn.s2$lb[, "Negativo"], ub = LP.tcn.s2$ub[, "Negativo"])
)
df_tcn_neg$spec <- factor(df_tcn_neg$spec, levels = c("Spec. 1 (interactuada)", "Spec. 2 (compartida)"))

p_tcn_neg <- ggplot(df_tcn_neg, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = spec), alpha = 0.35) +
  geom_line(aes(y = pe, color = spec), linewidth = 0.9) +
  scale_color_manual(values = colores_spec, name = NULL) +
  scale_fill_manual(values = colores_banda_spec, name = NULL) +
  labs(x = "Meses", y = "Respuesta acumulada (nivel)",
       title = "TCN, régimen negativo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10))
ggsave(here("output", "anexo_tcn_neg_spec.pdf"), p_tcn_neg, width = ancho_panel, height = alto_panel, units = "in")

# --- Gráfico A3: IPC, régimen Positivo, Spec 1 vs. Spec 2 ------------------- #
df_ipc_pos <- bind_rows(
  data.frame(horizonte = 0:H, spec = "Spec. 1 (interactuada)",
             pe = LP.ipc.s1$pe[, "Positivo"], lb = LP.ipc.s1$lb[, "Positivo"], ub = LP.ipc.s1$ub[, "Positivo"]),
  data.frame(horizonte = 0:H, spec = "Spec. 2 (compartida)",
             pe = LP.ipc.s2$pe[, "Positivo"], lb = LP.ipc.s2$lb[, "Positivo"], ub = LP.ipc.s2$ub[, "Positivo"])
)
df_ipc_pos$spec <- factor(df_ipc_pos$spec, levels = c("Spec. 1 (interactuada)", "Spec. 2 (compartida)"))

p_ipc_pos <- ggplot(df_ipc_pos, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = spec), alpha = 0.35) +
  geom_line(aes(y = pe, color = spec), linewidth = 0.9) +
  scale_color_manual(values = colores_spec, name = NULL) +
  scale_fill_manual(values = colores_banda_spec, name = NULL) +
  labs(x = "Meses", y = "Respuesta acumulada (nivel)",
       title = "IPC, régimen positivo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10))
ggsave(here("output", "anexo_ipc_pos_spec.pdf"), p_ipc_pos, width = ancho_panel, height = alto_panel, units = "in")

# --- Gráfico A4: IPC, régimen Negativo, Spec 1 vs. Spec 2 ------------------- #
df_ipc_neg <- bind_rows(
  data.frame(horizonte = 0:H, spec = "Spec. 1 (interactuada)",
             pe = LP.ipc.s1$pe[, "Negativo"], lb = LP.ipc.s1$lb[, "Negativo"], ub = LP.ipc.s1$ub[, "Negativo"]),
  data.frame(horizonte = 0:H, spec = "Spec. 2 (compartida)",
             pe = LP.ipc.s2$pe[, "Negativo"], lb = LP.ipc.s2$lb[, "Negativo"], ub = LP.ipc.s2$ub[, "Negativo"])
)
df_ipc_neg$spec <- factor(df_ipc_neg$spec, levels = c("Spec. 1 (interactuada)", "Spec. 2 (compartida)"))

p_ipc_neg <- ggplot(df_ipc_neg, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = spec), alpha = 0.35) +
  geom_line(aes(y = pe, color = spec), linewidth = 0.9) +
  scale_color_manual(values = colores_spec, name = NULL) +
  scale_fill_manual(values = colores_banda_spec, name = NULL) +
  labs(x = "Meses", y = "Respuesta acumulada (nivel)",
       title = "IPC, régimen negativo") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 10))
ggsave(here("output", "anexo_ipc_neg_spec.pdf"), p_ipc_neg, width = ancho_panel, height = alto_panel, units = "in")