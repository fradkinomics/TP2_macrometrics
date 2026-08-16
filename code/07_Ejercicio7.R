#------------------------------------------------------------------------------#
# Universidad de San Andrés - Macroeconometría 2026
# TP2: Exchange Rate Pass-Through - Colombia
#
# Inciso 7: ERPT asimétrico por signo del shock, para los 6 shocks del VAR7
#   (inciso 6). Extensión piecewise del LP-IV: para cada shock estructural se
#   separa el efecto según el signo (positivo vs. negativo) del shock.
#
# Este script NO reestima el VAR7 del inciso 6: LEE los objetos ya guardados
# desde 'data/objetos_06.rds' (struct.shocks, Y, p).
#
# Herramientas de la cátedra (viven en 'tools/'):
#   - PS3_LP_Tools.R : get.controls(), cumsum.h(), lp.multiplier.pw.shock(), HAC
#------------------------------------------------------------------------------#

remove(list = ls(all.names = TRUE))
gc()

library(here)
library(sandwich)   # vcovHAC, usado internamente por las funciones LP
library(ggplot2)
library(dplyr)
library(patchwork)

#------------------------------------------------------------------------------#
# 0. Carga de objetos del inciso 6 (NO se recalcula el VAR7)
#------------------------------------------------------------------------------#

obj06 <- readRDS(here("data", "objetos_06.rds"))
for (nm in names(obj06)) assign(nm, obj06[[nm]])
remove(obj06)

p <- unname(p)   # p venía como named int; lo dejamos como escalar limpio

#------------------------------------------------------------------------------#
# 1. Herramientas LP del curso (tools) y parámetros del inciso
#------------------------------------------------------------------------------#

source(here("tools", "PS3_LP_Tools.R"))

H.ERPT <- 36
gamma  <- 0.95
shocks.nombres <- c("Shock Fed (JK)", "Shock Petróleo (BH)", "Excess Bond Premium",
                    "Términos de Intercambio", "EMBI", "TCN")


#------------------------------------------------------------------------------#
# ERPT condicional por signo del shock, para los 6 shocks
#------------------------------------------------------------------------------#

H <- H.ERPT

ERPT.pw <- vector("list", 6)
names(ERPT.pw) <- shocks.nombres

for (s in 1:6) {
  shock.s <- struct.shocks[, s]
  D.s <- ifelse(shock.s > 0, 1, 0)
  
  ERPT.pw[[s]] <- lp.multiplier.pw.shock(Y, shock.s, D.s, p,
                                         idx.rl = 7, idx.rr = 6,
                                         H, gamma)
}

#------------------------------------------------------------------------------#
# Diagnóstico: composición de la muestra por signo del shock
#------------------------------------------------------------------------------#

cat("--- Composición de la muestra por signo del shock ---\n")
for (s in 1:6) {
  shock.s <- struct.shocks[, s]
  n.pos <- sum(shock.s > 0)
  n.neg <- sum(shock.s <= 0)
  cat(shocks.nombres[s], ": positivos =", n.pos, "| negativos =", n.neg, "\n")
}

#------------------------------------------------------------------------------#
# Gráfico: panel de ERPT por shock, con las dos líneas (positivo/negativo),
# eje y acotado por panel (zoom, sin descartar datos) via patchwork
#------------------------------------------------------------------------------#

library(ggplot2)
library(dplyr)
library(patchwork)

h <- 0:H

df_erpt7 <- lapply(1:6, function(s) {
  bind_rows(
    data.frame(horizonte = h, shock = shocks.nombres[s], regimen = "Positivo",
               pe = ERPT.pw[[s]]$pe[, "Positivo"],
               lb = ERPT.pw[[s]]$lb[, "Positivo"],
               ub = ERPT.pw[[s]]$ub[, "Positivo"]),
    data.frame(horizonte = h, shock = shocks.nombres[s], regimen = "Negativo",
               pe = ERPT.pw[[s]]$pe[, "Negativo"],
               lb = ERPT.pw[[s]]$lb[, "Negativo"],
               ub = ERPT.pw[[s]]$ub[, "Negativo"])
  )
}) %>% bind_rows()

df_erpt7$shock <- factor(df_erpt7$shock, levels = shocks.nombres)
df_erpt7$regimen <- factor(df_erpt7$regimen, levels = c("Positivo", "Negativo"))

colores_regimen <- c("Positivo" = "darkgreen", "Negativo" = "darkorange")
colores_banda_regimen <- c("Positivo" = "palegreen3", "Negativo" = "navajowhite2")

# Límites de zoom por panel: mediana +/- 4*IQR combinando ambos regímenes,
# para que las dos líneas del mismo shock compartan escala
limites7 <- df_erpt7 %>%
  group_by(shock) %>%
  summarise(
    centro = median(pe, na.rm = TRUE),
    escala = IQR(c(pe, lb, ub), na.rm = TRUE),
    ymin = centro - 4 * escala,
    ymax = centro + 4 * escala,
    .groups = "drop"
  )

hacer_panel7 <- function(shock.nombre) {
  lims <- limites7 %>% dplyr::filter(shock == shock.nombre)
  
  ggplot(df_erpt7 %>% dplyr::filter(shock == shock.nombre), aes(x = horizonte)) +
    geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = regimen), alpha = 0.30) +
    geom_line(aes(y = pe, color = regimen), linewidth = 0.9) +
    scale_color_manual(values = colores_regimen, name = NULL) +
    scale_fill_manual(values = colores_banda_regimen, name = NULL) +
    coord_cartesian(ylim = c(lims$ymin, lims$ymax)) +
    labs(x = NULL, y = NULL, title = shock.nombre) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          panel.grid.minor = element_blank(),
          legend.position = "none")
}

paneles7 <- lapply(shocks.nombres, hacer_panel7)

# Leyenda compartida (se extrae de un panel auxiliar con legend.position="bottom")
panel_legenda <- hacer_panel7(shocks.nombres[1]) +
  theme(legend.position = "bottom")
legenda <- cowplot::get_legend(panel_legenda)

p_erpt7 <- wrap_plots(paneles7, nrow = 2) /
  wrap_elements(legenda) +
  plot_layout(heights = c(10, 1)) +
  plot_annotation(caption = "Meses (eje x) | ERPT condicional (%) (eje y)")

p_erpt7

if (!dir.exists(here("output"))) dir.create(here("output"))
ggsave(here("output", "erpt_punto7_asimetrico.pdf"), p_erpt7,
       width = 10, height = 7, units = "in")

#------------------------------------------------------------------------------#
# Tabla resumen: ERPT de largo plazo (h=36) por shock y régimen
#------------------------------------------------------------------------------#

tabla.erpt7 <- data.frame(
  Shock          = shocks.nombres,
  ERPT_pos_pe = sapply(ERPT.pw, function(x) round(x$pe[H + 1, "Positivo"], 2)),
  ERPT_pos_lb = sapply(ERPT.pw, function(x) round(x$lb[H + 1, "Positivo"], 2)),
  ERPT_pos_ub = sapply(ERPT.pw, function(x) round(x$ub[H + 1, "Positivo"], 2)),
  ERPT_neg_pe = sapply(ERPT.pw, function(x) round(x$pe[H + 1, "Negativo"], 2)),
  ERPT_neg_lb = sapply(ERPT.pw, function(x) round(x$lb[H + 1, "Negativo"], 2)),
  ERPT_neg_ub = sapply(ERPT.pw, function(x) round(x$ub[H + 1, "Negativo"], 2))
)
print(tabla.erpt7, row.names = FALSE)

write.csv(tabla.erpt7, here("output", "tabla_erpt_punto7.csv"), row.names = FALSE)

Fstat.pw <- function(shock, D, Y, p, idx.rr, H) {
  m <- ncol(Y); Tn <- nrow(Y)
  Wt.full <- get.controls(Y, p, m)
  rsp.r <- Y[(p + 1):Tn, idx.rr]
  
  Fs <- matrix(NA, H + 1, 2)
  colnames(Fs) <- c("Positivo", "Negativo")
  
  for (h in 0:H) {
    n <- length(shock) - h
    Wt <- Wt.full[1:n, ]
    st <- shock[1:n]
    Dt <- D[1:n]
    
    yh.r <- cumsum.h(rsp.r, h)
    yh.r <- resid(lm(yh.r ~ -1 + Wt))
    
    for (reg in c(1, 0)) {
      idx <- Dt == reg
      fs <- summary(lm(yh.r[idx] ~ -1 + st[idx]))$fstatistic[1]
      Fs[h + 1, ifelse(reg == 1, 1, 2)] <- fs
    }
  }
  Fs
}

Fstat7.resumen <- sapply(1:6, function(s) {
  shock.s <- struct.shocks[, s]
  D.s <- ifelse(shock.s > 0, 1, 0)
  Fs <- Fstat.pw(shock.s, D.s, Y, p, idx.rr = 6, H)
  round(Fs[H + 1, ], 2)   # F en h=36, para Positivo y Negativo
})
colnames(Fstat7.resumen) <- shocks.nombres
print(t(Fstat7.resumen))


