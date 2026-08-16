#------------------------------------------------------------------------------#
# TP2 Macroeconometría 2026 — Traspaso cambiario (ERPT), Colombia
# INCISO 9: ERPT esperado vs. estimado (realizado forward), un solo shock (TCN)
# No recalcula: carga objetos del inciso 8 (04_LP_Expectativas.R)
#------------------------------------------------------------------------------#

remove(list = ls()); gc()

library(here)
library(sandwich)
library(ggplot2)

source(here("tools", "PS3_LP_Tools.R"))

# Objetos del inciso 8 (ajustá el nombre del .rds si guardaste con otro)
objetos_04 <- readRDS(here("data", "objetos_04.rds"))
list2env(objetos_04, envir = .GlobalEnv)
# esperados: Y.e, u.e, einf.e, etcn.e, ipc.fwd.e, tcn.fwd.e (y las lp.* del 8)



#---( Parámetros — mismos que inciso 8 )---#
p     <- 1
H     <- 24
gamma <- 0.95

#===============================================================================
#---( ERPT esperado: pi_e / dep_e )---#
erpt.esp <- lp.multiplier.shock.y(
  Y     = Y.e,
  shock = u.e,
  y.num = einf.e,   # numerador: inflación esperada 12m
  y.den = etcn.e,   # denominador: depreciación esperada 12m
  p = p, H = H, gamma = gamma
)

#---( ERPT estimado: pi_fwd / dep_fwd )---#
erpt.real <- lp.multiplier.shock.y(
  Y     = Y.e,
  shock = u.e,
  y.num = ipc.fwd.e,  # numerador: inflación realizada forward 12m
  y.den = tcn.fwd.e,  # denominador: depreciación realizada forward 12m
  p = p, H = H, gamma = gamma
)

#---( Inspección rápida en consola )---#
cat("ERPT esperado (h=0,6,12,24):",
    round(erpt.esp$pe[c(1,7,13,25)], 1), "\n")
cat("ERPT estimado (h=0,6,12,24):",
    round(erpt.real$pe[c(1,7,13,25)], 1), "\n")


#==============================================================================
#---( F de primera etapa por horizonte: denominador ~ shock )---#
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
    Fh[h + 1] <- summary(fit)$fstatistic[1]   # F del shock en la 1ra etapa
  }
  Fh
}

F.esp  <- first.stage.F(Y.e, u.e, etcn.e,    p, H)   # denominador esperado
F.real <- first.stage.F(Y.e, u.e, tcn.fwd.e, p, H)   # denominador realizado fwd

cat("F 1ra etapa — esperado (mediana):", round(median(F.esp), 1),
    " | realizado (mediana):", round(median(F.real), 1), "\n")
cat("F esperado (h=0,6,12,24):", round(F.esp[c(1,7,13,25)], 1), "\n")
cat("F realizado (h=0,6,12,24):", round(F.real[c(1,7,13,25)], 1), "\n")




#===============================================================================
#---( Gráfico ERPT esperado vs. estimado, eje recortado )---#
ylim_erpt <- c(-50, 50)   # ventana legible; valores fuera quedan clipeados, no borrados

df_erpt <- bind_rows(
  data.frame(horizonte = 0:H, panel = "ERPT esperado",
             pe = erpt.esp$pe,  lb = erpt.esp$lb,  ub = erpt.esp$ub,
             serie = "Esperada"),
  data.frame(horizonte = 0:H, panel = "ERPT estimado (realizado fwd)",
             pe = erpt.real$pe, lb = erpt.real$lb, ub = erpt.real$ub,
             serie = "Realizada")
)
df_erpt$panel <- factor(df_erpt$panel,
                        levels = c("ERPT esperado", "ERPT estimado (realizado fwd)"))

colores_serie <- c("Esperada" = "#2C5F8A", "Realizada" = "#C2410C")

p_erpt <- ggplot(df_erpt, aes(x = horizonte)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.4, linetype = "dashed") +
  geom_ribbon(aes(ymin = lb, ymax = ub, fill = serie), alpha = 0.35) +
  geom_line(aes(y = pe, color = serie), linewidth = 0.9) +
  scale_color_manual(values = colores_serie, name = NULL) +
  scale_fill_manual(values = colores_serie, name = NULL) +
  labs(x = "Meses", y = "ERPT acumulado (%)") +
  coord_cartesian(ylim = ylim_erpt) +          # <- recorta la VISTA, no los datos
  facet_wrap(~ panel) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 11))
p_erpt

#---( Guardado )---#
ggsave(file.path("output", "09_erpt_esp_vs_real.pdf"),
       p_erpt, width = 8, height = 5, units = "in")

saveRDS(list(erpt.esp = erpt.esp, erpt.real = erpt.real),
        here("data", "objetos_05.rds"))