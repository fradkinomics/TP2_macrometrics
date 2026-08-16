#------------------------------------------------------------------------------#
# Maestría en Economía
# Macroeconometría
# 2023, 3er trimestre 
# Profesor: Javier García-Cicco
# Tutor: Franco Nuñez

# Material basado en código de Luis Libonatti (usado en versiones anteriores de 
# la materia)
#------------------------------------------------------------------------------#

library(here)
library(sandwich)

cumsum.h <-  function(x, h) {
  diff(cumsum(c(0, x)), h + 1)
}

get.shock <- function(Y, p, m, idx.s) {
  
  Xc <- embed(Y, p + 1)
  
  W1 <- Xc[, (m + 1):(m * (p + 1))]
  if (idx.s == 1) {
    W0 <- NULL
  } else {
    W0 <- Xc[, -c(idx.s:m)]
  }
  Wt <- cbind(1, W0, W1)
  
  st <- Xc[, idx.s]
  st <- unname(resid(lm(st ~ -1 + Wt))) # The -1 term removes the intercept from the regression equation
  
  list(st = st, Wt = Wt)
  
}

get.shock.pw <- function(Y, D, p, m, idx.s) {
  
  Dt <- D[(p + 1):length(D)]
  
  Xc <- embed(Y, p + 1)
  
  W1 <- Xc[, (m + 1):(m * (p + 1))]
  if (idx.s == 1) {
    W0 <- NULL
  } else {
    W0 <- Xc[, -c(idx.s:m)]
  }
  Wt <- cbind(1, W0, W1)
  Wt <- cbind(Dt * Wt, (1 - Dt) * Wt)
  
  St <- Xc[, idx.s]
  St <- cbind(Dt * St, (1 - Dt) * St)
  St <- unname(resid(lm(St ~ -1 + Wt)))
  
  list(St = St, Wt = Wt, Dt = Dt)
  
}

lp <- function(Y, p, idx.s, idx.r, H, gamma, cumulative = FALSE) {
  
  m <- ncol(Y)
  T <- nrow(Y)
  
  reg <- get.shock(Y, p, m, idx.s)
  
  rsp <- Y[(p + 1):T, idx.r]
  imp <- reg$st
  ctr <- reg$Wt
  
  irf.pe <- rep(NA, H + 1)
  irf.lb <- rep(NA, length(irf.pe))
  irf.ub <- rep(NA, length(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    if (isTRUE(cumulative)) {
      yh <- cumsum.h(rsp, h)
    } else {
      yh <- rsp[(h + 1):(T - p)]
    }
    
    Wt <- ctr[1:(T - p - h), ]
    st <- imp[1:(T - p - h)]
    
    yh <- resid(lm(yh ~ -1 + Wt)) # The -1 term removes the intercept from the regression equation
    
    projection <- lm(yh ~ -1 + st)
    
    b.pe <- unname(coef(projection))
    b.se <- unname(sqrt(diag(vcovHAC(projection))))
    
    irf.pe[h + 1] <- b.pe
    irf.lb[h + 1] <- b.pe - z * b.se
    irf.ub[h + 1] <- b.pe + z * b.se
    
  }
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
  
}

lp.pw <- function(Y, D, p, idx.s, idx.r, H, gamma, cumulative = FALSE) {
  
  m <- ncol(Y)
  T <- nrow(Y)
  
  reg <- get.shock.pw(Y, D, p, m, idx.s)
  
  rsp <- Y[(p + 1):T, idx.r]
  imp <- reg$St
  ctr <- reg$Wt
  
  irf.pe <- array(NA, c(H + 1, ncol(imp)))
  irf.lb <- array(NA, dim(irf.pe))
  irf.ub <- array(NA, dim(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    if (isTRUE(cumulative)) {
      yh <- cumsum.h(rsp, h)
    } else {
      yh <- rsp[(h + 1):(T - p)]
    }
    
    Wt <- ctr[1:(T - p - h), ]
    St <- imp[1:(T - p - h), ]
    
    yh <- resid(lm(yh ~ -1 + Wt)) # The -1 term removes the intercept from the regression equation
    
    projection <- lm(yh ~ -1 + St)
    
    b.pe <- unname(coef(projection))
    b.se <- unname(sqrt(diag(vcovHAC(projection))))
    
    irf.pe[h + 1, ] <- b.pe
    irf.lb[h + 1, ] <- b.pe - z * b.se
    irf.ub[h + 1, ] <- b.pe + z * b.se
    
  }
  
  dimnames(irf.pe) <- list(0:H, c("D", "1-D"))
  dimnames(irf.lb) <- dimnames(irf.pe)
  dimnames(irf.ub) <- dimnames(irf.pe)
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
  
}

lp.multiplier <- function(Y, p, idx.s, idx.rl, idx.rr, H, gamma) {
  
  m <- ncol(Y)
  T <- nrow(Y)
  
  reg <- get.shock(Y, p, m, idx.s)
  
  rsp.l <- Y[(p + 1):T, idx.rl]
  rsp.r <- Y[(p + 1):T, idx.rr]
  imp <- reg$st
  ctr <- reg$Wt
  
  irf.pe <- rep(NA, H + 1)
  irf.lb <- rep(NA, length(irf.pe))
  irf.ub <- rep(NA, length(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    yh.l <- cumsum.h(rsp.l, h)
    yh.r <- cumsum.h(rsp.r, h)
    
    Wt <- ctr[1:(T - p - h), ]
    st <- imp[1:(T - p - h)]
    
    yh.l <- resid(lm(yh.l ~ -1 + Wt)) # The -1 term removes the intercept from the regression equation
    yh.r <- resid(lm(yh.r ~ -1 + Wt))
    
    yh.r <- fitted(lm(yh.r ~ -1 + st))
    projection <- lm(yh.l ~ -1 + yh.r)
    
    b.pe <- unname(coef(projection))
    b.se <- unname(sqrt(diag(vcovHAC(projection))))
    
    irf.pe[h + 1] <- b.pe
    irf.lb[h + 1] <- b.pe - z * b.se
    irf.ub[h + 1] <- b.pe + z * b.se
    
  }
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
  
}

lp.multiplier.pw <- function(Y, D, p, idx.s, idx.rl, idx.rr, H, gamma) {
  
  m <- ncol(Y)
  T <- nrow(Y)
  
  reg <- get.shock.pw(Y, D, p, m, idx.s)
  
  rsp.l <- Y[(p + 1):T, idx.rl]
  rsp.r <- Y[(p + 1):T, idx.rr]
  imp <- reg$St
  ctr <- reg$Wt
  idx <- reg$Dt
  
  irf.pe <- array(NA, c(H + 1, ncol(imp)))
  irf.lb <- array(NA, dim(irf.pe))
  irf.ub <- array(NA, dim(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    yh.l <- cumsum.h(rsp.l, h)
    Yh.r <- cumsum.h(rsp.r, h)
    
    Dt <- idx[1:(T - p - h)]
    Wt <- ctr[1:(T - p - h), ]
    St <- imp[1:(T - p - h), ]
    
    yh.l <- resid(lm(yh.l ~ -1 + Wt)) # The -1 term removes the intercept from the regression equation
    Yh.r <- cbind(Dt * Yh.r, (1 - Dt) * Yh.r)
    Yh.r <- resid(lm(Yh.r ~ -1 + Wt))
    
    Yh.r <- cbind(fitted(lm(Yh.r[, 1] ~ -1 + St[, 1])), fitted(lm(Yh.r[, 2] ~ -1 + St[, 2])))
    projection <- lm(yh.l ~ -1 + Yh.r)
    
    b.pe <- unname(coef(projection))
    b.se <- unname(sqrt(diag(vcovHAC(projection))))
    
    irf.pe[h + 1, ] <- b.pe
    irf.lb[h + 1, ] <- b.pe - z * b.se
    irf.ub[h + 1, ] <- b.pe + z * b.se
    
  }
  
  dimnames(irf.pe) <- list(0:H, c("D", "1-D"))
  dimnames(irf.lb) <- dimnames(irf.pe)
  dimnames(irf.ub) <- dimnames(irf.pe)
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
  
}

#EXTRA: funciones para shocks ya observados

# Controles: intercepto + rezagos 1:p de todo Y (sin contemporáneos, 
# porque el shock ya viene identificado desde el SVAR)
get.controls <- function(Y, p, m) {
  Xc <- embed(Y, p + 1)
  W1 <- Xc[, (m + 1):(m * (p + 1))]
  cbind(1, W1)
}

# LP para una sola respuesta (IRF), shock observado
lp.shock <- function(Y, shock, p, idx.r, H, gamma, cumulative = FALSE) {
  
  m <- ncol(Y)
  T <- nrow(Y)
  
  Wt.full <- get.controls(Y, p, m)
  rsp <- Y[(p + 1):T, idx.r]
  
  irf.pe <- rep(NA, H + 1)
  irf.lb <- rep(NA, length(irf.pe))
  irf.ub <- rep(NA, length(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    if (isTRUE(cumulative)) {
      yh <- cumsum.h(rsp, h)
    } else {
      yh <- rsp[(h + 1):(T - p)]
    }
    
    Wt <- Wt.full[1:(T - p - h), ]
    st <- shock[1:(T - p - h)]
    
    yh <- resid(lm(yh ~ -1 + Wt))
    projection <- lm(yh ~ -1 + st)
    
    b.pe <- unname(coef(projection))
    b.se <- unname(sqrt(diag(vcovHAC(projection))))
    
    irf.pe[h + 1] <- b.pe
    irf.lb[h + 1] <- b.pe - z * b.se
    irf.ub[h + 1] <- b.pe + z * b.se
    
  }
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
  
}

# LP para el RATIO (ERPT), shock observado usado como "instrumento"
lp.multiplier.shock <- function(Y, shock, p, idx.rl, idx.rr, H, gamma) {
  
  m <- ncol(Y)
  T <- nrow(Y)
  
  Wt.full <- get.controls(Y, p, m)   # intercepto + rezagos 1:p de Y
  rsp.l <- Y[(p + 1):T, idx.rl]      # numerador (ej: dipc)
  rsp.r <- Y[(p + 1):T, idx.rr]      # denominador (ej: dtcn)
  
  irf.pe <- rep(NA, H + 1)
  irf.lb <- rep(NA, length(irf.pe))
  irf.ub <- rep(NA, length(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    # Respuestas acumuladas hasta el horizonte h (nivel, no diferencias)
    yh.l <- cumsum.h(rsp.l, h)
    yh.r <- cumsum.h(rsp.r, h)
    
    Wt <- Wt.full[1:(T - p - h), ]
    st <- shock[1:(T - p - h)]
    
    # Purga de controles (rezagos) en ambas variables respuesta
    yh.l <- resid(lm(yh.l ~ -1 + Wt))
    yh.r <- resid(lm(yh.r ~ -1 + Wt))
    
    # Proyección del denominador sobre el shock (equivalente a un "primer
    # paso" de 2SLS, usando el shock como instrumento)
    yh.r <- fitted(lm(yh.r ~ -1 + st))
    
    # Regresión final: numerador contra el fitted del denominador.
    # El coeficiente ES el ratio de IRFs acumuladas, es decir, el ERPT_h
    projection <- lm(yh.l ~ -1 + yh.r)
    
    b.pe <- 100 * unname(coef(projection))                  # <- % de traspaso
    b.se <- 100 * unname(sqrt(diag(vcovHAC(projection))))   # <- SE en la misma escala
    
    irf.pe[h + 1] <- b.pe
    irf.lb[h + 1] <- b.pe - z * b.se
    irf.ub[h + 1] <- b.pe + z * b.se
    
  }
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
  
}



#===============================================================================
#INCISO 8
#===============================================================================

# LP con LHS provisto externamente (inciso 8): la respuesta NO sale de Y,
# se pasa como vector y.lhs (largo T, alineado con las filas de Y).
# Los controles siguen siendo los rezagos 1:p de Y (idénticos al inciso 3).
lp.shock.y <- function(Y, shock, y.lhs, p, H, gamma) {
  
  m <- ncol(Y)
  T <- nrow(Y)
  
  Wt.full <- get.controls(Y, p, m)     # intercepto + rezagos 1:p de Y
  rsp <- y.lhs[(p + 1):T]              # LHS alineado a partir de la fila p+1
  
  irf.pe <- rep(NA, H + 1)
  irf.lb <- rep(NA, length(irf.pe))
  irf.ub <- rep(NA, length(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    yh <- rsp[(h + 1):(T - p)]         # LHS en t+h (ya construido, sin cumsum)
    Wt <- Wt.full[1:(T - p - h), ]
    st <- shock[1:(T - p - h)]
    
    yh <- resid(lm(yh ~ -1 + Wt))
    projection <- lm(yh ~ -1 + st)
    
    b.pe <- unname(coef(projection))
    b.se <- unname(sqrt(diag(vcovHAC(projection))))
    
    irf.pe[h + 1] <- b.pe
    irf.lb[h + 1] <- b.pe - z * b.se
    irf.ub[h + 1] <- b.pe + z * b.se
  }
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
}




#===============================================================================
# INCISO 9
#===============================================================================

# ERPT vía LP-IV con AMBAS respuestas provistas como vectores externos.
# y.num, y.den: LHS ya construidas (largo T, alineadas con las filas de Y).
# NO se acumula: y.num e y.den ya son tasas/niveles de 12m (esperadas o forward).
# Controles = rezagos 1:p de Y (idénticos al inciso 3). 2SLS a mano: el shock
# instrumenta el denominador; el coef. de 2da etapa ES el ERPT_h (×100).
lp.multiplier.shock.y <- function(Y, shock, y.num, y.den, p, H, gamma) {
  
  m <- ncol(Y)
  T <- nrow(Y)
  
  Wt.full <- get.controls(Y, p, m)   # intercepto + rezagos 1:p de Y
  rsp.l <- y.num[(p + 1):T]          # numerador  (ej: pi esperado / pi forward)
  rsp.r <- y.den[(p + 1):T]          # denominador (ej: dep esperada / dep forward)
  
  irf.pe <- rep(NA, H + 1)
  irf.lb <- rep(NA, length(irf.pe))
  irf.ub <- rep(NA, length(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    yh.l <- rsp.l[(h + 1):(T - p)]   # SIN cumsum: la LHS ya es 12m
    yh.r <- rsp.r[(h + 1):(T - p)]
    
    Wt <- Wt.full[1:(T - p - h), ]
    st <- shock[1:(T - p - h)]
    
    yh.l <- resid(lm(yh.l ~ -1 + Wt))
    yh.r <- resid(lm(yh.r ~ -1 + Wt))
    
    yh.r <- fitted(lm(yh.r ~ -1 + st))     # 1ra etapa: denominador ~ shock
    projection <- lm(yh.l ~ -1 + yh.r)     # 2da etapa: numerador ~ den.ajustado
    
    b.pe <- 100 * unname(coef(projection))
    b.se <- 100 * unname(sqrt(diag(vcovHAC(projection))))
    
    irf.pe[h + 1] <- b.pe
    irf.lb[h + 1] <- b.pe - z * b.se
    irf.ub[h + 1] <- b.pe + z * b.se
    
  }
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
  
}


#------------------------------------------------------------------------------#
# AGREGADO POR LOS AUTORES DEL TP2 (no forma parte del PS3_LP_Tools.R original de
# la catedra): variantes piecewise-por-signo que reciben el shock como vector.
# Usadas en los incisos 4 y 7. Fusionadas desde la version del companero de grupo.
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Inciso 4: ERPT asimétrico (ratio directo vía IV, por régimen de signo)
#------------------------------------------------------------------------------#
# --- Función: LP con shock observado, IRF asimétrica por signo -------------- #
# interact.controls = TRUE  -> Spec 1 (controles interactuados con D)
# interact.controls = FALSE -> Spec 2 (controles compartidos, solo el shock 
#                                        interactuado con D)
lp.pw.shock <- function(Y, shock, D, p, idx.r, H, gamma, 
                        cumulative = FALSE, interact.controls = TRUE) {
  
  m <- ncol(Y); T <- nrow(Y)
  Wt.full <- get.controls(Y, p, m)
  rsp <- Y[(p + 1):T, idx.r]
  
  irf.pe <- array(NA, c(H + 1, 2))
  irf.lb <- array(NA, dim(irf.pe))
  irf.ub <- array(NA, dim(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    yh <- if (isTRUE(cumulative)) cumsum.h(rsp, h) else rsp[(h + 1):(T - p)]
    Wt <- Wt.full[1:(T - p - h), ]
    st <- shock[1:(T - p - h)]
    Dt <- D[1:(T - p - h)]
    
    if (isTRUE(interact.controls)) {
      # Spec 1: todo interactuado (intercepto, shock, y rezagos) x régimen
      X.pos <- Dt * cbind(1, st, Wt)
      X.neg <- (1 - Dt) * cbind(1, st, Wt)
      reg <- lm(yh ~ -1 + cbind(X.pos, X.neg))
      k <- ncol(Wt) + 2
      idx.pos <- 2       # posición del shock dentro del bloque "positivo"
      idx.neg <- k + 2   # posición del shock dentro del bloque "negativo"
    } else {
      # Spec 2: solo el shock interactuado, rezagos compartidos
      X <- cbind(Dt * st, (1 - Dt) * st, Wt)
      reg <- lm(yh ~ -1 + X)
      idx.pos <- 1
      idx.neg <- 2
    }
    
    b <- coef(reg)
    se <- sqrt(diag(vcovHAC(reg)))
    
    b.pos <- b[idx.pos]; se.pos <- se[idx.pos]
    b.neg <- -b[idx.neg]; se.neg <- se[idx.neg]  # <- normalización: signo invertido
    
    irf.pe[h + 1, ] <- c(b.pos, b.neg)
    irf.lb[h + 1, ] <- c(b.pos - z * se.pos, b.neg - z * se.neg)
    irf.ub[h + 1, ] <- c(b.pos + z * se.pos, b.neg + z * se.neg)
  }
  
  dimnames(irf.pe) <- list(0:H, c("Positivo", "Negativo"))
  dimnames(irf.lb) <- dimnames(irf.pe)
  dimnames(irf.ub) <- dimnames(irf.pe)
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
}

lp.multiplier.pw.shock <- function(Y, shock, D, p, idx.rl, idx.rr, H, gamma) {
  
  m <- ncol(Y); T <- nrow(Y)
  Wt.full <- get.controls(Y, p, m)
  rsp.l <- Y[(p + 1):T, idx.rl]  # numerador (IPC)
  rsp.r <- Y[(p + 1):T, idx.rr]  # denominador (TCN)
  
  irf.pe <- array(NA, c(H + 1, 2))
  irf.lb <- array(NA, dim(irf.pe))
  irf.ub <- array(NA, dim(irf.pe))
  
  z <- qnorm(1 - (1 - gamma) / 2)
  
  for (h in 0:H) {
    
    yh.l <- cumsum.h(rsp.l, h)
    yh.r <- cumsum.h(rsp.r, h)
    
    Wt <- Wt.full[1:(T - p - h), ]
    st <- shock[1:(T - p - h)]
    Dt <- D[1:(T - p - h)]
    
    # Controles totalmente interactuados por régimen (Spec 1)
    X.pos <- Dt * cbind(1, Wt)
    X.neg <- (1 - Dt) * cbind(1, Wt)
    Xc <- cbind(X.pos, X.neg)
    
    yh.l <- resid(lm(yh.l ~ -1 + Xc))
    yh.r <- resid(lm(yh.r ~ -1 + Xc))
    
    # Proyección del TCN (residual) sobre el shock, también interactuado por régimen
    st.pos <- Dt * st
    st.neg <- (1 - Dt) * st
    yh.r.fit <- fitted(lm(yh.r ~ -1 + st.pos + st.neg))
    
    # Regresión final: coeficientes separados por régimen (el ratio ERPT no
    # requiere invertir signo en el régimen negativo, porque numerador y
    # denominador escalan con el mismo signo del shock y el cociente se cancela)
    proj <- lm(yh.l ~ -1 + I(Dt * yh.r.fit) + I((1 - Dt) * yh.r.fit))
    
    b <- 100 * unname(coef(proj))
    se <- 100 * unname(sqrt(diag(vcovHAC(proj))))
    
    irf.pe[h + 1, ] <- b
    irf.lb[h + 1, ] <- b - z * se
    irf.ub[h + 1, ] <- b + z * se
  }
  
  dimnames(irf.pe) <- list(0:H, c("Positivo", "Negativo"))
  dimnames(irf.lb) <- dimnames(irf.pe)
  dimnames(irf.ub) <- dimnames(irf.pe)
  
  list(lb = irf.lb, pe = irf.pe, ub = irf.ub)
}


