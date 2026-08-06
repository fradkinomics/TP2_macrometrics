#------------------------------------------------------------------------------#
# Maestría en Economía
# Macroeconometría
# 2025, 2do trimestre 
# Profesor: Javier Garcia-Cicco
# Tutor: Franco Nuñez

# Material basado en código de Luis Libonatti (usado en versiones anteriores de 
# la materia)
#------------------------------------------------------------------------------#

wash <- function(X) {
  unclass(as.matrix(unname(X)))[]
}

VAR.coefficients <- function(VAR) {
  B <- wash(Bcoef(VAR))
  list(Pi = B[, m * p + 1], Phi = B[, 1:(m * p)])
}