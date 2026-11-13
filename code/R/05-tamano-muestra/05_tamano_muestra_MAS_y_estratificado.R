
# Determinación del tamaño de muestra
# MAS sin reemplazo y MAS estratificado con asignación proporcional
# Valores de planificación supuestos de una encuesta anterior

# Librerías necesarias
library(data.table)

# Información previa y requisitos de precisión
N = 17872
z = qnorm(0.975)
media = 400
S2 = 40000
P = c(0.50, 0.01)

parametros = c("Media", "Proporción 0.50", "Proporción 0.01")
valores = c(media, P)
d = c(25, 0.05, 0.005)
c0 = c(0.05, 0.10, 0.20)
c_log = 0.10

# Muestreo aleatorio simple sin reemplazo
varianzas_mas = c(S2, N / (N - 1) * P * (1 - P))

# Tamaños según el margen de error absoluto
n_margen_mas = ceiling(
  N * z^2 * varianzas_mas / (N * d^2 + z^2 * varianzas_mas)
)

# Tamaños según el CV ordinario
n_cv_mas = ceiling(
  N * varianzas_mas / (N * c0^2 * valores^2 + varianzas_mas)
)

# Tamaño según el CV logarítmico para P = 0.01
p = min(P[2], 1 - P[2])
limite_ee_log = c_log * p * (-log(p))

n_log_mas = ceiling(
  N * varianzas_mas[3] / (N * limite_ee_log^2 + varianzas_mas[3])
)

# Tamaño final y precisión prevista
n_final_mas = max(n_margen_mas, n_cv_mas, n_log_mas)

EE_mas = sqrt((1 - n_final_mas / N) * varianzas_mas / n_final_mas)
margen_mas = z * EE_mas
CV_mas = EE_mas / abs(valores)
CV_log_mas = EE_mas[3] / (p * (-log(p)))

resultados_mas = data.table(
  Parametro = parametros,
  n_margen = n_margen_mas,
  n_CV = n_cv_mas,
  Margen_previsto = margen_mas,
  CV_previsto = CV_mas
)

resultados_mas
n_log_mas
CV_log_mas
n_final_mas

# Muestreo aleatorio simple estratificado
# Valores supuestos de una encuesta anterior
estratos = data.table(
  Estrato = c("Caribe", "Oriental", "Central", "Pacífica",
              "Orinoquía, Amazonía e Insular", "Bogotá D. C."),
  N_h = c(3100, 2800, 4200, 3600, 1172, 3000),
  S2_h = c(25000, 30000, 40000, 30000, 50000, 35000),
  media_h = c(320, 360, 400, 450, 280, 480),
  P50_h = c(0.38, 0.47, 0.55, 0.62, 0.42, 0.58),
  P01_h = c(0.004, 0.007, 0.012, 0.015, 0.006, 0.016)
)

estratos[, W_h := N_h / N]
estratos[, media_h := media_h + (media - sum(W_h * media_h))]
estratos[, P50_h := P50_h + (P[1] - sum(W_h * P50_h))]
estratos[, P01_h := P01_h * (P[2] / sum(W_h * P01_h))]

# Varianzas previstas y asignación proporcional
varianzas_h = cbind(
  estratos$S2_h,
  estratos$N_h / (estratos$N_h - 1) * estratos$P50_h * (1 - estratos$P50_h),
  estratos$N_h / (estratos$N_h - 1) * estratos$P01_h * (1 - estratos$P01_h)
)

A = colSums(estratos$W_h * varianzas_h)

# Tamaños según el margen de error absoluto
n_margen_est = ceiling(N * z^2 * A / (N * d^2 + z^2 * A))

# Tamaños según el CV ordinario
n_cv_est = ceiling(N * A / (N * c0^2 * valores^2 + A))

# Tamaño según el CV logarítmico para P = 0.01
n_log_est = ceiling(N * A[3] / (N * limite_ee_log^2 + A[3]))

# Función de asignación proporcional con tamaños enteros
asignar_proporcional = function(n, N_h) {
  objetivo = n * N_h / sum(N_h)
  asignacion = floor(objetivo)
  faltantes = n - sum(asignacion)
  if (faltantes > 0) {
    posiciones = order(objetivo - asignacion, decreasing = TRUE)[seq_len(faltantes)]
    asignacion[posiciones] = asignacion[posiciones] + 1L
  }
  as.integer(asignacion)
}

# Tamaño final y verificación de la precisión
n_final_est = max(n_margen_est, n_cv_est, n_log_est)

repeat {
  estratos[, n_h := asignar_proporcional(n_final_est, N_h)]

  V_est = colSums(
    (estratos$W_h^2 * (1 - estratos$n_h / estratos$N_h) /
       estratos$n_h) * varianzas_h
  )

  EE_est = sqrt(V_est)
  margen_est = z * EE_est
  CV_est = EE_est / abs(valores)
  CV_log_est = EE_est[3] / (p * (-log(p)))

  if (all(margen_est <= d) && all(CV_est <= c0) && CV_log_est <= c_log) break
  n_final_est = n_final_est + 1L
}

resultados_est = data.table(
  Parametro = parametros,
  n_margen = n_margen_est,
  n_CV = n_cv_est,
  Margen_previsto = margen_est,
  CV_previsto = CV_est
)

resultados_est
n_log_est
CV_log_est
n_final_est
estratos[, .(Estrato, N_h, n_h)]

# Comparación de tamaños de muestra
comparacion = data.table(
  Requisito = c("Media: margen", "P = 0.50: margen", "P = 0.01: margen",
                "Media: CV", "P = 0.50: CV", "P = 0.01: CV",
                "P = 0.01: CV logarítmico", "Tamaño final"),
  MAS = c(n_margen_mas, n_cv_mas, n_log_mas, n_final_mas),
  MAS_estratificado = c(n_margen_est, n_cv_est, n_log_est, n_final_est)
)

comparacion
