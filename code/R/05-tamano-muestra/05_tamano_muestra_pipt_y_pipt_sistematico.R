
# Determinación del tamaño de muestra
# Muestreo πPT (Brewer) y πPT sistemático sin reemplazo
# Simulación con selección paralela de muestras

# Configuración de salida
options(width = 150)

# Marco muestral y ordenación de las unidades
# Cargar los datos y establecer un orden fijo para el diseño sistemático
library(data.table)
library(sampling)
library(gt)
library(parallel)

poblacion = readRDS("datos/poblacion_educacion_2025.rds")

marco = copy(poblacion[, .(
  id_unidad, departamento, municipio, sector, cantidad_sedes
)])

setorder(marco, departamento, municipio, sector, cantidad_sedes, id_unidad)
N = nrow(marco) ; N

x = marco$cantidad_sedes

marco[1:10]

# Simulación de la población de planificación
# Generar matrícula relacionada con la cantidad de sedes y dos indicadores binarios
# Valores supuestos de una encuesta anterior
media = 400
S2 = 40000
P = c(0.50, 0.01)

# Simular las variables de interés
set.seed(123)
matricula = rgamma(N, shape = media^2 / S2, scale = S2 / media) * (x / mean(x))

indicador_50 = integer(N)
indicador_01 = integer(N)
indicador_50[sample.int(N, round(N * P[1]))] = 1L
indicador_01[sample.int(N, round(N * P[2]))] = 1L

poblacion_simulada = copy(marco)
poblacion_simulada[, `:=`(
  matricula = matricula,
  indicador_50 = indicador_50,
  indicador_01 = indicador_01
)]

parametros = poblacion_simulada[, .(
  Media = mean(matricula),
  Varianza = var(matricula),
  P50 = mean(indicador_50),
  P01 = mean(indicador_01)
)]

parametros

# Requisitos de precisión y simulación paralela
# Comparar la precisión de ambos diseños para cada tamaño candidato
# Requisitos de precisión
z = qnorm(0.975)
d = c(25, 0.05, 0.005)
c0 = c(0.05, 0.10, 0.20)
c_log = 0.10
B = 100 # 500

Y = as.matrix(poblacion_simulada[, .(matricula, indicador_50, indicador_01)])
valores = colMeans(Y)
p = min(valores[3], 1 - valores[3])

# Tamaños candidatos
n_candidatos = c(100, 200, 400, 800, 1200, 1500)
evaluacion = data.table()

# Detectar los núcleos disponibles y reservar uno para el sistema
nucleos_disponibles = detectCores(logical = FALSE)
n_nucleos = max(1L, nucleos_disponibles - 1L)

cat("Núcleos detectados:", nucleos_disponibles,
    "\nNúcleos utilizados:", n_nucleos, "\n")

cluster = makeCluster(n_nucleos)
invisible(clusterEvalQ(cluster, library(sampling)))
clusterExport(cluster, c("Y", "N"))
clusterSetRNGStream(cluster, 456)

for (n in n_candidatos) {
  pik = inclusionprobabilities(x, n)
  clusterExport(cluster, "pik")

  estimaciones_pipt = do.call(cbind, parLapply(cluster, seq_len(B), function(b) {
    seleccion = UPbrewer(pik) == 1
    colSums(Y[seleccion, , drop = FALSE] / pik[seleccion]) / N
  }))

  estimaciones_pipt_sis = do.call(cbind, parLapply(cluster, seq_len(B), function(b) {
    seleccion = UPsystematic(pik) == 1
    colSums(Y[seleccion, , drop = FALSE] / pik[seleccion]) / N
  }))

  for (diseno in c("πPT", "πPT sistemático")) {
    estimaciones = if (diseno == "πPT") estimaciones_pipt else estimaciones_pipt_sis

    EE = apply(estimaciones, 1, sd)
    margen = z * EE
    CV = EE / abs(valores)
    CV_log = EE[3] / (p * (-log(p)))

    evaluacion = rbind(evaluacion, data.table(
      Diseno = diseno,
      n = n,
      Margen_media = margen[1],
      Margen_P50 = margen[2],
      Margen_P01 = margen[3],
      CV_media = CV[1],
      CV_P50 = CV[2],
      CV_P01 = CV[3],
      CV_log_P01 = CV_log
    ))
  }
}

stopCluster(cluster)

# Verificar los requisitos de precisión
evaluacion[, Cumple :=
  Margen_media <= d[1] &
  Margen_P50 <= d[2] &
  #Margen_P01 <= d[3] &
  CV_media <= c0[1] &
  CV_P50 <= c0[2] &
  #CV_P01 <= c0[3] &
  CV_log_P01 <= c_log
]

# Resultados: muestreo πPT (Brewer)
# Presentar los márgenes y coeficientes de variación para πPT
evaluacion[Diseno == "πPT"][, Diseno := NULL] |>
  gt() |>
  fmt_number(columns = c(Margen_media, Margen_P50, Margen_P01,
                         CV_media, CV_P50, CV_P01, CV_log_P01), decimals = 4)

# Resultados: muestreo πPT sistemático
# Presentar los márgenes y coeficientes de variación para πPT sistemático
evaluacion[Diseno == "πPT sistemático"][, Diseno := NULL] |>
  gt() |>
  fmt_number(columns = c(Margen_media, Margen_P50, Margen_P01,
                         CV_media, CV_P50, CV_P01, CV_log_P01), decimals = 4)

# Determinación del tamaño de muestra
# Seleccionar el menor tamaño candidato que cumple los requisitos
n_final = evaluacion[, .(
  n_final = if (any(Cumple)) min(n[Cumple]) else NA_integer_
), by = Diseno]

n_final

# Nota: en la verificación original están desactivados los requisitos
# de margen absoluto y CV ordinario para la proporción de 0.01.
