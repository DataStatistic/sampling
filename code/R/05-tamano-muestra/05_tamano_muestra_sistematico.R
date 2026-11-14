
# Determinación del tamaño de muestra
# Muestreo sistemático sin reemplazo mediante simulación

# Librerías necesarias
library(data.table)

# Marco muestral y ordenación
poblacion = readRDS("datos/poblacion_educacion_2025.rds")

marco = copy(poblacion[, .(
  id_unidad, departamento, municipio, sector, cantidad_sedes
)])

setorder(marco, departamento, municipio, sector, cantidad_sedes, id_unidad)
N = nrow(marco) ; N

# Valores supuestos de una encuesta anterior
media = 400
S2 = 40000
P = c(0.50, 0.01)

# Simular las variables de interés
set.seed(123)
matricula = rgamma(N, shape = media^2 / S2, scale = S2 / media)

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

# Parámetros de la población simulada
parametros = poblacion_simulada[, .(
  Media = mean(matricula),
  Varianza = var(matricula),
  P50 = mean(indicador_50),
  P01 = mean(indicador_01)
)]

parametros

# Requisitos de precisión
z = qnorm(0.975)
d = c(25, 0.05, 0.005)
c0 = c(0.05, 0.10, 0.20)
c_log = 0.10
B = 500

Y = as.matrix(poblacion_simulada[, .(matricula, indicador_50, indicador_01)])
valores = colMeans(Y)
p = min(valores[3], 1 - valores[3])

# Tamaños candidatos
n_candidatos = c(50, 100, 200, 400, 800, 1200, 1500, 2000, 2500, 3000)
evaluacion = data.table()

# Evaluar la precisión para cada tamaño mediante simulación
set.seed(456)
for (n in n_candidatos) {
  estimaciones = replicate(B, {
    inicio = runif(1)
    posiciones = floor(((seq_len(n) - 1) + inicio) * N / n) + 1L
    colMeans(Y[posiciones, , drop = FALSE])
  })

  EE = apply(estimaciones, 1, sd)
  margen = z * EE
  CV = EE / abs(valores)
  CV_log = EE[3] / (p * (-log(p)))

  evaluacion = rbind(evaluacion, data.table(
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

# Verificar los requisitos de precisión
evaluacion[, Cumple :=
  Margen_media <= d[1] &
  Margen_P50 <= d[2] &
  Margen_P01 <= d[3] &
  CV_media <= c0[1] &
  CV_P50 <= c0[2] &
  CV_P01 <= c0[3] &
  CV_log_P01 <= c_log
]

# Resultados para todos los tamaños candidatos
evaluacion

# Menor tamaño candidato que cumple los requisitos
n_final = evaluacion[Cumple == TRUE, min(n)]
n_final
