
# Determinación del tamaño de muestra: muestreo en dos etapas πPT + MAS
# Población de referencia: viviendas del rural disperso del RMW 2024 (DANE).
# Primera etapa: TRI seleccionados mediante πPT (Brewer).
# Segunda etapa: viviendas seleccionadas mediante MAS sin reemplazo.
# Variables de interés simuladas: Y uniforme y proporciones de 0.50 y 0.01.
# Se evalúan las combinaciones candidatas mediante 500 replicaciones.

options(width = 150)

# Carga de paquetes y construcción del marco de TRI
library(data.table)
library(sampling)
library(gt)
library(parallel)

# Estructura de TRI y viviendas
poblacion = fread(
  "datos/rmw_vivienda.csv",
  colClasses = list(character = c("COD_CUESTIONARIO", "COD_TRI"))
)
poblacion = poblacion[U_CLASE == 3 & !is.na(COD_TRI) & COD_TRI != ""]

# Marco de TRI
marco_tri = poblacion[, .(
  M_i = .N,
  municipio = first(U_MPIO),
  nombre_tri = first(NOM_TRI)
), by = COD_TRI]

K = nrow(marco_tri)
N = nrow(poblacion)
M = marco_tri$M_i

# Correspondencia de cada vivienda con su TRI
id_tri = match(poblacion$COD_TRI, marco_tri$COD_TRI)
indices_tri = split(seq_len(N), factor(id_tri, levels = seq_len(K)))

if (!identical(as.integer(lengths(indices_tri)), as.integer(M))) {
  stop("Los tamaños de TRI no coinciden con el marco.")
}

data.table(Concepto = c("TRI", "Viviendas"), Total = c(K, N))

# Generación de las variables de interés
# Generación de la población de interés
set.seed(123)
P = c(0.50, 0.01)

y = runif(N, min = 0, max = 100)

indicador_50 = integer(N)
indicador_50[sample.int(N, round(N * P[1]))] = 1L

indicador_01 = integer(N)
indicador_01[sample.int(N, round(N * P[2]))] = 1L

Y = cbind(y = y, indicador_50 = indicador_50, indicador_01 = indicador_01)
valores = colMeans(Y)
n_variables = ncol(Y)

# Parámetros de la población simulada
resumen_poblacion = data.table(
  Parametro = c("Media de Y", "Proporción P = 0.50", "Proporción P = 0.01"),
  Valor = as.numeric(valores)
)

resumen_poblacion |>
  gt() |>
  fmt_number(columns = Valor, decimals = 4)

# Requisitos de precisión y tamaños candidatos
# Requisitos de precisión
z = qnorm(0.975)
d = c(y = 5, indicador_50 = 0.05, indicador_01 = 0.005)
c0 = c(y = 0.10, indicador_50 = 0.10)
c_log = 0.10

# Combinaciones de tamaños de primera y segunda etapa
n1_candidatos = c(10L, 20L, 30L, 40L, 60L)
n1_candidatos = n1_candidatos[n1_candidatos <= K]
n2_candidatos = c(200L, 400L, 600L, 800L, 1200L,
                  1400L, 1600L, 1800L, 2000L)

# Límites de viviendas seleccionadas en cada TRI
minimo = 10L
maximo = 150L

# Número de replicaciones por combinación
B = 500L

# Función de asignación proporcional de segunda etapa
# Asignación aproximadamente proporcional con límites por TRI
asignar_segunda_etapa = function(N2, n2_total, minimo = 10L, maximo = 150L) {
  N2 = as.integer(N2)
  limite_inferior = pmin(minimo, N2)
  limite_superior = pmin(maximo, N2)

  # El total es un objetivo, no una restricción de igualdad
  n2_total = min(max(n2_total, sum(limite_inferior)), sum(limite_superior))

  objetivo_total = function(lambda) {
    sum(pmin(limite_superior, pmax(limite_inferior, lambda * N2)))
  }

  if (n2_total == sum(limite_inferior)) return(limite_inferior)
  if (n2_total == sum(limite_superior)) return(limite_superior)

  lambda = uniroot(
    function(x) objetivo_total(x) - n2_total,
    interval = c(0, max(limite_superior / N2) + 1)
  )$root

  objetivo = pmin(limite_superior, pmax(limite_inferior, lambda * N2))
  m = as.integer(floor(objetivo))
  faltan = n2_total - sum(m)

  if (faltan > 0) {
    candidatos = which(m < limite_superior)
    orden = candidatos[order(objetivo[candidatos] - m[candidatos], decreasing = TRUE)]
    m[head(orden, faltan)] = m[head(orden, faltan)] + 1L
  }
  m
}

# Organización de las variables por TRI
# Valores de las variables de interés organizados por TRI
Y_tri = lapply(indices_tri, function(pos) Y[pos, , drop = FALSE])

# Función de simulación πPT + MAS
# Simular varias replicaciones sobre un mismo proceso de R
simular_bloque = function(replicaciones) {
  J = length(n2_candidatos)
  suma = matrix(0, nrow = n_variables, ncol = J)
  suma_cuadrados = matrix(0, nrow = n_variables, ncol = J)
  suma_tamanos = numeric(J)

  for (b in replicaciones) {
    # Primera etapa: selección πPT de TRI mediante Brewer
    certeza = which(pik >= 1 - 1e-10)
    restantes = which(pik < 1 - 1e-10)
    seleccionados = certeza
    n_restantes = as.integer(round(sum(pik[restantes])))

    if (n_restantes == length(restantes)) {
      seleccionados = c(seleccionados, restantes)
    } else if (n_restantes == 1L) {
      seleccionados = c(seleccionados, restantes[
        sample.int(length(restantes), 1L, prob = pik[restantes])
      ])
    } else if (n_restantes > 1L) {
      seleccionados = c(seleccionados,
                        restantes[UPbrewer(pik[restantes]) == 1])
    }

    N2 = M[seleccionados]
    # Asignar cada tamaño objetivo y registrar el tamaño efectivo
    asignaciones = lapply(n2_candidatos, function(n2) {
      asignar_segunda_etapa(N2, n2, minimo, maximo)
    })
    suma_tamanos = suma_tamanos + vapply(asignaciones, function(a) as.numeric(sum(a)), numeric(1))
    tamano_maximo = Reduce(pmax, asignaciones)
    estimaciones = matrix(0, nrow = n_variables, ncol = J)

    # Segunda etapa: una permutación MAS por TRI permite evaluar varios n2
    for (h in seq_along(seleccionados)) {
      k = seleccionados[h]
      posiciones = sample.int(M[k], tamano_maximo[h], replace = FALSE)
      acumulados = apply(Y_tri[[k]][posiciones, , drop = FALSE], 2, cumsum)
      dim(acumulados) = c(tamano_maximo[h], n_variables)

      for (a in seq_len(J)) {
        m_i = asignaciones[[a]][h]
        estimaciones[, a] = estimaciones[, a] +
          acumulados[m_i, ] * M[k] / (pik[k] * m_i * N)
      }
    }

    suma = suma + estimaciones
    suma_cuadrados = suma_cuadrados + estimaciones^2
  }

  list(suma = suma, suma_cuadrados = suma_cuadrados,
       suma_tamanos = suma_tamanos)
}

# Ejecución paralela y evaluación de la precisión
# Núcleos disponibles y bloques de replicaciones
nucleos_disponibles = detectCores(logical = FALSE)
if (is.na(nucleos_disponibles)) nucleos_disponibles = 1L
n_nucleos = min(B, max(1L, nucleos_disponibles - 1L))

bloques = split(seq_len(B), rep(seq_len(n_nucleos), length.out = B))

# Configuración del procesamiento paralelo
if (n_nucleos > 1L) {
  cluster = makeCluster(n_nucleos)
  invisible(clusterEvalQ(cluster, library(sampling)))
  clusterExport(cluster, c(
    "M", "N", "Y_tri", "n_variables", "n2_candidatos",
    "minimo", "maximo", "asignar_segunda_etapa", "simular_bloque"
  ))
  clusterSetRNGStream(cluster, 456)
} else {
  set.seed(456)
}

# Se conservan los totales de viviendas de segunda etapa como candidatos
resumen = data.table()
evaluacion = data.table()

for (n1 in n1_candidatos) {
  pik = inclusionprobabilities(M, n1)

  if (n_nucleos > 1L) {
    clusterExport(cluster, "pik")
    resultados = parLapply(cluster, bloques, simular_bloque)
  } else {
    resultados = lapply(bloques, simular_bloque)
  }

  suma = Reduce(`+`, lapply(resultados, function(x) x$suma))
  suma_cuadrados = Reduce(`+`, lapply(resultados, function(x) x$suma_cuadrados))
  suma_tamanos = Reduce(`+`, lapply(resultados, function(x) x$suma_tamanos))

  for (j in seq_along(n2_candidatos)) {
    n2 = n2_candidatos[j]
    viviendas_promedio = suma_tamanos[j] / B

    EE = sqrt(pmax((suma_cuadrados[, j] - suma[, j]^2 / B) / (B - 1L), 0))
    margen = z * EE

    # CV ordinario para la media y la proporción del 50 %
    CV = EE[1:2] / abs(valores[1:2])

    # CV logarítmico para la proporción del 1 %
    p = min(valores[3], 1 - valores[3])
    CV_log = EE[3] / (p * (-log(p)))

    CV_evaluado = c(CV, CV_log)
    tipo_cv = c("CV", "CV", "CV logarítmico")
    limites_cv = c(unname(c0), c_log)
    cumple_variable = margen <= d & CV_evaluado <= limites_cv
    cumple = all(cumple_variable)

    evaluacion = rbind(evaluacion, data.table(
      n1 = n1,
      n2_objetivo = n2,
      Viviendas_promedio = viviendas_promedio,
      Variable = colnames(Y),
      Valor = as.numeric(valores),
      EE = EE,
      Margen = margen,
      Medida = tipo_cv,
      CV_evaluado = CV_evaluado,
      Cumple_variable = cumple_variable
    ))

    resumen = rbind(resumen, data.table(
      n1 = n1,
      n2_objetivo = n2,
      Viviendas_promedio = viviendas_promedio,
      Cumple = cumple
    ))
  }
}

# Cierre de los procesos paralelos
if (n_nucleos > 1L) stopCluster(cluster)

setorder(resumen, Viviendas_promedio, n1)
setorder(evaluacion, Viviendas_promedio, n1, Variable)

# Resultados de precisión por variable
evaluacion[, .(n1, n2_objetivo, Viviendas_promedio, Variable, Valor, EE, Margen,
               Medida, CV_evaluado, Cumple_variable)] |>
  gt() |>
  cols_label(
    n1 = "TRI",
    n2_objetivo = "Viviendas objetivo",
    Viviendas_promedio = "Viviendas promedio",
    Variable = "Variable",
    Valor = "Valor poblacional",
    EE = "EE",
    Margen = "Margen",
    Medida = "Medida",
    CV_evaluado = "CV evaluado",
    Cumple_variable = "Cumple"
  ) |>
  fmt_number(columns = c(Valor, EE, Margen, CV_evaluado), decimals = 4) |>
  fmt_number(columns = Viviendas_promedio, decimals = 1) |>
  tab_options(table.width = pct(100))

# Resumen de combinaciones candidatas
resumen[, .(n1, n2_objetivo, Viviendas_promedio, Cumple)] |>
  gt() |>
  cols_label(
    n1 = "TRI seleccionados",
    n2_objetivo = "Viviendas objetivo",
    Viviendas_promedio = "Viviendas seleccionadas (promedio)",
    Cumple = "Cumple precisión"
  ) |>
  fmt_number(columns = Viviendas_promedio, decimals = 1) |>
  tab_options(table.width = pct(100))

# Determinación del tamaño de muestra
# Tamaños de ambas etapas que satisfacen la precisión requerida
n_final = resumen[Cumple == TRUE][order(Viviendas_promedio, n1)]

if (nrow(n_final) > 0L) {
  n_final = n_final[1, .(n1, n2_objetivo, Viviendas_promedio)]
} else {
  n_final = data.table(Resultado = "Ninguna combinación candidata cumple los requisitos de precisión")
}

n_final

