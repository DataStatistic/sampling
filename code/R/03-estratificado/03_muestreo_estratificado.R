
# Ejemplo Muestreo Estratificado

library(data.table)
library(sampling)
library(survey)
library(svrep)
library(gt)

# Marco muestral y valores de referencia

# Población
poblacion = readRDS("datos/poblacion_educacion_2025.rds")

# Variable de estratificación: "region" o "departamento"
variable_estratificacion = "region"

# Agrupación territorial de los departamentos en seis regiones
codigo_departamento = substr(sprintf("%02d", as.integer(poblacion$cod_dane_departamento)), 1, 2)
poblacion[, region := fcase(
  codigo_departamento %in% c("08", "13", "20", "23", "44", "47", "70"), "Caribe",
  codigo_departamento %in% c("15", "25", "50", "54", "68"), "Oriental",
  codigo_departamento %in% c("05", "17", "18", "41", "63", "66", "73"), "Central",
  codigo_departamento %in% c("19", "27", "52", "76"), "Pacífica",
  codigo_departamento %in% c("81", "85", "88", "91", "94", "95", "86", "97", "99"), "Orinoquía, A., I.",
  codigo_departamento == "11", "Bogotá D. C."
)]

poblacion[, estrato := if (variable_estratificacion == "region")
  region else as.character(cod_dane_departamento)]

# Conservar el orden original para el muestreo sistemático
poblacion[, orden_marco := .I]

# Tamaño poblacional y número de estratos
N = nrow(poblacion)
H = uniqueN(poblacion$estrato)

N

# Valores poblacionales para contrastar las estimaciones
parametros = poblacion[, .(
  total_matricula = sum(total_matricula),
  promedio_matricula = mean(total_matricula)
)]

parametros

# Tamaño de cada estrato y total de la medida auxiliar
estratos = poblacion[, .(
  N_h = .N,
  t_x_h = sum(cantidad_sedes)
), by = estrato]

setorder(estratos, estrato)
estratos

# Asignación del tamaño de muestra
# Asignación proporcional entera, con mínimo por estrato
n = 1200L
minimo = 3L

asignacion = copy(estratos)
asignacion[, objetivo := n * N_h / N]
asignacion[, n_h := pmin(N_h, minimo)]

# Distribuir las unidades restantes hacia los mayores déficits
while (sum(asignacion$n_h) < n) {
  disponibles = which(asignacion$n_h < asignacion$N_h)
  deficit = asignacion$objetivo[disponibles] - asignacion$n_h[disponibles]
  h = disponibles[which.max(deficit)]
  asignacion$n_h[h] = asignacion$n_h[h] + 1L
}

asignacion[, f_h := n_h / N_h]

# Comprobar el tamaño total y consultar la asignación
sum(asignacion$n_h)
asignacion[, .(estrato, N_h, n_h, f_h)]

# Muestreo aleatorio simple estratificado
# Seleccionar MAS de forma independiente dentro de cada estrato
set.seed(123)
muestras_mas = lapply(seq_len(nrow(asignacion)), function(h) {
  info = asignacion[h]
  marco_h = poblacion[estrato == info$estrato]
  indicador = srswor(info$n_h, info$N_h)
  muestra_h = copy(marco_h[indicador == 1])
  muestra_h[, `:=`(
    N_h = info$N_h,
    n_h = info$n_h,
    pi = info$n_h / info$N_h
  )]
  muestra_h[, peso := 1 / pi]
  muestra_h
})

muestra_mas = rbindlist(muestras_mas)

# Muestra y pesos
muestra_mas[, .(id_unidad, estrato, total_matricula, pi, peso)]

# Especificación del diseño muestral
# Diseño MAS estratificado sin reemplazo
diseno_mas = svydesign(
  ids = ~1,
  strata = ~estrato,
  weights = ~peso,
  fpc = ~N_h,
  data = muestra_mas
)

diseno_mas

# Estimación de los parámetros
# Estimar el total y el promedio nacionales
total_mas = svytotal(~total_matricula, design = diseno_mas)
promedio_mas = as.numeric(coef(total_mas)) / N
coef(total_mas)
promedio_mas

# Estimación de las varianzas y errores estándar
# Varianzas estimadas y errores estándar
varianza_total_mas = vcov(total_mas)[1, 1]
varianza_promedio_mas = varianza_total_mas / N^2

ee_total_mas = sqrt(varianza_total_mas)
ee_promedio_mas = sqrt(varianza_promedio_mas)

varianza_total_mas
varianza_promedio_mas

ee_total_mas
ee_promedio_mas

# Intervalos de confianza del 95 %
ic_total_mas = confint(total_mas, level = 0.95, df = degf(diseno_mas))
ic_promedio_mas = ic_total_mas / N

ic_total_mas
ic_promedio_mas

# Coeficientes de variación
cv_total_mas = 100 * cv(total_mas)
cv_promedio_mas = cv_total_mas

cv_total_mas
cv_promedio_mas

# Resultados por estrato

# Preparar las estimaciones y medidas de precisión por estrato
tabla_estratos = function(totales_h, varianzas_h, asignacion, critico) {
  N_h = asignacion$N_h
  ee_total = sqrt(varianzas_h)
  estimaciones = c(totales_h, totales_h / N_h)
  varianzas = c(varianzas_h, varianzas_h / N_h^2)
  errores = c(ee_total, ee_total / N_h)
  limites = rep(critico, 2)
  
  tabla = data.table(
    Estrato = rep(asignacion$estrato, 2),
    Parametro = rep(c("Total", "Promedio"), each = nrow(asignacion)),
    N_h = rep(N_h, 2),
    n_h = rep(asignacion$n_h, 2),
    Estimacion = round(estimaciones, 1),
    Varianza = round(varianzas, 1),
    EE = round(errores, 1),
    IC_inf = round(estimaciones - limites * errores, 1),
    IC_sup = round(estimaciones + limites * errores, 1),
    CV = round(100 * errores / abs(estimaciones), 3)
  )
  setorder(tabla, Estrato, Parametro)
  tabla
}

# Estimaciones y varianzas por estrato bajo MAS
estimaciones_h_mas = svyby(
  ~total_matricula, ~estrato, diseno_mas,
  svytotal, vartype = "se"
)
posiciones_mas = match(as.character(asignacion$estrato),
                       as.character(estimaciones_h_mas$estrato))

resultados_estratos_mas = tabla_estratos(
  as.numeric(coef(estimaciones_h_mas))[posiciones_mas],
  as.numeric(SE(estimaciones_h_mas))[posiciones_mas]^2,
  asignacion, qt(0.975, pmax(1, asignacion$n_h - 1))
)
resultados_estratos_mas

# Muestreo sistemático estratificado

# Seleccionar sistemáticamente dentro de cada estrato
set.seed(456)
muestras_sis = lapply(seq_len(nrow(asignacion)), function(h) {
  info = asignacion[h]
  marco_h = poblacion[estrato == info$estrato]
  setorder(marco_h, orden_marco)
  indicador = UPsystematic(rep(info$n_h / info$N_h, info$N_h))
  muestra_h = copy(marco_h[indicador == 1])
  muestra_h[, `:=`(
    N_h = info$N_h,
    n_h = info$n_h,
    pi = info$n_h / info$N_h
  )]
  muestra_h[, peso := 1 / pi]
  muestra_h
})

muestra_sis = rbindlist(muestras_sis)
setorder(muestra_sis, estrato, orden_marco)
muestra_sis[, .(id_unidad, estrato, total_matricula, pi, peso)]

# Estimación de los parámetros
# Total y promedio nacionales
total_sis = muestra_sis[, sum(peso * total_matricula)]
promedio_sis = total_sis / N

total_sis
promedio_sis

# Estimación de las varianzas y errores estándar
# Calcular SD1 por estrato y sumar sus contribuciones
varianzas_sis = muestra_sis[, {
  N_h_actual = N_h[1]
  n_h_actual = .N
  y_expandida = peso * total_matricula
  
  varianza_h = if (n_h_actual == N_h_actual) {
    0
  } else {
    matriz_sd1 = make_quad_form_matrix(
      variance_estimator = "SD1",
      cluster_ids = data.frame(id = seq_len(n_h_actual)),
      sort_order = seq_len(n_h_actual),
      strata_pop_sizes = data.frame(N = rep(N_h_actual, n_h_actual))
    )
    as.numeric(crossprod(y_expandida, matriz_sd1 %*% y_expandida))
  }
  .(varianza = varianza_h)
}, by = estrato]

varianza_total_sis = sum(varianzas_sis$varianza)
varianza_promedio_sis = varianza_total_sis / N^2

ee_total_sis = sqrt(varianza_total_sis)
ee_promedio_sis = sqrt(varianza_promedio_sis)

varianza_total_sis
varianza_promedio_sis

ee_total_sis
ee_promedio_sis

# Intervalos de confianza aproximados del 95 %
z = qnorm(0.975)
ic_total_sis = c(total_sis - z * ee_total_sis,
                 total_sis + z * ee_total_sis)
ic_promedio_sis = ic_total_sis / N
ic_total_sis
ic_promedio_sis

# Coeficientes de variación
cv_total_sis = 100 * ee_total_sis / abs(total_sis)
cv_promedio_sis = 100 * ee_promedio_sis / abs(promedio_sis)

cv_total_sis
cv_promedio_sis

# Resultados por estrato
# Estimaciones por estrato y varianzas SD1 ya calculadas
estimaciones_h_sis = muestra_sis[, .(total = sum(peso * total_matricula)),
                                 by = estrato]
posiciones_sis = match(asignacion$estrato,
                       estimaciones_h_sis$estrato)
posiciones_var_sis = match(asignacion$estrato,
                           varianzas_sis$estrato)

resultados_estratos_sis = tabla_estratos(
  estimaciones_h_sis$total[posiciones_sis],
  varianzas_sis$varianza[posiciones_var_sis],
  asignacion, z
)
resultados_estratos_sis

# Muestreo estratificado piPT

# Selección de la muestra
# Selección πPT de Brewer independiente en cada estrato
set.seed(789)
muestras_pipt = lapply(seq_len(nrow(asignacion)), function(h) {
  info = asignacion[h]
  marco_h = poblacion[estrato == info$estrato]
  pi_h = inclusionprobabilities(marco_h$cantidad_sedes, info$n_h)
  
  # Incorporar directamente los elementos de inclusión cierta
  certeza = pi_h >= 1 - 1e-10
  seleccion = certeza
  if (any(!certeza)) {
    seleccion[!certeza] = UPbrewer(pi_h[!certeza]) == 1
  }
  
  muestra_h = copy(marco_h[seleccion])
  muestra_h[, `:=`(
    N_h = info$N_h,
    n_h = info$n_h,
    pi = pi_h[seleccion]
  )]
  muestra_h[, peso := 1 / pi]
  muestra_h
})

muestra_pipt = rbindlist(muestras_pipt)
muestra_pipt[, .(id_unidad, estrato, cantidad_sedes,
                 total_matricula, pi, peso)]

# Especificación del diseño muestral
# Diseño estratificado πPT: aproximación de Brewer
diseno_pipt = svydesign(
  ids = ~1,
  strata = ~estrato,
  fpc = ~pi,
  pps = "brewer",
  variance = "HT",
  data = muestra_pipt
)

diseno_pipt

# Estimación de los parámetros
# Estimar el total con survey y obtener el promedio
total_pipt = svytotal(~total_matricula, design = diseno_pipt)
estimacion_total_pipt = as.numeric(coef(total_pipt))
estimacion_promedio_pipt = estimacion_total_pipt / N

estimacion_total_pipt
estimacion_promedio_pipt

# Estimación de las varianzas y errores estándar
# Varianzas estimadas y errores estándar
varianza_total_pipt = vcov(total_pipt)[1, 1]
varianza_promedio_pipt = varianza_total_pipt / N^2

ee_total_pipt = sqrt(varianza_total_pipt)
ee_promedio_pipt = sqrt(varianza_promedio_pipt)

varianza_total_pipt
varianza_promedio_pipt

ee_total_pipt
ee_promedio_pipt

# Intervalos de confianza del 95 % usando los grados de libertad del diseño
ic_total_pipt = confint(total_pipt, level = 0.95, df = degf(diseno_pipt))
ic_promedio_pipt = ic_total_pipt / N

ic_total_pipt
ic_promedio_pipt

# Coeficientes de variación

# Coeficientes de variación en porcentaje
cv_total_pipt = 100 * ee_total_pipt / abs(estimacion_total_pipt)
cv_promedio_pipt = 100 * ee_promedio_pipt / abs(estimacion_promedio_pipt)

cv_total_pipt
cv_promedio_pipt

# Resultados por estrato
# Estimaciones y varianzas por estrato con Brewer
estimaciones_h_pipt = svyby(
  ~total_matricula, ~estrato, diseno_pipt,
  svytotal, vartype = "se"
)
posiciones_pipt = match(as.character(asignacion$estrato),
                        as.character(estimaciones_h_pipt$estrato))

resultados_estratos_pipt = tabla_estratos(
  as.numeric(coef(estimaciones_h_pipt))[posiciones_pipt],
  as.numeric(SE(estimaciones_h_pipt))[posiciones_pipt]^2,
  asignacion, qt(0.975, pmax(1, asignacion$n_h - 1))
)
resultados_estratos_pipt

# Muestreo sistemático estratificado piPT

# Selección sistemática πPT independiente por estrato
set.seed(987)
muestras_pipt_sis = lapply(seq_len(nrow(asignacion)), function(h) {
  info = asignacion[h]
  marco_h = poblacion[estrato == info$estrato]
  setorder(marco_h, orden_marco)
  pi_h = inclusionprobabilities(marco_h$cantidad_sedes, info$n_h)
  indicador = UPsystematic(pi_h)
  muestra_h = copy(marco_h[indicador == 1])
  muestra_h[, `:=`(
    N_h = info$N_h,
    n_h = info$n_h,
    pi = pi_h[indicador == 1]
  )]
  muestra_h[, peso := 1 / pi]
  muestra_h
})

muestra_pipt_sis = rbindlist(muestras_pipt_sis)
setorder(muestra_pipt_sis, estrato, orden_marco)
muestra_pipt_sis[, .(id_unidad, estrato, cantidad_sedes,
                     total_matricula, pi, peso)]

# Estimación de los parámetros
# Total y promedio nacionales
total_pipt_sis = muestra_pipt_sis[, sum(peso * total_matricula)]
promedio_pipt_sis = total_pipt_sis / N

total_pipt_sis
promedio_pipt_sis

# Estimación de las varianzas y errores estándar

# SD1 por estrato con los valores de matrícula expandidos
varianzas_pipt_sis = muestra_pipt_sis[, {
  N_h_actual = N_h[1]
  n_h_actual = .N
  y_expandida = peso * total_matricula
  
  varianza_h = if (n_h_actual == N_h_actual) {
    0
  } else {
    matriz_sd1 = make_quad_form_matrix(
      variance_estimator = "SD1",
      cluster_ids = data.frame(id = seq_len(n_h_actual)),
      sort_order = seq_len(n_h_actual),
      strata_pop_sizes = data.frame(N = rep(N_h_actual, n_h_actual))
    )
    as.numeric(crossprod(y_expandida, matriz_sd1 %*% y_expandida))
  }
  .(varianza = varianza_h)
}, by = estrato]

varianza_total_pipt_sis = sum(varianzas_pipt_sis$varianza)
varianza_promedio_pipt_sis = varianza_total_pipt_sis / N^2

ee_total_pipt_sis = sqrt(varianza_total_pipt_sis)
ee_promedio_pipt_sis = sqrt(varianza_promedio_pipt_sis)

varianza_total_pipt_sis
varianza_promedio_pipt_sis

ee_total_pipt_sis
ee_promedio_pipt_sis

# Intervalos de confianza
ic_total_pipt_sis = c(total_pipt_sis - z * ee_total_pipt_sis,
                      total_pipt_sis + z * ee_total_pipt_sis)
ic_promedio_pipt_sis = ic_total_pipt_sis / N

ic_total_pipt_sis
ic_promedio_pipt_sis

# Coeficientes de variación en porcentaje
cv_total_pipt_sis = 100 * ee_total_pipt_sis / abs(total_pipt_sis)
cv_promedio_pipt_sis = 100 * ee_promedio_pipt_sis / abs(promedio_pipt_sis)

cv_total_pipt_sis
cv_promedio_pipt_sis

# Resultados por estrato
# Estimaciones por estrato y varianzas SD1 ya calculadas
estimaciones_h_pipt_sis = muestra_pipt_sis[, .(total = sum(peso * total_matricula)),
                                           by = estrato]
posiciones_pipt_sis = match(asignacion$estrato,
                            estimaciones_h_pipt_sis$estrato)
posiciones_var_pipt_sis = match(asignacion$estrato,
                                varianzas_pipt_sis$estrato)

resultados_estratos_pipt_sis = tabla_estratos(
  estimaciones_h_pipt_sis$total[posiciones_pipt_sis],
  varianzas_pipt_sis$varianza[posiciones_var_pipt_sis],
  asignacion, z
)
resultados_estratos_pipt_sis

# Comparación de los diseños muestrales
# Comparación de las estimaciones

comparacion = rbindlist(list(
  copy(resultados_estratos_mas)[, Diseno := "ESTMAS"],
  copy(resultados_estratos_sis)[, Diseno := "ESTSIS"],
  copy(resultados_estratos_pipt)[, Diseno := "ESTπPT Brewer"],
  copy(resultados_estratos_pipt_sis)[, Diseno := "ESTπPT SIS"]
), use.names = TRUE)

setcolorder(comparacion, c("Diseno", "Parametro", "Estimacion",
                           "Varianza", "EE", "IC_inf", "IC_sup", "CV"))

comparacion






