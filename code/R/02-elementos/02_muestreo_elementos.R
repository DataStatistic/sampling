
# Ejemplo: 
# Diseños de muestreo de elementos Muestreo aleatorio 
# simple sin reemplazo, muestreo sistemático y muestreo
# con probabilidades proporcionales a una medida de tamaño

# Población: 
# Registros de 17.872 establecimientos educativos
# distribuidos en 33 departamentos y 1.118 municipios

# Variables:
# - id_unidad: identificador único de la unidad de muestreo
# - codigo_dane: código administrativo del establecimiento
# - nombre_establecimiento: nombre del establecimiento educativo
# - cod_dane_departamento: código del departamento
# - departamento`: nombre del departamento
# - cod_dane_municipio`: código del municipio
# - municipio`: nombre del municipio
# - sector`: clasificación del establecimiento como oficial o no oficial
# - calendario`: calendario académico
# - cantidad_sedes`: número de sedes del establecimiento
# - total_matricula: número de estudiantes matriculados

# Librerias necesarias
library(data.table)
library(sampling)
library(survey)

# Población y valores de referencia
poblacion = readRDS("datos/poblacion_educacion_2025.rds")
N = nrow(poblacion) ; N # tamaño de la población
View(poblacion) # población

# Parámetros
parametros = poblacion[, .(
  total_matricula = sum(total_matricula),
  promedio_matricula = mean(total_matricula)
)]

parametros

# Muestreo aleatorio simple sin reemplazo

# Definir los tamaños de población y muestra
N = nrow(poblacion)
n = 400 # se debe buscar

# Seleccionar una muestra aleatoria simple sin reemplazo
set.seed(123)
seleccion = srswor(n, N)
muestra_mas = copy(poblacion[seleccion == 1])

# Incorporar las probabilidades de inclusión y factores de expansión
muestra_mas[, pi := n / N]
muestra_mas[, peso := 1 / pi]

# Muestra
muestra_mas[, .(id_unidad, pi, peso, total_matricula)]

# Especificación del diseño muestral
# Incorporar el tamaño de la población finita
muestra_mas[, N := N]

# Especificar el diseño MAS sin reemplazo
diseno_mas = svydesign(
  ids = ~ 1,
  weights = ~ peso,
  fpc = ~ N,
  data = muestra_mas
)

diseno_mas

# Estimación de los parámetros

# Estimar el total poblacional
total_mas = svytotal(
  ~total_matricula,
  design = diseno_mas
)

# Estimar el promedio poblacional
promedio_mas = svymean(
  ~total_matricula,
  design = diseno_mas
)

# Estimaciones
coef(total_mas)
coef(promedio_mas)

# Estimación de la varianza

# Estimar las varianzas
varianza_total_mas = vcov(total_mas)[1, 1]
varianza_promedio_mas = vcov(promedio_mas)[1, 1]

# Obtener los errores estándar
ee_total_mas = SE(total_mas)
ee_promedio_mas = SE(promedio_mas)

# Varianzas
varianza_total_mas
varianza_promedio_mas

# Errores estándar
ee_total_mas
ee_promedio_mas

# Intervalos de confianza

# Intervalo de confianza para el total poblacional
ic_total_mas = confint(total_mas, level = 0.95, df = degf(diseno_mas))

# Intervalo de confianza para el promedio poblacional
ic_promedio_mas = confint(promedio_mas, level = 0.95, df = degf(diseno_mas))

ic_total_mas
ic_promedio_mas

# Coeficientes de variación
cv_total_mas = 100 * cv(total_mas)
cv_promedio_mas = 100 * cv(promedio_mas)

cv_total_mas
cv_promedio_mas

# Resumen de resultados
resultados_mas = data.table(
  Parametro = c("Total", "Promedio"),
  Estimacion = c(coef(total_mas), coef(promedio_mas)),
  Varianza = c(varianza_total_mas, varianza_promedio_mas),
  EE = c(ee_total_mas, ee_promedio_mas),
  IC_inf = c(ic_total_mas[1, 1], ic_promedio_mas[1, 1]),
  IC_sup = c(ic_total_mas[1, 2], ic_promedio_mas[1, 2]),
  CV = c(cv_total_mas, cv_promedio_mas)
)

resultados_mas

# Muestreo Sistemático sin reemplazo

# Definir los tamaños de población y muestra
N = nrow(poblacion)
n = 400

# Seleccionar la muestra sistemática
set.seed(456)
seleccion_sis = UPsystematic(rep(n / N, N))
muestra_sis = copy(poblacion[seleccion_sis == 1])

# Conservar el orden del marco muestral
setorder(muestra_sis, id_unidad)

# Incorporar probabilidades de inclusión y factores de expansión
muestra_sis[, pi := n / N]
muestra_sis[, peso := 1 / pi]

# Muestra
muestra_sis[, .(id_unidad, pi, peso, total_matricula)]

# Estimar el total y el promedio poblacionales
total_sis = muestra_sis[
  , sum(peso * total_matricula)
]

promedio_sis = total_sis / N

total_sis

promedio_sis

library(svrep)

# Construir la matriz del estimador de varianza SD1
matriz_sd1 = make_quad_form_matrix(
  variance_estimator = "SD1",
  cluster_ids = data.frame(id = seq_len(n)),
  sort_order = seq_len(n),
  strata_pop_sizes = data.frame(N = rep(N, n))
)

# Valores de matrícula expandidos
y_expandida = muestra_sis$peso * muestra_sis$total_matricula

# Estimar la varianza del total mediante SD1
varianza_total_sis = as.numeric(
  crossprod(y_expandida, matriz_sd1 %*% y_expandida)
)

# Obtener la varianza del promedio
varianza_promedio_sis = varianza_total_sis / N^2

# Calcular los errores estándar
ee_total_sis = sqrt(varianza_total_sis)
ee_promedio_sis = sqrt(varianza_promedio_sis)

# Varianzas
varianza_total_sis

varianza_promedio_sis

# Errores estándar
ee_total_sis

ee_promedio_sis

# Valor crítico para un nivel de confianza del 95 %
z = qnorm(0.975)

# Intervalos de confianza para el total y el promedio
ic_total_sis = c(
  total_sis - z * ee_total_sis,
  total_sis + z * ee_total_sis
)

ic_promedio_sis = c(
  promedio_sis - z * ee_promedio_sis,
  promedio_sis + z * ee_promedio_sis
)

ic_total_sis

ic_promedio_sis

# Calcular los coeficientes de variación en porcentaje
cv_total_sis = 100 * ee_total_sis / abs(total_sis)
cv_promedio_sis = 100 * ee_promedio_sis / abs(promedio_sis)

cv_total_sis

cv_promedio_sis

resultados_sis = data.table(
  Parametro = c("Total", "Promedio"),
  Estimacion = c(total_sis, promedio_sis),
  Varianza = c(varianza_total_sis, varianza_promedio_sis),
  Error_estandar = c(ee_total_sis, ee_promedio_sis),
  IC_inferior = c(ic_total_sis[1], ic_promedio_sis[1]),
  IC_superior = c(ic_total_sis[2], ic_promedio_sis[2]),
  CV = c(cv_total_sis, cv_promedio_sis)
)

resultados_sis

# Muestreo piPT sin reemplazo

# Definir los tamaños de población y muestra
N = nrow(poblacion)
n = 400

# Calcular las probabilidades de inclusión
pik = inclusionprobabilities(
  poblacion$cantidad_sedes,
  n
)

# Seleccionar una muestra πPT mediante el método de Brewer
set.seed(789)
seleccion_pipt = UPbrewer(pik)
muestra_pipt = copy(poblacion[seleccion_pipt == 1])

# Incorporar las probabilidades y los factores de expansión
muestra_pipt[, pi := pik[seleccion_pipt == 1]]
muestra_pipt[, peso := 1 / pi]

# Muestra
muestra_pipt[ , .(id_unidad, cantidad_sedes, total_matricula, pi, peso)]

# Especificar el diseño πPT sin reemplazo
diseno_pipt = svydesign(
  ids = ~1,
  fpc = ~pi,
  pps = "brewer",
  variance = "HT",
  data = muestra_pipt
)

# Estimar el total poblacional
total_pipt = svytotal(
  ~total_matricula,
  design = diseno_pipt
)

# Obtener las estimaciones del total y el promedio
estimacion_total_pipt = as.numeric(coef(total_pipt))
estimacion_promedio_pipt = estimacion_total_pipt / N

estimacion_total_pipt

estimacion_promedio_pipt

# Obtener las varianzas estimadas
varianza_total_pipt = vcov(total_pipt)[1, 1]
varianza_promedio_pipt = varianza_total_pipt / N^2

# Obtener los errores estándar
ee_total_pipt = SE(total_pipt)
ee_promedio_pipt = ee_total_pipt / N

# Vianzas
varianza_total_pipt

varianza_promedio_pipt

# Errores estándar
ee_total_pipt

ee_promedio_pipt

# Intervalo de confianza del 95 % para el total
ic_total_pipt = confint(
  total_pipt,
  level = 0.95,
  df = degf(diseno_pipt)
)

# Intervalo de confianza del 95 % para el promedio
ic_promedio_pipt = ic_total_pipt / N

ic_total_pipt

ic_promedio_pipt

# Calcular los coeficientes de variación en porcentaje
cv_total_pipt = 100 * cv(total_pipt)
cv_promedio_pipt = cv_total_pipt

cv_total_pipt

cv_promedio_pipt

resultados_pipt = data.table(
  Parametro = c("Total", "Promedio"),
  Estimacion = c(estimacion_total_pipt, estimacion_promedio_pipt),
  Varianza = c(varianza_total_pipt, varianza_promedio_pipt),
  Error_estandar = c(ee_total_pipt, ee_promedio_pipt),
  IC_inferior = c(ic_total_pipt[1, 1], ic_promedio_pipt[1, 1]),
  IC_superior = c(ic_total_pipt[1, 2], ic_promedio_pipt[1, 2]),
  CV = c(cv_total_pipt, cv_promedio_pipt)
)

resultados_pipt

# Muestreo Sistemático piPT

# Definir los tamaños de población y muestra
N = nrow(poblacion)
n = 400

# Calcular las probabilidades de inclusión proporcionales al tamaño
pik_pipt_sis = inclusionprobabilities(
  poblacion$cantidad_sedes,
  n
)

# Seleccionar una muestra sistemática πPT
set.seed(987)
seleccion_pipt_sis = UPsystematic(pik_pipt_sis)

muestra_pipt_sis = copy(poblacion[seleccion_pipt_sis == 1])

# Conservar el orden del marco muestral
setorder(muestra_pipt_sis, id_unidad)

# Incorporar probabilidades de inclusión y factores de expansión
muestra_pipt_sis[
  , pi := pik_pipt_sis[seleccion_pipt_sis == 1]
]
muestra_pipt_sis[, peso := 1 / pi]

# Muestra
muestra_pipt_sis[ , .(id_unidad, cantidad_sedes, pi, peso, total_matricula)]

# Estimar el total y el promedio poblacionales
total_pipt_sis = muestra_pipt_sis[
  , sum(peso * total_matricula)
]

promedio_pipt_sis = total_pipt_sis / N

total_pipt_sis

promedio_pipt_sis

# Construir la matriz del estimador de varianza SD1
matriz_sd1_pipt_sis = make_quad_form_matrix(
  variance_estimator = "SD1",
  cluster_ids = data.frame(id = seq_len(n)),
  sort_order = seq_len(n),
  strata_pop_sizes = data.frame(N = rep(N, n))
)

# Valores expandidos de matrícula
y_expandida_pipt_sis = muestra_pipt_sis$peso *
  muestra_pipt_sis$total_matricula

# Estimar la varianza del total mediante SD1
varianza_total_pipt_sis = as.numeric(
  crossprod(
    y_expandida_pipt_sis,
    matriz_sd1_pipt_sis %*% y_expandida_pipt_sis
  )
)

# Obtener la varianza del promedio
varianza_promedio_pipt_sis = varianza_total_pipt_sis / N^2

# Calcular los errores estándar
ee_total_pipt_sis = sqrt(varianza_total_pipt_sis)
ee_promedio_pipt_sis = sqrt(varianza_promedio_pipt_sis)

# Varianzas
varianza_total_pipt_sis

varianza_promedio_pipt_sis

# Errores estándar
ee_total_pipt_sis

ee_promedio_pipt_sis

# Valor crítico
z = qnorm(0.975)

# Intervalos de confianza para el total y el promedio
ic_total_pipt_sis = c(
  total_pipt_sis - z * ee_total_pipt_sis,
  total_pipt_sis + z * ee_total_pipt_sis
)

ic_promedio_pipt_sis = c(
  promedio_pipt_sis - z * ee_promedio_pipt_sis,
  promedio_pipt_sis + z * ee_promedio_pipt_sis
)

ic_total_pipt_sis

ic_promedio_pipt_sis

# Calcular los coeficientes de variación en porcentaje
cv_total_pipt_sis = 100 * ee_total_pipt_sis /
  abs(total_pipt_sis)

cv_promedio_pipt_sis = 100 * ee_promedio_pipt_sis /
  abs(promedio_pipt_sis)

cv_total_pipt_sis

cv_promedio_pipt_sis

resultados_pipt_sis = data.table(
  Parametro = c("Total", "Promedio"),
  Estimacion = c(total_pipt_sis, promedio_pipt_sis),
  Varianza = c(varianza_total_pipt_sis, varianza_promedio_pipt_sis),
  Error_estandar = c(ee_total_pipt_sis, ee_promedio_pipt_sis),
  IC_inferior = c(ic_total_pipt_sis[1], ic_promedio_pipt_sis[1]),
  IC_superior = c(ic_total_pipt_sis[2], ic_promedio_pipt_sis[2]),
  CV = c(cv_total_pipt_sis, cv_promedio_pipt_sis)
)

resultados_pipt_sis

comparacion = rbindlist(
  list(
    copy(resultados_mas)[
      , .(Diseno = "MAS", Parametro, Estimacion, Varianza,
          EE, IC_inf, IC_sup, CV)
    ],
    copy(resultados_sis)[
      , .(Diseno = "SIS", Parametro, Estimacion, Varianza,
          EE = Error_estandar,
          IC_inf = IC_inferior,
          IC_sup = IC_superior, CV)
    ],
    copy(resultados_pipt)[
      , .(Diseno = "πPT", Parametro, Estimacion, Varianza,
          EE = Error_estandar,
          IC_inf = IC_inferior,
          IC_sup = IC_superior, CV)
    ],
    copy(resultados_pipt_sis)[
      , .(Diseno = "πPT SIS", Parametro, Estimacion, Varianza,
          EE = Error_estandar,
          IC_inf = IC_inferior,
          IC_sup = IC_superior, CV)
    ]
  )
)

comparacion




