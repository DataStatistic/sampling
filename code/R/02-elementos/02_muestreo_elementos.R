
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


