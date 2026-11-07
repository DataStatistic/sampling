
# Ejemplo muestreo en dos etapas

# Carga de paquetes, lectura de datos, 
# construcción de la población y marco de TRI

library(data.table)
library(sampling)
library(survey)
library(gt)

poblacion = fread(
  "datos/rmw_vivienda.csv",
  colClasses = list(character = c("COD_CUESTIONARIO", "COD_TRI"))
)

poblacion = poblacion[U_CLASE == 3 & !is.na(COD_TRI) & COD_TRI != ""]

poblacion[, `:=`(
  acueducto = as.integer(VIV_SERV_ACUE == 1),
  energia = as.integer(VIV_SERV_ENER == 1),
  internet = as.integer(VIV_SERV_INT == 1)
)]

# Marco de primera etapa: un registro por TRI
marco_tri = poblacion[, .(
  M_i = .N,
  municipio = first(U_MPIO),
  nombre_tri = first(NOM_TRI)
), by = COD_TRI]

K = nrow(marco_tri)
N = nrow(poblacion)
c(TRI = K, Viviendas = N)

# Visualización del marco de primera etapa
gt(head(marco_tri[order(COD_TRI)], 10)) |>
  cols_label(COD_TRI = "Código TRI", M_i = "Viviendas",
             municipio = "Municipio", nombre_tri = "Nombre del TRI") |>
  tab_options(table.width = pct(95), table.font.size = px(12))

# Proporciones poblacionales de referencia
parametros = data.table(
  Variable = c("Acueducto", "Energía eléctrica", "Internet"),
  Proporcion = c(mean(poblacion$acueducto),
                 mean(poblacion$energia),
                 mean(poblacion$internet))
)
parametros

# Tamaños muestrales y función de asignación proporcional restringida
n = 20L
n2_total = 1200L
minimo = 10L
maximo = 150L

asignar_segunda_etapa = function(N2, n2_total, minimo = 10L, maximo = 150L) {
  N2 = as.integer(N2)
  limite_inferior = pmin(minimo, N2)
  limite_superior = pmin(maximo, N2)

  if (n2_total < sum(limite_inferior) || n2_total > sum(limite_superior)) {
    stop("El tamaño de segunda etapa es incompatible con los límites establecidos.")
  }

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

# MAS–MAS: selección de TRI y asignación de viviendas
set.seed(123)
tri_mas = copy(marco_tri[COD_TRI %in%
  sample(marco_tri$COD_TRI, size = n, replace = FALSE)])

tri_mas[, m_i := asignar_segunda_etapa(M_i, n2_total, minimo, maximo)]
tri_mas[, `:=`(pi1 = n / K, pi2 = m_i / M_i)]

tri_mas[, .(COD_TRI, M_i, m_i, pi1, pi2)][order(COD_TRI)]

# MAS–MAS: marco de viviendas de un TRI seleccionado
tri_ejemplo = tri_mas$COD_TRI[1]
marco_segunda_ejemplo = poblacion[COD_TRI == tri_ejemplo,
  .(COD_TRI, COD_CUESTIONARIO)]

gt(head(marco_segunda_ejemplo, 10)) |>
  cols_label(COD_TRI = "Código TRI", COD_CUESTIONARIO = "Código vivienda") |>
  tab_options(table.width = pct(70), table.font.size = px(12))

# MAS–MAS: selección de viviendas y cálculo de probabilidades y pesos
muestra_mas = merge(poblacion, tri_mas[, .(COD_TRI, M_i, m_i, pi1, pi2)],
                    by = "COD_TRI", sort = FALSE)

muestra_mas = muestra_mas[, .SD[sample.int(.N, m_i[1], replace = FALSE)],
                          by = COD_TRI]

muestra_mas[, `:=`(
  K = K,
  pi = pi1 * pi2,
  peso = 1 / (pi1 * pi2)
)]

# MAS–MAS: visualización de las viviendas seleccionadas
gt(head(muestra_mas[order(COD_TRI, COD_CUESTIONARIO),
                    .(COD_TRI, COD_CUESTIONARIO, M_i, m_i,
                      acueducto, energia, internet, pi, peso)], 10)) |>
  fmt_number(columns = c(pi, peso), decimals = 4) |>
  tab_options(table.width = pct(100), table.font.size = px(11))

# MAS–MAS: especificación del diseño
diseno_mas = svydesign(
  ids = ~ COD_TRI + COD_CUESTIONARIO,
  fpc = ~ K + M_i,
  data = muestra_mas
)

diseno_mas

# MAS–MAS: estimación de proporciones
estimacion_mas = svytotal(
  ~ acueducto + energia + internet,
  design = diseno_mas
)

coef(estimacion_mas) / N

# MAS–MAS: varianzas, errores estándar, CV e intervalos de confianza
varianzas_mas = diag(vcov(estimacion_mas)) / N^2
errores_mas = SE(estimacion_mas) / N
cv_mas = 100 * errores_mas / (as.numeric(coef(estimacion_mas)) / N)
critico_mas = qt(0.975, df = degf(diseno_mas))
ic_mas = cbind(
  Inferior = as.numeric(coef(estimacion_mas)) / N - critico_mas * errores_mas,
  Superior = as.numeric(coef(estimacion_mas)) / N + critico_mas * errores_mas
)

varianzas_mas
errores_mas
cv_mas
ic_mas

# Función de resultados: CV convencional y CV logarítmico
# Función de organización de resultados de los dos diseños
resumir_estimaciones = function(objeto, diseno, N) {
  proporcion = as.numeric(coef(objeto)) / N
  varianza = as.numeric(diag(vcov(objeto))) / N^2
  error = as.numeric(SE(objeto)) / N
  critico = qt(0.975, df = degf(diseno))
  cv = rep(NA_real_, length(proporcion))
  cv_log = rep(NA_real_, length(proporcion))
  validas_cv = is.finite(proporcion) & proporcion > 0
  cv[validas_cv] = 100 * error[validas_cv] / proporcion[validas_cv]
  validas_log = is.finite(proporcion) & proporcion > 0 & proporcion < 1
  p_menor = pmin(proporcion[validas_log], 1 - proporcion[validas_log])
  cv_log[validas_log] = 100 * error[validas_log] / (-p_menor * log(p_menor))

  data.table(
    Variable = c("Acueducto", "Energía eléctrica", "Internet"),
    Estimacion = proporcion,
    Varianza = varianza,
    EE = error,
    IC_inf = proporcion - critico * error,
    IC_sup = proporcion + critico * error,
    CV = cv,
    CVLog = cv_log
  )
}
resultados_mas = resumir_estimaciones(estimacion_mas, diseno_mas, N)

resultados_mas[, .(
  Variable, Estimacion = round(Estimacion, 4),
  Varianza = signif(Varianza, 4), EE = round(EE, 5),
  IC_inf = round(IC_inf, 4), IC_sup = round(IC_sup, 4),
  CV = round(CV, 2), CVLog = round(CVLog, 2)
)]

# πPT–MAS: probabilidades, selección de TRI y asignación de viviendas
pik = inclusionprobabilities(marco_tri$M_i, n)

set.seed(789)
seleccion_pipt = integer(K)
certeza = which(pik >= 1 - 1e-10)
restantes = which(pik < 1 - 1e-10)
seleccion_pipt[certeza] = 1L
seleccion_pipt[restantes] = UPbrewer(pik[restantes])

tri_pipt = copy(marco_tri[seleccion_pipt == 1])
tri_pipt[, pi1 := pik[seleccion_pipt == 1]]
tri_pipt[, m_i := asignar_segunda_etapa(M_i, n2_total, minimo, maximo)]
tri_pipt[, pi2 := m_i / M_i]

tri_pipt[, .(COD_TRI, M_i, m_i, pi1, pi2)][order(COD_TRI)]

# πPT–MAS: selección de viviendas y cálculo de probabilidades y pesos
muestra_pipt = merge(poblacion,
                      tri_pipt[, .(COD_TRI, M_i, m_i, pi1, pi2)],
                      by = "COD_TRI", sort = FALSE)

muestra_pipt = muestra_pipt[, .SD[sample.int(.N, m_i[1], replace = FALSE)],
                            by = COD_TRI]

muestra_pipt[, `:=`(
  pi = pi1 * pi2,
  peso = 1 / (pi1 * pi2)
)]

# πPT–MAS: visualización de las viviendas seleccionadas
gt(head(muestra_pipt[order(COD_TRI, COD_CUESTIONARIO),
                     .(COD_TRI, COD_CUESTIONARIO, M_i, m_i,
                       acueducto, energia, internet, pi, peso)], 10)) |>
  fmt_number(columns = c(pi, peso), decimals = 4) |>
  tab_options(table.width = pct(100), table.font.size = px(11))

# πPT–MAS: especificación del diseño
diseno_pipt = svydesign(
  ids = ~ COD_TRI + COD_CUESTIONARIO,
  fpc = ~ pi1 + pi2,
  pps = "brewer",
  variance = "HT",
  data = muestra_pipt
)

diseno_pipt

# πPT–MAS: estimación de proporciones
estimacion_pipt = svytotal(
  ~ acueducto + energia + internet,
  design = diseno_pipt
)

coef(estimacion_pipt) / N

# πPT–MAS: varianzas, errores estándar, CV e intervalos de confianza
varianzas_pipt = diag(vcov(estimacion_pipt)) / N^2
errores_pipt = SE(estimacion_pipt) / N
cv_pipt = 100 * errores_pipt / (as.numeric(coef(estimacion_pipt)) / N)
critico_pipt = qt(0.975, df = degf(diseno_pipt))
ic_pipt = cbind(
  Inferior = as.numeric(coef(estimacion_pipt)) / N - critico_pipt * errores_pipt,
  Superior = as.numeric(coef(estimacion_pipt)) / N + critico_pipt * errores_pipt
)

varianzas_pipt
errores_pipt
cv_pipt
ic_pipt

# πPT–MAS: resultados resumidos
resultados_pipt = resumir_estimaciones(estimacion_pipt, diseno_pipt, N)

resultados_pipt[, .(
  Variable, Estimacion = round(Estimacion, 4),
  Varianza = signif(Varianza, 4), EE = round(EE, 5),
  IC_inf = round(IC_inf, 4), IC_sup = round(IC_sup, 4),
  CV = round(CV, 2), CVLog = round(CVLog, 2)
)]

# Comparación de los dos diseños
comparacion = rbindlist(list(
  copy(resultados_mas)[, Diseno := "MAS–MAS"],
  copy(resultados_pipt)[, Diseno := "πPT–MAS"]
), use.names = TRUE)

setcolorder(comparacion,
            c("Diseno", "Variable", "Estimacion", "Varianza",
              "EE", "IC_inf", "IC_sup", "CV", "CVLog"))

comparacion[, .(
  Diseno, Variable, Estimacion = round(Estimacion, 4),
  Varianza = signif(Varianza, 4), EE = round(EE, 5),
  IC_inf = round(IC_inf, 4), IC_sup = round(IC_sup, 4),
  CV = round(CV, 2), CVLog = round(CVLog, 2)
)]

# Comparación con los parámetros poblacionales
referencia = parametros[, .(Variable, Referencia = Proporcion)]

comparacion = merge(comparacion, referencia,
                    by = "Variable", sort = FALSE)
comparacion[, Error := Estimacion - Referencia]

comparacion[, .(Diseno, Variable,
                Referencia = round(Referencia, 4),
                Estimacion = round(Estimacion, 4),
                Error = round(Error, 4),
                CV = round(CV, 3),
                CVLog = round(CVLog, 3))]

