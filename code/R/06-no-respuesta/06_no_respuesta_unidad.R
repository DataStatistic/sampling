# Ejemplo No respuesta de unidad: ajustes por clases y propensión

# Construcción de la población

library(data.table)
library(survey)
library(svrep)
library(gt)

options(width = 150, scipen = 6)

poblacion = fread(
  "datos/rmw_vivienda.csv",
  colClasses = list(character = c("COD_CUESTIONARIO", "COD_TRI"))
)

poblacion = poblacion[
  U_CLASE == 3 & !is.na(COD_TRI) & COD_TRI != ""
]

servicios = c("VIV_SERV_ACUE", "VIV_SERV_ENER", "VIV_SERV_INT")
for (v in servicios) set(poblacion, j = v, value = as.integer(poblacion[[v]]))

poblacion = poblacion[
  VIV_SERV_ACUE %in% 1:2 & VIV_SERV_ENER %in% 1:2 &
    VIV_SERV_INT %in% 1:2 & !is.na(U_MPIO)
]

poblacion[, `:=`(
  acueducto = as.integer(VIV_SERV_ACUE == 1L),
  energia = as.integer(VIV_SERV_ENER == 1L),
  internet = as.integer(VIV_SERV_INT == 1L),
  municipio = factor(U_MPIO)
)]

marco_tri = poblacion[, .(M_i = .N), by = COD_TRI]
K = nrow(marco_tri)
N = nrow(poblacion)

c(TRI = K, Viviendas = N)

# Valor poblacional de referencia

P = mean(poblacion$internet)
data.table(Parametro = "Proporción con internet", Referencia = P)

# Diseño muestral

asignar_segunda_etapa = function(N2, n2_total, minimo = 10L, maximo = 150L) {
  N2 = as.integer(N2)
  a = pmin(minimo, N2)
  b = pmin(maximo, N2)
  if (n2_total < sum(a) || n2_total > sum(b)) {
    stop("El tamaño de segunda etapa no satisface los límites.")
  }
  if (n2_total == sum(a)) return(a)
  if (n2_total == sum(b)) return(b)
  objetivo_total = function(lambda) sum(pmin(b, pmax(a, lambda * N2)))
  lambda = uniroot(
    function(x) objetivo_total(x) - n2_total,
    interval = c(0, max(b / N2) + 1)
  )$root
  objetivo = pmin(b, pmax(a, lambda * N2))
  m = as.integer(floor(objetivo))
  faltan = n2_total - sum(m)
  if (faltan > 0) {
    candidatos = which(m < b)
    orden = candidatos[order(objetivo[candidatos] - m[candidatos],
                              decreasing = TRUE)]
    m[head(orden, faltan)] = m[head(orden, faltan)] + 1L
  }
  m
}

set.seed(123)
n = 20L
n2_total = 1200L
tri_muestra = copy(marco_tri[COD_TRI %in%
  sample(marco_tri$COD_TRI, n, replace = FALSE)])
tri_muestra[, m_i := asignar_segunda_etapa(M_i, n2_total)]

muestra = merge(
  poblacion,
  tri_muestra,
  by = "COD_TRI",
  sort = FALSE
)
muestra = muestra[, .SD[sample.int(.N, m_i[1])], by = COD_TRI]
muestra[, `:=`(
  K = K,
  pi1 = n / K,
  pi2 = m_i / M_i,
  peso = (K / n) * (M_i / m_i)
)]
setorder(muestra, COD_TRI, COD_CUESTIONARIO)

stopifnot(nrow(muestra) == n2_total, uniqueN(muestra$COD_TRI) == n)

gt(head(muestra[, .(COD_TRI, COD_CUESTIONARIO, M_i, m_i,
                    internet, peso)], 10)) |>
  fmt_number(columns = peso, decimals = 3) |>
  tab_options(table.width = pct(95), table.font.size = px(12))

# Especificación del diseño

diseno = svydesign(
  ids = ~ COD_TRI + COD_CUESTIONARIO,
  fpc = ~ K + M_i,
  data = muestra
)

set.seed(2026)
diseno_rep = as_bootstrap_design(
  diseno,
  type = "Rao-Wu-Yue-Beaumont",
  replicates = 500,
  mse = TRUE
)

# Estimación con respuesta completa

estimacion_completa = svytotal(~internet, diseno)
P_completa = as.numeric(coef(estimacion_completa)) / N
V_completa = as.numeric(vcov(estimacion_completa)) / N^2
c(Estimacion = P_completa, Varianza = V_completa)

# Generación de no respuesta de unidad

set.seed(404)

muestra[, municipio_id := as.integer(factor(municipio))]
muestra[, log_M := log(M_i)]
muestra[, log_M_c := log_M - mean(log_M)]

muestra[, rho_generadora := plogis(
  qlogis(0.80) - 0.28 * (municipio_id - mean(municipio_id)) +
    0.27 * log_M_c
)]
muestra[, respuesta := rbinom(.N, 1L, rho_generadora)]

# La variable de interés queda oculta para las viviendas no respondientes.
muestra[, internet_obs := fifelse(respuesta == 1L, internet, NA_integer_)]

c(
  Muestra = nrow(muestra),
  Respondientes = sum(muestra$respuesta),
  No_respondientes = sum(muestra$respuesta == 0L),
  Tasa_respuesta = mean(muestra$respuesta)
)

# Generación de no respuesta de unidad

gt(head(muestra[, .(COD_TRI, COD_CUESTIONARIO, municipio,
                    M_i, peso, respuesta, internet_obs)], 10)) |>
  fmt_number(columns = peso, decimals = 3) |>
  tab_options(table.width = pct(100), table.font.size = px(12))

# Estimación sin ajuste por no respuesta

total_sin_ajuste = sum(muestra$peso[muestra$respuesta == 1] *
                         muestra$internet_obs[muestra$respuesta == 1])
P_sin_ajuste = total_sin_ajuste / N
P_sin_ajuste

# Ajuste mediante clases de respuesta

limite = median(marco_tri$M_i)
muestra[, clase_tamano := factor(
  fifelse(M_i <= limite, "TRI pequeño", "TRI grande")
)]
muestra[, clase := interaction(municipio, clase_tamano, drop = TRUE)]

control_clases = muestra[, .(
  Seleccionadas = .N,
  Respondientes = sum(respuesta),
  Tasa_ponderada = weighted.mean(respuesta, peso)
), by = clase]

if (any(control_clases$Respondientes == 0L)) {
  stop("Existe una clase sin respondientes; reagrupe las clases.")
}

gt(control_clases) |>
  fmt_number(columns = Tasa_ponderada, decimals = 3) |>
  tab_options(table.width = pct(75), table.font.size = px(12))

# Ajuste mediante clases de respuesta

diseno_rep = update(
  diseno_rep,
  respuesta = muestra$respuesta,
  clase = muestra$clase,
  internet_obs = muestra$internet_obs
)

diseno_clases = redistribute_weights(
  design = diseno_rep,
  reduce_if = respuesta == 0,
  increase_if = respuesta == 1,
  by = "clase"
)

estimacion_clases = svytotal(~internet_obs,
                              subset(diseno_clases, respuesta == 1))

resultados_clases = data.frame(
  Parametro = "Total de viviendas con internet",
  Estimacion = as.numeric(coef(estimacion_clases)),
  Varianza = as.numeric(vcov(estimacion_clases)),
  Error_estandar = as.numeric(SE(estimacion_clases)),
  CV = as.numeric(cv(estimacion_clases)) * 100,
  IC_inferior = confint(estimacion_clases)[1, 1],
  IC_superior = confint(estimacion_clases)[1, 2]
)

resultados_clases |>
  gt() |>
  fmt_number(
    columns = c(Estimacion, Varianza, Error_estandar,
                CV, IC_inferior, IC_superior),
    decimals = 2
  ) |>
  cols_label(
    Parametro = "Parámetro",
    Estimacion = "Estimación",
    Varianza = "Varianza",
    Error_estandar = "Error estándar",
    CV = "CV (%)",
    IC_inferior = "IC inferior",
    IC_superior = "IC superior"
  )

# Ajuste mediante propensión de respuesta

estimar_propension = function(w, datos) {
  positivos = w > 0
  modelo = glm(
    respuesta ~ municipio + log_M_c,
    data = datos,
    weights = w / mean(w[positivos]),
    family = quasibinomial(link = "logit")
  )
  rho = predict(modelo, newdata = datos, type = "response")
  if (any(!is.finite(rho))) stop("No se pudieron estimar las propensiones.")
  r = datos$respuesta == 1L
  sum((w[r] / rho[r]) * datos$internet_obs[r])
}

set.seed(405)
estimacion_propension = withReplicates(
  diseno_rep,
  theta = estimar_propension
)

resultados_propension = data.frame(
  Parametro = "Total de viviendas con internet",
  Estimacion = as.numeric(coef(estimacion_propension)),
  Varianza = as.numeric(vcov(estimacion_propension)),
  Error_estandar = as.numeric(SE(estimacion_propension)),
  CV = as.numeric(cv(estimacion_propension)) * 100,
  IC_inferior = confint(estimacion_propension)[1, 1],
  IC_superior = confint(estimacion_propension)[1, 2]
)

resultados_propension |>
  gt() |>
  fmt_number(
    columns = c(Estimacion, Varianza, Error_estandar,
                CV, IC_inferior, IC_superior),
    decimals = 2
  )

# Resultados de los ajustes

resumir = function(nombre, estimacion, diseno_replicado = diseno_rep) {
  p = as.numeric(coef(estimacion)) / N
  v = as.numeric(vcov(estimacion)) / N^2
  ee = sqrt(v)
  tcrit = qt(0.975, df = degf(diseno_replicado))
  data.table(
    Metodo = nombre,
    Estimacion = p,
    Varianza = v,
    EE = ee,
    CV = 100 * ee / abs(p),
    IC_inf = p - tcrit * ee,
    IC_sup = p + tcrit * ee,
    Error = p - P
  )
}

resultados = rbindlist(list(
  resumir("Clases de respuesta", estimacion_clases),
  resumir("Propensión logística", estimacion_propension)
))

gt(resultados[, .(Metodo, Estimacion, Varianza, EE, CV, IC_inf, IC_sup)]) |>
  fmt_number(columns = c(Estimacion, EE, IC_inf, IC_sup), decimals = 4) |>
  fmt_number(columns = Varianza, decimals = 7) |>
  fmt_number(columns = CV, decimals = 2) |>
  cols_label(IC_inf = "IC inferior", IC_sup = "IC superior") |>
  tab_options(table.width = pct(100), table.font.size = px(12))

# Comparación con los valores de referencia

comparacion = rbindlist(list(
  data.table(Metodo = "Población", Estimacion = P),
  data.table(Metodo = "Respuesta completa", Estimacion = P_completa),
  data.table(Metodo = "Respondientes sin ajuste", Estimacion = P_sin_ajuste),
  resultados[, .(Metodo, Estimacion)]
))
comparacion[, Error := Estimacion - P]

gt(comparacion) |>
  fmt_number(columns = c(Estimacion, Error), decimals = 4) |>
  tab_options(table.width = pct(75), table.font.size = px(12))
