# Ejemplo No respuesta de ítem: imputación hot-deck ponderada

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

# Generación de no respuesta de ítem

set.seed(505)
muestra[, rho_item := plogis(
  qlogis(0.82) + 0.34 * energia -
    0.18 * (as.integer(municipio) - mean(as.integer(municipio)))
)]
muestra[, respuesta_item := rbinom(.N, 1L, rho_item)]
muestra[, internet_obs := fifelse(
  respuesta_item == 1L, internet, NA_integer_
)]

c(
  Muestra = nrow(muestra),
  Items_observados = sum(muestra$respuesta_item),
  Items_faltantes = sum(muestra$respuesta_item == 0L),
  Tasa_respuesta_item = mean(muestra$respuesta_item)
)

# Generación de no respuesta de ítem

gt(head(muestra[, .(COD_TRI, COD_CUESTIONARIO, municipio,
                    acueducto, energia, internet_obs,
                    respuesta_item, peso)], 10)) |>
  fmt_number(columns = peso, decimals = 3) |>
  tab_options(table.width = pct(100), table.font.size = px(12))

# Celdas de imputación

muestra[, celda := interaction(municipio, energia, drop = TRUE)]

control_celdas = muestra[, .(
  Viviendas = .N,
  Donantes = sum(respuesta_item),
  Receptores = sum(respuesta_item == 0L)
), by = celda]

if (any(control_celdas$Donantes[control_celdas$Receptores > 0] < 2L)) {
  stop("Alguna celda tiene menos de dos donantes; se deben reagrupar las celdas.")
}

gt(control_celdas) |>
  tab_options(table.width = pct(70), table.font.size = px(12))

# Imputación hot-deck aleatoria ponderada

set.seed(506)

donante = seq_len(nrow(muestra))
faltantes = which(muestra$respuesta_item == 0L)
observados = which(muestra$respuesta_item == 1L)

for (c in unique(as.character(muestra$celda[faltantes]))) {
  receptores = faltantes[as.character(muestra$celda[faltantes]) == c]
  elegibles = observados[as.character(muestra$celda[observados]) == c]
  seleccion = sample.int(
    length(elegibles),
    length(receptores),
    replace = TRUE,
    prob = muestra$peso[elegibles]
  )
  donante[receptores] = elegibles[seleccion]
}

muestra[, `:=`(
  donante = donante,
  internet_imp = internet_obs[donante]
)]

stopifnot(!anyNA(muestra$internet_imp),
          all(muestra$celda[faltantes] == muestra$celda[donante[faltantes]]))

gt(head(muestra[respuesta_item == 0L,
                .(COD_CUESTIONARIO, celda, internet_obs,
                  internet_imp, peso, donante)], 10)) |>
  fmt_number(columns = peso, decimals = 3) |>
  tab_options(table.width = pct(95), table.font.size = px(12))

# Estimación con valores imputados

P_imputada = sum(muestra$peso * muestra$internet_imp) / N
c(Referencia = P, Respuesta_completa = P_completa,
  Imputacion = P_imputada)

# Información requerida para la varianza

varianza_celda = muestra[respuesta_item == 1L, {
  w = peso
  y = internet_obs
  mu = weighted.mean(y, w)
  .(sigma2 = sum(w * (y - mu)^2) / (sum(w) - sum(w^2) / sum(w)))
}, by = celda]

if (any(!is.finite(varianza_celda$sigma2))) {
  stop("No es posible estimar la varianza residual en alguna celda.")
}

gt(varianza_celda) |>
  fmt_number(columns = sigma2, decimals = 5) |>
  tab_options(table.width = pct(60), table.font.size = px(12))

# Corrección de varianza de Brick, Kalton y Kim (2004)

varianza_hotdeck = function(datos, diseno_rep, varianzas) {
  w = datos$peso
  W = as.matrix(weights(diseno_rep, type = "analysis"))
  c = diseno_rep$scale * diseno_rep$rscales
  D = sweep(W, 1, w, "-")

  y = datos$internet_imp
  m = which(datos$respuesta_item == 0L)
  don = datos$donante

  # Varianza del diseño aplicada a los valores imputados.
  V0 = sum(c * colSums(D * y)^2)

  # Diagonal de los coeficientes Omega de la forma cuadrática.
  omega_diag = rowSums(sweep(D^2, 2, c, "*"))

  B0 = 0
  Vimp = 0
  Vmix = 0
  grupos_donante = split(m, don[m])

  for (idx in grupos_donante) {
    j = don[idx[1]]
    celda_j = as.character(datos$celda[j])
    sigma2 = varianzas$sigma2[match(celda_j, as.character(varianzas$celda))]

    z = colSums(D[idx, , drop = FALSE])
    a = sum(w[idx])
    s2w = sum(w[idx]^2)

    # Lema 1: sesgo del estimador V0.
    B0 = B0 + sigma2 * (
      2 * sum(c * D[j, ] * z) +
        sum(c * z^2) - sum(omega_diag[idx])
    )

    # Lemas 2 y 3: componentes de imputación y mixto.
    Vimp = Vimp + sigma2 * (a^2 + s2w)
    Vmix = Vmix + sigma2 * (w[j] * a - s2w)
  }

  Vsam = V0 - B0
  Vtotal = Vsam + Vimp + 2 * Vmix

  list(
    V0 = V0,
    Correccion = B0,
    V_muestreo = Vsam,
    V_imputacion = Vimp,
    V_mixto = Vmix,
    V_total = Vtotal
  )
}

componentes = varianza_hotdeck(muestra, diseno_rep, varianza_celda)
if (!is.finite(componentes$V_total) || componentes$V_total <= 0) {
  stop("La varianza total no es positiva; revisar las condiciones de imputación.")
}

as.data.table(componentes)

# Precisión de la proporción imputada

V_P = componentes$V_total / N^2
EE = sqrt(V_P)
CV = 100 * EE / abs(P_imputada)
tcrit = qt(0.975, df = degf(diseno))

resultados = data.table(
  Parametro = "Viviendas con internet",
  Estimacion = P_imputada,
  Varianza = V_P,
  EE = EE,
  CV = CV,
  IC_inf = P_imputada - tcrit * EE,
  IC_sup = P_imputada + tcrit * EE
)

gt(resultados) |>
  fmt_number(columns = c(Estimacion, EE, IC_inf, IC_sup), decimals = 4) |>
  fmt_number(columns = Varianza, decimals = 7) |>
  fmt_number(columns = CV, decimals = 2) |>
  cols_label(IC_inf = "IC inferior", IC_sup = "IC superior") |>
  tab_options(table.width = pct(100), table.font.size = px(12))

# Comparación con la población

comparacion = data.table(
  Metodo = c("Población", "Respuesta completa", "Hot-deck ponderado"),
  Estimacion = c(P, P_completa, P_imputada)
)
comparacion[, Error := Estimacion - P]

gt(comparacion) |>
  fmt_number(columns = c(Estimacion, Error), decimals = 4) |>
  tab_options(table.width = pct(75), table.font.size = px(12))
