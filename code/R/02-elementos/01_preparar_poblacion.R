
# Base de datos: establecimientos educativos de Colombia
#
# Fuente: Datos Abiertos Colombia.
# Conjunto de datos: Sedes Educativas.
#
# La base contiene información administrativa sobre establecimientos
# educativos oficiales y no oficiales de Colombia, incluyendo su
# ubicación geográfica, sector, calendario académico, matrícula total
# y número de sedes.
#
# Población de estudio:
# Registros de establecimientos educativos de Colombia correspondientes
# al año 2025.
#
# Unidad de muestreo:
# Cada registro de establecimiento educativo en la base de 2025.
#
# Variables cuantitativas:
# - total_matricula: número de estudiantes matriculados.
# - cantidad_sedes: número de sedes del establecimiento.
#
# Variables de clasificación:
# - departamento y municipio.
# - sector educativo.
# - calendario académico.

library(data.table)

# Descargar los registros correspondientes a 2025

url = paste0(
  "https://www.datos.gov.co/resource/cfw5-qzt5.csv",
  "?$where=a_o%3D2025&$limit=50000"
)

educacion = fread(url, encoding = "UTF-8")


# Seleccionar las variables necesarias para los ejercicios

variables = c(
  "codigo_dane",
  "cod_dane_departamento",
  "departamento",
  "cod_dane_municipio",
  "municipio",
  "nombre_establecimiento",
  "sector",
  "calendario",
  "total_matricula",
  "cantidad_sedes"
)

poblacion = educacion[, ..variables]


# Organizar los establecimientos por ubicación e identificador

setorder(
  poblacion,
  cod_dane_departamento,
  cod_dane_municipio,
  codigo_dane
)


# Asignar un identificador único a cada registro

poblacion[, id_unidad := .I]


# Dimensiones

print(dim(poblacion))


# Guardar la población en el directorio de datos del proyecto

saveRDS(
  poblacion,
  "datos/poblacion_educacion_2025.rds"
)
