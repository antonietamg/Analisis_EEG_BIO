############################################################
# ANALISIS DE SENALES EN BASE R
# Ruidos de colores, espectro de potencia, matriz canal-tiempo,
# beta en log-log y DFA
############################################################


############################################################
# FUNCIONES GENERALES
############################################################

beta_del_color <- function(color) {
  color <- tolower(color)
  
  if (color %in% c("blanco", "white")) {
    return(0)
  }
  
  if (color %in% c("rosa", "pink")) {
    return(1)
  }
  
  if (color %in% c("cafe", "café", "marron", "marrón", "brown", "rojo", "red")) {
    return(2)
  }
  
  if (color %in% c("azul", "blue")) {
    return(-1)
  }
  
  if (color %in% c("violeta", "morado", "purple", "violet")) {
    return(-2)
  }
  
  stop("Color no reconocido. Use: blanco, rosa, cafe, azul o violeta.")
}


crear_ruido_color <- function(n, color) {
  beta <- beta_del_color(color)
  
  x <- rnorm(n)
  X <- fft(x)
  
  k <- 0:(n - 1)
  k <- pmin(k, n - k)
  k[1] <- 1
  
  factor <- 1 / (k^(beta / 2))
  factor[1] <- 0
  
  y <- Re(fft(X * factor, inverse = TRUE)) / n
  
  y <- y - mean(y)
  y <- y / sd(y)
  
  return(y)
}


calcular_espectro <- function(x, fs) {
  x <- x[!is.na(x)]
  n <- length(x)
  
  x <- x - mean(x)
  
  ventana <- 0.5 - 0.5 * cos(2 * pi * (0:(n - 1)) / (n - 1))
  xw <- x * ventana
  
  X <- fft(xw)
  
  potencia <- (Mod(X)^2) / (fs * sum(ventana^2))
  
  m <- floor(n / 2) + 1
  frecuencia <- (0:(m - 1)) * fs / n
  
  potencia <- potencia[1:m]
  
  if (m > 2) {
    if (n %% 2 == 0) {
      indices <- 2:(m - 1)
    } else {
      indices <- 2:m
    }
    potencia[indices] <- 2 * potencia[indices]
  }
  
  potencia_db <- 10 * log10(potencia + .Machine$double.eps)
  
  resultado <- data.frame(
    frecuencia = frecuencia,
    potencia = potencia,
    potencia_db = potencia_db
  )
  
  return(resultado)
}


graficar_espectros_db <- function(datos, columnas, fs, rango_frec, colores_linea) {
  if (missing(rango_frec)) {
    rango_frec <- c(0, fs / 2)
  }
  
  if (missing(colores_linea)) {
    colores_linea <- c("black", "pink", "orange", "blue", 
                       "purple", "yellow", "brown", "cyan4")
  }
  
  colores_linea <- rep(colores_linea, length.out = length(columnas))
  
  espectros <- vector("list", length(columnas))
  
  minimo <- Inf
  maximo <- -Inf
  
  for (i in 1:length(columnas)) {
    esp <- calcular_espectro(datos[[columnas[i]]], fs)
    espectros[[i]] <- esp
    
    ok <- esp$frecuencia >= rango_frec[1] & esp$frecuencia <= rango_frec[2]
    
    minimo <- min(minimo, min(esp$potencia_db[ok]))
    maximo <- max(maximo, max(esp$potencia_db[ok]))
  }
  
  if (minimo == maximo) {
    minimo <- minimo - 1
    maximo <- maximo + 1
  }
  
  plot(rango_frec, c(minimo, maximo),
       type = "n",
       xlab = "Frecuencia (Hz)",
       ylab = "Potencia (dB)",
       main = "Espectro de potencia",
       ylim = c(minimo, maximo))
  
  for (i in 1:length(columnas)) {
    esp <- espectros[[i]]
    ok <- esp$frecuencia >= rango_frec[1] & esp$frecuencia <= rango_frec[2]
    
    lines(esp$frecuencia[ok], esp$potencia_db[ok],
          col = colores_linea[i],
          lwd = 2)
  }
  
  legend("topright",
         legend = columnas,
         col = colores_linea,
         lwd = 2,
         bty = "n")
  
  return(espectros)
}


############################################################
# FUNCION: MATRIZ BANDA x TIEMPO PARA UN SOLO CANAL
############################################################

matriz_bandas_tiempo <- function(datos, canal, fs, ventana_seg, paso_seg, bandas) {
  
  if (length(canal) != 1) {
    stop("Debe elegir un solo canal.")
  }
  
  if (!(canal %in% names(datos))) {
    stop("El canal elegido no existe en el dataframe.")
  }
  
  x <- datos[[canal]]
  n <- length(x)
  
  ventana <- round(ventana_seg * fs)
  paso <- round(paso_seg * fs)
  
  if (ventana > n) {
    stop("La ventana es más larga que la señal.")
  }
  
  if (paso < 1) {
    paso <- 1
  }
  
  inicios <- seq(1, n - ventana + 1, by = paso)
  
  puntos_centro <- inicios + (ventana - 1) / 2
  tiempos_centro_seg <- (puntos_centro - 1) / fs
  
  n_bandas <- nrow(bandas)
  
  nombres_bandas <- rownames(bandas)
  
  if (is.null(nombres_bandas)) {
    nombres_bandas <- paste("Banda", 1:n_bandas)
  }
  
  matriz_db <- matrix(NA,
                      nrow = n_bandas,
                      ncol = length(inicios))
  
  rownames(matriz_db) <- nombres_bandas
  colnames(matriz_db) <- round(tiempos_centro_seg, 3)
  
  for (j in 1:length(inicios)) {
    
    ind <- inicios[j]:(inicios[j] + ventana - 1)
    
    x_ventana <- x[ind]
    
    esp <- calcular_espectro(x_ventana, fs)
    
    for (i in 1:n_bandas) {
      
      f1 <- bandas[i, 1]
      f2 <- bandas[i, 2]
      
      if (f2 > fs / 2) {
        f2 <- fs / 2
      }
      
      ok <- esp$frecuencia >= f1 & esp$frecuencia <= f2
      
      if (sum(ok) > 0) {
        potencia_promedio <- mean(esp$potencia[ok])
        matriz_db[i, j] <- 10 * log10(potencia_promedio + .Machine$double.eps)
      }
    }
  }
  
  resultado <- list(
    matriz_db = matriz_db,
    tiempos_seg = tiempos_centro_seg,
    puntos = puntos_centro,
    canal = canal,
    fs = fs,
    bandas = bandas,
    ventana_seg = ventana_seg,
    paso_seg = paso_seg
  )
  
  return(resultado)
}



############################################################
# FUNCION: GRAFICAR MATRIZ BANDA x TIEMPO
############################################################

graficar_matriz_bandas <- function(resultado,
                                   eje_x = "segundos",
                                   titulo = "",
                                   onsets_seg = NULL,
                                   limites_db = NULL) {
  
  matriz_db <- resultado$matriz_db
  
  eje_x <- tolower(eje_x)
  
  if (eje_x %in% c("segundos", "segundo", "s", "seg")) {
    
    x <- resultado$tiempos_seg
    etiqueta_x <- "Tiempo (s)"
    
    if (is.null(onsets_seg)) {
      onsets_x <- NULL
    } else {
      onsets_x <- onsets_seg
    }
  }
  
  else if (eje_x %in% c("puntos", "punto", "muestras", "muestra")) {
    
    x <- resultado$puntos
    etiqueta_x <- "Puntos / muestras"
    
    if (is.null(onsets_seg)) {
      onsets_x <- NULL
    } else {
      onsets_x <- round(onsets_seg * resultado$fs) + 1
    }
  }
  
  else if (eje_x %in% c("puntos0", "muestras0")) {
    
    x <- resultado$puntos - 1
    etiqueta_x <- "Puntos / muestras desde 0"
    
    if (is.null(onsets_seg)) {
      onsets_x <- NULL
    } else {
      onsets_x <- round(onsets_seg * resultado$fs)
    }
  }
  
  else {
    stop("eje_x debe ser: 'segundos', 'puntos' o 'puntos0'.")
  }
  
  if (titulo == "") {
    titulo <- paste("Canal", resultado$canal, "- potencia por bandas")
  }
  
  z <- t(matriz_db)
  
  if (is.null(limites_db)) {
    zmin <- min(z, na.rm = TRUE)
    zmax <- max(z, na.rm = TRUE)
  } else {
    zmin <- limites_db[1]
    zmax <- limites_db[2]
  }
  
  if (zmin == zmax) {
    zmin <- zmin - 1
    zmax <- zmax + 1
  }
  
  colores <- colorRampPalette(c("blue", "cyan", "yellow", "red"))(100)
  
  par_anterior <- par(no.readonly = TRUE)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(4, 0.6))
  
  par(mar = c(5, 8, 4, 2))
  
  image(x = x,
        y = 1:nrow(matriz_db),
        z = z,
        col = colores,
        zlim = c(zmin, zmax),
        xlab = etiqueta_x,
        ylab = "",
        yaxt = "n",
        main = titulo)
  
  axis(2,
       at = 1:nrow(matriz_db),
       labels = rownames(matriz_db),
       las = 1)
  
  if (!is.null(onsets_x)) {
    
    abline(v = onsets_x,
           col = "black",
           lwd = 2,
           lty = 2)
    
    text(x = onsets_x,
         y = rep(nrow(matriz_db), length(onsets_x)),
         labels = paste("A", 1:length(onsets_x), sep = ""),
         pos = 3,
         cex = 0.8,
         xpd = TRUE)
  }
  
  par(mar = c(5, 1, 4, 4))
  
  y_leyenda <- seq(zmin, zmax, length.out = 100)
  z_leyenda <- matrix(rep(y_leyenda, each = 2), nrow = 2)
  
  image(x = c(0, 1),
        y = y_leyenda,
        z = z_leyenda,
        col = colores,
        zlim = c(zmin, zmax),
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = "",
        main = "dB")
  
  axis(4, las = 1)
  
  layout(1)
  par(par_anterior)
}


graficar_loglog_beta <- function(datos, columnas, fs, rango_beta, colores_linea) {
  if (missing(rango_beta)) {
    rango_beta <- c(1, fs / 2)
  }
  
  if (missing(colores_linea)) {
    colores_linea <- c("black", "pink", "orange", "blue", 
                       "purple", "yellow", "brown", "cyan4")
  }
  
  colores_linea <- rep(colores_linea, length.out = length(columnas))
  
  espectros <- vector("list", length(columnas))
  
  fmin <- Inf
  fmax <- -Inf
  pmin <- Inf
  pmax <- -Inf
  
  for (i in 1:length(columnas)) {
    esp <- calcular_espectro(datos[[columnas[i]]], fs)
    
    ok <- esp$frecuencia >= rango_beta[1] &
      esp$frecuencia <= rango_beta[2] &
      esp$potencia > 0
    
    esp <- esp[ok, ]
    
    espectros[[i]] <- esp
    
    fmin <- min(fmin, min(esp$frecuencia))
    fmax <- max(fmax, max(esp$frecuencia))
    pmin <- min(pmin, min(esp$potencia))
    pmax <- max(pmax, max(esp$potencia))
  }
  
  plot(c(fmin, fmax), c(pmin, pmax),
       type = "n",
       log = "xy",
       xlab = "Frecuencia (Hz)",
       ylab = "Potencia",
       main = "Espectro log-log y beta")
  
  resultados <- data.frame(
    senal = columnas,
    pendiente_beta = NA,
    beta_positivo_1_sobre_f = NA,
    pendiente_db_por_decada = NA
  )
  
  leyenda <- character(length(columnas))
  
  for (i in 1:length(columnas)) {
    esp <- espectros[[i]]
    
    modelo <- lm(log10(esp$potencia) ~ log10(esp$frecuencia))
    
    intercepto <- coef(modelo)[1]
    pendiente <- coef(modelo)[2]
    
    resultados$pendiente_beta[i] <- pendiente
    resultados$beta_positivo_1_sobre_f[i] <- -pendiente
    resultados$pendiente_db_por_decada[i] <- 10 * pendiente
    
    lines(esp$frecuencia,
          esp$potencia,
          col = colores_linea[i],
          lwd = 2)
    
    xx <- exp(seq(log(min(esp$frecuencia)),
                  log(max(esp$frecuencia)),
                  length.out = 100))
    
    yy <- 10^(intercepto + pendiente * log10(xx))
    
    lines(xx,
          yy,
          col = colores_linea[i],
          lwd = 2,
          lty = 2)
    
    leyenda[i] <- paste(columnas[i],
                        "pendiente =",
                        round(pendiente, 3))
  }
  
  legend("topright",
         legend = leyenda,
         col = colores_linea,
         lwd = 2,
         bty = "n")
  
  return(resultados)
}


calcular_dfa <- function(x, escalas) {
  x <- x[!is.na(x)]
  n <- length(x)
  
  perfil <- cumsum(x - mean(x))
  
  F <- rep(NA, length(escalas))
  
  for (i in 1:length(escalas)) {
    s <- escalas[i]
    n_bloques <- floor(n / s)
    
    if (n_bloques < 2) {
      next
    }
    
    rms <- rep(NA, n_bloques)
    
    for (j in 1:n_bloques) {
      ind <- ((j - 1) * s + 1):(j * s)
      tt <- 1:s
      
      modelo <- lm(perfil[ind] ~ tt)
      tendencia <- predict(modelo)
      
      residuo <- perfil[ind] - tendencia
      
      rms[j] <- sqrt(mean(residuo^2))
    }
    
    F[i] <- sqrt(mean(rms^2))
  }
  
  ok <- !is.na(F) & F > 0
  
  escalas_ok <- escalas[ok]
  F_ok <- F[ok]
  
  modelo_final <- lm(log10(F_ok) ~ log10(escalas_ok))
  
  alpha <- coef(modelo_final)[2]
  intercepto <- coef(modelo_final)[1]
  
  resultado <- list(
    escalas = escalas_ok,
    F = F_ok,
    alpha = alpha,
    intercepto = intercepto
  )
  
  return(resultado)
}


graficar_dfa <- function(datos, columnas, escalas, colores_linea) {
  if (missing(colores_linea)) {
    colores_linea <- c("black", "pink", "orange", "blue", 
                       "purple", "yellow", "brown", "cyan4")
  }
  
  colores_linea <- rep(colores_linea, length.out = length(columnas))
  
  resultados <- data.frame(
    senal = columnas,
    alpha = NA
  )
  
  par_anterior <- par(no.readonly = TRUE)
  
  if (length(columnas) == 1) {
    par(mfrow = c(1, 1))
  } else {
    par(mfrow = c(ceiling(length(columnas) / 2), 2))
  }
  
  for (i in 1:length(columnas)) {
    x <- datos[[columnas[i]]]
    
    dfa <- calcular_dfa(x, escalas)
    
    resultados$alpha[i] <- dfa$alpha
    
    plot(dfa$escalas,
         dfa$F,
         log = "xy",
         pch = 16,
         col = colores_linea[i],
         xlab = "Escala",
         ylab = "F(s)",
         main = paste("DFA:", columnas[i]))
    
    xx <- exp(seq(log(min(dfa$escalas)),
                  log(max(dfa$escalas)),
                  length.out = 100))
    
    yy <- 10^(dfa$intercepto + dfa$alpha * log10(xx))
    
    lines(xx,
          yy,
          col = colores_linea[i],
          lwd = 2)
    
    legend("topleft",
           legend = paste("alpha =", round(dfa$alpha, 3)),
           bty = "n")
  }
  
  par(par_anterior)
  
  return(resultados)
}



############################################################
# 1) CREAR UNO O VARIOS RUIDOS DE COLORES
############################################################

set.seed(123)

# Frecuencia de muestreo en Hz
fs <- 250

# Longitud de la señal
duracion_seg <- 30
n_muestras <- fs * duracion_seg

# También se puede fijar directamente:
# n_muestras <- 5000

# Tipos posibles:
# "blanco", "rosa", "cafe", "azul", "violeta"

tipos_ruido <- c("blanco", "rosa", "cafe", "azul", "violeta")

nombres_senales <- c("ruido_blanco",
                     "ruido_rosa",
                     "ruido_cafe",
                     "ruido_azul",
                     "ruido_violeta")

if (length(tipos_ruido) != length(nombres_senales)) {
  stop("tipos_ruido y nombres_senales deben tener la misma longitud.")
}

tiempo <- (0:(n_muestras - 1)) / fs

datos <- data.frame(tiempo = tiempo)

for (i in 1:length(tipos_ruido)) {
  datos[[nombres_senales[i]]] <- crear_ruido_color(n_muestras, tipos_ruido[i])
}

# Ver las primeras filas
head(datos)



############################################################
# 2) CALCULAR Y GRAFICAR ESPECTRO DE POTENCIA EN dB
############################################################

# Para una sola señal:
# columnas_espectro <- c("ruido_rosa")

# ESPECIFICAR NOMBRE DEL DATAFRAME 
# datos <- whole_brain_p11_3_RAW_ds100Hz_bp_1.30Hz

# ESTABLECER FRECUENCIA DE MUESTREO
# fs = fs

# Para varias señales:
columnas_espectro <- c("ruido_blanco", "ruido_rosa", "ruido_cafe")

rango_frecuencia_grafica <- c(0, fs / 2)

espectros <- graficar_espectros_db(datos,
                                   columnas_espectro,
                                   fs,
                                   rango_frecuencia_grafica)



############################################################
# 3) MATRIZ BANDA x TIEMPO PARA UN SOLO CANAL
############################################################

# Canal que quieres analizar.
canal_matriz <- "ruido_rosa"

# Elegir el eje X:
# "segundos" = tiempo en segundos
# "puntos"   = número de muestra empezando en 1
# "puntos0"  = número de muestra empezando en 0
eje_x_matriz <- "segundos"

# Ventana y paso para el cálculo de potencia.
# Ventana más larga da mejor resolución en frecuencia.
# Paso más pequeño da más resolución temporal.
ventana_seg <- 2
paso_seg <- 0.25

# Bandas de frecuencia.
# Como antes estabas usando 1.6 Hz como límite inferior,
# dejo delta desde 1.6 hasta 4 Hz.
bandas_eeg <- matrix(c(1.6, 4,
                       4,   8,
                       8,   13,
                       13,  30),
                     ncol = 2,
                     byrow = TRUE)

rownames(bandas_eeg) <- c("Delta 1.6-4 Hz",
                          "Theta 4-8 Hz",
                          "Alpha 8-13 Hz",
                          "Beta 13-30 Hz")

colnames(bandas_eeg) <- c("f_min", "f_max")

# Onsets de los ataques en segundos.
# Si no quieres marcarlos, usa:
 onsets_seg <- NULL

# onsets_seg <- c(1517, 4397, 7201, 9243)

resultado_matriz <- matriz_bandas_tiempo(datos,
                                         canal_matriz,
                                         fs,
                                         ventana_seg,
                                         paso_seg,
                                         bandas_eeg)

graficar_matriz_bandas(resultado_matriz,
                       eje_x_matriz,
                       paste("Canal", canal_matriz, "- potencia por bandas"),
                       onsets_seg)



############################################################
# 4) ESPECTRO LOG-LOG Y COEFICIENTE BETA
############################################################

columnas_beta <- c("ruido_blanco", "ruido_rosa", "ruido_cafe",
                   "ruido_azul", "ruido_violeta")

# No usar 0 Hz para log-log.
rango_beta <- c(1, fs / 2)

resultado_beta <- graficar_loglog_beta(datos,
                                       columnas_beta,
                                       fs,
                                       rango_beta)

print(resultado_beta)



############################################################
# 5) DFA DE CADA SEÑAL Y VALOR ALPHA
############################################################

columnas_dfa <- c("ruido_blanco", "ruido_rosa", "ruido_cafe",
                  "ruido_azul", "ruido_violeta")

# Escalas usadas para DFA.
# Se usan valores entre 4 muestras y n/4.
escalas_dfa <- unique(round(exp(seq(log(4),
                                    log(floor(n_muestras / 4)),
                                    length.out = 20))))

resultado_dfa <- graficar_dfa(datos,
                              columnas_dfa,
                              escalas_dfa)

print(resultado_dfa)
