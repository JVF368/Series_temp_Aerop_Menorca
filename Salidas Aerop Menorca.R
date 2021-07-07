# CARGA Y MANIPULACIÓN PREVIA -------------------------------------------------
require(readxl)
require(tidyverse)
require(lubridate)
require(stringr)
require(forecast)
require(tseries)
require(TSA)
require(KFAS)

## 1 ## DATOS, LECTURA Y MANIPULACIÓN -----------------------------------------
# Lectura
ser <- read_excel("Vanrell_jordi_dat.xlsx", sheet = "Hoja1")

# Manipulación para homogeneización
# 1a parte (hasta 31/12/2005)
ser_1 <- ser[1:3287,] %>%
  mutate(t = gsub(" ", "", ...1)) %>%
  separate(t, into = c("resto", "dia_semana"), sep = -3, convert = T) %>%
  separate(resto, into = c("resto", "anyo"), sep = -4, convert = T) %>%
  separate(resto, into = c("dia", "mes"), sep = 2, convert = T)

ser_1$mes <- factor(ser_1$mes)
levels(ser_1$mes) <- c("04", "08", "12", "01", "02", "07", 
                       "06", "03", "05", "11", "10", "09")

ser_1 <- ser_1 %>%
  mutate(t = dmy(paste0(dia, "/", mes, "/", anyo))) %>%
  select(t, `salidas menorca`)

# 2a parte (01/01/2006 - 31/12/2007)
ser_2 <- ser[3288:4017,] %>%
  mutate(t = dmy(...1)) %>%
  select(t, `salidas menorca`)

# 3a parte (01/01/2008 - 31/08/2009)
ser_3 <- ser[4018:4626,] %>%
  mutate(t = seq(dmy('01-01-2008'), dmy('31-08-2009'), by = '1 day')) %>%
  select(t, `salidas menorca`)

# Agregación
ser <- rbind(ser_1, ser_2, ser_3)
rm(ser_1, ser_2, ser_3)

plot(x = ser$t, y = ser$`salidas menorca`, ylab = 'Pasajeros', xlab = 'Fecha', 
     type = "l", main = 'Salidas del aeropuerto de Menorca')

pascua = c(dmy('31-03-1997'), dmy('13-04-1998'), dmy('05-04-1999'),
           dmy('24-04-2000'), dmy('16-04-2001'), dmy('01-04-2002'), 
           dmy('21-04-2003'), dmy('12-04-2004'), dmy('28-03-2005'), 
           dmy('17-04-2006'), dmy('09-04-2007'), dmy('24-03-2008'), 
           dmy('13-04-2009'))

ser <- ser %>%
  mutate(festivo = ifelse(t %in% pascua, 1, 0)) %>%
  select(t, `salidas menorca`, festivo)
rm(pascua)

plot(x = ser$t, y = log(ser$`salidas menorca`), ylab = 'Pasajeros', xlab = 'Fecha', 
     type = "l", main = 'Salidas del aeropuerto de Menorca')


## 2 ## TRANSF. EN SERIE TEMPORAL, TRATAMIENTO DE OUTLIERS Y ESTACIONALIDAD ---

# Eliminación de outliers:
ser$`salidas menorca`[ser$t == "2004-12-29" | ser$t == "2005-01-26"] <- NA
# Construcción de la serie temporal completa
ser_ts <- ts(as.matrix(ser[,2], ncol = 1), start = c(1997, 1), frequency = 365)
# Rellena outliers por interpolación
ser_ts <- na.interp(ser_ts)
ts.plot(log(ser_ts), ylab = 'log(Pasajeros)', xlab = 'Fecha', 
        main = 'Salidas del aeropuerto de Menorca (outlier-free)')

# Separación en train y test
test <- ser_ts[3896:4626] # obs 3896 es el 01/09/2007, hasta el final; 731 observaciones
train <- ser_ts[1:4261] # obs 4261 es el 31 de agosto de 2008
# Entre las observaciones 3896 y 4261 hay 366 días (2008 fue bisiesto)

# Puntos clave:
# obs 3896 es el 01/09/2007; c(2007, 244)
# obs 4261 es el 31/08/2008; c(2008, 244) * 2008 fue bisiesto
# obs 4262 es el 01/09/2008; c(2008, 245) * 2008 fue bisiesto
# obs 4626 es el 31/08/2009; c(2009, 243)

# Construcción de la serie temporal de entrenamiento
train <- ts(as.matrix(train, ncol = 1), start = c(1997, 1), frequency = 365) # datos diarios: frequency = 365.

# Espectros
par(mfrow = c(1, 3))
spectrum(log(train))
spectrum(log(train), method = "ar")
spectrum(log(train), kernel("fejer", 100, r = 7)) # estacionalidad semanal

# Construcción de la parte determinista de la serie de entrenamiento
t <- seq(1, length(train))

det <- cos(0*t)
for (j in 1:182) {
  det = cbind(det, cos(2*pi*t*j/365), sin(2*pi*t*j/365))
  rm(j)
  }
end1 <- log(train)
exo1 <- cbind(det, t, t^2)
rho <- solve(t(exo1) %*% exo1) %*% (t(exo1) %*% end1) # MCO

rho_f <- rho
train_e <- log(train) - exo1 %*% rho # parte estacionaria
train_f <- exo1 %*% rho # parte fija estacional (train)
train_e <- ts(train_e, start = c(1997, 1), frequency = 365) # se convierte en ts
train_f <- ts(train_f, start = c(1997, 1), frequency = 365) # se convierte en ts

par(mfrow = c(2, 1))
ts.plot(train_e, main = "Parte estacionaria")
ts.plot(log(train), train_f, col = 1:2, 
        main = "Parte determinista sobre serie train")

# Espectros
par(mfrow = c(1, 3))
spectrum(train_e)
spectrum(train_e, method = "ar")
spectrum(train_e, kernel("fejer", 100, r = 7))

# Correlogramas
par(mfrow = c(1, 2))
acf(train_e, lag.max = 28)
pacf(train_e, lag.max = 28)

# Tomamos diferencias en el retardo 7
par(mfrow = c(1, 2))
acf(diff(train_e, 7), lag.max = 28)
pacf(diff(train_e, 7), lag.max = 28)


## 3 ## ARIMA + SARIMA --------------------------------------------------------

# Interpretación, dados los ACF y PACF
# En la ACF se aprecia que existe un AR (primeros retardos descienden exp.)
# Además, hay un coeficiente muy significativo en el retardo 7 y efectos satélite.
# El 7º retardo en la PACF (y el descenso exponencial cada 7 retardos) sugiere 
# que hay un MA estacional.

## 3.1 ## Modelo 1
# Modelización de ARIMA(1,0,0) con SARIMA (0, 1, 1) y dummies
mod1 <- arima(train_e, xreg = ser$festivo[1:4261], order = c(1, 0, 0), 
              seasonal = list(order = c(0, 1, 1), period = 7))
mod1

par(mfrow = c(2, 2))
acf(mod1$residuals, lag.max = 28, main = "Residuos")
pacf(mod1$residuals, lag.max = 28, main = "Residuos")
spectrum(mod1$residuals, main = "Periodograma")
spectrum(mod1$residuals, kernel("fejer", 100, r = 7), 
         main = "Periodograma suavizado")


## 3.2 ## Modelo 2
# Modelización de ARMA(1,0,1) con primeras diferencias estacionales + media móvil estacional
mod2 <- arima(train_e, order = c(1, 0, 0), 
              seasonal = list(order = c(0, 1, 1), period = 7))
mod2 # No funciona tan bien (aic)

par(mfrow = c(2, 2))
acf(mod2$residuals, lag.max = 28, main = "Residuos")
pacf(mod2$residuals, lag.max = 28, main = "Residuos")
spectrum(mod2$residuals, main = "Periodograma")
spectrum(mod2$residuals, kernel("fejer", 100, r = 7), 
         main = "Periodograma suavizado")

## 3.3 ## Modelo 3
# Modelización de ARIMA(1,0,0) con SARIMA (0, 1, 1) y dummies
mod3 <- arima(train_e, xreg = ser$festivo[1:4261], order = c(1, 0, 0), 
              seasonal = list(order = c(1, 1, 2), period = 7))
mod3

par(mfrow = c(2, 2))
acf(mod3$residuals, lag.max = 28, main = "Residuos")
pacf(mod3$residuals, lag.max = 28, main = "Residuos")
spectrum(mod3$residuals, main = "Periodograma")
spectrum(mod3$residuals, kernel("fejer", 100, r = 7), main = "Periodograma suavizado")

## 3.4 ## Valoración de los modelos 1, 2 y 3
# Función de valoración común para modelos ARIMA y representación gráfica

arimafunct <- function(modelo){
  # Función para la estimación del error en el ajuste de modelos ARIMA + SARIMA
  # y la extracción de series para su representación
  # INPUTS:
  # modelo: modelo ARIMA + SARIMA propuesto
  # OUTPUTS:
  # Objeto tipo lista en el que se almacenan las series temporales para 
  # su representación y medidas de error del modelo.
  fitmod <- fitted(modelo) # Valores modelizados (parte estacionaria) hasta obs. 4261
  ffitmod <- fitmod + train_f # suma de la parte determinista y los valores del modelo de la parte estacionaria
  ffitmod_e <- log(train) - ffitmod # errores que no recoge el modelo
  
  t <- seq(1, length(ser_ts))
  det_part = cos(0*t)
  for (j in 1:182) { # parte determinista
    det_part <- cbind(det_part, cos(2*pi*t*j/365), sin(2*pi*t*j/365))
    rm(j)
  }
  exo1p <- cbind(det_part, t, t^2)
  for_det <- exo1p %*% rho_f # prediccion parte determinista
  for_det2 <- for_det[4262:4626] # parte determinista del final de la muestra test
  
  # La función escoge la forma de predecir de acuerdo con si hay o no dummies en 
  # el modelo
  if (sum(str_detect(as.character(modelo$call), "festivo")) == 1)
    {predmod <- predict(modelo, newxreg = ser$festivo[4262:4626], n.ahead = 365, se.fit = T)}
  else 
    {predmod <- predict(modelo, n.ahead = 365, se.fit = T)}
  
  # Series para representación de las predicciones con intervalos de confianza
  predict_st <- ts(c(rep(NA, length(train_e)), predmod$pred), start = c(1997, 1), frequency = 365)
  upper <- ts(c(rep(NA, length(train_e)), predmod$pred+2*predmod$se), start = c(1997, 1), frequency = 365) 
  lower <- ts(c(rep(NA, length(train_e)), predmod$pred-2*predmod$se), start = c(1997, 1), frequency = 365) 
  observed_st <- ts(c(train_e[3896:4261], rep(NA, length(predmod$pred))), frequency = 365, start = c(2007, 244), end = c(2009, 244))
  
  for_tot <- predmod$pred + for_det2 
  # valores predichos en la serie en niveles entre 4262 y 4626
  
  # Series para la representación de las series predichas y observadas en log.
  predict <- ts(c(rep(NA, 366), for_tot), start = c(2007, 244), frequency = 365)
  observed <- ts(log(test), start = c(2007, 244), frequency = 365)
  
  # Error en la muestra de train
  err_tr <- as.vector(ffitmod) - as.vector(log(train))
  ECM_tr <- (t(err_tr) %*% err_tr)/length(log(train)) #ECM 
  EAM_tr <- (sum(abs(err_tr)))/length(log(train)) #EAM
  
  # Error en la muestra de test
  err_tst <- na.omit(log(test) - predict)
  ECM_tst <- (t(err_tst) %*% err_tst)/length(na.omit(log(test) - predict)) #ECM 
  EAM_tst <- (sum(abs(err_tst)))/length(na.omit(log(test) - predict)) #EAM
  
  statistics <- c(ECM_tr = ECM_tr, EAM_tr = EAM_tr, ECM_tst = ECM_tst, EAM_tst = EAM_tst)
  grafseries <- list(fitmod = fitmod, ffitmod = ffitmod, ffitmod_e = ffitmod_e, 
                     for_det = for_det, predict_st = predict_st, upper = upper, 
                     lower = lower, observed_st = observed_st, predict = predict, 
                     observed = observed, predmod = predmod)
  results <- list(accuracy = statistics, grafseries = grafseries)
  return(results)
}

## 3.4.1 ## Resultados del Modelo 1
results_mod1 <- arimafunct(mod1) # Llama a la función anterior

# Representación gráfica de la predicción sobre la muestra de entrenamiento
par(mfrow = c(2, 2))
ts.plot(log(ser_ts), results_mod1$grafseries$for_det, col = 1:2, main = "Ajuste de la parte determinista (serie completa)")
ts.plot(train_e, results_mod1$grafseries$fitmod, col = 1:2, main = "Ajuste de la parte estacionaria (entrenamiento)")
ts.plot(log(train), results_mod1$grafseries$ffitmod, col = 1:2, main = "Ajuste agregado (entrenamiento)") # gráfico de los valores predichos sobre la serie original
ts.plot(results_mod1$grafseries$ffitmod_e, ylab = "", main = "Errores del modelo (entrenamiento)")

# Plot of actual and forecasted values (parte estacionaria)
par(mfrow = c(1, 2))
plot(results_mod1$grafseries$observed_st, type = "l", xlab = "", ylab = "",
     main = "Serie y predicción (parte estacionaria)", ylim = c(min(results_mod1$grafseries$predmod$pred-2*results_mod1$grafseries$predmod$se), max(results_mod1$grafseries$predmod$pred+2*results_mod1$grafseries$predmod$se))) 
lines(results_mod1$grafseries$predict_st, col = "blue", lty = 2)
lines(results_mod1$grafseries$lower, col = "red", lty = 5)
lines(results_mod1$grafseries$upper, col = "red",lty = 5)
ts.plot(results_mod1$grafseries$observed, results_mod1$grafseries$predict, 
        col = 1:2, main = "Ajuste sobre la serie test")

# ECM_tr
results_mod1$accuracy[[1]]
# EAM_tr
results_mod1$accuracy[[2]]
# ECM_tst
results_mod1$accuracy[[3]]
# EAM_tst
results_mod1$accuracy[[4]]


## 3.4.2 ## Resultados del Modelo 2
results_mod2 <- arimafunct(mod2)  # Llama a la función anterior

# Representación gráfica de la predicción sobre la muestra de entrenamiento
par(mfrow = c(2, 2))
ts.plot(log(ser_ts), results_mod2$grafseries$for_det, col = 1:2, main = "Ajuste de la parte determinista (serie completa)")
ts.plot(train_e, results_mod2$grafseries$fitmod, col = 1:2, main = "Ajuste de la parte estacionaria (entrenamiento)")
ts.plot(log(train), results_mod2$grafseries$ffitmod, col = 1:2, main = "Ajuste agregado (entrenamiento)") # gráfico de los valores predichos sobre la serie original
ts.plot(results_mod2$grafseries$ffitmod_e, ylab = "", main = "Errores del modelo (entrenamiento)")

# Plot of actual and forecasted values (parte estacionaria)
par(mfrow = c(1, 2))
plot(results_mod2$grafseries$observed_st, type = "l", xlab = "", ylab = "",
     main = "Serie y predicción (parte estacionaria)", ylim = c(min(results_mod2$grafseries$predmod$pred-2*results_mod2$grafseries$predmod$se), max(results_mod2$grafseries$predmod$pred+2*results_mod2$grafseries$predmod$se))) 
lines(results_mod2$grafseries$predict_st, col = "blue", lty = 2)
lines(results_mod2$grafseries$lower, col = "red", lty = 5)
lines(results_mod2$grafseries$upper, col = "red",lty = 5)
ts.plot(results_mod2$grafseries$observed, results_mod2$grafseries$predict, 
        col = 1:2, main = "Ajuste sobre la serie test")

# ECM_tr
results_mod2$accuracy[[1]]
# EAM_tr
results_mod2$accuracy[[2]]
# ECM_tst
results_mod2$accuracy[[3]]
# EAM_tst
results_mod2$accuracy[[4]]


## 3.4.3 ## Resultados del Modelo 3
results_mod3 <- arimafunct(mod3) # Llama a la función anterior

# Representación gráfica de la predicción sobre la muestra de entrenamiento
par(mfrow = c(2, 2))
ts.plot(log(ser_ts), results_mod3$grafseries$for_det, col = 1:2, main = "Ajuste de la parte determinista (serie completa)")
ts.plot(train_e, results_mod3$grafseries$fitmod, col = 1:2, main = "Ajuste de la parte estacionaria (entrenamiento)")
ts.plot(log(train), results_mod3$grafseries$ffitmod, col = 1:2, main = "Ajuste agregado (entrenamiento)") # gráfico de los valores predichos sobre la serie original
ts.plot(results_mod3$grafseries$ffitmod_e, ylab = "", main = "Errores del modelo (entrenamiento)")

# Plot of actual and forecasted values (parte estacionaria)
par(mfrow = c(1, 2))
plot(results_mod3$grafseries$observed_st, type = "l", xlab = "", ylab = "",
     main = "Serie y predicción (parte estacionaria)", ylim = c(min(results_mod3$grafseries$predmod$pred-2*results_mod3$grafseries$predmod$se), max(results_mod3$grafseries$predmod$pred+2*results_mod3$grafseries$predmod$se))) 
lines(results_mod3$grafseries$predict_st, col = "blue", lty = 2)
lines(results_mod3$grafseries$lower, col = "red", lty = 5)
lines(results_mod3$grafseries$upper, col = "red",lty = 5)
ts.plot(results_mod3$grafseries$observed, results_mod3$grafseries$predict, 
        col = 1:2, main = "Ajuste sobre la serie test")

# ECM_tr
results_mod3$accuracy[[1]]
# EAM_tr
results_mod3$accuracy[[2]]
# ECM_tst
results_mod3$accuracy[[3]]
# EAM_tst
results_mod3$accuracy[[4]]


## 4 ## ANÁLISIS DE VOLATILIDAD ------------------------------------------------

# Se toman diferencias en el séptimo retardo
dtrain_e = diff(train_e, 7) # train_e es la serie sin la parte determinista
# Se representa la volatilidad:
par(mfrow = c(1, 3))
acf(dtrain_e^2, lag.max = 28)
pacf(dtrain_e^2, lag.max = 28) 
McLeod.Li.test(y = dtrain_e) # p-valores por debajo de 0.05.
# El contraste de McLeod-Li:
# Bajo H0 n^2 es un ruido blanco. Aquí se rechaza que sea un ruido blanco.

# Ahora se modelizan las diferencias al cuadrado como un ARMA(2,1)
gmod1 = arima(dtrain_e^2, order = c(2, 0, 1))
gmod11 <- garch(dtrain_e, order = c(2, 1))
# Representación gráfica
par(mfrow = c(1, 3))
acf(gmod1$residuals, lag.max = 36)
pacf(gmod1$residuals, lag.max = 36)
McLeod.Li.test(y=gmod11$residuals)
# H0: No hay efectos ARCH antes del Lag k (7)


## 5 ## MODELOS ESTRUCTURALES DE SERIES TEMPORALES -----------------------------

# Se rescata la parte determinista calculada antes de los modelos ARIMA
for_det <- results_mod1$grafseries$for_det

structfunct <- function(seatype){
  # Función para la estimación del error en el ajuste de modelos esructurales
  # y la extracción de series para su representación
  # INPUTS:
  # seatype es la especificacición: "dummy" or "trigonometric" debe especificarse
  # entre comillas.
  # OUTPUTS:
  # Objeto tipo lista con series para representación y medidas de error
  mod <- SSModel(train_e ~ SSMtrend(degree = 2, Q = list(matrix(NA), matrix(NA))) +
                    SSMseasonal(period = 7, sea.type = seatype, Q = NA), H = NA)
  fitmod <- fitSSM(mod, inits = c(0, 0, 0, 0, 0, 0, 0, 0, 0), method = "BFGS")
  outmod <- KFS(fitmod$model)
  
  h = fitmod$model$H
  if (seatype == "dummy")
  {QQ = matrix(fitmod$model$Q,3,3)}
  else
  {QQ = matrix(fitmod$model$Q,8,8)}
  
  # Residuos:
  resmod <- ts(residuals(outmod), start = c(1997, 1), frequency = 365)
  
  # Prediccion con BSM (de la parte estacionaria)
  modpred <- SSModel(train_e ~ SSMtrend(degree = 2, Q = list(matrix(QQ[1,1]), matrix(QQ[2,2]))) + 
                        SSMseasonal(period = 7, sea.type = seatype, Q = QQ[3,3]), H = h)
  predmod <- predict(modpred, n.ahead = 365, se.fit = T)
  
  # Series para la representación de las predicciones con intervalos de confianza
  predict_st <- ts(c(rep(NA, length(train_e)), predmod[,1]), start = c(1997, 1), frequency = 365)
  upper <- ts(c(rep(NA, length(train_e)), predmod[,1]+2*predmod[,2]), start = c(1997, 1), frequency = 365) 
  lower <- ts(c(rep(NA, length(train_e)), predmod[,1]-2*predmod[,2]), start = c(1997, 1), frequency = 365) 
  observed_st <- ts(c(train_e[3896:4261], rep(NA, length(predmod[,1]))), frequency = 365, start = c(2007, 244), end = c(2009, 244))
  
  # Prediccion on Fourrier (parte determinista y agregación con estacionaria)
  # agrega la parte estacionaria (descontados errores) observada y la predicción 365 obs en adelante
  train_eplus <- rbind(train_e - resmod, matrix(predmod[,1], ncol = 1, byrow = F))
  predtot <- for_det + train_eplus
  
  predtot_tr <- predtot[1:4261] # predicciones de la parte train
  predtot_tst <- predtot[4262:4626] # predicciones exclusivas de la parte test
  
  # Series para la representación de las series predichas y bservada en log
  predict <- ts(c(rep(NA, 366), predtot_tst), start = c(2007, 244), frequency = 365)
  observed <- ts(log(test), start = c(2007, 244), frequency = 365)
  
  # Error train
  err_tr <- log(train) - predtot_tr
  ECM_tr <- (t(err_tr) %*% err_tr)/length(err_tr)
  EAM_tr <- sum(abs(err_tr))/length(err_tr)
  # Error test
  err_tst <- log(test[367:length(test)]) - predtot_tst
  ECM_tst <- (t(err_tst) %*% err_tst)/length(predmod[,1])
  EAM_tst <- sum(abs(err_tst))/length(predmod[,1])
  
  statistics <- c(ECM_tr = ECM_tr, EAM_tr = EAM_tr, ECM_tst = ECM_tst, EAM_tst = EAM_tst)
  grafseries <- list(predict_st = predict_st, upper = upper, lower = lower, 
                     observed_st = observed_st, predict = predict, 
                     observed = observed, predmod = predmod, resmod = resmod, 
                     predtot = predtot)
  results <- list(accuracy = statistics, grafseries = grafseries)
}


## 5.1 ## Resultados del Modelo 4
# Por variables instrumentales
results_mod4 <- structfunct("dummy")

# Errores y periodograma
par(mfrow = c(2, 2))
acf(results_mod4$grafseries$resmod, lag.max = 28, main = "Residuos")
pacf(results_mod4$grafseries$resmod, lag.max = 28, main = "Residuos")
spectrum(results_mod4$grafseries$resmod, main = "Periodograma")
spectrum(results_mod4$grafseries$resmod, kernel("fejer", 100, r = 7), 
         main = "Periodograma suavizado")

# Representación gráfica de la predicción sobre la muestra de entrenamiento
par(mfrow = c(2, 2))
ts.plot(log(ser_ts), for_det, col = 1:2, 
        main = "Ajuste de la parte determinista (serie completa)")
ts.plot(log(ser_ts), results_mod4$grafseries$predtot, col = 1:2, 
        main = "Ajuste agregado (serie completa)") # gráfico de los valores
# predichos sobre la serie original
ts.plot(results_mod4$grafseries$resmod, ylab = "", 
        main = "Errores del modelo (entrenamiento)")

# Plot of actual and forecasted values
par(mfrow = c(1, 2))
plot(results_mod4$grafseries$observed_st, type = "l", xlab = "", ylab = "", 
     main = "Serie y predicción (parte estacionaria)", ylim = c(min(results_mod4$grafseries$predmod[,1]-2*results_mod4$grafseries$predmod[,2]), max(results_mod4$grafseries$predmod[,1]+2*results_mod4$grafseries$predmod[,2]))) 
lines(results_mod4$grafseries$predict_st, col = "blue", lty = 2)
lines(results_mod4$grafseries$lower, col = "red", lty = 5)
lines(results_mod4$grafseries$upper, col = "red",lty = 5)
ts.plot(results_mod4$grafseries$observed, results_mod4$grafseries$predict, 
        col = 1:2, main = "Ajuste sobre la serie test")

# ECM_tr
results_mod4$accuracy[[1]]
# EAM_tr
results_mod4$accuracy[[2]]
# ECM_tst
results_mod4$accuracy[[3]]
# EAM_tst
results_mod4$accuracy[[4]]


## 5.2 ## Resultados del Modelo 5
# Especificación trigonométrica
results_mod5 <- structfunct("trigonometric")

# Errores y periodograma
par(mfrow = c(2, 2))
acf(results_mod5$grafseries$resmod, lag.max = 28, main = "Residuos")
pacf(results_mod5$grafseries$resmod, lag.max = 28, main = "Residuos")
spectrum(results_mod5$grafseries$resmod, main = "Periodograma")
spectrum(results_mod5$grafseries$resmod, kernel("fejer", 100, r = 7), 
         main = "Periodograma suavizado")

# Representación gráfica de la predicción sobre la muestra de entrenamiento
par(mfrow = c(2, 2))
ts.plot(log(ser_ts), for_det, col = 1:2, 
        main = "Ajuste de la parte determinista (serie completa)")
ts.plot(log(ser_ts), results_mod5$grafseries$predtot, col = 1:2, 
        main = "Ajuste agregado (serie completa)") # gráfico de los valores 
# predichos sobre la serie original
ts.plot(results_mod5$grafseries$resmod, ylab = "", 
        main = "Errores del modelo (entrenamiento)")

# Plot of actual and forecasted values (parte estacionaria)
par(mfrow = c(1, 2))
plot(results_mod5$grafseries$observed_st, type = "l", xlab = "", ylab = "", 
     main = "Serie y predicción (parte estacionaria)", ylim = c(min(results_mod5$grafseries$predmod[,1]-2*results_mod5$grafseries$predmod[,2]), max(results_mod5$grafseries$predmod[,1]+2*results_mod5$grafseries$predmod[,2]))) 
lines(results_mod5$grafseries$predict_st, col = "blue", lty = 2)
lines(results_mod5$grafseries$lower, col = "red", lty = 5)
lines(results_mod5$grafseries$upper, col = "red", lty = 5)
ts.plot(results_mod5$grafseries$observed, results_mod5$grafseries$predict, 
        col = 1:2, main = "Ajuste sobre la serie test")

# ECM_tr
results_mod5$accuracy[[1]]
# EAM_tr
results_mod5$accuracy[[2]]
# ECM_tst
results_mod5$accuracy[[3]]
# EAM_tst
results_mod5$accuracy[[4]]


## 6 ## ALISADO EXPONENCIAL ---------------------------------------------------
# Alisado exponencial por State Space

# Es absolutamente necesario bajar la frecuencia:
train_e <- ts(train_e, start = c(1997, 1), frequency = 7)

etsfunct <- function(model, opt){
  # Función para la estimación del error en el ajuste de modelos de alisado
  # exponencial con State Space y 
  # INPUTS:
  # model: un string de longitud tres en donde 
  # - la 1a letra denota el tipo de error;
  # - la 2a denota el tipo de tendencia
  # - la 3a denota el tipo de estación
  # Los valores posibles de las 3 letras son: 
  # "N" = none, "A"= additive, "M" = multiplicative and "Z" = automatically selected
  # opt: criterio de optimización
  # posibles valores: 
  # "lik" (log-Verosimilitud) "mse" (Error cuadrático medio) 
  # "sigma" (desviación típica de los residuos) "mae" (Error absoluto medio)
  # OUTPUTS:
  # Objeto tipo lista con series para la representación y medidas de error.
  mod <- ets(train_e, model = model, opt.crit = opt)
  residuals <- mod$residuals
  
  modpred <- forecast(mod, h = 365)
  
  train_eplus <- rbind(modpred$fitted, as.matrix(modpred$mean, ncol = 1, byrow = F))
  for_tot <- as.matrix(for_det, ncol = 1, byrow = F) + train_eplus
  
  predict_st <- ts(c(rep(NA, length(train_e)), modpred$mean), start = c(1997, 1), frequency = 365)
  upper <- ts(c(rep(NA, length(train_e)), modpred$upper[,"95%"]), start = c(1997, 1), frequency = 365) 
  lower <- ts(c(rep(NA, length(train_e)), modpred$lower[,"95%"]), start = c(1997, 1), frequency = 365) 
  observed_st <- ts(c(train_e[3896:4261], rep(NA, length(modpred$mean))), frequency = 365, start = c(2007, 244), end = c(2009, 244))
  
  predict <- ts(c(rep(NA, 366), for_tot[4262:4626]), start = c(2007, 244), frequency = 365)
  observed <- ts(log(test), start = c(2007, 244), frequency = 365)
  
  # Error train
  ECM_tr = mod$mse
  EAM_tr = sum(abs(modpred$residuals))/length(modpred$residuals)
  
  # Error test
  err_tst <- log(test[367:length(test)]) - for_tot[4262:4626]
  ECM_tst <- (t(err_tst) %*% err_tst)/length(err_tst)
  EAM_tst <- sum(abs(err_tst))/length(err_tst)
  
  statistics <- c(ECM_tr = ECM_tr, EAM_tr = EAM_tr, ECM_tst = ECM_tst, EAM_tst = EAM_tst)
  grafseries <- list(predict_st = predict_st, upper = upper, lower = lower, 
                     observed_st = observed_st, predict = predict, for_tot = for_tot,
                     observed = observed, modpred = modpred, residuals = residuals)
  results <- list(accuracy = statistics, grafseries = grafseries)
}


## 6.1 ## Resultados del Modelo 6
# Con coeficientes aditivos
results_mod6 <- etsfunct("AAA", "mse")

# Errores y periodograma
par(mfrow = c(2, 2))
acf(results_mod6$grafseries$residuals, lag.max = 28, main = "Residuos")
pacf(results_mod6$grafseries$residuals, lag.max = 28, main = "Residuos")
spectrum(results_mod6$grafseries$residuals, main = "Periodograma")
spectrum(results_mod6$grafseries$residuals, kernel("fejer", 100, r = 7), 
         main = "Periodograma suavizado")

# Representación gráfica de la predicción sobre la muestra de entrenamiento
par(mfrow = c(2, 2))
ts.plot(log(ser_ts), for_det, col = 1:2, 
        main = "Ajuste de la parte determinista (serie completa)")
ts.plot(train_e, results_mod6$grafseries$modpred$fitted, col = 1:2, 
        main = "Ajuste de la parte estacionaria (entrenamiento)")
ts.plot(log(ser_ts), results_mod6$grafseries$for_tot, col=1:2, 
        main = "Ajuste agregado (serie completa)")

# Plot of actual and forecasted values (parte estacionaria)
par(mfrow = c(1, 2))
plot(results_mod6$grafseries$observed_st, type = "l", xlab = "", ylab = "", 
     main = "Serie y predicción (parte estacionaria)", 
     ylim = c(min(results_mod6$grafseries$modpred$lower[,"95%"]), max(results_mod6$grafseries$modpred$upper[,"95%"]))) 
lines(results_mod6$grafseries$predict_st, col = "blue", lty = 2)
lines(results_mod6$grafseries$lower, col = "red", lty = 5)
lines(results_mod6$grafseries$upper, col = "red", lty = 5)
ts.plot(results_mod6$grafseries$observed, results_mod6$grafseries$predict, 
        col = 1:2, main = "Ajuste sobre la serie test")

# ECM_tr
results_mod6$accuracy[[1]]
# EAM_tr
results_mod6$accuracy[[2]]
# ECM_tst
results_mod6$accuracy[[3]]
# EAM_tst
results_mod6$accuracy[[4]]


## 6.2 ## Resultados del Modelo 7
# Con coeficientes automáticos
results_mod7 <- etsfunct("ZZZ", "mse")

# Errores y periodograma
par(mfrow = c(2, 2))
acf(results_mod7$grafseries$residuals, lag.max = 28, main = "Residuos")
pacf(results_mod7$grafseries$residuals, lag.max = 28, main = "Residuos")
spectrum(results_mod7$grafseries$residuals, main = "Periodograma")
spectrum(results_mod7$grafseries$residuals, kernel("fejer", 100, r = 7), main = "Periodograma suavizado")

# Representación gráfica de la predicción sobre la muestra de entrenamiento
par(mfrow = c(2, 2))
ts.plot(log(ser_ts), for_det, col = 1:2, 
        main = "Ajuste de la parte determinista (serie completa)")
ts.plot(train_e, results_mod7$grafseries$modpred$fitted, col = 1:2, 
        main = "Ajuste de la parte estacionaria (entrenamiento)")
ts.plot(log(ser_ts), results_mod7$grafseries$for_tot, col=1:2, 
        main = "Ajuste agregado (serie completa)")

# Plot of actual and forecasted values (parte estacionaria)
par(mfrow = c(1, 2))
plot(results_mod7$grafseries$observed_st, type = "l", xlab = "", ylab = "", 
     main = "Serie y predicción (parte estacionaria)", 
     ylim = c(min(results_mod7$grafseries$modpred$lower[,"95%"]), max(results_mod7$grafseries$modpred$upper[,"95%"]))) 
lines(results_mod7$grafseries$predict_st, col = "blue", lty = 2)
lines(results_mod7$grafseries$lower, col = "red", lty = 5)
lines(results_mod7$grafseries$upper, col = "red", lty = 5)
ts.plot(results_mod7$grafseries$observed, results_mod7$grafseries$predict, 
        col = 1:2, main = "Ajuste sobre la serie test")

# ECM_tr
results_mod7$accuracy[[1]]
# EAM_tr
results_mod7$accuracy[[2]]
# ECM_tst
results_mod7$accuracy[[3]]
# EAM_tst
results_mod7$accuracy[[4]]


## 6.3 ## Resultados del Modelo 8
# Doble Alisado exponencial
# No sirve la función definida en el punto 6 !!!

# Se convierte la serie de entrenamiento en una serie con dos frecuencias.
train2 <- msts(train, seasonal.periods = c(7, 365.25))
mod8 <- bats(log(train2), use.parallel = T, num.cores = 8) 
# el parámetro use.parallel = T permite la ejecución en paralelo para ahorrar 
# tiempo. En num.cores se especifica el número de cores del PC que se desean usar.

predicted_v <- forecast(mod8, h = 365)
# contiene los valores predichos en $mean y los extremos de los 
# intervalos de confianza

# Series de los valores predichos y observados en log
predict <- ts(c(rep(NA, 366), predicted_v$mean), start = c(2007, 244), frequency = 365)
observed <- ts(log(test), start = c(2007, 244), frequency = 365)

# Errores y periodograma
par(mfrow = c(2, 2))
acf(residuals(mod8), lag.max = 28, main = "Residuos")
pacf(residuals(mod8), lag.max = 28, main = "Residuos")
spectrum(residuals(mod8), main = "Periodograma")
spectrum(residuals(mod8), kernel("fejer", 100, r = 7), main = "Periodograma suavizado")

# Representación gráfica de la predicción sobre la muestra de entrenamiento
par(mfrow = c(2, 1))
ts.plot(log(train), fitted.values(mod8), col=1:2, 
        main = "Ajuste agregado (entrenamiento)")
plot(residuals(mod8), main = "Errores (entrenamiento)")

# Plot of actual and forecasted values (parte estacionaria)
par(mfrow = c(1, 1))
plot(predicted_v, xlim = c(2007.7, 2009.7), main = "Serie y predicción")
ts.plot(observed, predict, col = 1:2, main = "Ajuste sobre la serie test")

# Error train
ECM_tr_mod8 <- (t(residuals(mod8)) %*% residuals(mod8))/length(residuals(mod8))
EAM_tr_mod8 <- sum(abs(residuals(mod8)))/length(residuals(mod8))
ECM_tr_mod8
EAM_tr_mod8

# Error test
err_tst_mod8 <- log(test[367:length(test)]) - predicted_v$mean
ECM_tst_mod8 <- (t(err_tst_mod8) %*% err_tst_mod8)/length(err_tst_mod8)
EAM_tst_mod8 <- sum(abs(err_tst_mod8))/length(err_tst_mod8)
ECM_tst_mod8
EAM_tst_mod8