#WD

# Data setting ------------------------------------------------------------
 
a.min <- 0 # età minima dell'intervallo di fit
a.max <- 100 # età massima dell'intervallo di fit
omega <-110 # età estrema delle tavole
A.fit <-c(a.min:a.max) 
y.fit.min <- 1950 # anno minimo dell'intervallo di fit
y.fit.max <- 2021 # anno massimo dell'intervallo di fit
Y.fit <- c(y.fit.min:y.fit.max)
y.pred <- 67 # orizzonte d proiezione
Y.pred <- c((y.fit.max+1):(y.fit.max+y.pred))
n.sim <- 100 # numero di simulazioni
n.boot <- 20 # numero di simulazioni per il bootstrap
int_vect <- c(rep(0.02,omega-a.min)) # vettore di tassi d'interesse per scadenza
v_vect <- cumprod((1+int_vect)^-1) # vettore di fattori di sconto per scadenza
#################################### fine data setting


# Importazione dati demografici HMD ---------------------------------------

library(demography)

SPAdata <- read.demogdata(file = "Mx_1x1_spagna.txt",
                          popfile = "Exposures_1x1_spagna.txt",
                          type = "mortality",
                          label = "Spain")
names(SPAdata)


# Il fit dei modelli ------------------------------------------------------

library("StMoMo")
source("Functions.r") # caricamento di funzioni da un altro script

# Trasfomazione dei dati e definizione dell'intervallo di stima
SPAmStMoMoC <- StMoMoData(SPAdata, series = "male")     # esposizioni centrali
SPAfStMoMoC <- StMoMoData(SPAdata, series = "female")   
SPAmStMoMoI <- central2initial(SPAmStMoMoC)             # esposizioni iniziali
SPAfStMoMoI <- central2initial(SPAfStMoMoC)
wxt <- genWeightMat(ages = A.fit, years = Y.fit, clip = 4) # pesi
SPAmRates <- SPAmStMoMoC$Dxt/SPAmStMoMoC$Ext
SPAmRates <- SPAmRates[A.fit+1,tail(SPAmStMoMoC$years+1,
                                    length(Y.fit))-SPAmStMoMoC$years[1]]
SPAfRates <- SPAfStMoMoC$Dxt/SPAfStMoMoC$Ext
SPAfRates <- SPAfRates[A.fit+1,tail(SPAfStMoMoC$years+1,
                                    length(Y.fit))-SPAfStMoMoC$years[1]]
 
################### Models fitting
LCfit_SPAm <- fit(lc(link = "logit"), #se seguo l'approccio logit, non lavoro sulla forza di motalità ma sugli esposti per cui devo prendere gli esposti
                  data = SPAmStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt)
LCfit_SPAf <- fit(lc(link = "logit"),
                  data = SPAfStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt)
RHfit_SPAm <- fit(rh(link = "logit", cohortAgeFun="1"), ## anche in questo forzo con la logit, B0 è 1 in quanto uso rh semplificato 
                  data = SPAmStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt,
                  start.ax = LCfit_SPAm$ax,
                  start.bx = LCfit_SPAm$bx,
                  start.kt = LCfit_SPAm$kt)  #la lentezza nel run è dato da un problema di convergenza per semplificare possiamo farlo partire dai valori di partenza del lee carter con le 3 stringhe finali  
RHfit_SPAf <- fit(rh(link = "logit", cohortAgeFun="1"),
                  data = SPAfStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt,
                  start.ax = LCfit_SPAf$ax,
                  start.bx = LCfit_SPAf$bx,
                  start.kt = LCfit_SPAf$kt)
APCfit_SPAm <- fit(apc(link = "logit"),
                   data = SPAmStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max,
                   wxt = wxt)
APCfit_SPAf <- fit(apc(link = "logit"),
                   data = SPAfStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max,
                   wxt = wxt)
CBDfit_SPAm <- fit(cbd(),
                   data = SPAmStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max,
                   wxt = wxt)
CBDfit_SPAf <- fit(cbd(),
                   data = SPAfStMoMoI,
                   ages.fit = a.min:a.max,
                   years.fit=y.fit.min:y.fit.max,
                   wxt = wxt)
M7fit_SPAm <- fit(m7(link = "logit"),
                  data =SPAmStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt)
M7fit_SPAf <- fit(m7(link = "logit"),
                  data = SPAfStMoMoI,
                  ages.fit = a.min:a.max,
                  years.fit=y.fit.min:y.fit.max,
                  wxt = wxt)
PLATfit_SPAm <- fit(PLAT,
                    data = SPAmStMoMoI,
                    ages.fit = a.min:a.max,
                    years.fit=y.fit.min:y.fit.max,
                    wxt = wxt)
PLATfit_SPAf <- fit(PLAT,
                    data = SPAfStMoMoI,
                    ages.fit = a.min:a.max,
                    years.fit=y.fit.min:y.fit.max,
                    wxt = wxt)

##il lavoro di stima è completato, ora passiamo ala ## Rappresentazione dei parametri

dev.new(width = 550, height = 330, unit = "px")
plot(LCfit_SPAm, nCol = 2)
plot(LCfit_SPAf, nCol = 3)
#dev.new(width = 550, height = 330, unit = "px")
plot(RHfit_SPAm, nCol = 3) #l'effetto coorte non devrebbe essere lineare. in questo caso, la risposta del rh è che non avendo un effetto coorte abbiamo bisogno di due K 
plot(RHfit_SPAf, nCol = 3)
plot(APCfit_SPAm, parametricbx = F)
plot(APCfit_SPAf, parametricbx = F)
plot(CBDfit_SPAm, parametricbx = F)
plot(CBDfit_SPAf, parametricbx = F)
plot(M7fit_SPAm, parametricbx = F)
plot(M7fit_SPAf, parametricbx = F)
plot(PLATfit_SPAm, parametricbx = F)
plot(PLATfit_SPAf, parametricbx = F)

#### le linee seguenti si riferiscono al modello Lee-Carter
#### istruzioni simili possono essere scritte per gli altri modelli -> parto dal fittato e isolo i vettori ad esempio alpha, beta ecc per 
par(mfrow=c(2, 3))
plot(LCfit_SPAm$ax, main=c(expression(paste("Male population - ",alpha[x]))))
plot(LCfit_SPAm$bx, main=c(expression(paste("Male population - ",beta[x]))))
plot(t(LCfit_SPAm$kt), main=c(expression(paste("Male population - ",kappa[t]))))
plot(LCfit_SPAf$ax, main=c(expression(paste("Female population - ",alpha[x]))))
plot(LCfit_SPAf$bx, main=c(expression(paste("Female population - ",beta[x]))))
plot(t(LCfit_SPAf$kt), main=c(expression(paste("Female population - ",kappa[t]))))

View(LCfor_SPAm$rates)

# Analisi e confronto tra modelli -----------------------------------------

######### Confronto tra modelli con indicatore di bontà
### Log-like più è alta più funziona bene 
#Maschi
LCfit_SPAm$loglik 
CBDfit_SPAm$loglik
RHfit_SPAm$loglik
APCfit_SPAm$loglik
M7fit_SPAm$loglik
PLATfit_SPAm$loglik
#Femmine
LCfit_SPAf$loglik
CBDfit_SPAf$loglik
RHfit_SPAf$loglik
APCfit_SPAf$loglik
M7fit_SPAf$loglik
PLATfit_SPAf$loglik

###### BIC 
#Maschi
BIC(LCfit_SPAm)
BIC(CBDfit_SPAm)
BIC(RHfit_SPAm)
BIC(APCfit_SPAm)
BIC(M7fit_SPAm)
BIC(PLATfit_SPAm)
#Femmine
BIC(LCfit_SPAf)
BIC(CBDfit_SPAf)
BIC(RHfit_SPAf)
BIC(APCfit_SPAf)
BIC(M7fit_SPAf)
BIC(PLATfit_SPAf)

###### Residui
LCres_SPAm <- residuals(LCfit_SPAm)
LCres_SPAf <- residuals(LCfit_SPAf)
CBDres_SPAm <- residuals(CBDfit_SPAm)
CBDres_SPAf <- residuals(CBDfit_SPAf)
RHres_SPAm <- residuals(RHfit_SPAm)
RHres_SPAf <- residuals(RHfit_SPAf)
APCres_SPAm <- residuals(APCfit_SPAm)
APCres_SPAf <- residuals(APCfit_SPAf)
M7res_SPAm <- residuals(M7fit_SPAm)
M7res_SPAf <- residuals(M7fit_SPAf)
PLATres_SPAm <- residuals(PLATfit_SPAm)
PLATres_SPAf <- residuals(PLATfit_SPAf)

### rappresentazioni alternative dei residui
TYPE="scatter"#"scatter","signplot"
par(mfrow=c(1, 1))
plot(LCres_SPAm, type=TYPE)
plot(LCres_SPAf, type=TYPE)
plot(CBDres_SPAm, type=TYPE)
plot(CBDres_SPAf, type=TYPE)
plot(RHres_SPAm, type=TYPE)
plot(RHres_SPAf, type=TYPE)
plot(APCres_SPAm, type=TYPE)
plot(APCres_SPAf, type=TYPE)
plot(M7res_SPAm, type=TYPE)
plot(M7res_SPAf, type=TYPE)
plot(PLATres_SPAm, type=TYPE)
plot(PLATres_SPAf, type=TYPE)

dev.new()
plot(LCres_SPAm, type = "colourmap", reslim = c(-3.5, 3.5))
dev.new(width = 550, height = 330, unit = "px")
plot(RHres_SPAm, type = "colourmap", reslim = c(-3.5, 3.5))

# Proiezioni e simulazioni con il modello Lee-Carter ----------------------

LCfor_SPAm <- forecast(LCfit_SPAm, h=y.pred)
### istruzioni alternative
#LCfor_ITAmArima <- forecast(LCfit_ITAm, h=y.pred, kt.method = "iarima",  kt.order = c(1, 1, 2))  # ARIMA order specified
#LCfor_ITAmBArima <- forecast(LCfit_ITAm, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # auto.arima
LCfor_SPAf <- forecast(LCfit_SPAf, h=y.pred)

#####  rappresentazione grafica - evoluzione dell'indice temporale
plot(LCfor_SPAm, only.kt = TRUE)
plot(LCfor_SPAf, only.kt = TRUE)

#####  rappresentazione grafica - confronto tra osservazioni e proiezioni
dev.new()
par(mfrow=c(1, 2), oma=c(2,0,0,0))
plot(c(a.min:a.max),log(LCfor_SPAm$rates[,y.pred]), bty="l", xlab="Ages", ylab="Male death rates (log scale)",
     ylim=c(min(log(LCfor_SPAm$rates[,y.pred])),max(log((LCfit_SPAm$Dxt/LCfit_SPAm$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((LCfit_SPAm$Dxt/LCfit_SPAm$Ext)[,length(Y.fit)]), type="p")
plot(c(a.min:a.max),log(LCfor_SPAf$rates[,y.pred]), bty="l", xlab="Ages", ylab="Female death rates (log scale)",
     ylim=c(min(log(LCfor_SPAf$rates[,y.pred])),max(log((LCfit_SPAf$Dxt/LCfit_SPAf$Ext)[,length(Y.fit)]))), type="l")
lines(c(a.min:a.max),log((LCfit_SPAf$Dxt/LCfit_SPAf$Ext)[,length(Y.fit)]), type="p")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("last observed", "last projected"), pch=c(1,-1), lty=c(0,1), lwd=c(1,1), col=c(1:1), cex=0.9)

#####  generazione di probabilità dai tassi (osservati e proiettati)
rates_LCfit_SPAm <- fitted(LCfit_SPAm, type = "rates")
rates_LCfor_SPAm <- LCfor_SPAm$rates
rates_LC_SPAm <- cbind(SPAmRates,rates_LCfor_SPAm)
q_LC_SPAm <- 1- exp(-rates_LC_SPAm)

# estensione delle probabilità qxt ad età estreme
q_LC_SPAm.ext  <- extrapolation.fit(q_LC_SPAm)
rates_LCfit_SPAf <- fitted(LCfit_SPAf, type = "rates")
rates_LCfor_SPAf <- LCfor_SPAf$rates
rates_LC_SPAf <- cbind(SPAfRates,rates_LCfor_SPAf)
rates_LC_ITAf <- cbind(SPAfRates,rates_LCfor_SPAf)
q_LC_SPAf <- 1- exp(-rates_LC_SPAf)
#### estensione delle probabilità qxt ad età estreme
q_LC_SPAf.ext  <- extrapolation.fit(q_LC_SPAf)

q_LC.SPAm<-as.data.frame(q_LC_SPAm.ext)
q_LC.SPAf<-as.data.frame(q_LC_SPAf.ext)



### rappresentazione grafica delle probabilità proiettate
dev.new(width = 550, height = 330, unit = "px")
plot(extractCohort(q_LC_SPAm.ext, cohort = 1950), type = "p", log = "y", xlab = "age", ylab = "q(x)",
     main = "Mortality rates (log scale) for the cohort 1950")
lines(extractCohort(q_LC_SPAm, cohort = 1950))
dev.new(width = 550, height = 330, unit = "px")
plot(extractCohort(log(q_LC_SPAm.ext/(1-q_LC_SPAm.ext)), cohort = 1950), type = "p", xlab = "age", ylab = "logit(q(x))",
     main = "Logit transform of mortality rates for the cohort 1950")
lines(extractCohort(log(q_LC_SPAm/(1-q_LC_SPAm)), cohort = 1950))


# Probabilità di sopravvivenza e aspettative di vita proiettate 
p_LC_SPAm.ext <- 1-q_LC_SPAm.ext                    # Probabilità di sopravvivenza a 1 anno
p0n_LC_SPAm.ext <- apply(p_LC_SPAm.ext, 2, cumprod) # Probabilità di sopravvivenza a n anni
ex_LC_SPAm.ext <- life.exp(q_LC_SPAm.ext)           # Aspettativa di vita

p_LC_SPAf.ext <- 1-q_LC_SPAf.ext # 1 year survival probability
p0n_LC_SPAf.ext <- apply(p_LC_SPAf.ext, 2, cumprod) # n year survival probability
ex_LC_SPAf.ext <- life.exp(q_LC_SPAf.ext)

p_LC.SPAm<-as.data.frame(p_LC_SPAm.ext)
p0n_LC.SPAm<-as.data.frame(p0n_LC_SPAm.ext)
ex_LC.SPAm<-as.data.frame(ex_LC_SPAm.ext)
p_LC.SPAf<-as.data.frame(p_LC_SPAf.ext)
p0n_LC.SPAf<-as.data.frame(p0n_LC_SPAf.ext)
ex_LC.SPAf<-as.data.frame(ex_LC_SPAf.ext)


# Simulazioni (RANDOM WALK WITH DRIFT) ------------------------------------
library(fanplot)
PLATsim_SPAm.mrwd <- simulate(PLATfit_SPAm, nsim = n.sim, h=y.pred)
### Istruzioni alternative
#PLATsim_ITAmArima <- forecast(PLATfit_ITAm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = c(1, 1, 2))  # ARIMA order specified
#PLATsim_ITAmBArima <- forecast(PLATfit_ITAm, nsim = n.sim, h=y.pred, kt.method = "iarima",  kt.order = NULL)  # auto.arima
rates_PLAT_SPAm.st <- PLATsim_SPAm.mrwd$rates
q_PLAT_SPAm.st <- 1- exp(-rates_PLAT_SPAm.st)
q_PLAT_SPAm.st.ext <-  extrapolation.sim(q_PLAT_SPAm.st)

PLATsim_SPAf.mrwd <- simulate(PLATfit_SPAf, nsim = n.sim, h=y.pred)
rates_PLAT_SPAf.st <- PLATsim_SPAf.mrwd$rates
q_PLAT_SPAf.st <- 1- exp(-rates_PLAT_SPAf.st)
q_PLAT_SPAf.st.ext <-  extrapolation.sim(q_PLAT_SPAf.st)

### rappresentazione grafica
par(mfrow=c(1, 2))
plot(PLATfit_SPAm$years, PLATfit_SPAm$kt[1, ], xlim = range(PLATfit_SPAm$years, PLATsim_SPAm.mrwd$kt.s$years),
     ylim = range(PLATfit_SPAm$kt, PLATsim_SPAm.mrwd$kt.s$sim), type = "l", xlab = "year", ylab = "kt",
     main = "Lee-Carter: Period index (mrwd)")
matlines(PLATsim_SPAm.mrwd$kt.s$years, PLATsim_SPAm.mrwd$kt.s$sim[1, , ], type = "l", lty = 1)
plot(PLATfit_SPAm$years, (PLATfit_SPAm$Dxt / PLATfit_SPAm$Ext)["65", ], xlim = range(PLATfit_SPAm$years, PLATsim_SPAm.mrwd$years),
     ylim = range((PLATfit_SPAm$Dxt / PLATfit_SPAm$Ext)["65", ], PLATsim_SPAm.mrwd$rates["65", , ]), type = "l", xlab = "year", ylab = "rate",
     main = "Lee-Carter: Simulated mortality rates at age 65")
matlines(PLATsim_SPAm.mrwd$years, PLATsim_SPAm.mrwd$rates["65", , ], type = "l", lty = 1)

library(fanplot)
par(mfrow=c(1, 1))
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
par(mfrow=c(1, 2))
matplot(PLATfit_SPAm$years, t(q_PLAT_SPAm[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.max+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "male mortality rate (log scale)")
fan(t(PLATsim_SPAm.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(PLATsim_SPAm.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(PLATsim_SPAm.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_PLAT_SPAm[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))
matplot(PLATfit_SPAf$years, t(q_PLAT_SPAf[c("65", "75", "85"),c(1:length(Y.fit)) ]),
        xlim = c(y.fit.min, (y.fit.max+y.pred)), ylim = c(0.001, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "female mortality rate (log scale)")
fan(t(PLATsim_SPAf.mrwd$rates["65", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("black", "white")), ln = NULL)
fan(t(PLATsim_SPAf.mrwd$rates["75", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("red", "white")), ln = NULL)
fan(t(PLATsim_SPAf.mrwd$rates["85", , ]), start = y.fit.max+1, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
text(y.fit.min+4, q_PLAT_SPAf[c("65", "75", "85"), 30],
     labels = c("x = 65", "x = 75", "x = 85"))


# estrazione di probabilità di coorte
q_LC_SPAf.1955 <- extractCohort(q_LC_ITAf.ext, cohort = 1955)
q_LC_SPAf.st.1955 <- extractCohort(q_LC_ITAf.st.ext, cohort = 1955)
q_LC_SPA.1980 <- extractCohort(q_LC_ITAf.ext, cohort = 1980)
q_LC_SPAf.st.1980 <- extractCohort(q_LC_ITAf.st.ext, cohort = 1980)


# Bootstrap ---------------------------------------------------------------

LCboot_SPAm <- bootstrap(LCfit_SPAm, nBoot = n.boot, type = "residual")
par(mfrow=c(1, 1))
dev.new()
plot(LCboot_SPAm, nCol = 3)
LCsim_SPAm.boot <- simulate(LCboot_SPAm, nsim = n.sim/n.boot, h = y.pred)
rates_LC_SPAm.boot.st <- LCsim_SPAm.boot$rates
q_LC_SPAm.boot.st <- 1- exp(-rates_LC_SPAm.boot.st)
q_LC_SPAm.boot.st.ext <-  extrapolation.sim(q_LC_SPAm.boot.st)

LCboot_SPAf <- bootstrap(LCfit_SPAf, nBoot = n.boot, type = "semiparametric")
dev.new()
plot(LCboot_SPAf, nCol = 3)
LCsim_SPAf.boot <- simulate(LCboot_SPAf, nsim = n.sim/n.boot, h = y.pred)
rates_LC_SPAf.boot.st <- LCsim_SPAf.boot$rates
q_LC_SPAf.boot.st <- 1- exp(-rates_LC_SPAf.boot.st)
q_LC_SPAf.boot.st.ext <-  extrapolation.sim(q_LC_SPAf.boot.st)


# Intervalli di confidenza ------------------------------------------------

##### Probabilità di morte entro 1 anno
conf.lev <- 0.9
q_LC_SPAm.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_SPAm.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_SPAm.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_SPAm.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_LC_SPAm.q95[j,k] <- quantile(q_LC_SPAm.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_SPAm.q05[j,k] <- quantile(q_LC_SPAm.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_LC_SPAm.boot.q95[j,k] <- quantile(q_LC_SPAm.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_SPAm.boot.q05[j,k] <- quantile(q_LC_SPAm.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}

q_LC_SPAf.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_SPAf.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_SPAf.boot.q95 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
q_LC_SPAf.boot.q05 <- matrix(0, nrow=(omega-a.min),ncol=y.pred, dimnames=list(c(a.min:(omega-1)), c((tail(Y.fit,1)+1):(tail(Y.fit,1)+y.pred))))
for (j in 1:(omega-a.min)){
  for (k in 1:y.pred){
    q_LC_SPAf.q95[j,k] <- quantile(q_LC_SPAf.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_SPAf.q05[j,k] <- quantile(q_LC_SPAf.st.ext[j,k,], probs=(1-conf.lev)/2)
    q_LC_SPAf.boot.q95[j,k] <- quantile(q_LC_SPAf.boot.st.ext[j,k,], probs=0.5+conf.lev/2)
    q_LC_SPAf.boot.q05[j,k] <- quantile(q_LC_SPAf.boot.st.ext[j,k,], probs=(1-conf.lev)/2)
  }
}
q_LC.SPAm.q95<-as.data.frame(q_LC_SPAm.q95)
q_LC.SPAm.q05<-as.data.frame(q_LC_SPAm.q05)
q_LC.SPAm.boot.q95<-as.data.frame(q_LC_SPAm.boot.q95)
q_LC.SPAm.boot.q05<-as.data.frame(q_LC_SPAm.boot.q05)
q_LC.SPAf.q95<-as.data.frame(q_LC_SPAf.q95)
q_LC.SPAf.q05<-as.data.frame(q_LC_SPAf.q05)
q_LC.SPAf.boot.q95<-as.data.frame(q_LC_SPAf.boot.q95)
q_LC.SPAf.boot.q05<-as.data.frame(q_LC_SPAf.boot.q05)

plot(q_LC_SPAm.ext["65",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_SPAm.q95["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_SPAm.q05["65",], x=Y.pred, type="l", lty=2)
lines(q_LC_SPAm.boot.q95["65",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_SPAm.boot.q05["65",], x=Y.pred, type="l", lty=2, col=2)

plot(q_LC_SPAm.ext["90",],x=c(Y.fit,Y.pred), type="l", lty=1) # evolution with t (same age)
lines(q_LC_SPAm.q95["90",], x=Y.pred, type="l", lty=2)
lines(q_LC_SPAm.q05["90",], x=Y.pred, type="l", lty=2)
lines(q_LC_SPAm.boot.q95["90",], x=Y.pred, type="l", lty=2, col=2)
lines(q_LC_SPAm.boot.q05["90",], x=Y.pred, type="l", lty=2, col=2)

##### Probabilità di sopravvivenza a 1 anno,
##### Probabilità di sopravvivenza a n anni,
##### Aspettative di vita
p_LC_SPAm.q95 <- 1-q_LC_SPAm.q95# 1 year survival probability
p0n_LC_SPAm.q95 <- apply(p_LC_SPAm.q95, 2, cumprod) # n year survival probability
ex_LC_SPAm.q95 <- life.exp(q_LC_SPAm.q95) # life expectancy

p_LC_SPAm.q05 <- 1-q_LC_ITAm.q05# 1 year survival probability
p0n_LC_SPAm.q05 <- apply(p_LC_SPAm.q05, 2, cumprod) # n year survival probability
ex_LC_SPAm.q05 <- life.exp(q_LC_SPAm.q05) # life expectancy

p_LC_SPAm.boot.q95 <- 1-q_LC_SPAm.boot.q95# 1 year survival probability
p0n_LC_SPAm.boot.q95 <- apply(p_LC_SPAm.boot.q95, 2, cumprod) # n year survival probability
ex_LC_SPAm.boot.q95 <- life.exp(q_LC_SPAm.boot.q95) # life expectancy

p_LC_SPAm.boot.q05 <- 1-q_LC_SPAm.boot.q05# 1 year survival probability
p0n_LC_SPAm.boot.q05 <- apply(p_LC_SPAm.boot.q05, 2, cumprod) # n year survival probability
ex_LC_SPAm.boot.q05 <- life.exp(q_LC_SPAm.boot.q05) # life expectancy

p_LC_SPAf.q95 <- 1-q_LC_SPAf.q95# 1 year survival probability
p0n_LC_SPAf.q95 <- apply(p_LC_SPAf.q95, 2, cumprod) # n year survival probability
ex_LC_SPAf.q95 <- life.exp(q_LC_SPAf.q95) # life expectancy

p_LC_SPAf.q05 <- 1-q_LC_SPAf.q05# 1 year survival probability
p0n_LC_SPAf.q05 <- apply(p_LC_SPAf.q05, 2, cumprod) # n year survival probability
ex_LC_SPAf.q05 <- life.exp(q_LC_SPAf.q05) # life expectancy

p_LC_SPAf.boot.q95 <- 1-q_LC_SPAf.boot.q95# 1 year survival probability
p0n_LC_SPAf.boot.q95 <- apply(p_LC_SPAf.boot.q95, 2, cumprod) # n year survival probability
ex_LC_SPAf.boot.q95 <- life.exp(q_LC_SPAf.boot.q95) # life expectancy

p_LC_SPAf.boot.q05 <- 1-q_LC_SPAf.boot.q05# 1 year survival probability
p0n_LC_SPAf.boot.q05 <- apply(p_LC_SPAf.boot.q05, 2, cumprod) # n year survival probability
ex_LC_SPAf.boot.q05 <- life.exp(q_LC_SPAf.boot.q05) # life expectancy

ex_LC.SPAm.q05<-as.data.frame(ex_LC_SPAm.q05)
ex_LC.SPAm.q95<-as.data.frame(ex_LC_SPAm.q95)
ex_LC.SPAm.boot.q05<-as.data.frame(ex_LC_SPAm.boot.q05)
ex_LC.SPAm.boot.q95<-as.data.frame(ex_LC_SPAm.boot.q95)
ex_LC.SPAf.q05<-as.data.frame(ex_LC_SPAf.q05)
ex_LC.SPAf.q95<-as.data.frame(ex_LC_SPAf.q95)
ex_LC.SPAf.boot.q05<-as.data.frame(ex_LC_SPAf.boot.q05)
ex_LC.SPAf.boot.q95<-as.data.frame(ex_LC_SPAf.boot.q95)


## Valore attuale di rendite vitalizie 

p_LC_SPAm.st.ext <- 1-q_LC_SPAm.st.ext
ann_LC_SPAm.st <- annuity.st(q_LC_SPAm.st.ext,v_vect,0.99)
ann_LC_SPAm.mean <- ann_LC_SPAm.st[[1]]
ann_LC_SPAm.q95 <- ann_LC_SPAm.st[[2]]
ann_LC_SPAm.q05 <- ann_LC_SPAm.st[[3]]

dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 2))
plot(ann_LC_SPAm.mean["65",], type="l", lty=1) # evoluzione rispetto a t (stessa età x)
lines(ann_LC_SPAm.q95["65",], type="l", lty=2)
lines(ann_LC_SPAm.q05["65",], type="l", lty=2)
plot(ann_LC_SPAm.mean[,"2050"], type="l", lty=1) # evoluzione rispetto a x (stesso anno t)
lines(ann_LC_SPAm.q95[,"2050"], type="l", lty=2)
lines(ann_LC_SPAm.q05[,"2050"], type="l", lty=2)

p_LC_SPAf.st.ext <- 1-q_LC_SPAf.st.ext
ann_LC_SPAf.st <- annuity.st(q_LC_SPAf.st.ext,v_vect,0.9)
ann_LC_SPAf.mean <- ann_LC_SPAf.st[[1]]
ann_LC_SPAf.q95 <- ann_LC_SPAf.st[[2]]
ann_LC_SPAf.q05 <- ann_LC_SPAf.st[[3]]

dev.new(width = 550, height = 330, unit = "px")
par(mfrow=c(1, 2))
plot(ann_LC_SPAf.mean["65",], type="l", lty=1) # evolution with t (same age)
lines(ann_LC_SPAf.q95["65",], type="l", lty=2)
lines(ann_LC_SPAf.q05["65",], type="l", lty=2)
plot(ann_LC_SPAf.mean[,"2050"], type="l", lty=1) # evolution with x (same year)
lines(ann_LC_SPAf.q95[,"2050"], type="l", lty=2)
lines(ann_LC_SPAf.q05[,"2050"], type="l", lty=2)

ann_LC.SPAm.m <- as.data.frame(ann_LC_SPAm.mean)
ann_LC.SPAm.q95 <- as.data.frame(ann_LC_SPAm.q95)
ann_LC.SPAm.q05 <- as.data.frame(ann_LC_SPAm.q05)
ann_LC.SPAf.m <- as.data.frame(ann_LC_SPAf.mean)
ann_LC.SPAf.q95 <-as.data.frame(ann_LC_SPAf.q95)
ann_LC.SPAf.q05<- as.data.frame(ann_LC_SPAf.q05)
