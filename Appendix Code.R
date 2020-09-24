#### Analysis of distance sampling data: live observations and dung piles
#### By simulating a plausible set of forward distances

library(mvtnorm)

#### Load 2D distance functions
source("com_hfunctions.R")
source("com_pifunctions.R")
source("com_likelihoodutilities.R")

#### Load dataset
data <- read.csv("montiCORRIGE20122013.csv", sep=";")
colnames(data)[11] = "x" ## rename perpendicular distance to x
data$heure_deb = as.POSIXct(as.character(data$heure_deb),format = "%d/%m/%Y %H:%M:%S")
data$heure_fin = as.POSIXct(as.character(data$heure_fin),format = "%d/%m/%Y %H:%M:%S")
data$heure_obs = as.POSIXct(as.character(data$heure_obs),format = "%d/%m/%Y %H:%M:%S")

#### Transect lengths 
transect.names = c("A","I","J","K","L")
transect.lengths = c(8, 6, 5.8, 6, 5.4)   # A I J K L
names(transect.lengths) =  transect.names
data = data[data$transect %in% transect.names,]
dim(data)
table(data$transect, as.POSIXlt(data$heure_obs)$year)

#### Generate forward distance data
####    (scenario letters correspond to Fig. 3 in the main text)
if (scenario == "a")
  data$y = pmax(0.01, data$x + rnorm(length(data$x), -1+(data$contact=="Crotte"), 1))
if (scenario == "b") {
  mx=max(data$x)
  data$y = pmax(0.01, (mx/2)*exp(-data$x/(mx/2)) + rnorm(length(data$x), -1+(data$contact=="Crotte"), 1))
}
if (scenario == "c") {
  b = c(0.2,0.2)
  x<-runif(length(data$x), px(13-0.01,rev(b),h.RE,13,nint=100),1)
  f <- function(x, u) px(x,rev(b),h.RE,13,nint=100) - u
  my.uniroot <- function(x) uniroot(f, interval=c(0, 13), tol = 0.0001, u = x)$root
  data$y <- vapply(x, my.uniroot, numeric(1))
}


#### Divide into year- and cue-specific datasets
montiVU2012 = data[as.POSIXlt(data$heure_obs)$year==112 & data$contact=="Vu" & data$type_obs=="Nocturne",]
montiCR2012 = data[as.POSIXlt(data$heure_obs)$year==112 & data$contact=="Crotte" & data$type_obs=="Diurne",]
montiVU2013 = data[as.POSIXlt(data$heure_obs)$year==113 & data$contact=="Vu" & data$type_obs=="Nocturne",]
montiCR2013 = data[as.POSIXlt(data$heure_obs)$year==113 & data$contact=="Crotte" & data$type_obs=="Diurne",]
unlist(lapply(list(montiVU2012,montiCR2012,montiVU2013,montiCR2013), function(x) dim(x)[1]))


#### Fit 2D distance sampling model using multiple initial value to avoid local minima in the deviance
 y = c(montiVU2012$y , montiVU2013$y)
 x = c(montiVU2012$x , montiVU2013$x)
 hr = h.RE
 pi.x = pi.sigmo
 ystart = ceiling(max(y))
 w = ceiling(max(x))
 length.b = 2
 debug=FALSE

 FIT=list(); dev=NULL
 for (m in 1:10) {
   pars = rnorm(4, c(0.25,0.25,-4,-1), 3)
   tmp0 <- tryCatch.W.E (
     fityx(y,x,pars[1:length.b],hr,ystart,pi.x,pars[(length.b+1):length(pars)],w,control=list(),hessian=TRUE,corrFlag=0.7,debug=FALSE)
   )
   fit = NA
   if(! "error" %in% class(tmp0$value)) {
     fit <- tmp0$value
     fit$vcov <-  matrix(Matrix::nearPD(fit$vcov)$mat,4,4)
   }
   FIT[[m]] = fit
   if(is.na(fit[1])) dev=c(dev, 1e12) else dev = c(dev, fit$val)
 }
 fitVU = FIT[[which.min(dev)]]
 tabVU = matrix(NA,2,3)
 if(is.na(fitVU[1])) tabVU = matrix(NA,2,3) else {
 tmp1 <- tryCatch.W.E (boot(fitVU))
 if(! "error" %in% class(tmp1$value))  tabVU=tmp1$value
 }
 
 y = c(montiCR2012$y , montiCR2013$y) 
 x = c(montiCR2012$x , montiCR2013$x)
 hr = h.RE
 pi.x = pi.sigmo
 ystart = ceiling(max(y))
 w = ceiling(max(x))
 length.b = 2
 debug=FALSE

 FIT=list(); dev=NULL
 for (m in 1:10) {
   pars = rnorm(4, 0, 3)
   tmp0 <- tryCatch.W.E (
     fityx(y,x,pars[1:length.b],hr,ystart,pi.x,pars[(length.b+1):length(pars)],w,control=list(),hessian=TRUE,corrFlag=0.7,debug=FALSE)
   )
   fit = NA
   if(! "error" %in% class(tmp0$value)) {
     fit <- tmp0$value
     fit$vcov <-  matrix(Matrix::nearPD(fit$vcov)$mat,4,4)
   }
   FIT[[m]] = fit
   if(is.na(fit[1])) dev=c(dev, 1e12) else dev = c(dev, fit$val)
 }
 fitCR = FIT[[which.min(dev)]]
 tabCR = matrix(NA,2,3)
 if(is.na(fitCR[1])) tabCR = matrix(NA,2,3) else {
 tmp1 <- tryCatch.W.E (boot(fitCR))
 if(! "error" %in% class(tmp1$value))  tabCR=tmp1$value
 }
 
### Conversion Ratio dung/live
tabCR[2,]*(fitVU$w/fitCR$w)/tabVU[2,]
### Population density
tabVU[2,]/sum(transect.lengths)/(fitVU$w/1000) /2 ## divided by two to get the average over 2 years
### Fit details (parameters values, fitted detection and avoidance functions) are in the objects fitVU and fitCR for live and dung data respectively

