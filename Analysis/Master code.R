
############################################################
#
# Master code to examine data 
# Futu Chen
# 1/26/2023
# Update Notes:
# dlnm + mixmeta (1.26)
# I LOVE PRESSURE
#############################################################
setwd("/Users/madaopt/Desktop/APPLY FOR JOB/HSPH BST assignment/wildfire_covid19_casestudy/Analysis")
dat <- read.csv("/Users/madaopt/Desktop/APPLY FOR JOB/HSPH BST assignment/wildfire_covid19_casestudy/Data/moddat_Feb2021.csv")
model_dir <- "/Users/madaopt/Desktop/APPLY FOR JOB/HSPH BST assignment/wildfire_covid19_casestudy/Data/"

library(lubridate)
library(mgcv)
library(splines)
library(dlnm)
library(mixmeta)
library(dplyr)
library(parameters)

#check dataset
str(dat)
summary(duplicated(dat$FIPS)) #92 unique FIPS
summary(dat$pm25) #current day pm, imputed
summary(dat$wildfire) #T/F whether wildfire event 
summary(dat$cases) #daily cases (contains NA)
summary(dat$pm_ambient) #df$pm_ambient = df$pm25 - df$pm_wildfire, PM from wildfire contribution to pm
summary(dat$relative_change_feb) #1662 NA, movement
summary(dat$ratio_travelers) #1662 NA, movement
summary(dat$tmmx) #climate
summary(dat$rmax) #climate 
dat$date_dt <- as.Date(dat$date, "%Y-%m-%d")
# dat$population: offset; exposure: pm25; outcome: cases/deaths ?cumu_deaths?
# how many days contains zero cases
checkzero <-  dat %>% group_by(FIPS) %>% 
  filter(cases== 0) %>% 
  summarise(pctmiss = n()/277 *100) %>%  
  ungroup()
checkzero <-  dat %>% group_by(FIPS) %>% 
  filter(deaths== 0) %>% 
  summarise(pctmiss = n()/277 *100) %>%  
  ungroup()
#wow greater than I thought. Have to consider the zeros.

table(dat$StateFIPS, dat$State)

# Figure 2 provided time series plots - refer to paper 
# seems to have seasonality. Cannot find a significant jump that cooccur with exposure - exclude ITS

#exposure window use paper: PM 0-1 avg 
lag <-function(x,n) {
  nn<-length(x)
  xn<-c(rep(NA,n),x[1:(nn-n)])
  return(xn)
}
dat$pm25_lag1 <- lag(dat$pm25,1)
dat$pm25_01avg <- (dat$pm25+dat$pm25_lag1)/2
check <- select(dat, c("pm25", "pm25_lag1","pm25_01avg"))

#covaraites: consider splines
#control for meteorological season?
dat$season <- ifelse(dat$month %in% c(12,1,2),"winter",ifelse(dat$month %in% c(3,4,5),"spring", ifelse(dat$month %in% c(6,7,8), "summer", ifelse(dat$month %in% c(9,10,11),"fall", NA))))
table(dat$season)
#or use a spline with a pattern going up and down as the seasonality - df=5 for X years
length(table(year(dat$date_dt))) #ah i forgot only one year

#tmmx and rmax, first check linearity with penalized spline
smd1 <- gam(cases ~ pm25_01avg + s(date_num,bs = "cr",fx=TRUE,k=5) + s(tmmx, bs="cr",fx=FALSE) + s(rmax, bs="cr",fx=FALSE) + as.factor(dayofweek) + s(relative_change_feb, bs="cr",fx=FALSE) + s(ratio_travelers, bs="cr",fx=FALSE), family=ziP(), offset=log(population), data = dat, na.action=na.omit)#zero-inflated 
summary(smd1) 
sum(residuals(smd1)^2) / smd1$df.residual #...overdispersion
plot(smd1,scale=0,pages=1)
# choose a df: keep df=5 for time, tmmx df=3? rmax  df= 5 natural spline
# ratio_travelers looks as if linear. df2 maybe. relative_change_feb hard to tell, keep df=3

#exclude 1662 NA
smd2 <- gam(cases ~ pm25_01avg + s(date_num,bs = "cr",fx=TRUE,k=5) + s(tmmx, bs="cr",fx=FALSE) + s(rmax, bs="cr",fx=FALSE) + as.factor(dayofweek) , family=ziP(), offset=log(population), data = dat, na.action=na.omit)#zero-inflated 
summary(smd2) 
smd2$aic
plot(smd2,scale=0,pages=1) #looks like 2-3 df for temp(df=2)  and rh (df=3) is fine

smd3 <- gam(cases ~ pm25_01avg + s(date_num,bs = "cr",fx=TRUE,k=5) + s(tmmx, bs="cr",fx=FALSE, k=2) + s(rmax, bs="cr",fx=FALSE, k=3) + as.factor(dayofweek) , family=ziP(), offset=log(population), data = dat, na.action=na.omit)#zero-inflated 
summary(smd3) 
plot(smd3,scale=0,pages=1)
smd3$aic # not too bad

mdlist <- list(smd1,smd2,smd3)
saveRDS(mdlist, paste0(model_dir,"cov_explore.rds"))


#additional outcome contributed by wildfire event - use proxiy (additional PM). interaction term in 2nd stage
#use same lag as main exposure 
dat$pm_ambient_lag1 <- lag(dat$pm_ambient,1)
dat$pm_ambient_01avg <- (dat$pm_ambient+dat$pm_ambient_lag1)/2
check <- select(dat, c("pm_ambient", "pm_ambient_lag1","pm_ambient_01avg"))
dat$FIPS.f <- factor(dat$FIPS)
fips.list <- unique(dat$FIPS.f)

# I can either do a 2-staged meta regression or gamm
# go dlnm and mixveta!
# create the crossbasis objects
library(dlnm)
# Here is the problem: if dlnm work with gam, then cannot load into mixmeta
# if dlnm with glm, then can use mixmeta pacakge for 2nd stage, but has to do a quai possion and no penalization
# choose the second option 

# start with simplest linear functions 
cb1 <- crossbasis(dat$pm25, lag=28, argvar=list("lin"), arglag=list("lin")) 
cb2 <- crossbasis(dat$pm25, lag=28, argvar=list("poly",degree=2), arglag=list("lin")) # exposure function a poly of degree 2
cb3 <- crossbasis(dat$pm25, lag=28, argvar=list("poly",degree=4), arglag=list("lin")) # exposure function a poly of degree 4
cb4 <- crossbasis(dat$pm25, lag=28, argvar=list("lin"), arglag=list(fun="strata"))  #lag function 
cb5 <- crossbasis(dat$pm25, lag=28, argvar=list("lin"), arglag=list(fun="ns",knots=logknots(28,4))) #change knots
cb5 <- crossbasis(dat$pm25, lag=28, argvar=list("lin"), arglag=list(fun="ns",knots=logknots(28,8)))
cb6 <- crossbasis(dat$pm25, lag=28, argvar=list("poly",degree=2), arglag=list(fun="ns",knots=logknots(28,4)))

md.cb1 <- glm(cases ~ cb1 + ns(date_num, df=5) + ns(tmmx, df=2) + ns(rmax, df=3) + as.factor(dayofweek), dat ,family=quasipoisson, offset=log(population), na.action="na.exclude") 
md.cb2 <- glm(cases ~ cb2 + ns(date_num, df=5) + ns(tmmx, df=2) + ns(rmax, df=3) + as.factor(dayofweek), dat ,family=quasipoisson, offset=log(population), na.action="na.exclude") 
md.cb3 <- glm(cases ~ cb3 + ns(date_num, df=5) + ns(tmmx, df=2) + ns(rmax, df=3) + as.factor(dayofweek), dat ,family=quasipoisson, offset=log(population), na.action="na.exclude") 
md.cb4 <- glm(cases ~ cb4 + ns(date_num, df=5) + ns(tmmx, df=2) + ns(rmax, df=3) + as.factor(dayofweek), dat ,family=quasipoisson, offset=log(population), na.action="na.exclude") 
md.cb5 <- glm(cases ~ cb5 + ns(date_num, df=5) + ns(tmmx, df=2) + ns(rmax, df=3) + as.factor(dayofweek), dat ,family=quasipoisson, offset=log(population), na.action="na.exclude") 
md.cb6 <- glm(cases ~ cb6 + ns(date_num, df=5) + ns(tmmx, df=2) + ns(rmax, df=3) + as.factor(dayofweek), dat ,family=quasipoisson, offset=log(population), na.action="na.exclude") 

mdlist <- list(md.cb1,md.cb2,md.cb3,md.cb4,md.cb5,md.cb6)
saveRDS(mdlist, paste0(model_dir,"cb_explore.rds"))
plot(md.cb1, which=1)
plot(md.cb2, which=1)
plot(md.cb3, which=1)
plot(md.cb4, which=1)
plot(md.cb5, which=1)
plot(md.cb6, which=1)

pm.cb <- cb1
# selected cb
md.try <- glm(cases ~ pm.cb + ns(date_num, df=5) + ns(tmmx, df=2) + ns(rmax, df=3) + as.factor(dayofweek), dat ,family=quasipoisson, offset=log(population), na.action="na.exclude") 
#md1 <- gam(cases ~ pm25_01avg + pm.cb + s(date_num,bs = "cr",fx=TRUE,k=5) + s(tmmx, bs="cr",fx=FALSE, k=2) + s(rmax, bs="cr",fx=FALSE, k=3) + as.factor(dayofweek) , family=ziP(), offset=log(population), data = dat, na.action=na.omit)
summary(md.try) 
#reduce to one dim
md.reduce <- crossreduce(pm.cb, md.try,cen=mean(dat$pm25,na.rm=T))
mmt <- md.reduce$predvar[which.min(md.reduce$fit)[[1]]]
md.try2 <- crossreduce(pm.cb, md.try,cen=mmt)
coef(md.try2)
vcov(md.try2)
plot(md.try2)
pred1.pm <- crosspred(pm.cb, md.try, at=0:28, bylag=1, cumul=TRUE)
plot(pred1.pm, "slices", var=10, col=3, ylab="RR", ci.arg=list(density=15,lwd=2), main="Association with a 10-unit increase in PM25") # longer lag seemed to be protective 
plot(pred1.pm, "slices", var=10, col=2, cumul=TRUE, ylab="Cumulative RR", main="Cumulative association with a 10-unit increase in PM25") # cumulatively seemed to max around 15 days and then decay over time 


############### expand this to each county #################
fips <- split(dat, dat$FIPS.f)
names(fips)

coefall <- matrix(NA,length(fips),2,dimnames=list(names(fips)))
vcovall <- vector("list",length(fips))
names(vcovall) <- names(fips)
mdall <- vector("list",length(fips))

for (i in 1:length(fips.list)) {
  data <- filter(dat, FIPS.f %in% fips.list[i])
  # cb for pm
  cb <-  crossbasis(data$pm25, lag=28, argvar=list("lin"), arglag=list(fun="lin")) 
  md <- glm(cases ~ cb + ns(date_num, df=5) + ns(tmmx, df=2) + ns(rmax, df=3) + as.factor(dayofweek), data , family=quasipoisson, offset=log(population), na.action="na.exclude") 
  # #reduce to one dim
  # md.reduce <- crossreduce(cb, md,cen=mean(data$pm25,na.rm=T))
  # mmt <- md.reduce$predvar[which.min(md.reduce$fit)[[1]]]
  # mdall[[i]] <- crossreduce(cb, md,cen=mmt)
  # coefall[i,] <- coef(mdall[[i]])
  # vcovall[[i]] <- vcov(mdall[[i]])
  predmean <- crosspred(cb,md,cen=mean(data$pm25,na.rm=T))
  mdall[[i]] <- predmean
   coefall[i,] <- coef(mdall[[i]])
   vcovall[[i]] <- vcov(mdall[[i]])
 }#FIPS
for(i in seq(length(mdall))) {
  plot(mdall[[i]], main= paste0("1st stage ",names(fips[i])))
}



# mixmeta
# 2nd stage
# null
md.null <- mixmeta(coefall ~ 1, vcovall, control=list(showiter=T, igls.inititer=10),method="reml")
summary(md.null)
# Multivariate Cochran Q-test for heterogeneity:
#   Q = 409.8073 (df = 182), p-value = 0.0000
# I-square statistic = 55.6%

#something is not correct...

blupall <- blup(md.null,vcov=T)

# I do not understand anymore. I cannot keep going. 


#####ignore####



# go, dlnm and gamm! 
# mgcv approach
# per 10 unit increase in PM and weather variables
varlist <- names(dat[c(33,37, 11,12)])
for (i in 1:length(varlist)) {
  dat[[paste0(varlist[i],"_10")]] <- dat[[paste0(varlist[i])]]/10
}

dat$wildfire.f <- factor(dat$wildfire)
dat$FIPS.f <- factor(dat$FIPS)
# mixed effects models with a random intercept for each county. Not including mobility!
pm.cb <- crossbasis(dat$pm25, lag=28, argvar=list("poly",degree=2,scale=1), # exposure function a poly of degree 2
                    arglag=list("lin")) # linear lag function

md2 <- gamm(cases ~ pm25_01avg + pm.cb + s(date_num,bs = "cr",fx=TRUE,k=5) + s(tmmx, bs="cr",fx=FALSE, k=2) + s(rmax, bs="cr",fx=FALSE, k=3) + as.factor(dayofweek), random=list(FIPS.f =~1) , family=ziP(), offset=log(population), data = dat, na.action=na.omit) #failed to run. family not supported 
#alternative, bam with random intercepts - adjust the height of other modelterms with a constant value
md2 <- bam(cases ~ pm25_01avg + pm.cb + s(date_num,bs = "cr",fx=TRUE,k=5) + s(tmmx, bs="cr",fx=FALSE, k=2) + s(rmax, bs="cr",fx=FALSE, k=3) + as.factor(dayofweek) +s(FIPS.f, bs="re") , family=ziP(), offset=log(population), data = dat, na.action=na.omit)
summary(md2)  
 


