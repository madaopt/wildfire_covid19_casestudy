
############################################################
#
# Master code to examine data 
# Futu Chen
# 1/25/2023
# Update Notes:
# 
# I LOVE PRESSURE
#############################################################
setwd("/Users/madaopt/Desktop/APPLY FOR JOB/HSPH BST assignment/wildfire_covid19_casestudy/Analysis")
dat <- read.csv("/Users/madaopt/Desktop/APPLY FOR JOB/HSPH BST assignment/wildfire_covid19_casestudy/Data/moddat_Feb2021.csv")

library(lubridate)
library(mgcv)
library(splines)
library(AER)
library(dplyr)
library(mixmeta)
library(dlnm)
library(mgwrsar)


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

# Figure 2 provided time series plots - refer to paper 
# seems to have seasonality. Cannot find a jump that occur with exposure - exclude ITS

#exposure window use paper: PM 0-1 avg 
lag <-function(x,n) {
  nn<-length(x)
  xn<-c(rep(NA,n),x[1:(nn-n)])
  return(xn)
}
dat$pm25_lag1 <- lagm(dat$pm25,1)
dat$pm25_01avg <- (dat$pm25+dat$pm25_lag1)/2
check <- select(dat, c("pm25", "pm25_lag1","pm25_01avg"))
#delayed outcome - distributed lag 28 days - LATER
#covaraites: consider splines - GCV

#control for meteorological season
dat$season <- ifelse(dat$month %in% c(12,1,2),"winter",ifelse(dat$month %in% c(3,4,5),"spring", ifelse(dat$month %in% c(6,7,8), "summer", ifelse(dat$month %in% c(9,10,11),"fall", NA))))
table(dat$season)

#or use a spline with a pattern going up and down as the seasonality - df=5 for X years
length(table(year(dat$date_dt))) #ah i forgot only one year
#tmmx and rmax, first check linearity with penalized spline
smd1 <- gam(cases ~ pm25_01avg + s(date_num,bs = "cr",fx=TRUE,k=5) + s(tmmx, bs="cr",fx=FALSE) + s(rmax, bs="cr",fx=FALSE) + as.factor(dayofweek) + s(relative_change_feb, bs="cr",fx=FALSE) + s(ratio_travelers, bs="cr",fx=FALSE), family=ziP(), offset=log(population), data = dat, na.action=na.omit)#zero-inflated 
summary(smd1) # I don't want to use df=8 or 9 (keep some wiggling in the bottom,)
sum(residuals(smd1)^2) / smd1$df.residual #...overdispersion
plot(smd1,scale=0,pages=1)
# choose a df: keep df=5 for time, tmmx df=3? rmax  df= 5 natural spline
# ratio_travelers looks as if linear. df2 maybe. relative_change_feb hard to tell, keep df=3

#include 1662 NA
smd2 <- gam(cases ~ pm25_01avg + s(date_num,bs = "cr",fx=TRUE,k=5) + s(tmmx, bs="cr",fx=FALSE) + s(rmax, bs="cr",fx=FALSE) + as.factor(dayofweek) , family=ziP(), offset=log(population), data = dat, na.action=na.omit)#zero-inflated 
summary(smd2) 
plot(smd2,scale=0,pages=1) #looks like 2-3 df for temp and rh is fine

#additional outcome contributed by wildfire event - use proxiy (additional PM). interaction term in 2nd stage
#use same lag as main exposure 
dat$pm_ambient_lag1 <- lagm(dat$pm_ambient,1)
dat$pm_ambient_01avg <- (dat$pm_ambient+dat$pm_ambient_lag1)/2
check <- select(dat, c("pm_ambient", "pm_ambient_lag1","pm_ambient_01avg"))

# I can either do a 2-staged meta regression or gamm
# sorry I don't know how to code stage 2..how can I carry standared errors into stage 2
# go, gamm! 


# per 10 unit increase in PM and weather variables
varlist <- names(dat[c(33,37, 11,12)])
for (i in 1:length(varlist)) {
  dat[[paste0(varlist[i],"_10")]] <- dat[[paste0(varlist[i])]]/10
}

#linear mixed effects models with a random intercept for each county.




#MGWR
#limited ability, collapse the deepth of the data 
dat2 <-dat
dat2 <- dat2 %>% group_by(FIPS)  %>%
  summarize(avg_pm=mean(pm25), avg_cases=mean(cases, na.rm=TRUE), avg_tmmx = mean(tmmx, na.rm=TRUE), avg_rmax = mean(rmax, na.rm=TRUE), avg_change = mean(relative_change_feb, na.rm=TRUE), avg_travel=mean(ratio_travelers, na.rm=TRUE), population=mean(population))

 #       mutate(avg_pm=mean(pm25), avg_cases=mean(cases, na.rm=TRUE), avg_tmmx = mean(tmmx, na.rm=TRUE), avg_rmax = mean(rmax, na.rm=TRUE), avg_change = mean(relative_change_feb, na.rm=TRUE), avg_travel=mean(ratio_travelers, na.rm=TRUE))

dat2$y_rate = dat2$avg_cases/dat2$population
coord <- select(dat, c("FIPS","Lat","Long"))
coord=as.matrix(dat[,c("Lat","Long")])
 ## Creating a spatial weight matrix ## of 4 nearest neighbors with 0 in diagonal 
W = kernel_matW(H=4,kernels='rectangle',coord_i=coord,NN=4,adaptive=TRUE, diagnull=TRUE,rowNorm=TRUE)

space.md <-MGWRSAR(formula = 'y_rate ~ avg_pm+ avg_tmmx + avg_rmax+ season + avg_change + avg_travel', data = dat2, coord=coord, fixed_vars=c("season","dayofweek"),kernels=c('gauss'),H=20, Model = 'MGWRSAR_0_kc_kv', control=list(SE=FALSE,adaptive=TRUE,W=W))  #model choice: MGWRSAR_1_kc_kv

summary_mgwrsar(mgwrsar_0_kc_kv)
  
