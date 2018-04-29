# The Fate of Haloacetonitriles in Drinking Waters
# Author: Yun Yu, University of Massachusetts Amherst - 2016
# Refined/updated by: William Raseman, University of Colorado Boulder - 2018

# Dependencies: before running this program, you must have JAGS-4.x.y.exe (for any x >=0, y>=0) installed. 
#   This can be found at the following website: https://sourceforge.net/projects/mcmc-jags/files/rjags/4/

# load packages
library(rjags)      # Bayesian graphical models using MCMC
library(mnormt)     # multivariate normal and t distributions
library(Matrix)     # sparse and dense matrix classes and methods
library(scales)     # scale functions for visualizations
library(tidyverse)  # ggplot2, dplyr, stringr 
library(reshape2)   # reshape data: melt() function

#------ Determine Chlorine Decay Logarithmic Model Coefficients ------#

Chlorine_Decay <- read.csv("./data/Chlorine Decay.csv")  # read in chlorine decay data
pH_Cl2_levels <- unique(Chlorine_Decay$pH_Cl2)

CT_decay_model <- data.frame()
pH_Cl2 <- vector()
count <-1
for(i in 1:length(pH_Cl2_levels))
{
  current_CT <- filter(Chlorine_Decay, pH_Cl2==pH_Cl2_levels[i])
  CT_log_coeff <- lm(formula=(current_CT$Ct)/71/1000 ~ log(current_CT$Time))
  CT_decay_model <- rbind(CT_decay_model,t(data.frame(CT_log_coeff$coefficients)))
  pH_Cl2[count] <- as.character(pH_Cl2_levels[i])
  count <- count+1
}
rownames(CT_decay_model)=NULL
CT_decay_model <- cbind(pH_Cl2,CT_decay_model)
colnames(CT_decay_model)=c("pH_Cl2","Intercept","Slope")

#----- Calculate Ct=fn(Time) Integral -----#
Data <- read.csv("./data/Hydrolysis & Chlorination Data.csv")  # read in hydrolysis and chlorination data

Data_Ct_Integral <- data.frame()
for(j in 1:length(pH_Cl2_levels))
{
  Current_Data <- filter(Data,pH_Cl2==as.character(CT_decay_model$pH_Cl2[j]))
  Current_Data$Ct_Integral <- CT_decay_model$Slope[j]*(Current_Data$Time*log(Current_Data$Time)-Current_Data$Time)+CT_decay_model$Intercept[j]*Current_Data$Time
  Data_Ct_Integral <- rbind(Data_Ct_Integral, Current_Data)
}
Data <- Data_Ct_Integral

# define variables
pH_Cl2_levels <- unique(Data$pH_Cl2)
x <- which(Data$Time==1 & Data$Replicate==1)
pH_levels <- Data$pH[x]
Ct_0_levels <- Data$Ct_0[x]
Time <- unique(Data$Time)
Cmpd_levels <- names(Data)[6:(ncol(Data)-1)]
num_pH_Cl2 <- length(pH_Cl2_levels)
num_Cmpd <- length(Cmpd_levels)
num_Time <- length(Time)
pKa <- c(7.536,7.536,7.536,7.536,7.533,7.533,7.533,7.533,7.533,7.522,7.522,7.522,7.522,7.522,7.510,7.510,7.507)
dose <- c(250,250,250,125,125,250,250)
MW <- c(75.5,119.95,109.94,154.39,198.84,144.39,188.84)
ln_dose <- vector()
for(i in 1:length(dose))
{
  ln_dose[i] <- log(dose[i]/MW[i])
}

# define number of replicates for each pH and check how to loop through full dataset
replicate_array <- array(0,(num_pH_Cl2+1))
for (p in 1:num_pH_Cl2) 
{
  Initial <- which(Data$Time==1 & Data$pH_Cl2==pH_Cl2_levels[p])
  replicate_array[p+1] <- length(Initial)
}

# alter the data to remove zeros and take logs of concentrations
melt_Data <- melt(Data, id.vars=c("pH_Cl2","Ct_Integral","Replicate","Time"), measure.vars=c("MCAN","MBAN","DCAN","BCAN","DBAN","TCAN","BDCAN"), variable.name="Compounds", value.name="Concentration")
melt_Data$Log_Conc <- log(melt_Data$Concentration)
melt_Data$Log_Conc[melt_Data$Log_Conc==-Inf] <- NA
melt_Ct_Integral <- select(melt_Data, pH_Cl2, Ct_Integral, Time, Replicate)


FINAL_DATA <- acast(melt_Data,Time ~ pH_Cl2+Replicate ~ Compounds, value.var="Log_Conc")
CT_INTEGRAL <- acast(melt_Ct_Integral,Time ~ pH_Cl2+Replicate, value.var="Ct_Integral",fun.aggregate=mean)

#----- Bayesian Calibration_HOMOSCEDASTIC -----#
win.data <- list(FINAL_DATA = FINAL_DATA, num_cmpd = num_Cmpd, num_pH_Cl2 = num_pH_Cl2, replicate_array = replicate_array,
                 pH_levels = pH_levels, ln_dose = ln_dose, pKa = pKa, CT_INTEGRAL = CT_INTEGRAL, Time = Time, num_Time = num_Time)

params <- c("pH","k_H2O","k_OH","k_HOCl","k_OCl","precission_e","lnC0")
source("./R/Bayesian Model_Dynamic CT.R")

##Setting up JAGS run
nc <- 3
ni <- 20000
nb <- 20000

nt <- 5
hymod <- jags.model(
  file = "model.txt",
  data = win.data,
  n.chains = nc,
  n.adapt = nb
)

(mid1 <- Sys.time())

out <- jags.samples(
  model = hymod,
  variable.names = params,
  n.iter = ni,
  thin = nt,
  progress.bar = 'text'
)
save(out, file="Output_4k.RData")


#----- Check for Convergence -----#

com <- 6
par(mfrow=c(2,2))
plot(out$k_H2O[com,,1],type="l",col="black")
lines(out$k_H2O[com,,2],type="l",col="red")
lines(out$k_H2O[com,,3],type="l",col="blue")

plot(out$k_OH[com,,1],type="l",col="black")
lines(out$k_OH[com,,2],type="l",col="red")
lines(out$k_OH[com,,3],type="l",col="blue")

plot(out$k_HOCl[com,,1],type="l",col="black")
lines(out$k_HOCl[com,,2],type="l",col="red")
lines(out$k_HOCl[com,,3],type="l",col="blue")

plot(out$k_OCl[com,,1],type="l",col="black")
lines(out$k_OCl[com,,2],type="l",col="red")
lines(out$k_OCl[com,,3],type="l",col="blue")

# plot posterior distribution of k_H2O, k_OH, k_HOCl, k_OCl
par(mfrow=c(7,4),mar=c(2,2,2,2))
for (com in 1:7) 
{
  hist(out$k_H2O[com,,],main="",xlab="",ylab="")
  lines(c(0,0),c(-10000,10000),col="red",lwd=3)
  hist(out$k_OH[com,,],main="",xlab="",ylab="")
  lines(c(0,0),c(-10000,10000),col="red",lwd=3)  
  hist(out$k_HOCl[com,,],main="",xlab="",ylab="")
  lines(c(0,0),c(-10000,10000),col="red",lwd=3)  
  hist(out$k_OCl[com,,],main="",xlab="",ylab="")
  lines(c(0,0),c(-10000,10000),col="red",lwd=3)  
}

# plot posterior distribution of pH
par(mfrow=c(5,4),mar=c(2,0,0,0))
for(p in 1:length(pH_levels))
{
  hist(out$pH[p,,], main="", xlab="", ylab="")
}

# plot posterior distribution of precission_e
par(mfrow=c(10,8), mar=c(2,0,0,0))
for(t in 1:74)
{
  hist(out$precission_e[,t,,], main="", xlab="", ylab="")
}

#----- Evaluate Model -----#

# posterior of k_H2O, k_OH, k_HOCl, k_OCl and pH
k_H2O_posterior <- as.data.frame(t(apply(out$k_H2O,FUN=quantile,1,c(0.025,0.25,0.5,0.75,0.975))))
k_H2O_posterior <- cbind(Cmpd_levels,k_H2O_posterior)

k_OH_posterior <- as.data.frame(t(apply(out$k_OH,FUN=quantile,1,c(0.025,0.25,0.5,0.75,0.975))))
k_OH_posterior <- cbind(Cmpd_levels,k_OH_posterior)

k_HOCl_posterior <- as.data.frame(t(apply(out$k_HOCl,FUN=quantile,1,c(0.025,0.25,0.5,0.75,0.975))))
k_HOCl_posterior <- cbind(Cmpd_levels,k_HOCl_posterior)

k_OCl_posterior <- as.data.frame(t(apply(out$k_OCl,FUN=quantile,1,c(0.025,0.25,0.5,0.75,0.975))))
k_OCl_posterior <- cbind(Cmpd_levels, k_OCl_posterior)

pH_posterior <- as.data.frame(t(apply(out$pH, FUN=quantile,1,c(0.025,0.25,0.5,0.75,0.975))))
pH_Estimates <- as.data.frame(apply(out$pH, FUN=quantile, 1, c(0.5)))
pH_Estimates <- cbind(pH_Cl2_levels, pH_levels, pH_Estimates)
colnames(pH_Estimates) <- c("pH_Cl2","pH_levels","pH_Estimates")

# k_median, k_upper and k_lower
k_H2O_median <- apply(out$k_H2O,FUN=quantile,1,0.5)
k_H2O_upper <- apply(out$k_H2O,FUN=quantile,1,0.975)
k_H2O_lower <- apply(out$k_H2O,FUN=quantile,1,0.025)

k_OH_median <- apply(out$k_OH,FUN=quantile,1,0.5)
k_OH_upper <- apply(out$k_OH,FUN=quantile,1,0.975)
k_OH_lower <- apply(out$k_OH,FUN=quantile,1,0.025)

k_HOCl_median <- apply(out$k_HOCl,FUN=quantile,1,0.5)
k_HOCl_upper <- apply(out$k_HOCl,FUN=quantile,1,0.975)
k_HOCl_lower <- apply(out$k_HOCl,FUN=quantile,1,0.025)

k_OCl_median <- apply(out$k_OCl,FUN=quantile,1,0.5)
k_OCl_upper <- apply(out$k_OCl,FUN=quantile,1,0.975)
k_OCl_lower <- apply(out$k_OCl,FUN=quantile,1,0.025)

k_Estimates <- as.data.frame(t(rbind(k_H2O_median,k_H2O_upper,k_H2O_lower,k_OH_median,k_OH_upper,k_OH_lower, k_HOCl_median, k_HOCl_upper,k_HOCl_lower,k_OCl_median,k_OCl_upper,k_OCl_lower)))
k_Estimates <- cbind(Cmpd_levels, k_Estimates)
names(k_Estimates)[names(k_Estimates)=="Cmpd_levels"] <- "Compounds"

# median lnC0 for individual trials
lnC0 <- as.data.frame(t(apply(out$lnC0, FUN=quantile,c(1,2),c(0.5))))
lnC0 <- cbind(colnames(FINAL_DATA), lnC0)
colnames(lnC0) <- c("pH_Cl2_Rep", "MCAN","MBAN","DCAN","BCAN","DBAN","TCAN","BDCAN")
Replicate <- str_sub(lnC0$pH_Cl2_Rep, start=-1, end=-1)
lnC0 <- cbind(Replicate, lnC0)
melt_lnC0 <- melt(lnC0, id.vars=c("pH_Cl2_Rep","Replicate"), measure.vars=c("MCAN","MBAN","DCAN","BCAN","DBAN","TCAN","BDCAN"), variable.name="Compounds", value.name="lnC0")


#----- Plot predicted ln[HAN] and actual ln[HAN] versus Time for all the compounds -----#

Data <- left_join(Data, pH_Estimates, by="pH_Cl2", copy=TRUE) %>% select(-pH_levels)

alpha <- data.frame(pH=unique(pH_Estimates$pH_Estimates),pKa=c(7.536,7.536,7.536,7.536,7.533,7.533,7.533,7.533,7.533,7.522,7.522,7.522,7.522,7.522,7.510,7.508,7.507))
alpha$alpha_0 <- 10^(-alpha$pH)/(10^(-alpha$pH)+10^(-alpha$pKa))
alpha$alpha_1 <- 10^(-alpha$pKa)/(10^(-alpha$pH)+10^(-alpha$pKa))

for(com in 1:num_Cmpd)
{
  Current_Cmpd <- cbind(Data[,1],Data[,14],Data[,2:5],Data[,13],Data[,5+com])
  colnames(Current_Cmpd)=c("pH","pH_Estimates","Ct_0","pH_Cl2","Replicate","Time","Ct_Integral","Concentration")
  Current_Cmpd$pH_Cl2_Rep <- paste(Current_Cmpd$pH_Cl2, Current_Cmpd$Replicate, sep="_", collapse=NULL)
  Current_Cmpd$Log_Conc <- log(Current_Cmpd$Concentration)
  Current_Cmpd$Log_Conc[Current_Cmpd$Log_Conc==-Inf] <- NA
  
  k_Cmpd <- k_Estimates[com,]
  
  lnC0_Cmpd <- filter(melt_lnC0, Compounds==Cmpd_levels[com])
  
  Model_Prediction <- data.frame()
  for(r in 1:length(lnC0_Cmpd$pH_Cl2_Rep))
  {
    Current_Data <- filter(Current_Cmpd, pH_Cl2_Rep==lnC0_Cmpd$pH_Cl2_Rep[r])
    Current_Data$Predict_median <- lnC0_Cmpd$lnC0[r]-(k_Cmpd$k_H2O_median + k_Cmpd$k_OH_median*10^(-14+Current_Data$pH_Estimates))*Current_Data$Time - (k_Cmpd$k_HOCl_median*alpha$alpha_0[alpha$pH==unique(Current_Data$pH_Estimates)] + k_Cmpd$k_OCl_median*alpha$alpha_1[alpha$pH==unique(Current_Data$pH_Estimates)])*Current_Data$Ct_Integral 
    Current_Data$Predict_upper <- lnC0_Cmpd$lnC0[r]-(k_Cmpd$k_H2O_upper + k_Cmpd$k_OH_upper*10^(-14+Current_Data$pH_Estimates))*Current_Data$Time - (k_Cmpd$k_HOCl_upper*alpha$alpha_0[alpha$pH==unique(Current_Data$pH_Estimates)] + k_Cmpd$k_OCl_upper*alpha$alpha_1[alpha$pH==unique(Current_Data$pH_Estimates)])*Current_Data$Ct_Integral 
    Current_Data$Predict_lower <- lnC0_Cmpd$lnC0[r]-(k_Cmpd$k_H2O_lower + k_Cmpd$k_OH_lower*10^(-14+Current_Data$pH_Estimates))*Current_Data$Time - (k_Cmpd$k_HOCl_lower*alpha$alpha_0[alpha$pH==unique(Current_Data$pH_Estimates)] + k_Cmpd$k_OCl_lower*alpha$alpha_1[alpha$pH==unique(Current_Data$pH_Estimates)])*Current_Data$Ct_Integral 
    Model_Prediction <- rbind(Model_Prediction, Current_Data)
    Model_Prediction <- na.omit(Model_Prediction)
  }
  
  Fit <- ggplot() + geom_point(data=Model_Prediction, aes(x=Time, y=Log_Conc), pch=16, size=4, na.rm=TRUE) +
                    geom_line(data=Model_Prediction, aes(x=Time, y=Predict_median), linetype=1) +
                    #geom_line(data=Model_Prediction, aes(x=Time, y=Predict_upper), linetype=4) +
                    #geom_line(data=Model_Prediction, aes(x=Time, y=Predict_lower), linetype=4) +
                    facet_wrap(~pH_Cl2, ncol=5)
}

#----- Lower-level Liner Regression of Log_Conc on Time -----#

Hydrolysis_Data <- filter(Data, Ct_0==0) %>% select(pH, Replicate, Time, MCAN, MBAN, DCAN, BCAN, DBAN, TCAN, BDCAN, pH_Estimates)
melt_Hydrolysis_Data <- melt(Hydrolysis_Data, id.vars=c("pH","pH_Estimates","Replicate","Time"), measure.vars=c("MCAN","MBAN","DCAN","BCAN","DBAN","TCAN","BDCAN"), variable.name="Compounds", value.name="Concentration")
melt_Hydrolysis_Data$Log_Conc <- log(melt_Hydrolysis_Data$Concentration)
melt_Hydrolysis_Data$Log_Conc[melt_Hydrolysis_Data$Log_Conc==-Inf] <- NA

Trial <- with(melt_Hydrolysis_Data, paste(pH,Replicate,Compounds,sep = "_"))
melt_Hydrolysis_Data$Trial <- Trial

Hydrolysis_split <- split(melt_Hydrolysis_Data, f=Trial)
Hydrolysis_split <- lapply(Hydrolysis_split, na.omit)
Hydrolysis_split <- Hydrolysis_split[sapply(Hydrolysis_split, nrow) > 0]

lm <- lapply(Hydrolysis_split, lm, formula=Log_Conc ~ Time)
slopes <- vapply(lm, function(mod) mod$coef[2], numeric(1))
intercepts <- vapply(lm, function(mod) mod$coef[1], numeric(1))

lower_level_reg <- data.frame(Trial=names(slopes), Kobs=-slopes, Intercept=intercepts) %>%
  left_join(melt_Hydrolysis_Data, by = "Trial", copy=FALSE) %>%
  select(-Time, -Concentration, -Log_Conc) %>%
  unique()
lower_level_reg$Trial <- as.factor(lower_level_reg$Trial)

Hydrolysis <- ggplot() + geom_point(data=lower_level_reg, aes(x=10^(-14+pH_Estimates), y=Kobs), pch=1, size=4, color="blue", na.rm=TRUE) +
                         geom_abline(data=k_Estimates, aes(slope=k_OH_median, intercept=k_H2O_median), linetype=1, color="black") +
                         #geom_abline(data=k_Estimates, aes(slope=k_OH_upper, intercept=k_H2O_upper), linetype=4, color="black") +
                         #geom_abline(data=k_Estimates, aes(slope=k_OH_lower, intercept=k_H2O_lower), linetype=4, color="black") +
                         facet_wrap(~Compounds,ncol=4,scales="free")
Hydrolysis +theme_bw() +
            theme(panel.border=element_rect(color="black")) +
            labs(title="HANs Hydrolysis Predictive Model") + xlab("pH") + ylab(expression(bold(paste("K"[obs],"(1/hr)")))) +         
            scale_x_continuous(breaks=c(10^(-14+6),10^(-14+7),10^(-14+8),10^(-14+8.5),10^(-14+9),10^(-14+9.2)),labels=c("6","","8","8.5","9","9.2")) +
            theme(plot.title=element_text(size=rel(1.5),face="bold", vjust=1.5),
                  axis.title.x=element_text(size=14,face="bold",vjust=0.1),
                  axis.title.y=element_text(size=14,face="bold",vjust=1.0),
                  axis.text.x=element_text(size=10,vjust=0.7,face="bold"),
                  axis.text.y=element_text(size=10,face="bold"),
                  strip.text.x=element_text(size=10, face="bold"), 
                  strip.background=element_rect(color="white",fill="white"))
            
