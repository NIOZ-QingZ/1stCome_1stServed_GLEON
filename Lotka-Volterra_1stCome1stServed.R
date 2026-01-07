# Title: "Lotka-Volterra_1stCome1stServed.R"
# Description: Model is based on Lotka-Volterra competition equations. 
# Version 1: two species and four parameters including initial biomass, growth rates, competition effects, carrying capacity, species introudcing day.
# The code includes:
#   section 1: the model; 
# notes: NA
# section 2; sensitivity analysis; 
# notes: NA
# section 3: visualization
# notes: the sensitivity analysis takes more than 1hr to run. The simulation results are saved in the csv file entitled "sen_tIntro_K_P.csv". If you just want to play with the simulation results and don't do new sensitivy analysis, you can skip this section and jump to the next section about visualization.


########################################
#### section 1: Lotka-Volterra model############
#############################
# Load required libraries
library(deSolve)
library(ggplot2)
library(reshape2)  # For reshaping the data for ggplot

# Define the phytoplankton competition model function
phytoplankton_model <- function(t, state, parameters) {
  # Unpack state variables (phytoplankton species populations)
  P1 <- state[1]  # Population of phytoplankton species 1
  P2 <- state[2]  # Population of phytoplankton species 2
  
  # Unpack parameters
  r1 <- parameters["r1"]  # Growth rate of species 1
  r2 <- parameters["r2"]  # Growth rate of species 2
  alpha12 <- parameters["alpha12"]  # Effect of species 2 on species 1
  alpha21 <- parameters["alpha21"]  # Effect of species 1 on species 2
  K1 <- parameters["K1"]  # Carrying capacity for species 1 (nutrient/light limitation)
  K2 <- parameters["K2"]  # Carrying capacity for species 2 (nutrient/light limitation)
  t_intro2 <- parameters["t_intro2"]  # Time when species 2 is introduced
  
  # Apply condition for the introduction of species 2
  if (t < t_intro2) {
    P2 <- 0  # Before introduction time, species 2 is absent
  }
  
  # Define phytoplankton-specific growth equations
  dP1_dt <- r1 * P1 * (1 - (P1 + alpha12 * P2) / K1)  # Logistic growth with competition for species 1
  dP2_dt <- r2 * P2 * (1 - (P2 + alpha21 * P1) / K2)  # Logistic growth with competition for species 2
  
  return(list(c(dP1_dt, dP2_dt)))
}

phytoplankton_model_2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{

    # Define phytoplankton-specific growth equations
    dP1_dt <- r1 * P1 * (1 - (alpha11*P1 + alpha12 * P2))  # Logistic growth with competition for species 1
    dP2_dt <- r2 * P2 * (1 - (alpha22*P2 + alpha21 * P1))
    
    return(list(c(dP1_dt, dP2_dt)))
  })
}

##### sampling until meeting conditions based on fitness ratio and niche overlap
set.seed(123)  # For reproducibility

# Function to check the conditions
priority_effect_conditions <- function(a11, a12, a22, a21) {
  f2_to_f1 = sqrt((a11 * a12) / (a22 * a21))
  rho =  sqrt((a12 * a21) / (a11 * a22))
  
  cond1 <- f2_to_f1 < rho & f2_to_f1 > 1 / rho
  cond2 <- (1-rho) < 0
  
  return(cond1 && cond2)
}
coexistence_conditions <- function(a11, a12, a22, a21) {
  f2_to_f1 = sqrt((a11 * a12) / (a22 * a21))
  rho =  sqrt((a12 * a21) / (a11 * a22))
  
  cond1 <- f2_to_f1 > rho & f2_to_f1 < 1 / rho
  cond2 <- (1-rho) > 0
  return(cond1 && cond2)
}
spec1_win_conditions <- function(a11, a12, a22, a21) {
  f2_to_f1 = sqrt((a11 * a12) / (a22 * a21))
  rho =  sqrt((a12 * a21) / (a11 * a22))
  
  cond1 <- f2_to_f1 > rho & f2_to_f1 < 1 / rho
  cond2 <- f2_to_f1 > 1 / rho & f2_to_f1 < rho
  cond3 <- f2_to_f1 > 1e0
  return(!cond1 && !cond2 && !cond3)
}
spec2_win_conditions <- function(a11, a12, a22, a21) {
  f2_to_f1 = sqrt((a11 * a12) / (a22 * a21))
  rho =  sqrt((a12 * a21) / (a11 * a22))
  
  cond1 <- f2_to_f1 > rho & f2_to_f1 < 1 / rho
  cond2 <- f2_to_f1 > 1 / rho & f2_to_f1 < rho
  cond3 <- f2_to_f1 > 1e0
  return(!cond1 && !cond2 && cond3)
}

priority_effect_conditions(a11,a12,a22,a21)
coexistence_conditions(a11,a12,a22,a21)
spec1_win_conditions(a11,a12,a22,a21)
spec2_win_conditions(a11,a12,a22,a21)


# Initialize variables
a11 <- a12 <- a22 <- a21 <- 0

# Loop until both priority effects conditions are met
repeat {
  a11 <- runif(1, min=0, max=0.01)
  a12 <- runif(1, min=0, max=1)
  a22 <- runif(1, min=0, max=0.01)
  a21 <- runif(1, min=0, max=1)
  
  if (coexistence_conditions(a11, a12, a22, a21)) {
    break
  }
}



# # Print the sampled values
# cat("a11:", a11, "\n")
# cat("a12:", a12, "\n")
# cat("a22:", a22, "\n")
# cat("a21:", a21, "\n")

####################
## dynamics of species 1 and 2
###################
 
# parameters <- c(
#     r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
#     r2 = 1.6,      # Growth rate of phytoplankton species 2
#     alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
#     alpha21 = 1, # Competition effect of species 1 on species 2
#     K1 = 100,      # Carrying capacity for species 1 (related to nutrient/light availability)
#     K2 = 100,      # Carrying capacity for species 2 (related to nutrient/light availability)
#     t_intro2 = t_intro2  # Introduce species 2 at time = 10
#   )

#######################  
parameters_2 <- c(
  r1 = 1,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
  r2 = 10,      # Growth rate of phytoplankton species 2
  alpha11 = a11, # Competition effect of species 2 on species 1 (higher competition)
  alpha22 = a22, # Competition effect of species 1 on species 2
  alpha12 = a12,      # Carrying capacity for species 1 (related to nutrient/light availability)
  alpha21 = a21,      # Carrying capacity for species 2 (related to nutrient/light availability)
  t_intro = 0  # Introduce species 2 at time = 10
)

# Define the time points for the simulation
times <- seq(0, 300, by = 1)

# Define initial conditions and parameters
initial_state <- c(P1 = 11, P2 = 10)  # Start with species 1 present, species 2 absent

t_intro2 = 0  
# Run the simulation using ode solver
# output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)

#### following chesson ####
output_2 <- as.data.frame(ode(y = initial_state, times = times, func = phytoplankton_model_2, parms = parameters_2),
                          method = "lsoda", events = NULL)

print(c(spec1_win_conditions=spec1_win_conditions(a11, a12, a22, a21),
        spec2_win_conditions=spec2_win_conditions(a11, a12, a22, a21),
        priority_effect_conditions=priority_effect_conditions(a11, a12, a22, a21),
        coexistence_conditions=coexistence_conditions(a11, a12, a22, a21)))
tail(output_2[,2:3],n=6L)
c(a11=a11,a12=a12,a22=a22,a21=a21)

output2_ylim=range(c(output_2[,2:3]),na.rm = T)

plot(output_2$time,output_2$P1,
     ylim=output2_ylim,col="red",t="l", 
     main="initial biomass P1>P2",
     ylab="Density",xlab="Time")
lines(output_2$time,output_2$P2,col="green")
legend("right",legend=c("species 1","species 2"), col=c("red","green"), lty=c(1,1),text.col = c("red","green"))

# Convert the output to a data frame
output_df <- as.data.frame(output)

plot(output_df$time, output_df$P1, col="red", ylab="species abundance",ylim=c(0,100),t="l")
lines(output_df$time, output_df$P2, col="green")
legend("topright",legend=c("species 1 (a=0.8)","species 2 (a=1.6)"), col=c("red","green"), lty=c(1,1),text.col = c("red","green"))

########################################
#### section 2";" sensitivity analysis############
#############################


################ 
### section 2.1. senstivity analysis of individual parameter
##############

#### 
## for introducing days of t
sen_t_intro <- NULL

t_seq = seq(0,10, by=0.5)

for (t_intro2 in t_seq) {
  
  # Define initial conditions and parameters
  initial_state <- c(P1 = 10, P2 = 10)  # Start with species 1 present, species 2 absent
  
  parameters <- c(
    r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
    r2 = 0.8,      # Growth rate of phytoplankton species 2
    alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
    alpha21 = 1, # Competition effect of species 1 on species 2
    K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
    K2 = 90,      # Carrying capacity for species 2 (related to nutrient/light availability)
    t_intro2 = t_intro2  # Introduce species 2 at time = 10
  )
  
  # Define the time points for the simulation
  times <- seq(0, 1000, by = 1)
  
  # Run the simulation using ode solver
  output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
  
  # Convert the output to a data frame
  output_df <- as.data.frame(output)
  
  # Melt the data frame for ggplot
  output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
  
  P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
  P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
  
  diff_P1_to_P2 = P1_end_biomass-P2_end_biomass
  P1_perc = P1_end_biomass
  
  sen_t_intro = rbind(sen_t_intro, data.frame(t_intro2=t_intro2, 
                                              diff_P1_to_P2=diff_P1_to_P2,
                                              P1_end_biomass=P1_end_biomass,
                                              P2_end_biomass=P2_end_biomass,
                                              sum_biomass=P1_end_biomass+P2_end_biomass))
}

# plot(sen_t_intro$t_intro2, sen_t_intro$diff_P1_to_P2)

plot(sen_t_intro$t_intro2, sen_t_intro$P1_end_biomass/sen_t_intro$sum_biomass, ylab="percentage of species 1", xlab="introducing day of species 2 (x days after species 1)")

######
## for carry capacity K
######

sen_K <- NULL

K2_to_K1_seq = seq(0.98,1.05, by=0.0005)

for (K2_to_K1 in K2_to_K1_seq) {
  
  # Define initial conditions and parameters
  initial_state <- c(P1 = 10, P2 = 10)  # Start with species 1 present, species 2 absent
  K1 = 90
  parameters <- c(
    r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
    r2 = 0.8,      # Growth rate of phytoplankton species 2
    alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
    alpha21 = 1, # Competition effect of species 1 on species 2
    K1 = 90,     # Carrying capacity for species 1 (related to nutrient/light availability)
    K2 = K2_to_K1*K1,      # Carrying capacity for species 2 (related to nutrient/light availability)
    t_intro2 = 0  # Introduce species 2 at time x
  )
  
  # Define the time points for the simulation
  times <- seq(0, 1000, by = 1)
  
  # Run the simulation using ode solver
  output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
  
  # Convert the output to a data frame
  output_df <- as.data.frame(output)
  
  # Melt the data frame for ggplot
  output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
  
  P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
  P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
  
  diff_P1_to_P2 = P1_end_biomass-P2_end_biomass
  P1_perc = P1_end_biomass
  
  sen_K = rbind(sen_K, data.frame(K2_to_K1=K2_to_K1, 
                                              diff_P1_to_P2=diff_P1_to_P2,
                                              P1_end_biomass=P1_end_biomass,
                                              P2_end_biomass=P2_end_biomass,
                                              sum_biomass=P1_end_biomass+P2_end_biomass))
}


plot(sen_K$K2_to_K1, sen_K$P1_end_biomass/sen_K$sum_biomass, ylab="percentage of species 1", xlab="carry capacity K2/K1")
abline(h=0.5, col="red")
abline(v=1, col="red")

######
## for initial states P
######

sen_P <- NULL

P2_to_P1_seq = seq(0.01,10, by=0.2)

for (P2_to_P1 in P2_to_P1_seq) {
  
  # Define initial conditions and parameters
  P1 = 10
  P2 = P1*P2_to_P1
  initial_state <- c(P1 = P1, P2 = P2)  # Start with species 1 present, species 2 absent
  
  parameters <- c(
    r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
    r2 = 0.8,      # Growth rate of phytoplankton species 2
    alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
    alpha21 = 1, # Competition effect of species 1 on species 2
    K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
    K2 = 90,      # Carrying capacity for species 2 (related to nutrient/light availability)
    t_intro2 = 0  # Introduce species 2 at time x
  )
  
  # Define the time points for the simulation
  times <- seq(0, 1000, by = 1)
  
  # Run the simulation using ode solver
  output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
  
  # Convert the output to a data frame
  output_df <- as.data.frame(output)
  
  # Melt the data frame for ggplot
  output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
  
  P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
  P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
  
  sen_P = rbind(sen_P, data.frame(P2_to_P1=P2_to_P1, 
                                  P1_end_biomass=P1_end_biomass,
                                  P2_end_biomass=P2_end_biomass,
                                  sum_biomass=P1_end_biomass+P2_end_biomass))
}


plot(sen_P$P2_to_P1, sen_P$P2_end_biomass/sen_P$sum_biomass, ylab="percentage of species 2", xlab="initial biomass P2/P1")
abline(h=0.5, col="red")
abline(v=1, col="red")


######
## for growth rate alpha
######

sen_growth_rate <- NULL

r2_to_r1_seq = seq(0.01,5, by=0.5)

for (r2_to_r1 in r2_to_r1_seq) {
  r1 = 0.8
  # Define initial conditions and parameters
  P1 = 10
  P2 = 10
  initial_state <- c(P1 = P1, P2 = P2)  # Start with species 1 present, species 2 absent
  
  parameters <- c(
    r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
    r2 = r1*r2_to_r1,      # Growth rate of phytoplankton species 2
    alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
    alpha21 = 1, # Competition effect of species 1 on species 2
    K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
    K2 = 90,      # Carrying capacity for species 2 (related to nutrient/light availability)
    t_intro2 = 0  # Introduce species 2 at time x
  )
  
  # Define the time points for the simulation
  times <- seq(0, 1000, by = 1)
  
  # Run the simulation using ode solver
  output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
  
  # Convert the output to a data frame
  output_df <- as.data.frame(output)
  
  # Melt the data frame for ggplot
  output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
  
  P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
  P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
  
  sen_growth_rate = rbind(sen_growth_rate, data.frame(r2_to_r1=r2_to_r1, 
                                  P1_end_biomass=P1_end_biomass,
                                  P2_end_biomass=P2_end_biomass,
                                  sum_biomass=P1_end_biomass+P2_end_biomass))
}


plot(sen_growth_rate$r2_to_r1, sen_growth_rate$P2_end_biomass/sen_growth_rate$sum_biomass, ylab="percentage of species 2", xlab="fitness (growth rate) r2/r1")
abline(h=0.5, col="red")
abline(v=1, col="red")


#####################
###### 2.2. sensitivity analysis of all parameters:
#####################
# !!!Attention!!!
# the sensitivity analysis takes more than 1hr to run. The simulation results are saved in the csv file entitled "sen_tIntro_K_P.csv". If you just want to play with the simulation results and don't do new sensitivy analysis, you can skip this section and jump to the next section about visualization.

sen_tIntro_K_P_r <- NULL

t_seq = seq(from=0,to=10, length.out=10)
# K2_to_K1_seq = seq(from=0.98,to=1.03, length.out=10)
K2_to_K1_seq = 1
P2_to_P1_seq = seq(from=0.01,to=10, length.out=10)
r2_to_r1_seq = seq(from=0.01,to=5, length.out=10)

for (r2_to_r1 in r2_to_r1_seq) {
  for (t_intro2 in t_seq) {
    for (K2_to_K1 in K2_to_K1_seq) {
      for (P2_to_P1 in P2_to_P1_seq) {
    
      # Define initial conditions and parameters
      P1 = 10
      P2 = P1*P2_to_P1
      initial_state <- c(P1 = P1, P2 = P2)  # Start with species 1 present, species 2 absent
      r1 = 0.8
      K1 = 90
      
      parameters <- c(
      r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
      r2 = r1*r2_to_r1,      # Growth rate of phytoplankton species 2
      alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
      alpha21 = 1, # Competition effect of species 1 on species 2
      K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
      K2 = K2_to_K1*K1,      # Carrying capacity for species 2 (related to nutrient/light availability)
      t_intro2 = t_intro2  # Introduce species 2 at time x
    )
    
    # Define the time points for the simulation
    times <- seq(0, 10000, by = 10)
    
    # Run the simulation using ode solver
    output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
    
    # Convert the output to a data frame
    output_df <- as.data.frame(output)
    
    # Melt the data frame for ggplot
    output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
    
    P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
    P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
    
    sen_tIntro_K_P_r = rbind(sen_tIntro_K_P_r, 
                         data.frame(r2_to_r1=r2_to_r1,
                                    t_intro2=t_intro2,
                                    K2_to_K1=K2_to_K1,
                                    P2_to_P1=P2_to_P1,
                                    P1_end_biomass=P1_end_biomass,
                                    P2_end_biomass=P2_end_biomass,
                                    sum_biomass=P1_end_biomass+P2_end_biomass,
                                    P1_Perc=P1_end_biomass/(P1_end_biomass+P2_end_biomass)))
  }
}
}
}

write.csv(sen_tIntro_K_P_r,"sen_tIntro_K_P_r.csv")


head(sen_tIntro_K_P_r)


######################
#### section 3: visualization of senstivity analysis
######################
require("plot3D")
## Loading required package: plot3D
## Warning: package 'plot3D' was built under R version 3.6.3
require(plot3Drgl) # In an extension package, plot3Drgl, a similar function, plotrgl, plots the graphs to the device, opened with rgl. This allows interactive zooming, rotating, etc... 

my_wd = readClipboard()
setwd(my_wd)

par(cex=1)
# the sensitivity analysis results.
sen_tIntro_K_P <- read.csv("sen_tIntro_K_P.csv")

x <- sen_tIntro_K_P$t_intro2
y <- sen_tIntro_K_P$K2_to_K1
z <- sen_tIntro_K_P$P2_to_P1
P1_Perc <- sen_tIntro_K_P$P1_Perc


png("sensitivity_analysis.png",width = 16,height = 12,units = "cm",res = 300)
scatter3D(x,y,z,
          colvar = P1_Perc,
          xlab="A: Spec 2 introducing d",
          ylab="B: Carrying cap K2/K1",
          zlab="C: Init biomass P2/P1",
          clab = c("D: Spec 1 percentage"),
          pch=20,phi=0,theta = 60,bty="g",main="",
          ticktype="detailed", # to have axes ticks
          cex=1,
          facets=NA,
          alpha=0.6, lighting=TRUE,
          colkey = list(tick=FALSE,labels=TRUE,cex.axis=0.8),
          col=(terrain.colors(100)))


dev.off()
jet.col (n = 100, alpha = 1)

# This opens a device window, which allows interactive zooming, rotating, etc...
scatter3Drgl(x,y,z,
             colvar = P1_Perc,
             xlab="A: Spec 2 introducing d",
             ylab="B: Carrying cap K2/K1",
             zlab="C: Init biomass P2/P1",
             clab = c("D: Spec 1 percentage"),
             pch=20,phi=0,theta = 60,bty="g",main="",
             ticktype="detailed", # to have axes ticks
             cex=1,
             facets=NA,
             alpha=0.6, lighting=TRUE,
             colkey = list(tick=FALSE,labels=TRUE,cex.axis=0.8),
             col=(terrain.colors(100)))

#######################################

# the sensitivity analysis results.
sen_tIntro_K_P_r <- read.csv("sen_tIntro_K_P_r.csv")

t_intro <- sen_tIntro_K_P_r$t_intro2
Carry_K <- sen_tIntro_K_P_r$K2_to_K1
Init_P <- sen_tIntro_K_P_r$P2_to_P1
Growth_r <- sen_tIntro_K_P_r$r2_to_r1

P1_Perc <- sen_tIntro_K_P_r$P1_Perc
P2_Perc <- 1-P1_Perc

png("sensitivity_analysis.png",width = 16,height = 12,units = "cm",res = 300)
scatter3D(t_intro,Carry_K,Init_P,
          colvar = P1_Perc,
          xlab="A: Spec 2 introducing d",
          ylab="B: Carrying cap K2/K1",
          zlab="C: Init biomass P2/P1",
          clab = c("D: Spec 1 percentage"),
          pch=20,phi=0,theta = 60,bty="g",main="",
          ticktype="detailed", # to have axes ticks
          cex=1,
          facets=NA,
          alpha=0.6, lighting=TRUE,
          colkey = list(tick=FALSE,labels=TRUE,cex.axis=0.8),
          col=(terrain.colors(100)))


dev.off()
jet.col (n = 100, alpha = 1)

# This opens a device window, which allows interactive zooming, rotating, etc...
scatter3Drgl(t_intro,Carry_K,Init_P,
             colvar = P2_Perc,
             xlab="A: Spec 2 introducing d",
             ylab="B: Carrying cap K2/K1",
             zlab="C: Init biomass P2/P1",
             clab = c("D: Spec 2 percentage"),
             pch=20,phi=0,theta = 60,bty="g",main="",
             ticktype="detailed", # to have axes ticks
             cex=1,
             facets=NA,
             alpha=0.6, lighting=TRUE,
             colkey = list(tick=FALSE,labels=TRUE,cex.axis=0.8),
             col=(terrain.colors(100)))

scatter3Drgl(Growth_r,Carry_K,Init_P,
             colvar = P2_Perc,
             xlab="A: growth r2/r1",
             ylab="B: Carrying cap K2/K1",
             zlab="C: Init biomass P2/P1",
             clab = c("D: Spec 2 percentage"),
             pch=20,phi=0,theta = 60,bty="g",main="",
             ticktype="detailed", # to have axes ticks
             cex=1,
             facets=NA,
             alpha=0.6, lighting=TRUE,
             colkey = list(tick=FALSE,labels=TRUE,cex.axis=0.8),
             col=(terrain.colors(100)))

scatter3Drgl(Growth_r,t_intro,Init_P,
             colvar = P2_Perc,
             xlab="A: growth r2/r1",
             ylab="B: Spec 2 introducing d",
             zlab="C: Init biomass P2/P1",
             clab = c("D: Spec 2 percentage"),
             pch=20,phi=0,theta = 60,bty="g",main="",
             ticktype="detailed", # to have axes ticks
             cex=1,
             facets=NA,
             expand=TRUE,
             rasterImage=TRUE,
             alpha=0.9, lighting=TRUE,
             colkey = list(tick=FALSE,labels=TRUE,cex.axis=0.8),
             col=rev(heat.colors(100)))


scatter3D(Growth_r,t_intro,Init_P,
          colvar = P2_Perc,
          xlab="A: growth r2/r1",
          ylab="B: Spec 2 introducing d",
          zlab="C: Init biomass P2/P1",
          clab = c("D: Spec 1 percentage"),
          pch=20,phi=0,theta = 60,bty="g",main="",
          ticktype="detailed", # to have axes ticks
          cex=1,
          facets=NA,
          alpha=0.8, lighting=TRUE,
          expand=TRUE,
          colkey = list(tick=FALSE,labels=TRUE,cex.axis=0.8),
          col=(terrain.colors(100)))
# scatter3D(x=df.regime2$Q_thres,y=df.regime2$DOC_thres,z=df.regime2$W_T*149.6,colvar=df.regime2$DOC_C*100,add = T,type="l",colkey = F,lwd=10)

# write.csv(sen_tIntro_K_P,".../1stCome1stServed/sen_tIntro_K_P.csv")

#####################
###### 2.3. sensitivity analysis of MORE parameters:
#####################
# !!!Attention!!!
# the sensitivity analysis takes more than 1hr to run. The simulation results are saved in the csv file entitled "sen_tIntro_K_P_r_alpha.csv". If you just want to play with the simulation results and don't do new sensitivy analysis, you can skip this section and jump to the next section about visualization.

sen_tIntro_K_P <- NULL

t_seq = seq(0,10, by=0.5)
K2_to_K1_seq = seq(0.98,1.03, by=0.001)
P2_to_P1_seq = seq(0.01,10, by=0.15)

for (t_intro2 in t_seq) {
  for (K2_to_K1 in K2_to_K1_seq) {
    for (P2_to_P1 in P2_to_P1_seq) {
      
      # Define initial conditions and parameters
      P1 = 10
      P2 = P1*P2_to_P1
      initial_state <- c(P1 = P1, P2 = P2)  # Start with species 1 present, species 2 absent
      
      parameters <- c(
        r1 = 0.8,      # Growth rate of phytoplankton species 1 (higher because phytoplankton can grow fast)
        r2 = 0.8,      # Growth rate of phytoplankton species 2
        alpha12 = 1, # Competition effect of species 2 on species 1 (higher competition)
        alpha21 = 1, # Competition effect of species 1 on species 2
        K1 = 90,      # Carrying capacity for species 1 (related to nutrient/light availability)
        K2 = K2_to_K1*K1,      # Carrying capacity for species 2 (related to nutrient/light availability)
        t_intro2 = t_intro2  # Introduce species 2 at time x
      )
      
      # Define the time points for the simulation
      times <- seq(0, 1000, by = 1)
      
      # Run the simulation using ode solver
      output <- ode(y = initial_state, times = times, func = phytoplankton_model, parms = parameters)
      
      # Convert the output to a data frame
      output_df <- as.data.frame(output)
      
      # Melt the data frame for ggplot
      output_melted <- melt(output_df, id.vars = "time", variable.name = "Species", value.name = "Population")
      
      P1_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P1")]
      P2_end_biomass = output_melted$Population[which(output_melted$time==max(times)&output_melted$Species=="P2")]
      
      sen_tIntro_K_P = rbind(sen_tIntro_K_P, 
                             data.frame(t_intro2=t_intro2,
                                        K2_to_K1=K2_to_K1,
                                        P2_to_P1=P2_to_P1,
                                        P1_end_biomass=P1_end_biomass,
                                        P2_end_biomass=P2_end_biomass,
                                        sum_biomass=P1_end_biomass+P2_end_biomass))
    }
  }
}

head(sen_tIntro_K_P)
sen_tIntro_K_P$P1_Perc = sen_tIntro_K_P$P1_end_biomass/sen_tIntro_K_P$sum_biomass

