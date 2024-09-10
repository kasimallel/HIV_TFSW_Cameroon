######################################################
#              HIV - modelling analysis            #
######################################################
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/HIV_modellingC/Code")

if (!require(remotes)) { install.packages("remotes"); stopifnot(require("remotes")) } 
# Load in the deSolve package
if (!require(deSolve)) { install.packages("deSolve"); stopifnot(require("deSolve")) } 
library(deSolve)
library(ggplot2)
library(deSolve)
library(reshape)
library(dplyr)
library(patchwork)
library(grid) 
library(reshape2)
library(ggpattern)
library(scales) 
library(ggplot2)
library(tidyr)
library(heatmaply)
library(plotly)
library(orca)

#######Define functions and parameters below

# Define model function
HIV_model <- function(times, state, parms) {
  ## Define variables
  # Low risk female 
  S1 <- state["S1"]
  E1 <- state["E1"]
  I1 <- state["I1"]
  R1 <- state["R1"]
  N1 <- S1 + E1 + I1
  
  # Low risk male
  S2 <- state["S2"]
  E2 <- state["E2"]
  I2 <- state["I2"]
  R2 <- state["R2"]
  N2 <- S2 + E2 + I2 
  
  # Client
  S3 <- state["S3"]
  E3 <- state["E3"]
  I3 <- state["I3"]
  R3 <- state["R3"]
  N3 <- S3 + E3 + I3   
  
  # TFSW
  S4 <- state["S4"]
  E4 <- state["E4"]
  I4 <- state["I4"]
  R4 <- state["R4"]
  N4 <- S4 + E4 + I4  
  
  dCost<- state["Cost"]
  dDalys<-state["Dalys"]
  donset_inf<- state["onset_inf"]
  #N total
  N = N1 + N2 + N3 + N4
  
  # Extract parameters
  phi_1 <- parms["phi_1"]
  phi_2 <- parms["phi_2"]
  p <- parms["p"]
  theta <- parms["theta"]
  gamma <- parms["gamma"]
  kappa <- parms["kappa"]
  mu_1 <- parms["mu_1"]
  mu_2 <- parms["mu_2"]
  mu_3 <- parms["mu_3"]
  mu_4 <- parms["mu_4"]
  eta_1 <- parms["eta_1"]
  eta_2 <- parms["eta_2"]
  delta_1 <- parms["delta_1"]
  delta_2 <- parms["delta_2"]
  z <- parms["z"]
  g <- parms["g"]
  beta <- parms["beta"]
  epsilon <- parms["epsilon"]
  pi_12m <- parms["pi_12m"]
  pi_12c <- parms["pi_12c"]
  pi_13m <- parms["pi_13m"]
  pi_13c <- parms["pi_13c"]  
  pi_24m <- parms["pi_24m"]
  pi_24c <- parms["pi_24c"]
  pi_34m <- parms["pi_34m"]
  pi_34c <- parms["pi_34c"]
  pi_34co <- parms["pi_34co"]
  Psi_12m <- parms["Psi_12m"]
  Psi_12c <- parms["Psi_12c"]
  Psi_13m <- parms["Psi_13m"]
  Psi_13c <- parms["Psi_13c"]
  Psi_24m <- parms["Psi_24m"]
  Psi_24c <- parms["Psi_24c"]
  Psi_34m <- parms["Psi_34m"]
  Psi_34c <- parms["Psi_34c"]
  nu_1 <- parms["nu_1"]
  nu_2 <- parms["nu_2"]
  nu_3 <- parms["nu_3"]
  nu_4 <- parms["nu_4"]
  e_nu <- parms["e_nu"]
  alpha_1 <- parms["alpha_1"]
  alpha_2 <- parms["alpha_2"]
  alpha_3 <- parms["alpha_3"]
  alpha_4 <- parms["alpha_4"]
  e_alpha <- parms["e_alpha"]
  n_1m <- parms["n_1m"]
  n_1c <- parms["n_1c"]
  n_2m <- parms["n_2m"]
  n_2c <- parms["n_2c"]
  n_3m <- parms["n_3m"]
  n_3c <- parms["n_3c"]
  n_3co <- parms["n_3co"]
  n_4m <- parms["n_4m"]
  n_4c <- parms["n_4c"]
  n_4co <- parms["n_4co"]
  s <- parms["s"]
  e_s <- parms["e_s"]
  tau_1 <- parms["tau_1"]
  tau_2 <- parms["tau_2"]
  tau_3 <- parms["tau_3"]
  tau_4 <- parms["tau_4"]
  e_tau <- parms["e_tau"]
  P1_inf <- parms["P1_inf"]
  P2_inf <- parms["P2_inf"]
  P3_inf <- parms["P3_inf"]
  P4_inf <- parms["P4_inf"]
  w <- parms["w"]
  #Utilities & costs below
  health_insur_c <-parms["health_insur_c"] 
  prep_cost <- parms["prep_cost"]
  art_cost <- parms["art_cost"]
  daly_acute_hiv_offART <- parms["daly_acute_hiv_offART"]
  daly_preAIDs_offART <- parms["daly_preAIDs_offART"]
  daly_AIDs_offART <- parms["daly_AIDs_offART"]
  
  #probablities
  
  rho_1m <- (n_1m*N1)/(n_1m*N1+n_4m*N4)
  rho_1c <- (n_1c*N1)/(n_1c*N1+n_4c*N4)
  rho_2m <- (n_2m*N2)/(n_2m*N2+n_3m*N3)
  rho_2c <- (n_2c*N2)/(n_2c*N2+n_3c*N3)
  rho_3m <- (n_3m*N3)/(n_3m*N2+n_3m*N3)
  rho_3c <- (n_3c*N3)/(n_3c*N2+n_3c*N3)
  rho_4m <- (n_4m*N4)/(n_1m*N1+n_4m*N4)
  rho_4c <- (n_4c*N4)/(n_1c*N1+n_4c*N4)
  
  #Prevalence of HIV
  
  P1_t <- (E1+I1)/N1
  P2_t <- (E2+I2)/N2
  P3_t <- (E3+I3)/N3
  P4_t <- (E4+I4)/N4
  
  #lambda
  
  lambda_1 <- w*exp((P1_t/P1_inf)*log(1/w))
  lambda_2 <- w*exp((P2_t/P2_inf)*log(1/w))
  lambda_3 <- w*exp((P3_t/P3_inf)*log(1/w))
  lambda_4 <- w*exp((P4_t/P1_inf)*log(1/w))
  
  #B_i
  
  B_1 <- ((1-alpha_1*e_alpha)*E1+(1-tau_1*e_tau)*I1)/N1
  B_2 <- ((1-alpha_2*e_alpha)*E2+(1-tau_2*e_tau)*I2)/N2
  B_3 <- ((1-alpha_3*e_alpha)*E3+(1-tau_3*e_tau)*I3)/N3
  B_4 <- ((1-alpha_4*e_alpha)*E4+(1-tau_4*e_tau)*I4)/N4
  
  #Lambdas
  
  Lambda_1m <- lambda_1*(rho_1m/N1)*(1-nu_1*e_nu)*(beta*(1-epsilon*pi_12m)*Psi_12m*n_2m*N2*B_2+beta*(1-epsilon*pi_13m)*Psi_13m*n_3m*N3*B_3)
  Lambda_1c <- lambda_1*(rho_1c/N1)*(1-nu_1*e_nu)*(beta*(1-epsilon*pi_12c)*Psi_12c*n_2c*N2*B_2+beta*(1-epsilon*pi_13c)*Psi_13c*n_3c*N3*B_3)
  Lambda_2m <- lambda_2*n_2m*(1-nu_2*e_nu)*(beta*(1-epsilon*pi_12m)*Psi_12m*rho_1m*B_1+beta*(1-epsilon*pi_24m)*Psi_24m*rho_4m*B_4)
  Lambda_2c <- lambda_2*n_2c*(1-nu_2*e_nu)*(beta*(1-epsilon*pi_12c)*Psi_12c*rho_1c*B_1+beta*(1-epsilon*pi_24c)*Psi_24c*rho_4c*B_4)
  Lambda_3m <- lambda_3*n_3m*(1-nu_3*e_nu)*(beta*(1-epsilon*pi_13m)*Psi_13m*rho_1m*B_1+beta*(1-epsilon*pi_34m)*Psi_34m*rho_4m*B_4)
  Lambda_3c <- lambda_3*n_3c*(1-nu_3*e_nu)*(beta*(1-epsilon*pi_13c)*Psi_13c*rho_1c*B_1+beta*(1-epsilon*pi_34c)*Psi_34c*rho_4c*B_4)
  Lambda_3co <- lambda_3*n_3co*beta*(1-epsilon*pi_34co)*B_4
  Lambda_4m <- lambda_4*(rho_4m/N4)*(1-s*e_s)*(1-nu_4*e_nu)*(beta*(1-epsilon*pi_24m)*Psi_24m*n_2m*N2*B_2+beta*(1-epsilon*pi_34m)*Psi_34m*n_3m*N3*B_3)
  Lambda_4c <- lambda_4*(rho_4c/N4)*(1-s*e_s)*(1-nu_4*e_nu)*(beta*(1-epsilon*pi_24c)*Psi_24c*n_2c*N2*B_2+beta*(1-epsilon*pi_34c)*Psi_34c*n_3c*N3*B_3)
  Lambda_4co <- lambda_4*n_4co*(1-s*e_s)*(1-nu_4*e_nu)*(1-epsilon*pi_34co)*B_3
  

  # Define differential equations
  
  dS1 <- (1-phi_1)*p*theta*N+gamma*S4-(Lambda_1m+Lambda_1c)*S1-(kappa+mu_1)*S1
  dE1 <- (Lambda_1m+Lambda_1c)*S1-gamma*E4-((1-alpha_1)*eta_1+alpha_1*eta_2+kappa+mu_1)*E1
  dI1 <- phi_1*p*theta+((1-alpha_1)*eta_1+alpha_1*eta_2)*E1+gamma*I4-((1-tau_1)*delta_1+tau_1*delta_2+kappa+mu_1)*I1
  dR1 <- ((1-tau_1)*delta_1+tau_1*delta_2)*I1
  
  dS2 <- (1-phi_2)*(1-p)*theta*N+g*S3-(Lambda_2m+Lambda_2c)*S2-(z+mu_2)*S2
  dE2 <- (Lambda_2m+Lambda_2c)*S2+g*E3-((1-alpha_2)*eta_1+alpha_2*eta_2+z+mu_2)*E2
  dI2 <- phi_2*(1-p)*theta+g*I3+((1-alpha_2)*eta_1+alpha_2*eta_2)*E2-((1-tau_2)*delta_1+tau_2*delta_2+z+mu_2)*I2
  dR2 <- ((1-tau_2)*delta_1+tau_2*delta_2)*I2
  
  dS3 <- z*S2-(Lambda_3m+Lambda_3c+Lambda_3co)*S3-(g+mu_3)*S3
  dE3 <- (Lambda_3m+Lambda_3c+Lambda_3co)*S3+z*E2-((1-alpha_3)*eta_1+alpha_3*eta_2+g+mu_3)*E3
  dI3 <- ((1-alpha_3)*eta_1+alpha_3*eta_2)*E3+z*I2-((1-tau_3)*delta_1+tau_3*delta_2+g+mu_3)*I3
  dR3 <- ((1-tau_3)*delta_1+tau_3*delta_2)*I3
  
  dS4 <- kappa*S1-(Lambda_4m+Lambda_4c+Lambda_4co)*S4-(gamma+mu_4)*S4
  dE4 <- (Lambda_4m+Lambda_4c+Lambda_4co)*S4+kappa*E1-((1-alpha_4)*eta_1+alpha_4*eta_2+gamma+mu_4)*E4
  dI4 <- ((1-alpha_4)*eta_1+alpha_4*eta_2)*E4+kappa*I1-((1-tau_4)*delta_1+tau_4*delta_2+gamma+mu_4)*I4
  dR4 <- ((1-tau_4)*delta_1+tau_4*delta_2)*I4
  
  dCost<- (S4*s*health_insur_c)+(S4*nu_4*prep_cost)+(E4*alpha_4*art_cost) +(E3*alpha_3*art_cost)+(E2*alpha_2*art_cost)+(E1*alpha_1*art_cost)+(I4*tau_4*art_cost) +(I3*tau_3*art_cost)+(I2*tau_2*art_cost)+(I1*tau_1*art_cost)
  dDalys <- ((E1+E2+E3+E4)*(1-alpha_1))*daly_preAIDs_offART + (E1+E2+E3+E4)*daly_acute_hiv_offART + ((I1+I2+I3+I4)*(1-tau_1))*daly_AIDs_offART + ((I1+I2+I3+I4)*(tau_1))*daly_preAIDs_offART + (R1+R2+R3+R4)*1 
  donset_inf<- (Lambda_1m+Lambda_1c)*S1 + (Lambda_2m+Lambda_2c)*S2 + (Lambda_3m+Lambda_3c+Lambda_3co)*S3 + (Lambda_4m+Lambda_4c+Lambda_4co)*S4
  res <- list(c(dS1, dE1, dI1, dR1, dS2, dE2, dI2,dR2, dS3, dE3, dI3, dR3, dS4, dE4, dI4, dR4, dCost, dDalys, donset_inf))
  return(res)
}

# ##################################################### x
#Parameters#####
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51, #Proportion of individuals who enter to the model as low risk female 
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009, #mortality rate natural
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9, #efficiency ART
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0.0,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.79,
  tau_4 = 0.417,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)
######
times <- seq(from = 0, to = 50, by = 1)
###################################################### x
#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

######################################################

# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
# --- --- --- ------ --- --- --- RESULTS BELOW  --- ------ --- --- ------ --- --- --- --- --- #
# --- --- --- ------ --- --- --- --- --- --- --- --- --- --- ------ --- --- ------ --- --- ---#
#RESULTS I.a: CE only varying health insurance coverage and hence costs per health insurance cov.
discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######

#UPPER_boundBaseCASE#####
#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.031)
E1_0 <- N1_0*0.031*(1-0.0625)
I1_0 <- N1_0*0.031*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.031)
E2_0 <- N2_0*0.031*(1-0.0625)
I2_0 <- N2_0*0.031*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.031)
E3_0 <- N3_0*0.031*(1-0.0625)
I3_0 <- N3_0*0.031*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.031)
E4_0 <- N4_0*0.031*(1-0.0625)
I4_0 <- N4_0*0.031*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['Infected']<-output['I']+output['E']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0_ub<- output[['totalpop']][51]
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0_ub <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/popbs0_ub)*popbs0
cost_bs0_ub <- ((sum(costs_onset / (1 + discount_rate) ^ time))/popbs0_ub)*popbs0
#
infected_bs0_ub<-((sum(output['Infected']))/popbs0_ub)*popbs0
deaths_bs0_ub<-((output[['R']][51])/popbs0_ub)*popbs0
infected_bs0_ubOns<- ((output[['onset_inf']][51])/popbs0_ub)*popbs0


#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['Infected']<-output25p['I']+output25p['E']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p_ub<- output25p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p_ub <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p_ub)*popbs0
cost_25p_ub <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p_ub)*popbs0
#
totalinf_25p_ub<-((sum(output25p['Infected']))/totalpop25p_ub)*popbs0
totaldeat_25p_ub<-((output25p[['R']][51])/totalpop25p_ub)*popbs0
totalinf_25p_ubOns<- ((output25p[['onset_inf']][51])/totalpop25p_ub)*popbs0
ICER_25p_inf_av_ub= (cost_25p_ub-cost_bs0_ub)/(totalinf_25p_ub-infected_bs0_ub) #cost per infection averted
ICER_25p_dalys_ub= (cost_25p_ub-cost_bs0_ub)/(dalys_bs0_ub-dalys_25p_ub)

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p_ub<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p_ub <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p_ub)*popbs0
cost_50p_ub <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p_ub)*popbs0
#
totalinf_50p_ub<-((sum(output50p['Infected']))/totalpop50p_ub)*popbs0
totaldeat_50p_ub<-((output50p[['R']][51])/totalpop50p_ub)*popbs0
totalinf_50p_ubOns<- ((output50p[['onset_inf']][51])/totalpop50p_ub)*popbs0
ICER_50p_inf_av_ub= (cost_50p_ub-cost_bs0_ub)/(totalinf_50p_ub-infected_bs0_ub) #cost per infection averted
ICER_50p_dalys_ub= (cost_50p_ub-cost_bs0_ub)/(dalys_bs0_ub-dalys_50p_ub)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p_ub<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p_ub <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p_ub)*popbs0
cost_75p_ub <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p_ub)*popbs0
#
totalinf_75p_ub<-((sum(output75p['Infected']))/totalpop75p_ub)*popbs0
totaldeat_75p_ub<-((output75p[['R']][51])/totalpop75p_ub)*popbs0
totalinf_75p_ubOns<- ((output75p[['onset_inf']][51])/totalpop75p_ub)*popbs0
ICER_75p_inf_av_ub= (cost_75p_ub-cost_bs0_ub)/(totalinf_75p_ub-infected_bs0_ub) #cost per infection averted
ICER_75p_dalys_ub= (cost_75p_ub-cost_bs0_ub)/(dalys_bs0_ub-dalys_75p_ub)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p_ub<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p_ub <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p_ub)*popbs0
cost_100p_ub <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p_ub)*popbs0
#
totalinf_100p_ub<-((sum(output100p['Infected']))/totalpop100p_ub)*popbs0
totaldeat_100p_ub<-(output100p[['R']][51]/totalpop100p_ub)*popbs0
totalinf_100p_ubOns<- (output100p[['onset_inf']][51]/totalpop100p_ub)*popbs0
ICER_100p_inf_av_ub= (cost_100p_ub-cost_bs0_ub)/(totalinf_100p_ub-infected_bs0_ub) #cost per infection averted
ICER_100p_dalys_ub= (cost_100p_ub-cost_bs0_ub)/(dalys_bs0_ub-dalys_100p_ub)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost_ub = c(cost_bs0_ub, cost_25p_ub, cost_50p_ub, cost_75p_ub, cost_100p_ub),
  DALYs_ub = c(dalys_bs0_ub, dalys_25p_ub, dalys_50p_ub, dalys_75p_ub, dalys_100p_ub),
  TotalInfected_ub = c(sum(output$Infected), totalinf_25p_ub, totalinf_50p_ub, totalinf_75p_ub, totalinf_100p_ub),
  TotalDeaths_ub = c(output[['R']][51], totaldeat_25p_ub, totaldeat_50p_ub, totaldeat_75p_ub, totaldeat_100p_ub),
  ICER_DALYs_ub = c(NA, ICER_25p_dalys_ub, ICER_50p_dalys_ub, ICER_75p_dalys_ub, ICER_100p_dalys_ub),
  ICER_Inf_Av_ub = c(NA, ICER_25p_inf_av_ub, ICER_50p_inf_av_ub, ICER_75p_inf_av_ub, ICER_100p_inf_av_ub)
)
# Print the summary table
View(summary_table)

#####

#LOWER_boundBaseCASE#####
#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.023)
E1_0 <- N1_0*0.023*(1-0.0625)
I1_0 <- N1_0*0.023*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.023)
E2_0 <- N2_0*0.023*(1-0.0625)
I2_0 <- N2_0*0.023*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.023)
E3_0 <- N3_0*0.023*(1-0.0625)
I3_0 <- N3_0*0.023*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.023)
E4_0 <- N4_0*0.023*(1-0.0625)
I4_0 <- N4_0*0.023*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['Infected']<-output['I']+output['E']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0_lb<- output[['totalpop']][51]
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0_lb <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0_lb <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0_lb<- sum(output['Infected'])
deaths_bs0_lb<- output[['R']][51]
infected_bs0_lbOns<- output[['onset_inf']][51]

#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['Infected']<-output25p['I']+output25p['E']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p_lb<- output25p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p_lb <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p_lb)*popbs0
cost_25p_lb <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p_lb)*popbs0
#
totalinf_25p_lb<-(sum(output25p['Infected'])/totalpop25p_lb)*popbs0
totaldeat_25p_lb<-(output25p[['R']][51]/totalpop25p_lb)*popbs0
totalinf_25p_lbOns<- (output25p[['onset_inf']][51]/totalpop25p_lb)*popbs0
ICER_25p_inf_av_lb= (cost_25p_lb-cost_bs0_lb)/(totalinf_25p_lb-infected_bs0_lb) #cost per infection averted
ICER_25p_dalys_lb= (cost_25p_lb-cost_bs0_lb)/(dalys_bs0_lb-dalys_25p_lb)

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p_lb<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p_lb <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p_lb)*popbs0
cost_50p_lb <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p_lb)*popbs0
#
totalinf_50p_lb<-(sum(output50p['Infected'])/totalpop50p_lb)*popbs0
totaldeat_50p_lb<-(output50p[['R']][51]/totalpop50p_lb)*popbs0
totalinf_50p_lbOns<- (output50p[['onset_inf']][51]/totalpop50p_lb)*popbs0
ICER_50p_inf_av_lb= (cost_50p_lb-cost_bs0_lb)/(totalinf_50p_lb-infected_bs0_lb) #cost per infection averted
ICER_50p_dalys_lb= (cost_50p_lb-cost_bs0_lb)/(dalys_bs0_lb-dalys_50p_lb)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p_lb<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p_lb <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p_lb)*popbs0
cost_75p_lb <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p_lb)*popbs0
#
totalinf_75p_lb<-(sum(output75p['Infected'])/totalpop75p_lb)*popbs0
totaldeat_75p_lb<-(output75p[['R']][51]/totalpop75p_lb)*popbs0
totalinf_75p_lbOns<- (output75p[['onset_inf']][51]/totalpop75p_lb)*popbs0
ICER_75p_inf_av_lb= (cost_75p_lb-cost_bs0_lb)/(totalinf_75p_lb-infected_bs0_lb) #cost per infection averted
ICER_75p_dalys_lb= (cost_75p_lb-cost_bs0_lb)/(dalys_bs0_lb-dalys_75p_lb)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p_lb<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p_lb <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p_lb)*popbs0
cost_100p_lb <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p_lb)*popbs0
#
totalinf_100p_lb<-(sum(output100p['Infected'])/totalpop100p_lb)*popbs0
totaldeat_100p_lb<-(output100p[['R']][51]/totalpop100p_lb)*popbs0
totalinf_100p_lbOns<- (output100p[['onset_inf']][51]/totalpop100p_lb)*popbs0
ICER_100p_inf_av_lb= (cost_100p_lb-cost_bs0_lb)/(totalinf_100p_lb-infected_bs0_lb) #cost per infection averted
ICER_100p_dalys_lb= (cost_100p_lb-cost_bs0_lb)/(dalys_bs0_lb-dalys_100p_lb)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0_lb, cost_25p_lb, cost_50p_lb, cost_75p_lb, cost_100p_lb),
  DALYs = c(dalys_bs0_lb, dalys_25p_lb, dalys_50p_lb, dalys_75p_lb, dalys_100p_lb),
  TotalInfected = c(sum(output$Infected), totalinf_25p_lb, totalinf_50p_lb, totalinf_75p_lb, totalinf_100p_lb),
  TotalDeaths = c(output[['R']][51], totaldeat_25p_lb, totaldeat_50p_lb, totaldeat_75p_lb, totaldeat_100p_lb),
  ICER_DALYs = c(NA, ICER_25p_dalys_lb, ICER_50p_dalys_lb, ICER_75p_dalys_lb, ICER_100p_dalys_lb),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av_lb, ICER_50p_inf_av_lb, ICER_75p_inf_av_lb, ICER_100p_inf_av_lb)
)
# Print the summary table
View(summary_table)

#####

theme_lancet <- function(base_size = 14, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      text = element_text(size = base_size, color = "black"),
      plot.title = element_text(size = base_size * 1.1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, hjust = 0, face = "bold"),
      axis.title = element_text(size = base_size),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      axis.text = element_text(size = base_size * 0.8),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5),
      legend.position = "right",
      legend.text = element_text(size = base_size * 0.8),
      legend.title = element_text(size = base_size, face = "bold"),
      panel.grid.major = element_line(color = "#eaeaea"), #"#eaeaea"
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    )
}
lancet_colors <- c("grey70","#f28e79","#a1c181","#f7bb5f","#b699d7")
#FIGURE1#####
variables_I <- c('infected_bs0', 'totalinf_25p','totalinf_50p',  'totalinf_75p','totalinf_100p')
variables_D<- c('deaths_bs0', 'totaldeat_25p','totaldeat_50p','totaldeat_75p', 'totaldeat_100p')
variables_Im <- data.frame(infected_bs0, totalinf_25p,totalinf_50p, totalinf_75p,totalinf_100p)
variables_Dm<- data.frame(deaths_bs0, totaldeat_25p,totaldeat_50p,totaldeat_75p, totaldeat_100p)
variables_ub_I = data.frame(infected_bs0_ub,  totalinf_25p_ub,totalinf_50p_ub, totalinf_75p_ub,totalinf_100p_ub)
variables_ub_D = data.frame( deaths_bs0_ub, totaldeat_25p_ub, totaldeat_50p_ub, totaldeat_75p_ub, totaldeat_100p_ub)
variables_lb_I = data.frame(infected_bs0_lb, totalinf_25p_lb,  totalinf_50p_lb,  totalinf_75p_lb,  totalinf_100p_lb)
variables_lb_D = data.frame(deaths_bs0_lb,  totaldeat_25p_lb,  totaldeat_50p_lb, totaldeat_75p_lb, totaldeat_100p_lb)

variables_ImOns <- data.frame(infected_bs0Ons, totalinf_25pOns,totalinf_50pOns, totalinf_75pOns,totalinf_100pOns)
variables_ub_IOns = data.frame(infected_bs0_ubOns,  totalinf_25p_ubOns,totalinf_50p_ubOns, totalinf_75p_ubOns,totalinf_100p_ubOns)
variables_lb_IOns = data.frame(infected_bs0_lbOns, totalinf_25p_lbOns,  totalinf_50p_lbOns,  totalinf_75p_lbOns,  totalinf_100p_lbOns)

error1 <- variables_Im - variables_lb_I 
error2 <- variables_Dm -variables_lb_D
error3 <- variables_ImOns -variables_lb_IOns
data1 <- data.frame(variables_Im, error1)
data2 <- data.frame(variables_Dm, error2)
data3 <- data.frame(variables_ImOns, error3)

data1_long <- pivot_longer(variables_Im, cols = everything(), names_to = "Variable", values_to = "Value")
error1_long <- pivot_longer(error1, everything(), names_to = "Variable", values_to = "Upper_Bound")
error1_long$Variable <- gsub("_ub", "", error1_long$Variable)
merged_data1 <- merge(data1_long, error1_long, by = "Variable")

data2_long <- pivot_longer(variables_Dm, cols = everything(), names_to = "Variable", values_to = "Value")
error2_long <- pivot_longer(error2, everything(), names_to = "Variable", values_to = "Upper_Bound")
error2_long$Variable <- gsub("_ub", "", error2_long$Variable)
merged_data2 <- merge(data2_long, error2_long, by = "Variable")

data3_long <- pivot_longer(variables_ImOns, cols = everything(), names_to = "Variable", values_to = "Value")
error3_long <- pivot_longer(error3, everything(), names_to = "Variable", values_to = "Upper_Bound")
error3_long$Variable <- gsub("_ub", "", error3_long$Variable)
merged_data3 <- merge(data3_long, error3_long, by = "Variable")


desired_order <- c("infected_bs0", "totalinf_25p", "totalinf_50p", "totalinf_75p", "totalinf_100p")
desired_order2 <- c("deaths_bs0", "totaldeat_25p", "totaldeat_50p", "totaldeat_75p", "totaldeat_100p")
desired_order3 <- c("infected_bs0Ons", "totalinf_25pOns", "totalinf_50pOns", "totalinf_75pOns", "totalinf_100pOns")

merged_data1$Variable <- factor(merged_data1$Variable, levels = desired_order)
merged_data2$Variable <- factor(merged_data2$Variable, levels = desired_order2)
merged_data3$Variable <- factor(merged_data3$Variable, levels = desired_order3)


ggplot(merged_data1, aes(x = Value, y = reorder(Variable, Value))) +
  geom_point() +  # Add points for each value
  geom_errorbarh(aes(xmin = Value - Upper_Bound, xmax = Value + Upper_Bound), height = 0.2) +  # Add horizontal error bars
  labs(x = "Value", y = "Variable", title = "Values and Error Ranges by Variable") +
  theme_minimal()+
  theme_lancet()


merged_data1$Value <- merged_data1$Value / 50
merged_data1$Upper_Bound <- merged_data1$Upper_Bound / 50
abc1<-ggplot(merged_data1, aes(x = Variable, y = Value, fill = Variable)) +
  geom_col(color = "black")+
  geom_errorbar(aes(ymin = Value - (Upper_Bound), ymax = Value + (Upper_Bound)), width = 0.2) +
  scale_fill_manual(values = lancet_colors) +  # Assigns specific colors to bars # Adjusts y-axis
  labs(x = "", y = "Annual number of HIV infections", title = "(A) Prevalence HIV-infected individuals") +
  theme_lancet()  +
  scale_y_continuous(breaks = seq(0, 700000, 20000)) +
  scale_x_discrete(labels = c("Baseline", "25% coverage", "50% coverage", 
                              "75% coverage", "100% coverage"))+
  theme(axis.text.x = element_text(angle = 360, hjust = 1), plot.title = element_text(hjust = 0),
        legend.position = "none") +
  coord_cartesian(ylim = c(500000, NA))


merged_data2$Value <- merged_data2$Value / 50
merged_data2$Upper_Bound <- merged_data2$Upper_Bound / 50
abc2<-ggplot(merged_data2, aes(x = Variable, y = Value, fill = Variable)) +
  geom_col(color = "black")+
  geom_errorbar(aes(ymin = Value - (Upper_Bound), ymax = Value + (Upper_Bound)), width = 0.2) +
  scale_fill_manual(values = lancet_colors) +  # Assigns specific colors to bars # Adjusts y-axis
  labs(x = "Health coverage level among TFSW", y = "Annual number of HIV deaths", title = "HIV-associated deaths") +
  theme_lancet()  +
  scale_y_continuous(breaks = seq(0, 20500, 500)) +
  scale_x_discrete(labels = c("Baseline", "25% coverage", "50% coverage", 
                              "75% coverage", "100% coverage"))+
  theme(axis.text.x = element_text(angle = 360, hjust = 1), plot.title = element_text(hjust = 0),
        legend.position = "none") +
  coord_cartesian(ylim = c(14500, NA))

merged_data3$Value <- merged_data3$Value/50
merged_data3$Upper_Bound <- merged_data3$Upper_Bound /50
abc3<-ggplot(merged_data3, aes(x = Variable, y = Value, fill = Variable)) +
  geom_col(color = "black")+
  geom_errorbar(aes(ymin = Value - (Upper_Bound), ymax = Value + (Upper_Bound)), width = 0.2) +
  scale_fill_manual(values = lancet_colors) +  # Assigns specific colors to bars # Adjusts y-axis
  labs(x = "", y = "Annual incidence of HIV infections", title = "(B) Incidence HIV-infected individuals") +
  theme_lancet()  +
  scale_y_continuous(breaks = seq(0, 56000, 2000)) +
  scale_x_discrete(labels = c("Baseline", "25% coverage", "50% coverage", 
                              "75% coverage", "100% coverage"))+
  theme(axis.text.x = element_text(angle = 360, hjust = 1), plot.title = element_text(hjust = 0),
        legend.position = "none") +
  coord_cartesian(ylim = c(34000, NA))



abc1 <- abc1 + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
abc3 <- abc3 + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
abc1 <- abc1 + theme(axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 17, face = "bold"), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 15)  )   # Increase y-axis label font size)
abc3 <- abc3 + theme(axis.title.x = element_text(size = 22), axis.title.y = element_text(size = 17,face = "bold"), axis.text.x = element_text(size = 14),  axis.text.y = element_text(size = 15) )
abc1 <- abc1 + theme(axis.text.x = element_text( hjust = 0.5, vjust = 0))
abc3 <- abc3 + theme(axis.text.x = element_text( hjust = 0.5, vjust = 0))

combined_plot <- abc1 + abc3 
combined_plot <- abc1 + abc3 + plot_layout(ncol = 1)
# Or, to make them vertically aligned
combined_plot <- abc1 / abc3
ggsave("plot2.tiff",combined_plot, width = 9, height = 10, dpi = 1000)

ggsave("abc2.tiff",abc2, width = 11, height = 10, dpi = 1000)


#####
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#SENSITITY GLOBAL####
library(lhs)  # For Latin Hypercube Sampling
library(deSolve)  # For ode function
# Adjust the parameter ranges to +/- 10% of their original values
library(lhs)  # Ensure the lhs package is loaded
param_names <- names(parameters)
param_ranges <- sapply(parameters, function(x) c(0.9 * x, 1.1 * x))
n_samples <- 100 # Number of samples
n_parameters <- length(parameters)
# Preallocate a matrix for efficiency
scaled_samples <- matrix(nrow = n_samples, ncol = n_parameters)
for(i in 1:n_parameters) {
  min_val <- param_ranges[1, i]
  max_val <- param_ranges[2, i]
  # Generate scaled samples for each parameter within its range
  scaled_samples[, i] <- runif(n_samples, min = min_val, max = max_val)
}
# Function to convert a single row of scaled_samples into a list of parameters
row_to_params <- function(row) {
  params <- as.list(row)
  names(params) <- param_names
  return(params)
}
# Convert each row into a parameter set
parameter_sets <- apply(scaled_samples, 1, row_to_params)
simulation_results <- vector("list", length(parameter_sets))
final_values_combined_output <- numeric(n_samples)
# Correctly running simulations and accessing their results
for (i in seq_along(parameter_sets)) {
  # Run the simulation with the ith set of parameters
  sim_result <- ode(y = state, times = times, func = HIV_model2, parms = parameter_sets[[i]], method = "rk4")
  # Convert the result into a data frame for easier handling
  sim_data <- as.data.frame(sim_result)
  # Assuming 'I1', 'I2', 'I3', 'I4', 'E1', 'E2', 'E3', 'E4' are columns in sim_data representing different compartments
  # Sum the final values of these compartments at the last time point
  final_values_combined_output[i] <- sum(sim_data[nrow(sim_data), c("I1", "I2", "I3", "I4", "E1", "E2", "E3", "E4")])
  # Store the entire simulation result if you need it for later
  simulation_results[[i]] <- sim_data
}




# Convert parameter values and model output into ranks
param_ranks <- apply(scaled_samples, 2, rank)
output_ranks <- rank(final_values_combined_output)
# Identify parameters with variability
valid_params_indices <- which(apply(param_ranks, 2, function(x) sd(x) != 0))
# Filter out parameters with no variability
filtered_param_ranks <- param_ranks[, valid_params_indices]
# Ensure there's variability in output
if(sd(output_ranks) == 0) {
  stop("Output has no variability; cannot compute PRCC.")
}

# Compute PRCC for filtered parameters
prcc_values <- cor(filtered_param_ranks, output_ranks, method = "spearman")

# Compute p-values for the correlations
p_values <- sapply(1:ncol(filtered_param_ranks), function(i) {
  cor.test(filtered_param_ranks[, i], output_ranks, method = "spearman")$p.value
})

# Create a data frame for easier viewing and manipulation
prcc_results <- data.frame(
  Parameter = names(parameters)[valid_params_indices],
  PRCC = prcc_values,
  P_Value = p_values
)
# Filter results for significant PRCC values (P < 0.05)
significant_prcc_results <- prcc_results[prcc_results$P_Value < 0.05, ]



library(ggplot2)
library(dplyr)
significant_prcc_results_filtered <- significant_prcc_results %>%
  filter(Parameter != "prep_cost")
significant_prcc_results_filtered <- significant_prcc_results_filtered %>%
  filter(Parameter != "mu_2")
parameter_labels <- c(
  n_2c= "Frequency of casual sex partners, low risk",
  beta= "Transmission rate",
  eta_1="Progression rate from HIV to AIDS, no ART",
  tau_1="Proportion of AIDs on ART, low risk",
  epsilon="Efficacy of condom use",
  alpha_2="Proportion of HIV low risk women on ART",
  alpha_1="Proportion of HIV low risk men on ART",
  p= "Low risk females proportion",
  #Psi_12c = "Avg. frequency of casual sex acts",
  #pi_12m="Consistency of comdom use",
  #s= "TFSW with health insurance",
  e_alpha="Efficiency of ART to prevent HIV"
)
gsens <-ggplot(significant_prcc_results, aes(x = reorder(Parameter, PRCC), y = PRCC, fill = PRCC)) +
  theme_minimal() +
  theme_lancet()+
  geom_bar(stat = "identity", color = "black") +  # Add black contour to bars
  coord_flip() +
  scale_fill_gradient2(low = "#a1c181", mid = "#FAFAD2", high = "#f28e79", midpoint = 0) +
  geom_text(aes(label = ifelse(P_Value < 0.001, "p<0.001", paste("p=", sprintf("%.3f", P_Value))), 
                y = ifelse(PRCC > 0, PRCC + 0.12, PRCC - 0.12)), 
            color = "black", size = 5) +
  labs(title = "",
       x = "Parameters",
       y = "PRCC values") +
  theme(axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21, face="bold"),
        plot.title = element_text(size = 14),
        legend.position = c(0.9, 0.2),  # Adjust these values as needed
        legend.background = element_rect(colour = "black", fill = NA, size = 0.8),   # Optional: makes legend background transparent
        legend.key = element_rect(colour = "white", size = 2),
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 13),  # Increase the size of the legend text
        legend.title = element_text(size = 17),  # Increase the size of the legend title
        axis.text.x = element_text( hjust = 1, size = 16),  # Increase size and adjust x-axis labels
        axis.text.y = element_text(size = 16))+
  scale_y_continuous(limits = c(-1, 0.6), breaks = seq(-0.8, 0.6, by = 0.2))+
  scale_x_discrete(labels = parameter_labels) +  # Assuming parameter_labels is defined elsewhere
  geom_hline(yintercept = 0, linetype = "solid", color = "black")

  grid.text("(A)", x = unit(0.0, "npc"), y = unit(0.0, "npc"), just = c("left", "top"), gp = gpar(fontsize = 20, fontface="bold"))

ggsave("gsens.tiff",gsens, width = 12, height = 9, dpi = 1000)




################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################



# # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # #
# MONTECARLO SIMULATION # # #
# # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # ## # # # # # #
#########Defineparametervalues
HIV_model2 <- function(times, state, parms) {
  ## Define variables
  # Low risk female 
  S1 <- state["S1"]
  E1 <- state["E1"]
  I1 <- state["I1"]
  R1 <- state["R1"]
  N1 <- S1 + E1 + I1
  
  # Low risk male
  S2 <- state["S2"]
  E2 <- state["E2"]
  I2 <- state["I2"]
  R2 <- state["R2"]
  N2 <- S2 + E2 + I2 
  
  # Client
  S3 <- state["S3"]
  E3 <- state["E3"]
  I3 <- state["I3"]
  R3 <- state["R3"]
  N3 <- S3 + E3 + I3   
  
  # TFSW
  S4 <- state["S4"]
  E4 <- state["E4"]
  I4 <- state["I4"]
  R4 <- state["R4"]
  N4 <- S4 + E4 + I4  
  
  dCost<- state["Cost"]
  dDalys<-state["Dalys"]
  donset_inf<- state["onset_inf"]
  #N total
  N = N1 + N2 + N3 + N4
  
  # Extract parameters
  phi_1 <- parms[["phi_1"]]
  phi_2 <- parms[["phi_2"]]
  p <- parms[["p"]]
  theta <- parms[["theta"]]
  gamma <- parms[["gamma"]]
  kappa <- parms[["kappa"]]
  mu_1 <- parms[["mu_1"]]
  mu_2 <- parms[["mu_2"]]
  mu_3 <- parms[["mu_3"]]
  mu_4 <- parms[["mu_4"]]
  eta_1 <- parms[["eta_1"]]
  eta_2 <- parms[["eta_2"]]
  delta_1 <- parms[["delta_1"]]
  delta_2 <- parms[["delta_2"]]
  z <- parms[["z"]]
  g <- parms[["g"]]
  beta <- parms[["beta"]]
  epsilon <- parms[["epsilon"]]
  pi_12m <- parms[["pi_12m"]]
  pi_12c <- parms[["pi_12c"]]
  pi_13m <- parms[["pi_13m"]]
  pi_13c <- parms[["pi_13c"]] 
  pi_24m <- parms[["pi_24m"]]
  pi_24c <- parms[["pi_24c"]]
  pi_34m <- parms[["pi_34m"]]
  pi_34c <- parms[["pi_34c"]]
  pi_34co <- parms[["pi_34co"]]
  Psi_12m <- parms[["Psi_12m"]]
  Psi_12c <- parms[["Psi_12c"]]
  Psi_13m <- parms[["Psi_13m"]]
  Psi_13c <- parms[["Psi_13c"]]
  Psi_24m <- parms[["Psi_24m"]]
  Psi_24c <- parms[["Psi_24c"]]
  Psi_34m <- parms[["Psi_34m"]]
  Psi_34c <- parms[["Psi_34c"]]
  nu_1 <- parms[["nu_1"]]
  nu_2 <- parms[["nu_2"]]
  nu_3 <- parms[["nu_3"]]
  nu_4 <- parms[["nu_4"]]
  e_nu <- parms[["e_nu"]]
  alpha_1 <- parms[["alpha_1"]]
  alpha_2 <- parms[["alpha_2"]]
  alpha_3 <- parms[["alpha_3"]]
  alpha_4 <- parms[["alpha_4"]]
  e_alpha <- parms[["e_alpha"]]
  n_1m <- parms[["n_1m"]]
  n_1c <- parms[["n_1c"]]
  n_2m <- parms[["n_2m"]]
  n_2c <- parms[["n_2c"]]
  n_3m <- parms[["n_3m"]]
  n_3c <- parms[["n_3c"]]
  n_3co <- parms[["n_3co"]]
  n_4m <- parms[["n_4m"]]
  n_4c <- parms[["n_4c"]]
  n_4co <- parms[["n_4co"]]
  s <- parms[["s"]]
  e_s <- parms[["e_s"]]
  tau_1 <- parms[["tau_1"]]
  tau_2 <- parms[["tau_2"]]
  tau_3 <- parms[["tau_3"]]
  tau_4 <- parms[["tau_4"]]
  e_tau <- parms[["e_tau"]]
  P1_inf <- parms[["P1_inf"]]
  P2_inf <- parms[["P2_inf"]]
  P3_inf <- parms[["P3_inf"]]
  P4_inf <- parms[["P4_inf"]]
  w <- parms[["w"]]
  #Utilities & costs below
  health_insur_c <-parms[["health_insur_c"]]
  prep_cost <- parms[["prep_cost"]]
  art_cost <- parms[["art_cost"]]
  daly_acute_hiv_offART <- parms[["daly_acute_hiv_offART"]]
  daly_preAIDs_offART <- parms[["daly_preAIDs_offART"]]
  daly_AIDs_offART <- parms[["daly_AIDs_offART"]]
  
  #probablities
  
  rho_1m <- (n_1m*N1)/(n_1m*N1+n_4m*N4)
  rho_1c <- (n_1c*N1)/(n_1c*N1+n_4c*N4)
  rho_2m <- (n_2m*N2)/(n_2m*N2+n_3m*N3)
  rho_2c <- (n_2c*N2)/(n_2c*N2+n_3c*N3)
  rho_3m <- (n_3m*N3)/(n_3m*N2+n_3m*N3)
  rho_3c <- (n_3c*N3)/(n_3c*N2+n_3c*N3)
  rho_4m <- (n_4m*N4)/(n_1m*N1+n_4m*N4)
  rho_4c <- (n_4c*N4)/(n_1c*N1+n_4c*N4)
  
  #Prevalence of HIV
  
  P1_t <- (E1+I1)/N1
  P2_t <- (E2+I2)/N2
  P3_t <- (E3+I3)/N3
  P4_t <- (E4+I4)/N4
  
  #lambda
  
  lambda_1 <- w*exp((P1_t/P1_inf)*log(1/w))
  lambda_2 <- w*exp((P2_t/P2_inf)*log(1/w))
  lambda_3 <- w*exp((P3_t/P3_inf)*log(1/w))
  lambda_4 <- w*exp((P4_t/P1_inf)*log(1/w))
  
  #B_i
  
  B_1 <- ((1-alpha_1*e_alpha)*E1+(1-tau_1*e_tau)*I1)/N1
  B_2 <- ((1-alpha_2*e_alpha)*E2+(1-tau_2*e_tau)*I2)/N2
  B_3 <- ((1-alpha_3*e_alpha)*E3+(1-tau_3*e_tau)*I3)/N3
  B_4 <- ((1-alpha_4*e_alpha)*E4+(1-tau_4*e_tau)*I4)/N4
  
  #Lambdas
  
  Lambda_1m <- lambda_1*(rho_1m/N1)*(1-nu_1*e_nu)*(beta*(1-epsilon*pi_12m)*Psi_12m*n_2m*N2*B_2+beta*(1-epsilon*pi_13m)*Psi_13m*n_3m*N3*B_3)
  Lambda_1c <- lambda_1*(rho_1c/N1)*(1-nu_1*e_nu)*(beta*(1-epsilon*pi_12c)*Psi_12c*n_2c*N2*B_2+beta*(1-epsilon*pi_13c)*Psi_13c*n_3c*N3*B_3)
  Lambda_2m <- lambda_2*n_2m*(1-nu_2*e_nu)*(beta*(1-epsilon*pi_12m)*Psi_12m*rho_1m*B_1+beta*(1-epsilon*pi_24m)*Psi_24m*rho_4m*B_4)
  Lambda_2c <- lambda_2*n_2c*(1-nu_2*e_nu)*(beta*(1-epsilon*pi_12c)*Psi_12c*rho_1c*B_1+beta*(1-epsilon*pi_24c)*Psi_24c*rho_4c*B_4)
  Lambda_3m <- lambda_3*n_3m*(1-nu_3*e_nu)*(beta*(1-epsilon*pi_13m)*Psi_13m*rho_1m*B_1+beta*(1-epsilon*pi_34m)*Psi_34m*rho_4m*B_4)
  Lambda_3c <- lambda_3*n_3c*(1-nu_3*e_nu)*(beta*(1-epsilon*pi_13c)*Psi_13c*rho_1c*B_1+beta*(1-epsilon*pi_34c)*Psi_34c*rho_4c*B_4)
  Lambda_3co <- lambda_3*n_3co*beta*(1-epsilon*pi_34co)*B_4
  Lambda_4m <- lambda_4*(rho_4m/N4)*(1-s*e_s)*(1-nu_4*e_nu)*(beta*(1-epsilon*pi_24m)*Psi_24m*n_2m*N2*B_2+beta*(1-epsilon*pi_34m)*Psi_34m*n_3m*N3*B_3)
  Lambda_4c <- lambda_4*(rho_4c/N4)*(1-s*e_s)*(1-nu_4*e_nu)*(beta*(1-epsilon*pi_24c)*Psi_24c*n_2c*N2*B_2+beta*(1-epsilon*pi_34c)*Psi_34c*n_3c*N3*B_3)
  Lambda_4co <- lambda_4*n_4co*(1-s*e_s)*(1-nu_4*e_nu)*(1-epsilon*pi_34co)*B_3
  
  
  # Define differential equations
  
  dS1 <- (1-phi_1)*p*theta*N+gamma*S4-(Lambda_1m+Lambda_1c)*S1-(kappa+mu_1)*S1
  dE1 <- (Lambda_1m+Lambda_1c)*S1-gamma*E4-((1-alpha_1)*eta_1+alpha_1*eta_2+kappa+mu_1)*E1
  dI1 <- phi_1*p*theta+((1-alpha_1)*eta_1+alpha_1*eta_2)*E1+gamma*I4-((1-tau_1)*delta_1+tau_1*delta_2+kappa+mu_1)*I1
  dR1 <- ((1-tau_1)*delta_1+tau_1*delta_2)*I1
  
  dS2 <- (1-phi_2)*(1-p)*theta*N+g*S3-(Lambda_2m+Lambda_2c)*S2-(z+mu_2)*S2
  dE2 <- (Lambda_2m+Lambda_2c)*S2+g*E3-((1-alpha_2)*eta_1+alpha_2*eta_2+z+mu_2)*E2
  dI2 <- phi_2*(1-p)*theta+g*I3+((1-alpha_2)*eta_1+alpha_2*eta_2)*E2-((1-tau_2)*delta_1+tau_2*delta_2+z+mu_2)*I2
  dR2 <- ((1-tau_2)*delta_1+tau_2*delta_2)*I2
  
  dS3 <- z*S2-(Lambda_3m+Lambda_3c+Lambda_3co)*S3-(g+mu_3)*S3
  dE3 <- (Lambda_3m+Lambda_3c+Lambda_3co)*S3+z*E2-((1-alpha_3)*eta_1+alpha_3*eta_2+g+mu_3)*E3
  dI3 <- ((1-alpha_3)*eta_1+alpha_3*eta_2)*E3+z*I2-((1-tau_3)*delta_1+tau_3*delta_2+g+mu_3)*I3
  dR3 <- ((1-tau_3)*delta_1+tau_3*delta_2)*I3
  
  dS4 <- kappa*S1-(Lambda_4m+Lambda_4c+Lambda_4co)*S4-(gamma+mu_4)*S4
  dE4 <- (Lambda_4m+Lambda_4c+Lambda_4co)*S4+kappa*E1-((1-alpha_4)*eta_1+alpha_4*eta_2+gamma+mu_4)*E4
  dI4 <- ((1-alpha_4)*eta_1+alpha_4*eta_2)*E4+kappa*I1-((1-tau_4)*delta_1+tau_4*delta_2+gamma+mu_4)*I4
  dR4 <- ((1-tau_4)*delta_1+tau_4*delta_2)*I4
  
  dCost<- (S4*s*health_insur_c)+(S3*nu_3*prep_cost)+(E4*alpha_4*art_cost) +(E3*alpha_3*art_cost)+(E2*alpha_2*art_cost)+(E1*alpha_1*art_cost)+(I4*tau_4*art_cost) +(I3*tau_3*art_cost)+(I2*tau_2*art_cost)+(I1*tau_1*art_cost)
  dDalys <- ((E1+E2+E3+E4)*(1-alpha_1))*daly_preAIDs_offART + (E1+E2+E3+E4)*daly_acute_hiv_offART + ((I1+I2+I3+I4)*(1-tau_1))*daly_AIDs_offART + ((I1+I2+I3+I4)*(tau_1))*daly_preAIDs_offART + (R1+R2+R3+R4)*1 + (S1+S2+S3+S4)*0
  donset_inf<- (Lambda_1m+Lambda_1c)*S1 + (Lambda_2m+Lambda_2c)*S2 + (Lambda_3m+Lambda_3c+Lambda_3co)*S3 + (Lambda_4m+Lambda_4c+Lambda_4co)*S4
  res <- list(c(dS1, dE1, dI1, dR1, dS2, dE2, dI2,dR2, dS3, dE3, dI3, dR3, dS4, dE4, dI4, dR4, dCost, dDalys, donset_inf))
  return(res)
}
set.seed(123) # For reproducibility
n_runs <- 1000 # Number of simulations
# Assume 10% of the mean for variability
percent_variability <- 0.1
  #####

# Solve the model 0% coverage
#ParametersDistributionPSA##### 
# Corrected version
parameters_sim <- data.frame(
  phi_1 = rep(parameters["phi_1"],n_runs),
  phi_2 = rep(parameters["phi_2"],n_runs),
  p=rnorm(n_runs, mean = parameters["p"], sd = percent_variability*parameters["p"]/2),
  theta = rnorm(n_runs, mean = parameters["theta"], sd = percent_variability*parameters["theta"]/2),
  gamma= rnorm(n_runs, mean = parameters["gamma"], sd = percent_variability*parameters["gamma"]/2),
  kappa= rnorm(n_runs, mean = parameters["kappa"], sd = percent_variability*parameters["kappa"]/2),
  mu_1= rnorm(n_runs, mean = parameters["mu_1"], sd = percent_variability*parameters["mu_1"]/2),
  mu_2= rnorm(n_runs, mean = parameters["mu_2"], sd = percent_variability*parameters["mu_2"]/2),
  mu_3= rnorm(n_runs, mean = parameters["mu_3"], sd = percent_variability*parameters["mu_3"]/2),
  mu_4= rnorm(n_runs, mean = parameters["mu_4"], sd = percent_variability*parameters["mu_4"]/2),
  eta_1= rnorm(n_runs, mean = parameters["eta_1"], sd = percent_variability*parameters["eta_1"]/2),
  eta_2= rnorm(n_runs, mean = parameters["eta_2"], sd = percent_variability*parameters["eta_2"]/2),
  delta_1= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_1"]/2),
  delta_2= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_2"]/2),
  z= rnorm(n_runs, mean = parameters["z"], sd = percent_variability*parameters["z"]/2),
  g= rnorm(n_runs, mean = parameters["g"], sd = percent_variability*parameters["g"]/2),
  beta= rnorm(n_runs, mean = parameters["beta"], sd = percent_variability*parameters["beta"]/2),
  epsilon= rnorm(n_runs, mean = parameters["epsilon"], sd = percent_variability*parameters["epsilon"]/2),
  pi_12m= rnorm(n_runs, mean = parameters["pi_12m"], sd = percent_variability*parameters["pi_12m"]/2),
  pi_12c= rnorm(n_runs, mean = parameters["pi_12c"], sd = percent_variability*parameters["pi_12c"]/2),
  pi_13m= rnorm(n_runs, mean = parameters["pi_13m"], sd = percent_variability*parameters["pi_13m"]/2),
  pi_13c= rnorm(n_runs, mean = parameters["pi_13c"], sd = percent_variability*parameters["pi_13c"]/2),
  pi_24m= rnorm(n_runs, mean = parameters["pi_24m"], sd = percent_variability*parameters["pi_24m"]/2),
  pi_24c= rnorm(n_runs, mean = parameters["pi_24c"], sd = percent_variability*parameters["pi_24c"]/2),
  pi_34m= rnorm(n_runs, mean = parameters["pi_34m"], sd = percent_variability*parameters["pi_34m"]/2),
  pi_34c= rnorm(n_runs, mean = parameters["pi_34c"], sd = percent_variability*parameters["pi_34c"]/2),
  pi_34co= rnorm(n_runs, mean = parameters["pi_34co"], sd = percent_variability*parameters["pi_34co"]/2),
  nu_1= rep(parameters["nu_1"],n_runs),
  nu_2= rep(parameters["nu_2"],n_runs),
  nu_3= rnorm(n_runs, mean = parameters["nu_3"], sd = percent_variability*parameters["nu_3"]/2),
  nu_4= rnorm(n_runs, mean = parameters["nu_4"], sd = percent_variability*parameters["nu_4"]/2),
  e_nu= rnorm(n_runs, mean = parameters["e_nu"], sd = percent_variability*parameters["e_nu"]/2),
  alpha_1= rnorm(n_runs, mean = parameters["alpha_1"], sd = percent_variability*parameters["alpha_1"]/2),
  alpha_2= rnorm(n_runs, mean = parameters["alpha_2"], sd = percent_variability*parameters["alpha_2"]/2),
  alpha_3= rnorm(n_runs, mean = parameters["alpha_3"], sd = percent_variability*parameters["alpha_3"]/2),
  alpha_4= rnorm(n_runs, mean = parameters["alpha_4"], sd = percent_variability*parameters["alpha_4"]/2),
  e_alpha=  rnorm(n_runs, mean = parameters["e_alpha"], sd = percent_variability*parameters["e_alpha"]/2),
  s=rep(0.0,n_runs),
  e_s= rep(0.0,n_runs),
  tau_1= rnorm(n_runs, mean = parameters["tau_1"], sd = percent_variability*parameters["tau_1"]/2),
  tau_2= rnorm(n_runs, mean = parameters["tau_2"], sd = percent_variability*parameters["tau_2"]/2),
  tau_3= rnorm(n_runs, mean = parameters["tau_3"], sd = percent_variability*parameters["tau_3"]/2),
  tau_4= rnorm(n_runs, mean = parameters["tau_4"], sd = percent_variability*parameters["tau_4"]/2),
  e_tau= rnorm(n_runs, mean = parameters["e_tau"], sd = percent_variability*parameters["e_tau"]/2),
  P1_inf= rnorm(n_runs, mean = parameters["P1_inf"], sd = percent_variability*parameters["P1_inf"]/2),
  P2_inf= rnorm(n_runs, mean = parameters["P2_inf"], sd = percent_variability*parameters["P2_inf"]/2),
  P3_inf= rnorm(n_runs, mean = parameters["P3_inf"], sd = percent_variability*parameters["P3_inf"]/2),
  P4_inf= rnorm(n_runs, mean = parameters["P4_inf"], sd = percent_variability*parameters["P4_inf"]/2),
  daly_acute_hiv_offART= rnorm(n_runs, mean = parameters["daly_acute_hiv_offART"], sd = percent_variability*parameters["daly_acute_hiv_offART"]/2),
  daly_preAIDs_offART= rnorm(n_runs, mean = parameters["daly_preAIDs_offART"], sd = percent_variability*parameters["daly_preAIDs_offART"]/2),
  daly_AIDs_offART= rnorm(n_runs, mean = parameters["daly_AIDs_offART"], sd = percent_variability*parameters["daly_AIDs_offART"]/2),
  Psi_12m= rnorm(n_runs, mean = parameters["Psi_12m"], sd = percent_variability*parameters["Psi_12m"]/2),
  Psi_12c= rnorm(n_runs, mean = parameters["Psi_12c"], sd = parameters["Psi_12c"] * percent_variability/2),
  Psi_13m= rnorm(n_runs, mean = parameters["Psi_13m"], sd = parameters["Psi_13m"] * percent_variability/2),
  Psi_13c= rnorm(n_runs, mean = parameters["Psi_13c"], sd = parameters["Psi_13c"] * percent_variability/2),
  Psi_24m= rnorm(n_runs, mean = parameters["Psi_24m"], sd = parameters["Psi_24m"] * percent_variability/2),
  Psi_24c= rnorm(n_runs, mean = parameters["Psi_24c"], sd = parameters["Psi_24c"] * percent_variability/2),
  Psi_34m= rnorm(n_runs, mean = parameters["Psi_34m"], sd = parameters["Psi_34m"] * percent_variability/2),
  Psi_34c= rnorm(n_runs, mean = parameters["Psi_34c"], sd = parameters["Psi_34c"] * percent_variability/2),
  n_1m= rnorm(n_runs, mean = parameters["n_1m"], sd = parameters["n_1m"] * percent_variability/2),
  n_1c= rnorm(n_runs, mean = parameters["n_1c"], sd = parameters["n_1c"] * percent_variability/2),
  n_2m= rnorm(n_runs, mean = parameters["n_2m"], sd = parameters["n_2m"]*percent_variability/2),
  n_2c= rnorm(n_runs, mean = parameters["n_2c"], sd = parameters["n_2c"] * percent_variability/2),
  n_3m= rnorm(n_runs, mean = parameters["n_3m"], sd = parameters["n_3m"]*percent_variability/2),
  n_3c= rnorm(n_runs, mean = parameters["n_3c"], sd = parameters["n_3c"] * percent_variability/2),
  n_3co= rnorm(n_runs, mean = parameters["n_3co"], sd = parameters["n_3co"] * percent_variability/2),
  n_4m= rnorm(n_runs, mean = parameters["n_4m"], sd = parameters["n_4m"]*percent_variability/2),
  n_4c= rnorm(n_runs, mean = parameters["n_4c"], sd = parameters["n_4c"] * percent_variability/2),
  n_4co= rnorm(n_runs, mean = parameters["n_4co"], sd = parameters["n_4co"] * percent_variability/2),
  w= rnorm(n_runs, mean = parameters["w"], sd = parameters["w"] * percent_variability/2),
  health_insur_c= rnorm(n_runs, mean = parameters["health_insur_c"], sd = parameters["health_insur_c"] * percent_variability/2),
  prep_cost= rnorm(n_runs, mean = parameters["prep_cost"], sd = parameters["prep_cost"] * percent_variability/2),
  art_cost= rnorm(n_runs, mean = parameters["art_cost"], sd = parameters["art_cost"] * percent_variability/2))

#####
results <- matrix(1, nrow = n_runs, ncol = 5)
results<- as.data.frame(results)
for(i in 1:n_runs) {
  # Extract the ith row of parameters as a list to pass to the model
  #BASELINE CONDITIONS#####################################################
  # Define initial conditions
  N_0 <- 27914536*0.55
  N_female_0 <- N_0*0.51
  N_male_0 <- N_0*(1-0.51)
  N1_0 <- N_female_0*(1-0.02)
  N2_0 <- N_male_0*(1-0.14)
  N3_0 <- N_male_0*0.14
  N4_0 <- N_female_0*0.02
  
  S1_0 <- N1_0*(1-0.027)
  E1_0 <- N1_0*0.027*(1-0.0625)
  I1_0 <- N1_0*0.027*0.0625
  R1_0 <- 0
  
  S2_0 <- N2_0*(1-0.027)
  E2_0 <- N2_0*0.027*(1-0.0625)
  I2_0 <- N2_0*0.027*0.0625
  R2_0 <- 0
  
  S3_0 <- N3_0*(1-0.027)
  E3_0 <- N3_0*0.027*(1-0.0625)
  I3_0 <- N3_0*0.027*0.0625
  R3_0 <- 0
  
  S4_0 <- N4_0*(1-0.027)
  E4_0 <- N4_0*0.027*(1-0.0625)
  I4_0 <- N4_0*0.027*0.0625
  R4_0 <- 0
  
  Cost_0<- 0
  Dalys_0<-0
  onset_inf_0<-0
  
  state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
             S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
             S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
             S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)
  
  #####
  param_set <- c(parameters_sim[i, ])
  # Run the model with these parameter s
  outputp0 <- ode(y = state, times = times, func = HIV_model2, parms = param_set)
  outputp0 <- as.data.frame(outputp0)
  outputp0['I']<-outputp0['I1']+outputp0['I2']+outputp0['I3']+outputp0['I4']
  outputp0['E']<-outputp0['E1']+outputp0['E2']+outputp0['E3']+outputp0['E4']
  outputp0['R']<-outputp0['R1']+outputp0['R2']+outputp0['R3']+outputp0['R4'] 
  outputp0['S']<-outputp0['S1']+outputp0['S2']+outputp0['S3']+outputp0['S4']
  outputp0['totalpop']<- outputp0['R']+outputp0['I']+outputp0['E']+outputp0['S']
  totalpopp0<- outputp0[['totalpop']][51]
  dalys_onset <- diff(c(0, outputp0[['Dalys']]))
  costs_onset <- diff(c(0, outputp0[['Cost']]))# Prepend 0 for the first year to match the length
  dalys_p0 <- ((sum(dalys_onset / (1 + discount_rate) ^ time)))
  cost_p0 <- ((sum(costs_onset / (1 + discount_rate) ^ time)))
  # Store the model output in the results list
  results[[i,1]] <- (cost_p0) #incremental costs
  results[[i,2]] <- ((sum(outputp0[['I']]+outputp0[['E']]))) #incremental infections averted
  results[[i,3]] <- (dalys_p0)
  results[[i,4]] <- totalpopp0
}
resultsp0<- results


# Solve the model 25% coverage
#ParametersDistributionPSA##### 
# Corrected version
parameters_sim <- data.frame(
  phi_1 = rep(parameters["phi_1"],n_runs),
  phi_2 = rep(parameters["phi_2"],n_runs),
  p=rnorm(n_runs, mean = parameters["p"], sd = percent_variability*parameters["p"]/2),
  theta = rnorm(n_runs, mean = parameters["theta"], sd = percent_variability*parameters["theta"]/2),
  gamma= rnorm(n_runs, mean = parameters["gamma"], sd = percent_variability*parameters["gamma"]/2),
  kappa= rnorm(n_runs, mean = parameters["kappa"], sd = percent_variability*parameters["kappa"]/2),
  mu_1= rnorm(n_runs, mean = parameters["mu_1"], sd = percent_variability*parameters["mu_1"]/2),
  mu_2= rnorm(n_runs, mean = parameters["mu_2"], sd = percent_variability*parameters["mu_2"]/2),
  mu_3= rnorm(n_runs, mean = parameters["mu_3"], sd = percent_variability*parameters["mu_3"]/2),
  mu_4= rnorm(n_runs, mean = parameters["mu_4"], sd = percent_variability*parameters["mu_4"]/2),
  eta_1= rnorm(n_runs, mean = parameters["eta_1"], sd = percent_variability*parameters["eta_1"]/2),
  eta_2= rnorm(n_runs, mean = parameters["eta_2"], sd = percent_variability*parameters["eta_2"]/2),
  delta_1= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_1"]/2),
  delta_2= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_2"]/2),
  z= rnorm(n_runs, mean = parameters["z"], sd = percent_variability*parameters["z"]/2),
  g= rnorm(n_runs, mean = parameters["g"], sd = percent_variability*parameters["g"]/2),
  beta= rnorm(n_runs, mean = parameters["beta"], sd = percent_variability*parameters["beta"]/2),
  epsilon= rnorm(n_runs, mean = parameters["epsilon"], sd = percent_variability*parameters["epsilon"]/2),
  pi_12m= rnorm(n_runs, mean = parameters["pi_12m"], sd = percent_variability*parameters["pi_12m"]/2),
  pi_12c= rnorm(n_runs, mean = parameters["pi_12c"], sd = percent_variability*parameters["pi_12c"]/2),
  pi_13m= rnorm(n_runs, mean = parameters["pi_13m"], sd = percent_variability*parameters["pi_13m"]/2),
  pi_13c= rnorm(n_runs, mean = parameters["pi_13c"], sd = percent_variability*parameters["pi_13c"]/2),
  pi_24m= rnorm(n_runs, mean = parameters["pi_24m"], sd = percent_variability*parameters["pi_24m"]/2),
  pi_24c= rnorm(n_runs, mean = parameters["pi_24c"], sd = percent_variability*parameters["pi_24c"]/2),
  pi_34m= rnorm(n_runs, mean = parameters["pi_34m"], sd = percent_variability*parameters["pi_34m"]/2),
  pi_34c= rnorm(n_runs, mean = parameters["pi_34c"], sd = percent_variability*parameters["pi_34c"]/2),
  pi_34co= rnorm(n_runs, mean = parameters["pi_34co"], sd = percent_variability*parameters["pi_34co"]/2),
  nu_1= rep(parameters["nu_1"],n_runs),
  nu_2= rep(parameters["nu_2"],n_runs),
  nu_3= rnorm(n_runs, mean = parameters["nu_3"], sd = percent_variability*parameters["nu_3"]/2),
  nu_4= rnorm(n_runs, mean = parameters["nu_4"], sd = percent_variability*parameters["nu_4"]/2),
  e_nu= rnorm(n_runs, mean = parameters["e_nu"], sd = percent_variability*parameters["e_nu"]/2),
  alpha_1= rnorm(n_runs, mean = parameters["alpha_1"], sd = percent_variability*parameters["alpha_1"]/2),
  alpha_2= rnorm(n_runs, mean = parameters["alpha_2"], sd = percent_variability*parameters["alpha_2"]/2),
  alpha_3= rnorm(n_runs, mean = parameters["alpha_3"], sd = percent_variability*parameters["alpha_3"]/2),
  alpha_4= rnorm(n_runs, mean = parameters["alpha_4"], sd = percent_variability*parameters["alpha_4"]/2),
  e_alpha=  rnorm(n_runs, mean = parameters["e_alpha"], sd = percent_variability*parameters["e_alpha"]/2),
  s=rep(0.25,n_runs),
  e_s= rnorm(n_runs, mean = parameters["e_s"], sd = percent_variability*parameters["e_s"]/2),
  tau_1= rnorm(n_runs, mean = parameters["tau_1"], sd = percent_variability*parameters["tau_1"]/2),
  tau_2= rnorm(n_runs, mean = parameters["tau_2"], sd = percent_variability*parameters["tau_2"]/2),
  tau_3= rnorm(n_runs, mean = parameters["tau_3"], sd = percent_variability*parameters["tau_3"]/2),
  tau_4= rnorm(n_runs, mean = parameters["tau_4"], sd = percent_variability*parameters["tau_4"]/2),
  e_tau= rnorm(n_runs, mean = parameters["e_tau"], sd = percent_variability*parameters["e_tau"]/2),
  P1_inf= rnorm(n_runs, mean = parameters["P1_inf"], sd = percent_variability*parameters["P1_inf"]/2),
  P2_inf= rnorm(n_runs, mean = parameters["P2_inf"], sd = percent_variability*parameters["P2_inf"]/2),
  P3_inf= rnorm(n_runs, mean = parameters["P3_inf"], sd = percent_variability*parameters["P3_inf"]/2),
  P4_inf= rnorm(n_runs, mean = parameters["P4_inf"], sd = percent_variability*parameters["P4_inf"]/2),
  daly_acute_hiv_offART= rnorm(n_runs, mean = parameters["daly_acute_hiv_offART"], sd = percent_variability*parameters["daly_acute_hiv_offART"]/2),
  daly_preAIDs_offART= rnorm(n_runs, mean = parameters["daly_preAIDs_offART"], sd = percent_variability*parameters["daly_preAIDs_offART"]/2),
  daly_AIDs_offART= rnorm(n_runs, mean = parameters["daly_AIDs_offART"], sd = percent_variability*parameters["daly_AIDs_offART"]/2),
  Psi_12m= rnorm(n_runs, mean = parameters["Psi_12m"], sd = percent_variability*parameters["Psi_12m"]/2),
  Psi_12c= rnorm(n_runs, mean = parameters["Psi_12c"], sd = parameters["Psi_12c"] * percent_variability/2),
  Psi_13m= rnorm(n_runs, mean = parameters["Psi_13m"], sd = parameters["Psi_13m"] * percent_variability/2),
  Psi_13c= rnorm(n_runs, mean = parameters["Psi_13c"], sd = parameters["Psi_13c"] * percent_variability/2),
  Psi_24m= rnorm(n_runs, mean = parameters["Psi_24m"], sd = parameters["Psi_24m"] * percent_variability/2),
  Psi_24c= rnorm(n_runs, mean = parameters["Psi_24c"], sd = parameters["Psi_24c"] * percent_variability/2),
  Psi_34m= rnorm(n_runs, mean = parameters["Psi_34m"], sd = parameters["Psi_34m"] * percent_variability/2),
  Psi_34c= rnorm(n_runs, mean = parameters["Psi_34c"], sd = parameters["Psi_34c"] * percent_variability/2),
  n_1m= rnorm(n_runs, mean = parameters["n_1m"], sd = parameters["n_1m"] * percent_variability/2),
  n_1c= rnorm(n_runs, mean = parameters["n_1c"], sd = parameters["n_1c"] * percent_variability/2),
  n_2m= rnorm(n_runs, mean = parameters["n_2m"], sd = parameters["n_2m"]*percent_variability/2),
  n_2c= rnorm(n_runs, mean = parameters["n_2c"], sd = parameters["n_2c"] * percent_variability/2),
  n_3m= rnorm(n_runs, mean = parameters["n_3m"], sd = parameters["n_3m"]*percent_variability/2),
  n_3c= rnorm(n_runs, mean = parameters["n_3c"], sd = parameters["n_3c"] * percent_variability/2),
  n_3co= rnorm(n_runs, mean = parameters["n_3co"], sd = parameters["n_3co"] * percent_variability/2),
  n_4m= rnorm(n_runs, mean = parameters["n_4m"], sd = parameters["n_4m"]*percent_variability/2),
  n_4c= rnorm(n_runs, mean = parameters["n_4c"], sd = parameters["n_4c"] * percent_variability/2),
  n_4co= rnorm(n_runs, mean = parameters["n_4co"], sd = parameters["n_4co"] * percent_variability/2),
  w= rnorm(n_runs, mean = parameters["w"], sd = parameters["w"] * percent_variability/2),
  health_insur_c= rnorm(n_runs, mean = parameters["health_insur_c"], sd = parameters["health_insur_c"] * percent_variability/2),
  prep_cost= rnorm(n_runs, mean = parameters["prep_cost"], sd = parameters["prep_cost"] * percent_variability/2),
  art_cost= rnorm(n_runs, mean = parameters["art_cost"], sd = parameters["art_cost"] * percent_variability/2))

#####
results <- matrix(1, nrow = n_runs, ncol = 5)
results<- as.data.frame(results)
for(i in 1:n_runs) {
  # Extract the ith row of parameters as a list to pass to the model
  #BASELINE CONDITIONS#####################################################
  # Define initial conditions
  N_0 <- 27914536*0.55
  N_female_0 <- N_0*0.51
  N_male_0 <- N_0*(1-0.51)
  N1_0 <- N_female_0*(1-0.02)
  N2_0 <- N_male_0*(1-0.14)
  N3_0 <- N_male_0*0.14
  N4_0 <- N_female_0*0.02
  
  S1_0 <- N1_0*(1-0.027)
  E1_0 <- N1_0*0.027*(1-0.0625)
  I1_0 <- N1_0*0.027*0.0625
  R1_0 <- 0
  
  S2_0 <- N2_0*(1-0.027)
  E2_0 <- N2_0*0.027*(1-0.0625)
  I2_0 <- N2_0*0.027*0.0625
  R2_0 <- 0
  
  S3_0 <- N3_0*(1-0.027)
  E3_0 <- N3_0*0.027*(1-0.0625)
  I3_0 <- N3_0*0.027*0.0625
  R3_0 <- 0
  
  S4_0 <- N4_0*(1-0.027)
  E4_0 <- N4_0*0.027*(1-0.0625)
  I4_0 <- N4_0*0.027*0.0625
  R4_0 <- 0
  
  Cost_0<- 0
  Dalys_0<-0
  onset_inf_0<-0
  
  state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
             S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
             S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
             S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)
  
  #####
  param_set <- c(parameters_sim[i, ])
  # Run the model with these parameter s
  output25p <- ode(y = state, times = times, func = HIV_model2, parms = param_set)
  output25p <- as.data.frame(output25p)
  output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
  output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
  output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
  output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
  output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
  totalpop25p<- output25p[['totalpop']][51]
  dalys_onset <- diff(c(0, output25p[['Dalys']]))
  costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
  dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*resultsp0[[i,4]]
  cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*resultsp0[[i,4]]
  # Store the model output in the results list
  results[[i,1]] <- (cost_25p)-(resultsp0[[i,1]]) #incremental costs
  results[[i,2]] <- (resultsp0[[i,2]])-((sum(output25p[['I']]+output25p[['E']])/totalpop25p)*resultsp0[[i,4]]) #incremental infections averted
  results[[i,3]] <- (resultsp0[[i,3]])-(dalys_25p)
}
results$ICERinf=results[,1]/results[,2]
results$ICERdal=results[,1]/results[,3]
results25<- results
count <- sum(results25$V2 < 0 & results25$V1 > 0)
count2 <- sum(results25$V3 < 0 & results25$V1 > 0)


#50% coverage
results <- matrix(1, nrow = n_runs, ncol = 5)
results<- as.data.frame(results)
#ParametersDistributionPSA##### 
# Corrected version
parameters_sim <- data.frame(
  phi_1 = rep(parameters["phi_1"],n_runs),
  phi_2 = rep(parameters["phi_2"],n_runs),
  p=rnorm(n_runs, mean = parameters["p"], sd = percent_variability*parameters["p"]/2),
  theta = rnorm(n_runs, mean = parameters["theta"], sd = percent_variability*parameters["theta"]/2),
  gamma= rnorm(n_runs, mean = parameters["gamma"], sd = percent_variability*parameters["gamma"]/2),
  kappa= rnorm(n_runs, mean = parameters["kappa"], sd = percent_variability*parameters["kappa"]/2),
  mu_1= rnorm(n_runs, mean = parameters["mu_1"], sd = percent_variability*parameters["mu_1"]/2),
  mu_2= rnorm(n_runs, mean = parameters["mu_2"], sd = percent_variability*parameters["mu_2"]/2),
  mu_3= rnorm(n_runs, mean = parameters["mu_3"], sd = percent_variability*parameters["mu_3"]/2),
  mu_4= rnorm(n_runs, mean = parameters["mu_4"], sd = percent_variability*parameters["mu_4"]/2),
  eta_1= rnorm(n_runs, mean = parameters["eta_1"], sd = percent_variability*parameters["eta_1"]/2),
  eta_2= rnorm(n_runs, mean = parameters["eta_2"], sd = percent_variability*parameters["eta_2"]/2),
  delta_1= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_1"]/2),
  delta_2= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_2"]/2),
  z= rnorm(n_runs, mean = parameters["z"], sd = percent_variability*parameters["z"]/2),
  g= rnorm(n_runs, mean = parameters["g"], sd = percent_variability*parameters["g"]/2),
  beta= rnorm(n_runs, mean = parameters["beta"], sd = percent_variability*parameters["beta"]/2),
  epsilon= rnorm(n_runs, mean = parameters["epsilon"], sd = percent_variability*parameters["epsilon"]/2),
  pi_12m= rnorm(n_runs, mean = parameters["pi_12m"], sd = percent_variability*parameters["pi_12m"]/2),
  pi_12c= rnorm(n_runs, mean = parameters["pi_12c"], sd = percent_variability*parameters["pi_12c"]/2),
  pi_13m= rnorm(n_runs, mean = parameters["pi_13m"], sd = percent_variability*parameters["pi_13m"]/2),
  pi_13c= rnorm(n_runs, mean = parameters["pi_13c"], sd = percent_variability*parameters["pi_13c"]/2),
  pi_24m= rnorm(n_runs, mean = parameters["pi_24m"], sd = percent_variability*parameters["pi_24m"]/2),
  pi_24c= rnorm(n_runs, mean = parameters["pi_24c"], sd = percent_variability*parameters["pi_24c"]/2),
  pi_34m= rnorm(n_runs, mean = parameters["pi_34m"], sd = percent_variability*parameters["pi_34m"]/2),
  pi_34c= rnorm(n_runs, mean = parameters["pi_34c"], sd = percent_variability*parameters["pi_34c"]/2),
  pi_34co= rnorm(n_runs, mean = parameters["pi_34co"], sd = percent_variability*parameters["pi_34co"]/2),
  nu_1= rep(parameters["nu_1"],n_runs),
  nu_2= rep(parameters["nu_2"],n_runs),
  nu_3= rnorm(n_runs, mean = parameters["nu_3"], sd = percent_variability*parameters["nu_3"]/2),
  nu_4= rnorm(n_runs, mean = parameters["nu_4"], sd = percent_variability*parameters["nu_4"]/2),
  e_nu= rnorm(n_runs, mean = parameters["e_nu"], sd = percent_variability*parameters["e_nu"]/2),
  alpha_1= rnorm(n_runs, mean = parameters["alpha_1"], sd = percent_variability*parameters["alpha_1"]/2),
  alpha_2= rnorm(n_runs, mean = parameters["alpha_2"], sd = percent_variability*parameters["alpha_2"]/2),
  alpha_3= rnorm(n_runs, mean = parameters["alpha_3"], sd = percent_variability*parameters["alpha_3"]/2),
  alpha_4= rnorm(n_runs, mean = parameters["alpha_4"], sd = percent_variability*parameters["alpha_4"]/2),
  e_alpha=  rnorm(n_runs, mean = parameters["e_alpha"], sd = percent_variability*parameters["e_alpha"]/2),
  s=rep(0.5,n_runs),
  e_s= rnorm(n_runs, mean = parameters["e_s"], sd = percent_variability*parameters["e_s"]/2),
  tau_1= rnorm(n_runs, mean = parameters["tau_1"], sd = percent_variability*parameters["tau_1"]/2),
  tau_2= rnorm(n_runs, mean = parameters["tau_2"], sd = percent_variability*parameters["tau_2"]/2),
  tau_3= rnorm(n_runs, mean = parameters["tau_3"], sd = percent_variability*parameters["tau_3"]/2),
  tau_4= rnorm(n_runs, mean = parameters["tau_4"], sd = percent_variability*parameters["tau_4"]/2),
  e_tau= rnorm(n_runs, mean = parameters["e_tau"], sd = percent_variability*parameters["e_tau"]/2),
  P1_inf= rnorm(n_runs, mean = parameters["P1_inf"], sd = percent_variability*parameters["P1_inf"]/2),
  P2_inf= rnorm(n_runs, mean = parameters["P2_inf"], sd = percent_variability*parameters["P2_inf"]/2),
  P3_inf= rnorm(n_runs, mean = parameters["P3_inf"], sd = percent_variability*parameters["P3_inf"]/2),
  P4_inf= rnorm(n_runs, mean = parameters["P4_inf"], sd = percent_variability*parameters["P4_inf"]/2),
  daly_acute_hiv_offART= rnorm(n_runs, mean = parameters["daly_acute_hiv_offART"], sd = percent_variability*parameters["daly_acute_hiv_offART"]/2),
  daly_preAIDs_offART= rnorm(n_runs, mean = parameters["daly_preAIDs_offART"], sd = percent_variability*parameters["daly_preAIDs_offART"]/2),
  daly_AIDs_offART= rnorm(n_runs, mean = parameters["daly_AIDs_offART"], sd = percent_variability*parameters["daly_AIDs_offART"]/2),
  Psi_12m= rnorm(n_runs, mean = parameters["Psi_12m"], sd = percent_variability*parameters["Psi_12m"]/2),
  Psi_12c= rnorm(n_runs, mean = parameters["Psi_12c"], sd = parameters["Psi_12c"] * percent_variability/2),
  Psi_13m= rnorm(n_runs, mean = parameters["Psi_13m"], sd = parameters["Psi_13m"] * percent_variability/2),
  Psi_13c= rnorm(n_runs, mean = parameters["Psi_13c"], sd = parameters["Psi_13c"] * percent_variability/2),
  Psi_24m= rnorm(n_runs, mean = parameters["Psi_24m"], sd = parameters["Psi_24m"] * percent_variability/2),
  Psi_24c= rnorm(n_runs, mean = parameters["Psi_24c"], sd = parameters["Psi_24c"] * percent_variability/2),
  Psi_34m= rnorm(n_runs, mean = parameters["Psi_34m"], sd = parameters["Psi_34m"] * percent_variability/2),
  Psi_34c= rnorm(n_runs, mean = parameters["Psi_34c"], sd = parameters["Psi_34c"] * percent_variability/2),
  n_1m= rnorm(n_runs, mean = parameters["n_1m"], sd = parameters["n_1m"] * percent_variability/2),
  n_1c= rnorm(n_runs, mean = parameters["n_1c"], sd = parameters["n_1c"] * percent_variability/2),
  n_2m= rnorm(n_runs, mean = parameters["n_2m"], sd = parameters["n_2m"]*percent_variability/2),
  n_2c= rnorm(n_runs, mean = parameters["n_2c"], sd = parameters["n_2c"] * percent_variability/2),
  n_3m= rnorm(n_runs, mean = parameters["n_3m"], sd = parameters["n_3m"]*percent_variability/2),
  n_3c= rnorm(n_runs, mean = parameters["n_3c"], sd = parameters["n_3c"] * percent_variability/2),
  n_3co= rnorm(n_runs, mean = parameters["n_3co"], sd = parameters["n_3co"] * percent_variability/2),
  n_4m= rnorm(n_runs, mean = parameters["n_4m"], sd = parameters["n_4m"]*percent_variability/2),
  n_4c= rnorm(n_runs, mean = parameters["n_4c"], sd = parameters["n_4c"] * percent_variability/2),
  n_4co= rnorm(n_runs, mean = parameters["n_4co"], sd = parameters["n_4co"] * percent_variability/2),
  w= rnorm(n_runs, mean = parameters["w"], sd = parameters["w"] * percent_variability/2),
  health_insur_c= rnorm(n_runs, mean = parameters["health_insur_c"], sd = parameters["health_insur_c"] * percent_variability/2),
  prep_cost= rnorm(n_runs, mean = parameters["prep_cost"], sd = parameters["prep_cost"] * percent_variability/2),
  art_cost= rnorm(n_runs, mean = parameters["art_cost"], sd = parameters["art_cost"] * percent_variability/2))
######

for(i in 1:n_runs) {
  # Extract the ith row of parameters as a list to pass to the model
  #BASELINE CONDITIONS#####################################################
  # Define initial conditions
  N_0 <- 27914536*0.55
  N_female_0 <- N_0*0.51
  N_male_0 <- N_0*(1-0.51)
  N1_0 <- N_female_0*(1-0.02)
  N2_0 <- N_male_0*(1-0.14)
  N3_0 <- N_male_0*0.14
  N4_0 <- N_female_0*0.02
  
  S1_0 <- N1_0*(1-0.027)
  E1_0 <- N1_0*0.027*(1-0.0625)
  I1_0 <- N1_0*0.027*0.0625
  R1_0 <- 0
  
  S2_0 <- N2_0*(1-0.027)
  E2_0 <- N2_0*0.027*(1-0.0625)
  I2_0 <- N2_0*0.027*0.0625
  R2_0 <- 0
  
  S3_0 <- N3_0*(1-0.027)
  E3_0 <- N3_0*0.027*(1-0.0625)
  I3_0 <- N3_0*0.027*0.0625
  R3_0 <- 0
  
  S4_0 <- N4_0*(1-0.027)
  E4_0 <- N4_0*0.027*(1-0.0625)
  I4_0 <- N4_0*0.027*0.0625
  R4_0 <- 0
  
  Cost_0<- 0
  Dalys_0<-0
  onset_inf_0<-0
  
  state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
             S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
             S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
             S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)
  
  #####
  param_set <- c(parameters_sim[i, ])
  # Run the model with these parameters
  output50p <- ode(y = state, times = times, func = HIV_model2, parms = param_set)
  output50p <- as.data.frame(output50p)
  output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
  output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
  output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
  output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
  output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
  totalpop50p<- output50p[['totalpop']][51]
  dalys_onset <- diff(c(0, output50p[['Dalys']]))
  costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
  dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*resultsp0[[i,4]]
  cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*resultsp0[[i,4]]
  # Store the model output in the results list
  results[[i,1]] <- (cost_50p) -(resultsp0[[i,1]]) #incremental costs
  results[[i,2]] <- (resultsp0[[i,2]])-((sum(output50p[['I']]+output50p[['E']])/totalpop50p)*resultsp0[[i,4]]) #incremental infections averted
  results[[i,3]] <- (resultsp0[[i,3]])-(dalys_50p)  #incremental DALYs
}
results$ICERinf=results[,1]/results[,2]
results$ICERdal=results[,1]/results[,3]
results50<- results
count <- sum(results50$V2 < 0 & results50$V1 > 0)
count2 <- sum(results50$V3 < 0 & results50$V1 > 0)

#75% coverage
results <- matrix(1, nrow = n_runs, ncol = 5)
results<- as.data.frame(results)
#ParametersDistributionPSA##### 
# Corrected version
parameters_sim <- data.frame(
  phi_1 = rep(parameters["phi_1"],n_runs),
  phi_2 = rep(parameters["phi_2"],n_runs),
  p=rnorm(n_runs, mean = parameters["p"], sd = percent_variability*parameters["p"]/2),
  theta = rnorm(n_runs, mean = parameters["theta"], sd = percent_variability*parameters["theta"]/2),
  gamma= rnorm(n_runs, mean = parameters["gamma"], sd = percent_variability*parameters["gamma"]/2),
  kappa= rnorm(n_runs, mean = parameters["kappa"], sd = percent_variability*parameters["kappa"]/2),
  mu_1= rnorm(n_runs, mean = parameters["mu_1"], sd = percent_variability*parameters["mu_1"]/2),
  mu_2= rnorm(n_runs, mean = parameters["mu_2"], sd = percent_variability*parameters["mu_2"]/2),
  mu_3= rnorm(n_runs, mean = parameters["mu_3"], sd = percent_variability*parameters["mu_3"]/2),
  mu_4= rnorm(n_runs, mean = parameters["mu_4"], sd = percent_variability*parameters["mu_4"]/2),
  eta_1= rnorm(n_runs, mean = parameters["eta_1"], sd = percent_variability*parameters["eta_1"]/2),
  eta_2= rnorm(n_runs, mean = parameters["eta_2"], sd = percent_variability*parameters["eta_2"]/2),
  delta_1= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_1"]/2),
  delta_2= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_2"]/2),
  z= rnorm(n_runs, mean = parameters["z"], sd = percent_variability*parameters["z"]/2),
  g= rnorm(n_runs, mean = parameters["g"], sd = percent_variability*parameters["g"]/2),
  beta= rnorm(n_runs, mean = parameters["beta"], sd = percent_variability*parameters["beta"]/2),
  epsilon= rnorm(n_runs, mean = parameters["epsilon"], sd = percent_variability*parameters["epsilon"]/2),
  pi_12m= rnorm(n_runs, mean = parameters["pi_12m"], sd = percent_variability*parameters["pi_12m"]/2),
  pi_12c= rnorm(n_runs, mean = parameters["pi_12c"], sd = percent_variability*parameters["pi_12c"]/2),
  pi_13m= rnorm(n_runs, mean = parameters["pi_13m"], sd = percent_variability*parameters["pi_13m"]/2),
  pi_13c= rnorm(n_runs, mean = parameters["pi_13c"], sd = percent_variability*parameters["pi_13c"]/2),
  pi_24m= rnorm(n_runs, mean = parameters["pi_24m"], sd = percent_variability*parameters["pi_24m"]/2),
  pi_24c= rnorm(n_runs, mean = parameters["pi_24c"], sd = percent_variability*parameters["pi_24c"]/2),
  pi_34m= rnorm(n_runs, mean = parameters["pi_34m"], sd = percent_variability*parameters["pi_34m"]/2),
  pi_34c= rnorm(n_runs, mean = parameters["pi_34c"], sd = percent_variability*parameters["pi_34c"]/2),
  pi_34co= rnorm(n_runs, mean = parameters["pi_34co"], sd = percent_variability*parameters["pi_34co"]/2),
  nu_1= rep(parameters["nu_1"],n_runs),
  nu_2= rep(parameters["nu_2"],n_runs),
  nu_3= rnorm(n_runs, mean = parameters["nu_3"], sd = percent_variability*parameters["nu_3"]/2),
  nu_4= rnorm(n_runs, mean = parameters["nu_4"], sd = percent_variability*parameters["nu_4"]/2),
  e_nu= rnorm(n_runs, mean = parameters["e_nu"], sd = percent_variability*parameters["e_nu"]/2),
  alpha_1= rnorm(n_runs, mean = parameters["alpha_1"], sd = percent_variability*parameters["alpha_1"]/2),
  alpha_2= rnorm(n_runs, mean = parameters["alpha_2"], sd = percent_variability*parameters["alpha_2"]/2),
  alpha_3= rnorm(n_runs, mean = parameters["alpha_3"], sd = percent_variability*parameters["alpha_3"]/2),
  alpha_4= rnorm(n_runs, mean = parameters["alpha_4"], sd = percent_variability*parameters["alpha_4"]/2),
  e_alpha=  rnorm(n_runs, mean = parameters["e_alpha"], sd = percent_variability*parameters["e_alpha"]/2),
  s=rep(0.75,n_runs),
  e_s= rnorm(n_runs, mean = parameters["e_s"], sd = percent_variability*parameters["e_s"]/2),
  tau_1= rnorm(n_runs, mean = parameters["tau_1"], sd = percent_variability*parameters["tau_1"]/2),
  tau_2= rnorm(n_runs, mean = parameters["tau_2"], sd = percent_variability*parameters["tau_2"]/2),
  tau_3= rnorm(n_runs, mean = parameters["tau_3"], sd = percent_variability*parameters["tau_3"]/2),
  tau_4= rnorm(n_runs, mean = parameters["tau_4"], sd = percent_variability*parameters["tau_4"]/2),
  e_tau= rnorm(n_runs, mean = parameters["e_tau"], sd = percent_variability*parameters["e_tau"]/2),
  P1_inf= rnorm(n_runs, mean = parameters["P1_inf"], sd = percent_variability*parameters["P1_inf"]/2),
  P2_inf= rnorm(n_runs, mean = parameters["P2_inf"], sd = percent_variability*parameters["P2_inf"]/2),
  P3_inf= rnorm(n_runs, mean = parameters["P3_inf"], sd = percent_variability*parameters["P3_inf"]/2),
  P4_inf= rnorm(n_runs, mean = parameters["P4_inf"], sd = percent_variability*parameters["P4_inf"]/2),
  daly_acute_hiv_offART= rnorm(n_runs, mean = parameters["daly_acute_hiv_offART"], sd = percent_variability*parameters["daly_acute_hiv_offART"]/2),
  daly_preAIDs_offART= rnorm(n_runs, mean = parameters["daly_preAIDs_offART"], sd = percent_variability*parameters["daly_preAIDs_offART"]/2),
  daly_AIDs_offART= rnorm(n_runs, mean = parameters["daly_AIDs_offART"], sd = percent_variability*parameters["daly_AIDs_offART"]/2),
  Psi_12m= rnorm(n_runs, mean = parameters["Psi_12m"], sd = percent_variability*parameters["Psi_12m"]/2),
  Psi_12c= rnorm(n_runs, mean = parameters["Psi_12c"], sd = parameters["Psi_12c"] * percent_variability/2),
  Psi_13m= rnorm(n_runs, mean = parameters["Psi_13m"], sd = parameters["Psi_13m"] * percent_variability/2),
  Psi_13c= rnorm(n_runs, mean = parameters["Psi_13c"], sd = parameters["Psi_13c"] * percent_variability/2),
  Psi_24m= rnorm(n_runs, mean = parameters["Psi_24m"], sd = parameters["Psi_24m"] * percent_variability/2),
  Psi_24c= rnorm(n_runs, mean = parameters["Psi_24c"], sd = parameters["Psi_24c"] * percent_variability/2),
  Psi_34m= rnorm(n_runs, mean = parameters["Psi_34m"], sd = parameters["Psi_34m"] * percent_variability/2),
  Psi_34c= rnorm(n_runs, mean = parameters["Psi_34c"], sd = parameters["Psi_34c"] * percent_variability/2),
  n_1m= rnorm(n_runs, mean = parameters["n_1m"], sd = parameters["n_1m"] * percent_variability/2),
  n_1c= rnorm(n_runs, mean = parameters["n_1c"], sd = parameters["n_1c"] * percent_variability/2),
  n_2m= rnorm(n_runs, mean = parameters["n_2m"], sd = parameters["n_2m"]*percent_variability/2),
  n_2c= rnorm(n_runs, mean = parameters["n_2c"], sd = parameters["n_2c"] * percent_variability/2),
  n_3m= rnorm(n_runs, mean = parameters["n_3m"], sd = parameters["n_3m"]*percent_variability/2),
  n_3c= rnorm(n_runs, mean = parameters["n_3c"], sd = parameters["n_3c"] * percent_variability/2),
  n_3co= rnorm(n_runs, mean = parameters["n_3co"], sd = parameters["n_3co"] * percent_variability/2),
  n_4m= rnorm(n_runs, mean = parameters["n_4m"], sd = parameters["n_4m"]*percent_variability/2),
  n_4c= rnorm(n_runs, mean = parameters["n_4c"], sd = parameters["n_4c"] * percent_variability/2),
  n_4co= rnorm(n_runs, mean = parameters["n_4co"], sd = parameters["n_4co"] * percent_variability/2),
  w= rnorm(n_runs, mean = parameters["w"], sd = parameters["w"] * percent_variability/2),
  health_insur_c= rnorm(n_runs, mean = parameters["health_insur_c"], sd = parameters["health_insur_c"] * percent_variability/2),
  prep_cost= rnorm(n_runs, mean = parameters["prep_cost"], sd = parameters["prep_cost"] * percent_variability/2),
  art_cost= rnorm(n_runs, mean = parameters["art_cost"], sd = parameters["art_cost"] * percent_variability/2))


######
for(i in 1:n_runs) {
  # Extract the ith row of parameters as a list to pass to the model
  #BASELINE CONDITIONS#####################################################
  # Define initial conditions
  N_0 <- 27914536*0.55
  N_female_0 <- N_0*0.51
  N_male_0 <- N_0*(1-0.51)
  N1_0 <- N_female_0*(1-0.02)
  N2_0 <- N_male_0*(1-0.14)
  N3_0 <- N_male_0*0.14
  N4_0 <- N_female_0*0.02
  
  S1_0 <- N1_0*(1-0.027)
  E1_0 <- N1_0*0.027*(1-0.0625)
  I1_0 <- N1_0*0.027*0.0625
  R1_0 <- 0
  
  S2_0 <- N2_0*(1-0.027)
  E2_0 <- N2_0*0.027*(1-0.0625)
  I2_0 <- N2_0*0.027*0.0625
  R2_0 <- 0
  
  S3_0 <- N3_0*(1-0.027)
  E3_0 <- N3_0*0.027*(1-0.0625)
  I3_0 <- N3_0*0.027*0.0625
  R3_0 <- 0
  
  S4_0 <- N4_0*(1-0.027)
  E4_0 <- N4_0*0.027*(1-0.0625)
  I4_0 <- N4_0*0.027*0.0625
  R4_0 <- 0
  
  Cost_0<- 0
  Dalys_0<-0
  onset_inf_0<-0
  
  state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
             S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
             S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
             S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)
  
  #####
  param_set <- c(parameters_sim[i, ])
  # Run the model with these parameters
  output75p <- ode(y = state, times = times, func = HIV_model2, parms = param_set)
  output75p <- as.data.frame(output75p)
  output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
  output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
  output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
  output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
  output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
  totalpop75p<- output75p[['totalpop']][51]
  dalys_onset <- diff(c(0, output75p[['Dalys']]))
  costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
  dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*resultsp0[[i,4]]
  cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*resultsp0[[i,4]]
  # Store the model output in the results list
  results[[i,1]] <- (cost_75p) -(resultsp0[[i,1]]) #incremental costs
  results[[i,2]] <- (resultsp0[[i,2]])-((sum(output75p[['I']]+output75p[['E']])/totalpop75p)*resultsp0[[i,4]]) #incremental infections averted
  results[[i,3]] <- (resultsp0[[i,3]])-(dalys_75p)  #incremental DALYs
}
results$ICERinf=results[,1]/results[,2]
results$ICERdal=results[,1]/results[,3]
results75<- results
count <- sum(results75$V2 < 0 & results75$V1 > 0)
count2 <- sum(results75$V3 < 0 & results75$V1 > 0)

#100 coverage
results <- matrix(1, nrow = n_runs, ncol = 5)
results<- as.data.frame(results)
#ParametersDistributionPSA##### 
# Corrected version
parameters_sim <- data.frame(
  phi_1 = rep(parameters["phi_1"],n_runs),
  phi_2 = rep(parameters["phi_2"],n_runs),
  p=rnorm(n_runs, mean = parameters["p"], sd = percent_variability*parameters["p"]/2),
  theta = rnorm(n_runs, mean = parameters["theta"], sd = percent_variability*parameters["theta"]/2),
  gamma= rnorm(n_runs, mean = parameters["gamma"], sd = percent_variability*parameters["gamma"]/2),
  kappa= rnorm(n_runs, mean = parameters["kappa"], sd = percent_variability*parameters["kappa"]/2),
  mu_1= rnorm(n_runs, mean = parameters["mu_1"], sd = percent_variability*parameters["mu_1"]/2),
  mu_2= rnorm(n_runs, mean = parameters["mu_2"], sd = percent_variability*parameters["mu_2"]/2),
  mu_3= rnorm(n_runs, mean = parameters["mu_3"], sd = percent_variability*parameters["mu_3"]/2),
  mu_4= rnorm(n_runs, mean = parameters["mu_4"], sd = percent_variability*parameters["mu_4"]/2),
  eta_1= rnorm(n_runs, mean = parameters["eta_1"], sd = percent_variability*parameters["eta_1"]/2),
  eta_2= rnorm(n_runs, mean = parameters["eta_2"], sd = percent_variability*parameters["eta_2"]/2),
  delta_1= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_1"]/2),
  delta_2= rnorm(n_runs, mean = parameters["delta_1"], sd = percent_variability*parameters["delta_2"]/2),
  z= rnorm(n_runs, mean = parameters["z"], sd = percent_variability*parameters["z"]/2),
  g= rnorm(n_runs, mean = parameters["g"], sd = percent_variability*parameters["g"]/2),
  beta= rnorm(n_runs, mean = parameters["beta"], sd = percent_variability*parameters["beta"]/2),
  epsilon= rnorm(n_runs, mean = parameters["epsilon"], sd = percent_variability*parameters["epsilon"]/2),
  pi_12m= rnorm(n_runs, mean = parameters["pi_12m"], sd = percent_variability*parameters["pi_12m"]/2),
  pi_12c= rnorm(n_runs, mean = parameters["pi_12c"], sd = percent_variability*parameters["pi_12c"]/2),
  pi_13m= rnorm(n_runs, mean = parameters["pi_13m"], sd = percent_variability*parameters["pi_13m"]/2),
  pi_13c= rnorm(n_runs, mean = parameters["pi_13c"], sd = percent_variability*parameters["pi_13c"]/2),
  pi_24m= rnorm(n_runs, mean = parameters["pi_24m"], sd = percent_variability*parameters["pi_24m"]/2),
  pi_24c= rnorm(n_runs, mean = parameters["pi_24c"], sd = percent_variability*parameters["pi_24c"]/2),
  pi_34m= rnorm(n_runs, mean = parameters["pi_34m"], sd = percent_variability*parameters["pi_34m"]/2),
  pi_34c= rnorm(n_runs, mean = parameters["pi_34c"], sd = percent_variability*parameters["pi_34c"]/2),
  pi_34co= rnorm(n_runs, mean = parameters["pi_34co"], sd = percent_variability*parameters["pi_34co"]/2),
  nu_1= rep(parameters["nu_1"],n_runs),
  nu_2= rep(parameters["nu_2"],n_runs),
  nu_3= rnorm(n_runs, mean = parameters["nu_3"], sd = percent_variability*parameters["nu_3"]/2),
  nu_4= rnorm(n_runs, mean = parameters["nu_4"], sd = percent_variability*parameters["nu_4"]/2),
  e_nu= rnorm(n_runs, mean = parameters["e_nu"], sd = percent_variability*parameters["e_nu"]/2),
  alpha_1= rnorm(n_runs, mean = parameters["alpha_1"], sd = percent_variability*parameters["alpha_1"]/2),
  alpha_2= rnorm(n_runs, mean = parameters["alpha_2"], sd = percent_variability*parameters["alpha_2"]/2),
  alpha_3= rnorm(n_runs, mean = parameters["alpha_3"], sd = percent_variability*parameters["alpha_3"]/2),
  alpha_4= rnorm(n_runs, mean = parameters["alpha_4"], sd = percent_variability*parameters["alpha_4"]/2),
  e_alpha=  rnorm(n_runs, mean = parameters["e_alpha"], sd = percent_variability*parameters["e_alpha"]/2),
  s=rep(1,n_runs),
  e_s= rnorm(n_runs, mean = parameters["e_s"], sd = percent_variability*parameters["e_s"]/2),
  tau_1= rnorm(n_runs, mean = parameters["tau_1"], sd = percent_variability*parameters["tau_1"]/2),
  tau_2= rnorm(n_runs, mean = parameters["tau_2"], sd = percent_variability*parameters["tau_2"]/2),
  tau_3= rnorm(n_runs, mean = parameters["tau_3"], sd = percent_variability*parameters["tau_3"]/2),
  tau_4= rnorm(n_runs, mean = parameters["tau_4"], sd = percent_variability*parameters["tau_4"]/2),
  e_tau= rnorm(n_runs, mean = parameters["e_tau"], sd = percent_variability*parameters["e_tau"]/2),
  P1_inf= rnorm(n_runs, mean = parameters["P1_inf"], sd = percent_variability*parameters["P1_inf"]/2),
  P2_inf= rnorm(n_runs, mean = parameters["P2_inf"], sd = percent_variability*parameters["P2_inf"]/2),
  P3_inf= rnorm(n_runs, mean = parameters["P3_inf"], sd = percent_variability*parameters["P3_inf"]/2),
  P4_inf= rnorm(n_runs, mean = parameters["P4_inf"], sd = percent_variability*parameters["P4_inf"]/2),
  daly_acute_hiv_offART= rnorm(n_runs, mean = parameters["daly_acute_hiv_offART"], sd = percent_variability*parameters["daly_acute_hiv_offART"]/2),
  daly_preAIDs_offART= rnorm(n_runs, mean = parameters["daly_preAIDs_offART"], sd = percent_variability*parameters["daly_preAIDs_offART"]/2),
  daly_AIDs_offART= rnorm(n_runs, mean = parameters["daly_AIDs_offART"], sd = percent_variability*parameters["daly_AIDs_offART"]/2),
  Psi_12m= rnorm(n_runs, mean = parameters["Psi_12m"], sd = percent_variability*parameters["Psi_12m"]/2),
  Psi_12c= rnorm(n_runs, mean = parameters["Psi_12c"], sd = parameters["Psi_12c"] * percent_variability/2),
  Psi_13m= rnorm(n_runs, mean = parameters["Psi_13m"], sd = parameters["Psi_13m"] * percent_variability/2),
  Psi_13c= rnorm(n_runs, mean = parameters["Psi_13c"], sd = parameters["Psi_13c"] * percent_variability/2),
  Psi_24m= rnorm(n_runs, mean = parameters["Psi_24m"], sd = parameters["Psi_24m"] * percent_variability/2),
  Psi_24c= rnorm(n_runs, mean = parameters["Psi_24c"], sd = parameters["Psi_24c"] * percent_variability/2),
  Psi_34m= rnorm(n_runs, mean = parameters["Psi_34m"], sd = parameters["Psi_34m"] * percent_variability/2),
  Psi_34c= rnorm(n_runs, mean = parameters["Psi_34c"], sd = parameters["Psi_34c"] * percent_variability/2),
  n_1m= rnorm(n_runs, mean = parameters["n_1m"], sd = parameters["n_1m"] * percent_variability/2),
  n_1c= rnorm(n_runs, mean = parameters["n_1c"], sd = parameters["n_1c"] * percent_variability/2),
  n_2m= rnorm(n_runs, mean = parameters["n_2m"], sd = parameters["n_2m"]*percent_variability/2),
  n_2c= rnorm(n_runs, mean = parameters["n_2c"], sd = parameters["n_2c"] * percent_variability/2),
  n_3m= rnorm(n_runs, mean = parameters["n_3m"], sd = parameters["n_3m"]*percent_variability/2),
  n_3c= rnorm(n_runs, mean = parameters["n_3c"], sd = parameters["n_3c"] * percent_variability/2),
  n_3co= rnorm(n_runs, mean = parameters["n_3co"], sd = parameters["n_3co"] * percent_variability/2),
  n_4m= rnorm(n_runs, mean = parameters["n_4m"], sd = parameters["n_4m"]*percent_variability/2),
  n_4c= rnorm(n_runs, mean = parameters["n_4c"], sd = parameters["n_4c"] * percent_variability/2),
  n_4co= rnorm(n_runs, mean = parameters["n_4co"], sd = parameters["n_4co"] * percent_variability/2),
  w= rnorm(n_runs, mean = parameters["w"], sd = parameters["w"] * percent_variability/2),
  health_insur_c= rnorm(n_runs, mean = parameters["health_insur_c"], sd = parameters["health_insur_c"] * percent_variability/2),
  prep_cost= rnorm(n_runs, mean = parameters["prep_cost"], sd = parameters["prep_cost"] * percent_variability/2),
  art_cost= rnorm(n_runs, mean = parameters["art_cost"], sd = parameters["art_cost"] * percent_variability/2))


#####
for(i in 1:n_runs) {
  # Extract the ith row of parameters as a list to pass to the model
  #BASELINE CONDITIONS#####################################################
  # Define initial conditions
  N_0 <- 27914536*0.55
  N_female_0 <- N_0*0.51
  N_male_0 <- N_0*(1-0.51)
  N1_0 <- N_female_0*(1-0.02)
  N2_0 <- N_male_0*(1-0.14)
  N3_0 <- N_male_0*0.14
  N4_0 <- N_female_0*0.02
  
  S1_0 <- N1_0*(1-0.027)
  E1_0 <- N1_0*0.027*(1-0.0625)
  I1_0 <- N1_0*0.027*0.0625
  R1_0 <- 0
  
  S2_0 <- N2_0*(1-0.027)
  E2_0 <- N2_0*0.027*(1-0.0625)
  I2_0 <- N2_0*0.027*0.0625
  R2_0 <- 0
  
  S3_0 <- N3_0*(1-0.027)
  E3_0 <- N3_0*0.027*(1-0.0625)
  I3_0 <- N3_0*0.027*0.0625
  R3_0 <- 0
  
  S4_0 <- N4_0*(1-0.027)
  E4_0 <- N4_0*0.027*(1-0.0625)
  I4_0 <- N4_0*0.027*0.0625
  R4_0 <- 0
  
  Cost_0<- 0
  Dalys_0<-0
  onset_inf_0<-0
  
  state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
             S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
             S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
             S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)
  
  #####
  param_set <- c(parameters_sim[i, ])
  # Run the model with these parameters
  output100p <- ode(y = state, times = times, func = HIV_model2, parms = param_set)
  output100p <- as.data.frame(output100p)
  output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
  output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
  output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
  output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
  output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
  totalpop100p<- output100p[['totalpop']][51]
  dalys_onset <- diff(c(0, output100p[['Dalys']]))
  costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
  dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*resultsp0[[i,4]]
  cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*resultsp0[[i,4]]
  # Store the model output in the results list
  results[[i,1]] <- (cost_100p) -(resultsp0[[i,1]]) #incremental costs
  results[[i,2]] <- (resultsp0[[i,2]])-((sum(output100p[['I']]+output100p[['E']])/totalpop100p)*resultsp0[[i,4]]) #incremental infections averted
  results[[i,3]] <- (resultsp0[[i,3]])-(dalys_100p)  #incremental DALYs
}
results$ICERinf=results[,1]/results[,2]
results$ICERdal=results[,1]/results[,3]
results100<- results
mean(results100$V1)
mean(results100$V2)
mean(results100$V3)
mean(results100$V4)
mean(results100$V5)
mean(results100$V1)/mean(results100$V2)
count <- sum(results100$V2 < 0 & results100$V1 > 0)
count2 <- sum(results100$V3 < 0 & results100$V1 > 0)


totalpop100p
dalys_100p
totalpop75p
dalys_75p
totalpop50p
dalys_50p
totalpop25p
dalys_25p



################################################################################################################
################################################################################################################
################################################################################################################
#GRAPHS:
lancet_colors <- c("grey70","#f28e79","#a1c181","#f7bb5f","#b699d7")
#Graph25c####
mean_xd1 <- mean(results25$V3 / 10000)
mean_xi1 <- median(results25$V2 / 10000)
mean_y1 <- mean(results25$V1 / 1000000)
p1_i <- ggplot(results25, aes(x = V2/10000, y = V1/1000000)) +
  geom_point(shape = 21, color = "black", fill = "#f28e79", size = 3) +  # Circle with border and fill
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 37, face = "bold"),
        plot.subtitle = element_text(size = 37),
        axis.title.x = element_text(size = 26, vjust = -0.2),
        axis.title.y = element_text(size = 26, vjust = 1.5)) +
  labs(x = "Incremental Infections averted (in 10,000s)",
       y = "Incremental Costs ($ in 1,000,000s)") +
  theme_lancet() +
  scale_y_continuous(labels = label_comma()) +
  scale_y_continuous(breaks = seq(-1000, 3000, by = 500)) +
  scale_x_continuous(breaks = seq(-2500, 2500, by = 500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  geom_point(aes(x = mean_xi1, y = mean_y1), shape = 23, color = "black", fill = "#add8e9", size = 5)+
  coord_cartesian(xlim = c(-2500, 2500)) +
  coord_cartesian(ylim = c(-2500, 3000))  # Shape 4 is a cross, change as needed

p1_d <- ggplot(results25, aes(x = V3/10000, y = V1/1000000)) +
  geom_point(shape = 21, color = "black", fill = "#f28e79", size = 3) +  # Circle with border and fill
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 37, face = "bold"),
        plot.subtitle = element_text(size = 37),
        axis.title.x = element_text(size = 26, vjust = -0.2),
        axis.title.y = element_text(size = 26, vjust = 1.5)) +
  labs(x = "Incremental DALYs (in 10,000s)",
       y = "Incremental Costs (£ in 10,000s)") +
  theme_lancet() +
  scale_y_continuous(labels = label_comma()) +
  scale_y_continuous(breaks = seq(-1500, 3000, by = 500)) +
  scale_x_continuous(breaks = seq(-500, 500, by = 100), limits = c(-500,500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  geom_point(aes(x = mean_xd1, y = mean_y1),  shape = 23, color = "black", fill = "#add8e9", size = 5)+
  coord_cartesian(ylim = c(-1501, 3001))# Shape 4 is a cross, change as needed


wtp_thresholds <- seq(0, 4500, by = 100)
# Calculate proportions and 95% CIs
ci_lower <- numeric(length(wtp_thresholds))
ci_upper <- numeric(length(wtp_thresholds))
proportion_cost_effective <- numeric(length(wtp_thresholds))
for (i in seq_along(wtp_thresholds)) {
  wtp = wtp_thresholds[i]
  successes = sum(results25$ICERinf <= wtp, na.rm = TRUE)
  trials = sum(!is.na(results25$ICERinf))
  ci <- binom.test(successes, trials)$conf.int
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
  proportion_cost_effective[i] <- successes / trials}
proportion_cost_effective<- proportion_cost_effective-0.42
ci_lower=ci_lower-0.42
ci_upper= ci_upper-0.42

data_for_plot <- data.frame(
  WTP = wtp_thresholds,
  Proportion = proportion_cost_effective,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
library(ggplot2)
k1i<-ggplot(data_for_plot, aes(x = WTP)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#f28e79", alpha = 0.4) + # Less saturated color with transparency for CIs
  geom_line(aes(y = Proportion), color = "#f28e79") +  # Less saturated color for line
  theme_minimal() +
  labs(x = "Willingness-To-Pay Threshold ($ per infection averted)",
       y = "Proportion of simulations being cost-effective",
       title = " ") +
  theme(plot.title = element_text(size = 14))+
  theme_lancet()+
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  scale_x_continuous(breaks = seq(0, 4500, by = 500)) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1), limits = c(0, 0.8)) 
  


wtp_thresholds <- seq(0, 4500, by = 100)
# Calculate proportions and 95% CIs
ci_lower <- numeric(length(wtp_thresholds))
ci_upper <- numeric(length(wtp_thresholds))
proportion_cost_effective <- numeric(length(wtp_thresholds))
for (i in seq_along(wtp_thresholds)) {
  wtp = wtp_thresholds[i]
  successes = sum(results25$ICERdal <= wtp, na.rm = TRUE)
  trials = sum(!is.na(results25$ICERdal))
  ci <- binom.test(successes, trials)$conf.int
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
  proportion_cost_effective[i] <- successes / trials}
proportion_cost_effective<- proportion_cost_effective-0.42
ci_lower=ci_lower-0.42
ci_upper= ci_upper-0.42

data_for_plot2 <- data.frame(
  WTP = wtp_thresholds,
  Proportion = proportion_cost_effective,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
library(ggplot2)
k1d<-ggplot(data_for_plot2, aes(x = WTP)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#f28e79", alpha = 0.4) + # Less saturated color with transparency for CIs
  geom_line(aes(y = Proportion), color = "#f28e79") +  # Less saturated color for line
  theme_minimal() +
  labs(x = "Willingness-To-Pay Threshold ($ per DALY averted)",
       y = "Probability of being cost-Effective",
       title = " ") +
  theme(plot.title = element_text(size = 14))+
  theme_lancet()+
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  scale_x_continuous(breaks = seq(0, 45000, by = 500), limits = c(0,4500)) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1), limits = c(0, 0.8)) 
#####
#Graph50c####
mean_xd2 <- mean(results50$V3 / 10000)
mean_xi2 <- median(results50$V2 / 10000)
mean_y2 <- mean(results50$V1 / 1000000)
p2_i <- ggplot(results50, aes(x = V2/10000, y = V1/1000000)) +
  geom_point(shape = 21, color = "black", fill = "#a1c181", size = 3) +  # Circle with border and fill
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 37, face = "bold"),
        plot.subtitle = element_text(size = 37),
        axis.title.x = element_text(size = 26, vjust = -0.2),
        axis.title.y = element_text(size = 26, vjust = 1.5)) +
  labs(x = "Incremental infections averted (in 10,000s)",
       y = "Incremental Costs (£ in 1,000,000s)") +
  theme_lancet() +
  scale_y_continuous(labels = label_comma()) +
  scale_y_continuous(breaks = seq(-800, 3600, by = 400)) +
  scale_x_continuous(breaks = seq(-2500, 2500, by = 500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  geom_point(aes(x = mean_xi2, y = mean_y2), shape = 23, color = "black", fill = "#add8e9", size = 5)+
  coord_cartesian(xlim = c(-2500, 2500))# Shape 4 is a cross, change as needed

p2_d <- ggplot(results50, aes(x = V3/10000, y = V1/1000000)) +
  geom_point(shape = 21, color = "black", fill = "#a1c181", size = 3) +  # Circle with border and fill
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 37, face = "bold"),
        plot.subtitle = element_text(size = 37),
        axis.title.x = element_text(size = 26, vjust = -0.2),
        axis.title.y = element_text(size = 26, vjust = 1.5)) +
  labs(x = "Incremental DALYs (in 10,000s)",
       y = "Incremental Costs (£ in 1,000,000s)") +
  theme_lancet() +
  scale_y_continuous(labels = label_comma()) +
  scale_y_continuous(breaks = seq(-1000, 3500, by = 500)) +
  scale_x_continuous(breaks = seq(-500, 500, by = 100), limits = c(-500,500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  geom_point(aes(x = mean_xd2, y = mean_y2),shape = 23, color = "black", fill = "#add8e9", size = 5)
#coord_cartesian(xlim = c(-100, 100))# Shape 4 is a cross, change as needed


wtp_thresholds <- seq(0, 4500, by = 100)
# Calculate proportions and 95% CIs
ci_lower <- numeric(length(wtp_thresholds))
ci_upper <- numeric(length(wtp_thresholds))
proportion_cost_effective <- numeric(length(wtp_thresholds))
for (i in seq_along(wtp_thresholds)) {
  wtp = wtp_thresholds[i]
  successes = sum(results50$ICERinf <= wtp, na.rm = TRUE)
  trials = sum(!is.na(results50$ICERinf))
  ci <- binom.test(successes, trials)$conf.int
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
  proportion_cost_effective[i] <- successes / trials}
proportion_cost_effective<- proportion_cost_effective-0.36
ci_lower=ci_lower-0.36
ci_upper= ci_upper-0.36

data_for_plot <- data.frame(
  WTP = wtp_thresholds,
  Proportion = proportion_cost_effective,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
library(ggplot2)
k2i<-ggplot(data_for_plot, aes(x = WTP)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#a1c181", alpha = 0.4) + # Less saturated color with transparency for CIs
  geom_line(aes(y = Proportion), color = "#a1c181") +  # Less saturated color for line
  theme_minimal() +
  labs(x = "Willingness-To-Pay Threshold ($ per infection averted)",
       y = "Proportion of simulations being cost-effective",
       title = " ") +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  theme(plot.title = element_text(size = 14))+
  theme_lancet()+
  scale_x_continuous(breaks = seq(0, 4500, by = 500)) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1), limits = c(0, 0.8)) 


wtp_thresholds <- seq(100, 4500, by = 100)
# Calculate proportions and 95% CIs
ci_lower <- numeric(length(wtp_thresholds))
ci_upper <- numeric(length(wtp_thresholds))
proportion_cost_effective <- numeric(length(wtp_thresholds))
for (i in seq_along(wtp_thresholds)) {
  wtp = wtp_thresholds[i]
  successes = sum(results50$ICERdal <= wtp, na.rm = TRUE)
  trials = sum(!is.na(results50$ICERdal))
  ci <- binom.test(successes, trials)$conf.int
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
  proportion_cost_effective[i] <- successes / trials}
proportion_cost_effective<- proportion_cost_effective-0.37
ci_lower=ci_lower-0.37
ci_upper= ci_upper-0.37

data_for_plot2 <- data.frame(
  WTP = wtp_thresholds,
  Proportion = proportion_cost_effective,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
library(ggplot2)
k2d<-ggplot(data_for_plot2, aes(x = WTP)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#a1c181", alpha = 0.4) + # Less saturated color with transparency for CIs
  geom_line(aes(y = Proportion), color = "#a1c181") +  # Less saturated color for line
  theme_minimal() +
  labs(x = "Willingness-To-Pay Threshold ($ per DALY averted)",
       y = "Proportion of simulations being cost-Effective",
       title = " ") +
  theme(plot.title = element_text(size = 14))+
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  theme_lancet()+
  scale_x_continuous(breaks = seq(0, 4500, by = 500), limits = c(0,4500)) +
  scale_y_continuous(breaks = seq(0.0, 0.8, by = 0.1), limits = c(-0.1, 0.8)) 
#####
#Graph75c####
mean_xd3 <- mean(results75$V3 / 10000)
mean_y3 <- mean(results75$V1 / 1000000)
mean_xi3 <- mean(results75$V2 / 10000)
p3_i <- ggplot(results75, aes(x = V2/10000, y = V1/1000000)) +
  geom_point(shape = 21, color = "black", fill = "#f7bb5f", size = 3) +  # Circle with border and fill
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 37, face = "bold"),
        plot.subtitle = element_text(size = 37),
        axis.title.x = element_text(size = 26, vjust = -0.2),
        axis.title.y = element_text(size = 26, vjust = 1.5)) +
  labs(x = "Incremental infections averted (in 10,000s)",
       y = "Incremental Costs (£ in 1,000,000s)") +
  theme_lancet() +
  scale_y_continuous(labels = label_comma()) +
  scale_y_continuous(breaks = seq(-1000, 3500, by = 500)) +
  scale_x_continuous(breaks = seq(-2500, 2500, by = 500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  geom_point(aes(x = mean_xi3, y = mean_y3), shape = 23, color = "black", fill = "#add8e9", size = 5)+
  coord_cartesian(xlim = c(-2500, 2500))+
  coord_cartesian(ylim = c(-1000, 3700))
  # Shape 4 is a cross, change as needed

p3_d <- ggplot(results75, aes(x = V3/10000, y = V1/1000000)) +
  geom_point(shape = 21, color = "black", fill = "#f7bb5f", size = 3) +  # Circle with border and fill
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 37, face = "bold"),
        plot.subtitle = element_text(size = 37),
        axis.title.x = element_text(size = 26, vjust = -0.2),
        axis.title.y = element_text(size = 26, vjust = 1.5)) +
  labs(x = "Incremental DALYs (in 10,000s)",
       y = "Incremental Costs (£ in 1,000,000s)") +
  theme_lancet() +
  scale_y_continuous(labels = label_comma()) +
  #scale_y_continuous(breaks = seq(-1000, 1000, by = 200)) +
  scale_y_continuous(breaks = seq(-500, 4500, by = 500)) +
  scale_x_continuous(breaks = seq(-500, 500, by = 100), limits = c(-500, 500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  geom_point(aes(x = mean_xd3, y = mean_y3), shape = 23, color = "black", fill = "#add8e9", size = 5)
#coord_cartesian(xlim = c(-100, 100))# Shape 4 is a cross, change as needed


wtp_thresholds <- seq(0, 4500, by = 100)
# Calculate proportions and 95% CIs
ci_lower <- numeric(length(wtp_thresholds))
ci_upper <- numeric(length(wtp_thresholds))
proportion_cost_effective <- numeric(length(wtp_thresholds))
for (i in seq_along(wtp_thresholds)) {
  wtp = wtp_thresholds[i]
  successes = sum(results75$ICERinf <= wtp, na.rm = TRUE)
  trials = sum(!is.na(results75$ICERinf))
  ci <- binom.test(successes, trials)$conf.int
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
  proportion_cost_effective[i] <- successes / trials}
proportion_cost_effective<- proportion_cost_effective-0.24
ci_lower=ci_lower-0.24
ci_upper= ci_upper-0.24

data_for_plot <- data.frame(
  WTP = wtp_thresholds,
  Proportion = proportion_cost_effective,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
library(ggplot2)
k3i<-ggplot(data_for_plot, aes(x = WTP)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#f7bb5f", alpha = 0.4) + # Less saturated color with transparency for CIs
  geom_line(aes(y = Proportion), color = "#f7bb5f") +  # Less saturated color for line
  theme_minimal() +
  labs(x = "Willingness-To-Pay Threshold ($ per infection averted)",
       y = "Proportion of simulations being cost-effective",
       title = " ") +
  theme(plot.title = element_text(size = 14))+
  theme_lancet()+
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  scale_x_continuous(breaks = seq(0, 4500, by = 500)) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1) ,limits = c(0, 0.8)) 


wtp_thresholds <- seq(100, 4500, by = 100)
# Calculate proportions and 95% CIs
ci_lower <- numeric(length(wtp_thresholds))
ci_upper <- numeric(length(wtp_thresholds))
proportion_cost_effective <- numeric(length(wtp_thresholds))
for (i in seq_along(wtp_thresholds)) {
  wtp = wtp_thresholds[i]
  successes = sum(results75$ICERdal <= wtp, na.rm = TRUE)
  trials = sum(!is.na(results75$ICERdal))
  ci <- binom.test(successes, trials)$conf.int
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
  proportion_cost_effective[i] <- successes / trials}
proportion_cost_effective<- proportion_cost_effective-0.26
ci_lower=ci_lower-0.26
ci_upper= ci_upper-0.26

data_for_plot2 <- data.frame(
  WTP = wtp_thresholds,
  Proportion = proportion_cost_effective,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
library(ggplot2)
k3d<-ggplot(data_for_plot2, aes(x = WTP)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#f7bb5f", alpha = 0.4) + # Less saturated color with transparency for CIs
  geom_line(aes(y = Proportion), color = "#f7bb5f") +  # Less saturated color for line
  theme_minimal() +
  labs(x = "Willingness-To-Pay Threshold ($ per DALY averted)",
       y = "Proportion of simulations being cost-Effective",
       title = " ") +
  theme(plot.title = element_text(size = 14))+
  theme_lancet()+
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  scale_x_continuous(breaks = seq(0, 4500, by = 500), limits = c(0,4500)) +
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1), limits=c(-0.05,0.8) ) 

#####
#Graph100c####
mean_xd4 <- mean(results100$V3 / 10000)
mean_xi4 <- mean(results100$V2 / 10000)
mean_y4 <- mean(results100$V1 / 1000000)
p4_i <- ggplot(results100, aes(x = V2/10000, y = V1/1000000)) +
  geom_point(shape = 21, color = "black", fill = "#b699d7", size = 3) +  # Circle with border and fill
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 37, face = "bold"),
        plot.subtitle = element_text(size = 37),
        axis.title.x = element_text(size = 26, vjust = -0.2),
        axis.title.y = element_text(size = 26, vjust = 1.5)) +
  labs(x = "Incremental infections averted (in 10,000s)",
       y = "Incremental Costs (£ in 1,000,000s)") +
  theme_lancet() +
  scale_y_continuous(labels = label_comma()) +
  scale_y_continuous(breaks = seq(-500, 4500, by = 500)) +
  scale_x_continuous(breaks = seq(-2500, 2500, by = 500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  geom_point(aes(x = mean_xi4, y = mean_y4), shape = 23, color = "black", fill = "#add8e9", size = 5)+
  coord_cartesian(xlim = c(-2500, 2500))# Shape 4 is a cross, change as needed

p4_d <- ggplot(results100, aes(x = V3/10000, y = V1/1000000)) +
  geom_point(shape = 21, color = "black", fill = "#b699d7", size = 3) +  # Circle with border and fill
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 26, vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size=26, angle = 0, vjust = 0.5, hjust = 1),
        plot.title = element_text(size = 37, face = "bold"),
        plot.subtitle = element_text(size = 37),
        axis.title.x = element_text(size = 26, vjust = -0.2),
        axis.title.y = element_text(size = 26, vjust = 1.5)) +
  labs(x = "Incremental DALYs (in 10,000s)",
       y = "Incremental Costs (£ in 1,000,000s)") +
  theme_lancet() +
  scale_y_continuous(labels = label_comma()) +
  #scale_y_continuous(breaks = seq(-1000, 1000, by = 200)) +
  scale_y_continuous(breaks = seq(0, 5000, by = 500)) +
  scale_x_continuous(breaks = seq(-500, 500, by = 100), limits = c(-500,500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  geom_point(aes(x = mean_xd4, y = mean_y4), shape = 23, color = "black", fill = "#add8e9", size = 5)
  #coord_cartesian(xlim = c(-100, 100))# Shape 4 is a cross, change as needed


wtp_thresholds <- seq(100, 4500, by = 100)
# Calculate proportions and 95% CIs
ci_lower <- numeric(length(wtp_thresholds))
ci_upper <- numeric(length(wtp_thresholds))
proportion_cost_effective <- numeric(length(wtp_thresholds))
for (i in seq_along(wtp_thresholds)) {
  wtp = wtp_thresholds[i]
  successes = sum(results100$ICERinf <= wtp, na.rm = TRUE)
  trials = sum(!is.na(results100$ICERinf))
  ci <- binom.test(successes, trials)$conf.int
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
  proportion_cost_effective[i] <- successes / trials}
proportion_cost_effective<- proportion_cost_effective-0.13
ci_lower=ci_lower-0.13
ci_upper= ci_upper-0.13

data_for_plot <- data.frame(
  WTP = wtp_thresholds,
  Proportion = proportion_cost_effective,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
library(ggplot2)
k4i<-ggplot(data_for_plot, aes(x = WTP)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#b699d7", alpha = 0.4) + # Less saturated color with transparency for CIs
  geom_line(aes(y = Proportion), color = "#b699d7") +  # Less saturated color for line
  theme_minimal() +
  labs(x = "Willingness-To-Pay Threshold (£ per infection averted)",
       y = "Proportion of simulations being cost-effective",
       title = " ") +
  theme(plot.title = element_text(size = 14))+
  theme_lancet()+
  scale_x_continuous(breaks = seq(0, 4500, by = 500))+
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1), limits = c(0, 0.86))  


wtp_thresholds <- seq(100, 4500, by = 100)
# Calculate proportions and 95% CIs
ci_lower <- numeric(length(wtp_thresholds))
ci_upper <- numeric(length(wtp_thresholds))
proportion_cost_effective <- numeric(length(wtp_thresholds))
for (i in seq_along(wtp_thresholds)) {
  wtp = wtp_thresholds[i]
  successes = sum(results100$ICERdal <= wtp, na.rm = TRUE)
  trials = sum(!is.na(results100$ICERdal))
  ci <- binom.test(successes, trials)$conf.int
  ci_lower[i] <- ci[1]
  ci_upper[i] <- ci[2]
  proportion_cost_effective[i] <- successes / trials}
proportion_cost_effective<- proportion_cost_effective-0.15
ci_lower=ci_lower-0.15
ci_upper= ci_upper-0.15

data_for_plot2 <- data.frame(
  WTP = wtp_thresholds,
  Proportion = proportion_cost_effective,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper)
library(ggplot2)
k4d<-ggplot(data_for_plot2, aes(x = WTP)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), fill = "#b699d7", alpha = 0.4) + # Less saturated color with transparency for CIs
  geom_line(aes(y = Proportion), color = "#b699d7") +  # Less saturated color for line
  theme_minimal() +
  labs(x = "Willingness-To-Pay Threshold ($ per DALY averted)",
       y = "Proportion of simulations being cost-Effective",
       title = " ") +
  theme(plot.title = element_text(size = 14))+
  theme_lancet()+
  scale_x_continuous(breaks = seq(0, 4500, by = 500), limits = c(0,4500)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey", size = 0.5)+
  scale_y_continuous(breaks = seq(0, 0.8, by = 0.1), limits=c(-0.05,0.8)) 
#####

library(patchwork)
# Correctly adjust each plot
p1_i <- p1_i + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
p2_i <- p2_i + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
p3_i <- p3_i + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
p4_i <- p4_i + theme(axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 10, l = 2, unit = "mm")) # Only bottom margin adjusted for last plot in column
k1i <- k1i + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
k2i <- k2i + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
k3i <- k3i + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
k4i <- k4i + theme(axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 10, l = 2, unit = "mm")) # Adjust for x-title in last plot
library(ggplot2)
library(patchwork)
library(cowplot)
# Assuming your combined columns plot is ready

p1_i <- p1_i + labs(subtitle = "25% coverage: PSA analysis")
p2_i <- p2_i + labs(subtitle = "50% coverage: PSA analysis")
p3_i <- p3_i + labs(subtitle = "75% coverage: PSA analysis")
p4_i <- p4_i + labs(subtitle = "100% coverage: PSA analysis")
k1i <- k1i + labs(subtitle = "25% coverage: WTP thresholds")
k2i <- k2i + labs(subtitle = "50% coverage: WTP thresholds")
k3i <- k3i + labs(subtitle = "75% coverage: WTP thresholds")
k4i <- k4i + labs(subtitle = "100% coverage: WTP thresholds")
# Now combine your plots
combined_columns = (p1_i / p2_i / p3_i / p4_i | k1i / k2i / k3i / k4i) + 
  plot_layout(guides = 'collect')
# Function to create a spacer plot with a y-title
create_y_title_plot <- function(title, y_pos = 0.5) {
  plot = ggdraw() +
    draw_label(title, angle = 90, hjust = 1, x = 0.5, y = y_pos, fontface = 'bold', size = 15) +
    theme_void() +
    theme(plot.margin = margin(0, 10, 0, 10)) # Adjust margins as needed
  return(plot)}
# Create y-title plots
y_title_1 = create_y_title_plot("Incremental costs (£ in 1,000,000s)",  y_pos = 0.7)
y_title_2 = create_y_title_plot("Proportion of simulations being cost-effective (ICER<WTP)", y_pos = 0.8)
# Combine everything
final_plot = (y_title_1 | combined_columns | y_title_2) + 
  plot_layout(widths = c(1, 20, 1)) # Adjust the widths as needed
# Print the final plot
print(final_plot)
ggsave("sens_anICERinf.tiff", plot = final_plot, width = 12, height = 13, dpi = 1000)

#2DALYS.
p1_d <- p1_d + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
p2_d <- p2_d + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
p3_d <- p3_d + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
p4_d <- p4_d + theme(axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 10, l = 2, unit = "mm")) # Only bottom margin adjusted for last plot in column
k1d <- k1d + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
k2d <- k2d + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
k3d <- k3d + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))
k4d <- k4d + theme(axis.title.y = element_blank(), plot.margin = margin(t = 2, r = 2, b = 10, l = 2, unit = "mm")) # Adjust for x-title in last plot


p1_d <- p1_d + labs(subtitle = "25% coverage: PSA analysis")
p2_d <- p2_d + labs(subtitle = "50% coverage: PSA analysis")
p3_d <- p3_d + labs(subtitle = "75% coverage: PSA analysis")
p4_d <- p4_d + labs(subtitle = "100% coverage: PSA analysis")
k1d <- k1d + labs(subtitle = "25% coverage: WTP thresholds")
k2d <- k2d + labs(subtitle = "50% coverage: WTP thresholds")
k3d <- k3d + labs(subtitle = "75% coverage: WTP thresholds")
k4d <- k4d + labs(subtitle = "100% coverage: WTP thresholds")


library(ggplot2)
library(patchwork)
library(cowplot)
# Assuming your combined columns plot is ready
combined_columns = (p1_d / p2_d / p3_d / p4_d | k1d / k2d / k3d / k4d) + plot_layout(guides = 'collect')
# Function to create a spacer plot with a y-title
create_y_title_plot <- function(title, y_pos = 0.5) {
  plot = ggdraw() +
    draw_label(title, angle = 90, hjust = 1, x = 0.5, y = y_pos, fontface = 'bold', size = 15) +
    theme_void() +
    theme(plot.margin = margin(0, 10, 0, 10)) # Adjust margins as needed
  return(plot)}
# Create y-title plots
y_title_1 = create_y_title_plot("Incremental costs (£ in 1,000,000s)",  y_pos = 0.7)
y_title_2 = create_y_title_plot("Proportion of simulations being cost-effective (ICER<WTP)", y_pos = 0.7)
# Combine everything
final_plot2 = (y_title_1 | combined_columns | y_title_2) + 
  plot_layout(widths = c(1, 20, 1)) # Adjust the widths as needed
# Print the final plot
print(final_plot2)
ggsave("sens_anICERinf_dalys.tiff", plot = final_plot2, width = 12, height = 13, dpi = 1000)





#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
#MATRIX of outputs
results <- function(s,e_s,parameters,times,state){
  
  #redefinimos parametros
  parameters['s'] <- s
  parameters['e_s'] <- e_s
  
  # Solve equations
  output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                    method = "rk4")
  
  # Convert to data frame for easy extraction of columns
  output <- as.data.frame(output_raw)
  
  # Total de la población
  output['N1'] <- output['S1']+output['E1']+output['I1']
  output['N2'] <- output['S2']+output['E2']+output['I2']
  output['N3'] <- output['S3']+output['E3']+output['I3']
  output['N4'] <- output['S4']+output['E4']+output['I4']
  output['N_men'] <- output['N2']+output['N3']
  output['N_women'] <- output['N1']+output['N4']
  output['N'] <- output['N_men']+output['N_women']
  
  #### DEFINICIÓN VARIABLES 
  output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] #muertes totales
  output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
  output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
  output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
  output['Infected']<-output['I']+output['E']
  
  
  #NEW INFECTED
  
  # Assuming your dataframe is named df
  output$Aux <- rep(0, nrow(output))
  
  # Define tau, delta, and mu vectors
  tau <- c(parameters['tau_1'],parameters['tau_2'],parameters['tau_3'], parameters['tau_4'])
  delta <- c(parameters['delta_1'], parameters['delta_2']) 
  mu <- c(parameters['mu_1'], parameters['mu_2'], parameters['mu_3'], parameters['mu_4'])
  
  # Perform the operation in a loop
  for (i in 1:4) {
    output$Aux <- output$Aux + ((1 - tau[i]) * delta[1] + tau[i] * delta[2]) * output[[paste0("I", i)]] +
      mu[i] * (output[[paste0("E", i)]] + output[[paste0("I", i)]])
  }
  
  #definimos nuevos infectados
  library(dplyr)
  output <- output %>%
    mutate(Infected_lag = lag(Infected))  # Asegúrate de que Infected_lag exista
  
  output<-output %>%
    mutate(newInfected=ifelse(row_number() == 1, 0, 
                              Infected - Infected_lag + Aux)) 
  
  
  #TOTAL COST
  
  s<-parameters['s']
  person_cost<-45
  
  output['totalCost']<-output['S4']*s*person_cost
  
  #Infectados totales
  totalInfected<-sum(output$newInfected)
  #Mortalidad
  totalDeath<-tail(output$R, 1)
  #Costo total
  totalCost<-sum(output$totalCost)
  
  result <- list(totalInfected = totalInfected, totalDeath = totalDeath, totalCost = totalCost)
  return(result)
}

#######VARIAMOS VALORES Y FORMAMOS MATRICES
# Definir las secuencias
lancet_colors <- c( "#a1c181", "#FAFAD2", "#f28e79")  
s_values <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,0.75, 0.8, 0.85, 0.9, 0.95,1)
e_s_values <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7,0.75, 0.8, 0.85, 0.9, 0.95,1)
infected_matrix <- matrix(NA, nrow = length(s_values), ncol = length(e_s_values))
death_matrix <- matrix(NA, nrow = length(s_values), ncol = length(e_s_values))

for (i in seq_along(s_values)) {
  for (j in seq_along(e_s_values)) {
    parameters["s"] <- s_values[i]
    parameters["e_s"] <- e_s_values[j]
    output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters, method = "rk4")
    output <- as.data.frame(output_raw)
    # Sum across compartments for infected and R (deaths)
    total_infected <- sum(output[, c("I1", "I2", "I3", "I4", "E1", "E2", "E3", "E4")])/51
    total_deaths <- sum(output[nrow(output), c("R1", "R2", "R3", "R4")])/51  # Store in matrices
    infected_matrix[i, j] <- total_infected
    death_matrix[i, j] <- total_deaths
  }
}
new_x_labels <- c("0","5","10","15", "20","25", "30","35", "40","45", "50","55", "60","65", "70","75", "80","85", "90","95", "100")
# Define new y-axis labels (for rows), assuming you want to change them as well
new_y_labels <- c("0","5","10","15", "20","25", "30","35", "40","45", "50","55", "60","65", "70","75", "80","85", "90","95", "100")
# Apply the new labels to your matrices
colnames(infected_matrix) <- new_x_labels
rownames(infected_matrix) <- new_y_labels
colnames(death_matrix) <- new_x_labels
rownames(death_matrix) <- new_y_labels
infected_matrix<- infected_matrix[nrow(infected_matrix):1, ]
death_matrix <- death_matrix[nrow(death_matrix):1, ]
#Infections:
min_value <- min(infected_matrix, na.rm = TRUE)
max_value <- max(infected_matrix, na.rm = TRUE)
# Create the heatmap
infected_heatmap <- heatmaply(
  infected_matrix,
  xlab = "", 
  ylab = "Health coverage (%)", 
  main = " ",  # Set main to empty since we'll use annotations for a customized title
  dendrogram = "none", 
  colors = colorRampPalette(lancet_colors)(1056)
)
infected_heatmap <- infected_heatmap  %>% layout(
  xaxis = list(tickangle = 0,  # Ensure x-axis labels are horizontal
               title = list(text = "Health insurance efficiency (%)",
                            font = list(size = 21))),
  yaxis = list(title = list(
    text = "Health coverage (%)",
    font = list(size = 21))),
  annotations = list(list(text = "",
                          x = 0, y = 1.055,
                          xref = "paper", yref = "paper",
                          showarrow = FALSE,
                          font = list(size = 21, color = "black"),  # Adjust title font size if needed
                          align = "left")))

# Display the heatmap with the left-aligned title
infected_heatmap


# Convert the matrix/data.frame to a long format suitable for ggplot2
infected_long <- melt(infected_matrix)
# Create the heatmap
p <- ggplot(infected_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRampPalette(lancet_colors)(256), name= 'Number of HIV \nInfections') +
  labs(x = "Health insurance efficiency (%)", y = "Health coverage (%)") +
  theme_minimal() +
  theme_lancet()+
  theme(axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21),
        plot.title = element_text(size = 21, hjust = 0.5))+
  scale_x_continuous(breaks = seq(0, 100, by = 10) ) +
  scale_y_continuous(breaks = seq(0, 100, by = 10))+
  theme(axis.title.x = element_text(size = 21),
        axis.title.y = element_text(size = 21, face="bold"),
        plot.title = element_text(size = 14),
        #legend.position = c(0.9, 0.2),  # Adjust these values as needed
        legend.background = element_rect(colour = "black", fill = NA, size = 0.8),   # Optional: makes legend background transparent
        legend.key = element_rect(colour = "white", size = 2),
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 13),  # Increase the size of the legend text
        legend.title = element_text(size = 17),  # Increase the size of the legend title
        axis.text.x = element_text( hjust = 1, size = 16),  # Increase size and adjust x-axis labels
        axis.text.y = element_text(size = 16))
  
# Save the plot as a TIFF image with 1000 DPI
ggsave("infected_heatmap.tiff", plot = p, width = 12, height = 9, dpi = 1000)






#FIGURE COMBINED FOR SENS ANALYSES.
library(patchwork)
# Add (A) and (B) labels as plot annotations
gsenss <- gsens + annotate("text", label = "(A)", x = 0, y = Inf, hjust = 1.1, vjust = 2, size = 6, fontface = "bold")
ps <- p + annotate("text", label = "(B)", x = 0, y = Inf, hjust = 1.1, vjust = 2, size = 6, fontface = "bold")
# Combine the plots vertically
combined_plotxo <- gsens / p
# Display the combined plot
combined_plotxo
# If you want to save the combined plot to a file
ggsave("combined_plot.tiff", combined_plotxo, width = 14, height = 15, dpi = 1000)




















#deaths:
library(heatmaply)

# Calculate min and max values for the death matrix
min_value <- min(death_matrix, na.rm = TRUE)
max_value <- max(death_matrix, na.rm = TRUE)
# Create the heatmap for death_matrix with adjustments
death_heatmap <- heatmaply(
  death_matrix,
  xlab = "Health insurance efficiency (%)", 
  ylab = "Health coverage (%)", 
  main = " ",  # Main title is set empty; custom title will be added as an annotation
  dendrogram = "none", 
  colors = colorRampPalette(lancet_colors)(1056) )
# Customize the plot layout
death_heatmap <- death_heatmap %>% layout(
  xaxis = list(tickangle = 0,  # Ensure x-axis labels are horizontal
    title = list(text = "Health insurance efficiency (%)",
      font = list(size = 21))),
  yaxis = list(title = list(
      text = "Health coverage (%)",
      font = list(size = 21))),
  annotations = list(list(text = "",
      x = 0, y = 1.055,
      xref = "paper", yref = "paper",
      showarrow = FALSE,
      font = list(size = 21, color = "black"),  # Adjust title font size if needed
      align = "left")))
death_heatmap

#COMBO:
combined_heatmap <- subplot(
  infected_heatmap,
  death_heatmap,
  nrows = 2, # Arrange in 2 rows
  shareX = TRUE, # Share the X axis if applicable
  shareY =TRUE, # Share the Y axis if applicable
  titleX = TRUE, # Show the x-axis titles
  titleY = TRUE, # Show the y-axis titles
  margin = 0.1 # Adjust margins to try and create more space
) %>% layout(
  margin = list(t = 40), # Increase top margin to allow space for titles
  annotations = list( # Reiterate the annotations for titles if necessary
    list(text = "(A) Total HIV infected individuals",
      x = 0.5, y = 1.05, xref = "paper", yref = "paper",
      xanchor = "center", showarrow = FALSE, font = list(size = 16)),
    list(text = "(B) Total HIV associated deaths",
      x = 0.5, y = 0.525, xref = "paper", yref = "paper",
      xanchor = "center", showarrow = FALSE, font = list(size = 16))))


# Display the combined figure
combined_heatmap

########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
###############SCENARIO ANALYSIS! & UNIVARIATE########################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

####################################################
#. GREAT AND LOWER POPULATION SIZES AMONG TFSW
###################################################
discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
###################################################### x
#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.0336)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.004

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######
#CHANGE THE POPULATION SIZE FOR WTGS


####################################################
#. Time horizon
###################################################
discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=10
time <- 0:10
times <- seq(from = 0, to = 10, by = 1)

#baselineCharact####
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.031)
E1_0 <- N1_0*0.031*(1-0.0625)
I1_0 <- N1_0*0.031*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.031)
E2_0 <- N2_0*0.031*(1-0.0625)
I2_0 <- N2_0*0.031*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.031)
E3_0 <- N3_0*0.031*(1-0.0625)
I3_0 <- N3_0*0.031*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.031)
E4_0 <- N4_0*0.031*(1-0.0625)
I4_0 <- N4_0*0.031*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#10 years
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][11]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][11]
infected_bs0Ons<- output[['onset_inf']][11]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][11]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][11]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][11]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][11]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][11]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][11]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][11]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][11]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][11]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][11]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][11]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][11]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][11], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######

#25 years
n=25
time <- 0:25
times <- seq(from = 0, to = 25, by = 1)
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][26]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][26]
infected_bs0Ons<- output[['onset_inf']][26]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][26]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][26]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][26]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][26]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][26]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][26]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][26]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][26]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][26]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][26]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][26]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][26]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][26], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)



####################################################
#. Discount rate
###################################################
discount_rate <- 0.06
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)


######


####################################################
#. Higher/Lower TFSW population
###################################################
discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#baselineCharact####
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.004)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.004

S1_0 <- N1_0*(1-0.031)
E1_0 <- N1_0*0.031*(1-0.0625)
I1_0 <- N1_0*0.031*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.031)
E2_0 <- N2_0*0.031*(1-0.0625)
I2_0 <- N2_0*0.031*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.031)
E3_0 <- N3_0*0.031*(1-0.0625)
I3_0 <- N3_0*0.031*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.031)
E4_0 <- N4_0*0.031*(1-0.0625)
I4_0 <- N4_0*0.031*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######



####################################################
#. Higher/Lower Health insurance cost per year
###################################################
discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#Parameters#####
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 =  0.79,
  tau_4 = 0.417,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 120,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)



#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)


######

discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#Parameters#####
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.79,
  tau_4 = 0.417,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 40,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)







#####
#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)


######


####################################################
#. Higher/Lower efficiency over health insurance
###################################################
discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#Parameters#####
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 1,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.79,
  tau_4 = 0.4179,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)



#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)


######


discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#Parameters#####
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.80,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.79,
  tau_4 = 0.417,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)







#####
#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)


######


####################################################
#. Higher ART coverage among TFSW, equal to populations: 79%
###################################################
discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#Parameters#####
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.79,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.79,
  tau_4 = 0.79,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)



#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)


######


discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#Parameters#####
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.354,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.79,
  tau_4 = 0.354,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)



#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)


######




####################################################
#. PreP coverage among TFSW.
###################################################
discount_rate <- 0.03
#https://bmjopen.bmj.com/content/10/11/e041346
n=50
time <- 0:50
times <- seq(from = 0, to = 50, by = 1)
#Parameters20%#######
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0.20,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.417,
  tau_4 = 0.79,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)













#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######

#Parameters40%#######
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0.40,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.417,
  tau_4 = 0.79,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)













#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######


#Parameters60%#######
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0.60,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.417,
  tau_4 = 0.79,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)













#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######



#Parameters80%#######
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 0.80,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.417,
  tau_4 = 0.79,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)













#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######


#Parameters100%#######
#########Defineparametervalues
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.51,
  theta = 0.035,
  gamma = 1/6,
  kappa = 0.02,
  mu_1 = 0.009,
  mu_2 = 0.009,
  mu_3 = 0.009,
  mu_4 = 0.009,
  eta_1 = 1/10,
  eta_2 = 1/40,
  delta_1 = 0.235,
  delta_2 = 0.163,
  z = 0.01,
  g = 1/3.2,
  beta = 0.002762108,
  epsilon = 0.9,
  pi_12m = 0.7/1.25,
  pi_12c = 0.7,
  pi_13m = 0.7/1.25,
  pi_13c = 0.7,
  pi_24m = 0.505,
  pi_24c = 0.505*1.25,
  pi_34m = 0.505,
  pi_34c = 0.505*1.25,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  nu_1 = 0.0,
  nu_2 = 0.0,
  nu_3 = 0.0,
  nu_4 = 1.0,
  e_nu = 0.51,
  alpha_1 = 0.79,
  alpha_2 = 0.79,
  alpha_3 = 0.79,
  alpha_4 = 0.417,
  e_alpha = 0.96,
  n_1m = 1.5,
  n_1c = 8.1,
  n_2m = 1.5,
  n_2c = 8.1,
  n_3m = 1.95,
  n_3c = 8.1,
  n_3co = 37.55,
  n_4m = 2.05,
  n_4c = 9.5,
  n_4co = 2.7,
  s = 0.0,
  e_s = 0.89,
  tau_1 = 0.79,
  tau_2 = 0.79,
  tau_3 = 0.417,
  tau_4 = 0.79,
  e_tau = 0.96,
  P1_inf = 0.027,
  P2_inf = 0.027,
  P3_inf = 0.03,
  P4_inf = 0.034,
  w = 2,
  #cost & utilities below
  health_insur_c = 80,
  prep_cost= 320,
  art_cost = 220,
  daly_acute_hiv_offART=0.078,
  daly_preAIDs_offART=0.274,
  daly_AIDs_offART=0.582
)













#BASELINE CONDITIONS#####################################################
# Define initial conditions
N_0 <- 27914536*0.55
N_female_0 <- N_0*0.51
N_male_0 <- N_0*(1-0.51)
N1_0 <- N_female_0*(1-0.02)
N2_0 <- N_male_0*(1-0.14)
N3_0 <- N_male_0*0.14
N4_0 <- N_female_0*0.02

S1_0 <- N1_0*(1-0.027)
E1_0 <- N1_0*0.027*(1-0.0625)
I1_0 <- N1_0*0.027*0.0625
R1_0 <- 0

S2_0 <- N2_0*(1-0.027)
E2_0 <- N2_0*0.027*(1-0.0625)
I2_0 <- N2_0*0.027*0.0625
R2_0 <- 0

S3_0 <- N3_0*(1-0.027)
E3_0 <- N3_0*0.027*(1-0.0625)
I3_0 <- N3_0*0.027*0.0625
R3_0 <- 0

S4_0 <- N4_0*(1-0.027)
E4_0 <- N4_0*0.027*(1-0.0625)
I4_0 <- N4_0*0.027*0.0625
R4_0 <- 0

Cost_0<- 0
Dalys_0<-0
onset_inf_0<-0

state <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0, Cost=Cost_0, Dalys= Dalys_0, onset_inf=onset_inf_0)

#####
#MAINresults_baseCase######
# 0 health coverage
parameters[["s"]]<-0
output_raw <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output <- as.data.frame(output_raw)
output['R']<-output['R1']+output['R2']+output['R3']+output['R4'] 
output['I']<-output['I1']+output['I2']+output['I3']+output['I4']
output['E']<-output['E1']+output['E2']+output['E3']+output['E4']
output['S']<-output['S1']+output['S2']+output['S3']+output['S4']
output['totalpop']<- output['R']+output['I']+output['E']+output['S']
popbs0<- output[['totalpop']][51]
output['Infected']<-output['I']+output['E']
#
dalys_onset <- diff(c(0, output[['Dalys']]))
costs_onset <- diff(c(0, output[['Cost']]))# Prepend 0 for the first year to match the length
dalys_bs0 <- sum(dalys_onset / (1 + discount_rate) ^ time)
cost_bs0 <- sum(costs_onset / (1 + discount_rate) ^ time)
#
infected_bs0<-(sum(output['Infected']))
deaths_bs0<-output[['R']][51]
infected_bs0Ons<- output[['onset_inf']][51]



#25% health coverage
parameters[["s"]]<-0.25
output25p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                 method = "rk4")
output25p <- as.data.frame(output25p)
output25p['R']<-output25p['R1']+output25p['R2']+output25p['R3']+output25p['R4'] 
output25p['I']<-output25p['I1']+output25p['I2']+output25p['I3']+output25p['I4']
output25p['E']<-output25p['E1']+output25p['E2']+output25p['E3']+output25p['E4']
output25p['S']<-output25p['S1']+output25p['S2']+output25p['S3']+output25p['S4']
output25p['totalpop']<- output25p['R']+output25p['I']+output25p['E']+output25p['S']
totalpop25p<- output25p[['totalpop']][51]
output25p['Infected']<-output25p['I']+output25p['E']
#
dalys_onset <- diff(c(0, output25p[['Dalys']]))
costs_onset <- diff(c(0, output25p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_25p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
cost_25p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop25p)*popbs0
#
totalinf_25p<-(sum(output25p['Infected'])/totalpop25p)*popbs0
totaldeat_25p<-output25p[['R']][51]
ICER_25p_inf_av= (cost_25p-cost_bs0)/(infected_bs0-totalinf_25p) #cost per infection averted
ICER_25p_dalys= (cost_25p-cost_bs0)/(dalys_bs0-(dalys_25p))
totalinf_25pOns<- (output25p[['onset_inf']][51]/totalpop25p)*popbs0

#50% health coverage
parameters[["s"]]<-0.5
output_50p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output50p <- as.data.frame(output_50p)
output50p['R']<-output50p['R1']+output50p['R2']+output50p['R3']+output50p['R4'] 
output50p['I']<-output50p['I1']+output50p['I2']+output50p['I3']+output50p['I4']
output50p['E']<-output50p['E1']+output50p['E2']+output50p['E3']+output50p['E4']
output50p['S']<-output50p['S1']+output50p['S2']+output50p['S3']+output50p['S4']
output50p['Infected']<-output50p['I']+output50p['E']
output50p['totalpop']<- output50p['R']+output50p['I']+output50p['E']+output50p['S']
totalpop50p<- output50p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output50p[['Dalys']]))
costs_onset <- diff(c(0, output50p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_50p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
cost_50p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop50p)*popbs0
#
totalinf_50p<-(sum(output50p['Infected'])/totalpop50p)*popbs0
totaldeat_50p<-(output50p[['R']][51]/totalpop25p)*popbs0
totalinf_50pOns<- (output50p[['onset_inf']][51]/totalpop50p)*popbs0
ICER_50p_inf_av= (cost_50p-cost_bs0)/(infected_bs0-totalinf_50p) #cost per infection averted
ICER_50p_dalys= (cost_50p-cost_bs0)/(dalys_bs0-dalys_50p)

#75% health coverage
parameters[["s"]]<-0.75
output_75p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                  method = "rk4")
output75p <- as.data.frame(output_75p)
output75p['R']<-output75p['R1']+output75p['R2']+output75p['R3']+output75p['R4'] 
output75p['I']<-output75p['I1']+output75p['I2']+output75p['I3']+output75p['I4']
output75p['E']<-output75p['E1']+output75p['E2']+output75p['E3']+output75p['E4']
output75p['S']<-output75p['S1']+output75p['S2']+output75p['S3']+output75p['S4']
output75p['Infected']<-output75p['I']+output75p['E']
output75p['totalpop']<- output75p['R']+output75p['I']+output75p['E']+output75p['S']
totalpop75p<- output75p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output75p[['Dalys']]))
costs_onset <- diff(c(0, output75p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_75p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
cost_75p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop75p)*popbs0
#
totalinf_75p<-(sum(output75p['Infected'])/totalpop75p)*popbs0
totaldeat_75p<-(output75p[['R']][51]/totalpop75p)*popbs0
totalinf_75pOns<- (output75p[['onset_inf']][51]/totalpop75p)*popbs0
ICER_75p_inf_av= (cost_75p-cost_bs0)/(infected_bs0-totalinf_75p) #cost per infection averted
ICER_75p_dalys= (cost_75p-cost_bs0)/(dalys_bs0-dalys_75p)


#100% health coverage
parameters[["s"]]<-1
output_100p <- ode(y = state, times = times, func = HIV_model, parms = parameters,
                   method = "rk4")
output100p <- as.data.frame(output_100p)
output100p['R']<-output100p['R1']+output100p['R2']+output100p['R3']+output100p['R4'] 
output100p['I']<-output100p['I1']+output100p['I2']+output100p['I3']+output100p['I4']
output100p['E']<-output100p['E1']+output100p['E2']+output100p['E3']+output100p['E4']
output100p['S']<-output100p['S1']+output100p['S2']+output100p['S3']+output100p['S4']
output100p['Infected']<-output100p['I']+output100p['E']
output100p['totalpop']<- output100p['R']+output100p['I']+output100p['E']+output100p['S']
totalpop100p<- output100p[['totalpop']][51]
#
dalys_onset <- diff(c(0, output100p[['Dalys']]))
costs_onset <- diff(c(0, output100p[['Cost']]))# Prepend 0 for the first year to match the length
dalys_100p <- ((sum(dalys_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
cost_100p <- ((sum(costs_onset / (1 + discount_rate) ^ time))/totalpop100p)*popbs0
#
totalinf_100p<-(sum(output100p['Infected'])/totalpop100p)*popbs0
totaldeat_100p<-(output100p[['R']][51]/totalpop100p)*popbs0
totalinf_100pOns<- (output100p[['onset_inf']][51]/totalpop100p)*popbs0
ICER_100p_inf_av= (cost_100p-cost_bs0)/(infected_bs0-totalinf_100p) #cost per infection averted
ICER_100p_dalys= (cost_100p-cost_bs0)/(dalys_bs0-dalys_100p)

#RESULTS I.b: CE only varying health insurance coverage and costs per ART & PrEP.
summary_table <- data.frame(
  `Health Coverage` = c(0, 25, 50, 75, 100),
  Cost = c(cost_bs0, cost_25p, cost_50p, cost_75p, cost_100p),
  DALYs = c(dalys_bs0, dalys_25p, dalys_50p, dalys_75p, dalys_100p),
  TotalInfected = c(sum(output$Infected), totalinf_25p, totalinf_50p, totalinf_75p, totalinf_100p),
  TotalDeaths = c(output[['R']][51], totaldeat_25p, totaldeat_50p, totaldeat_75p, totaldeat_100p),
  ICER_DALYs = c(NA, ICER_25p_dalys, ICER_50p_dalys, ICER_75p_dalys, ICER_100p_dalys),
  ICER_Inf_Av = c(NA, ICER_25p_inf_av, ICER_50p_inf_av, ICER_75p_inf_av, ICER_100p_inf_av))
# Print the summary table
View(summary_table)

######








###