######################################################
#              Beta Fitting        #
######################################################
setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/HIV_modellingC/Code")
if (!require(remotes)) { install.packages("remotes"); stopifnot(require("remotes")) } 

# Load in the deSolve package
if (!require(deSolve)) { install.packages("deSolve"); stopifnot(require("deSolve")) } 
library(deSolve)
library(stats)

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
  Psi_34co <- parms["Psi_34co"]
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
  Lambda_4co <- lambda_4*Psi_34co*(1-s*e_s)*(1-nu_4*e_nu)*(1-epsilon*pi_34co)*B_3
  
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
  
  res <- list(c(dS1, dE1, dI1, dR1, dS2, dE2, dI2,dR2, dS3, dE3, dI3, dR3, dS4, dE4, dI4, dR4))
  return(res)
}

#PARAMETERS####
parameters <- c(
  phi_1 = 0,
  phi_2 = 0,
  p = 0.501,
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
  beta = 0.009549963,
  epsilon = 0.9,
  pi_12m = 0.3355,
  pi_12c = 0.592,
  pi_13m = 0.3365,
  pi_13c = 0.592,
  pi_24m = 0.3365,
  pi_24c = 0.592,
  pi_34m = 0.3365,
  pi_34c = 0.592,
  pi_34co = 0.5,
  Psi_12m = 75,
  Psi_12c = 5.2,
  Psi_13m = 75,
  Psi_13c = 5.35,
  Psi_24m = 84,
  Psi_24c = 4.25,
  Psi_34m = 75,
  Psi_34c = 5.2,
  Psi_34co = 203,
  nu_1 = 0,
  nu_2 = 0,
  nu_3 = 0.44,
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
  n_3co = 37.75,
  n_4m = 2.05,
  n_4c = 9.5,
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
  w = 2
)

######
# Define time to solve equations
time <- seq(from = 0, to = 21, by = 1)

# Define initial conditions
#BASELINE_conditions####
N_0 <- 15493253*0.52 #initial pop at 1990
N_female_0 <- N_0*0.504
N_male_0 <- N_0*(1-0.504)
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
#####
initial_values <- c(S1 = S1_0, E1 = E1_0, I1 = I1_0, R1 = R1_0,
           S2 = S2_0, E2 = E2_0, I2 = I2_0, R2 = R2_0,
           S3 = S3_0, E3 = E3_0, I3 = I3_0, R3 = R3_0,
           S4 = S4_0, E4 = E4_0, I4 = I4_0, R4 = R4_0)

# Solving differential equations
hiv_model_output <- ode(y = initial_values, times = time, func = HIV_model, parms = parameters,
                  method = "rk4")


#Observed data
library(readxl)
#setwd("/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/HIV_modellingC/Code")
ruta_archivo <- "/Users/lsh1807578/CISS Dropbox/kasim allel henriquez/B_Projects/HIV_modellingC/Code/adult_prevalence.xlsx"
prevalence_data <- read_excel(ruta_archivo)
observed_data <- prevalence_data[,'adult prevalence']/100

#cost function
cost_function <- function(beta,observed_data,hiv_model_output){
  parameters['beta']<-beta
  model_output <- ode(y=initial_values, times = time, func=HIV_model, parms=parameters)
  
  # Total de la población
  N1 <- model_output[,'S1']+model_output[,'E1']+model_output[,'I1']
  N2 <- model_output[,'S2']+model_output[,'E2']+model_output[,'I2']
  N3 <- model_output[,'S3']+model_output[,'E3']+model_output[,'I3']
  N4 <- model_output[,'S4']+model_output[,'E4']+model_output[,'I4']
  N <-N1+N2+N3+N4
  
  #Infectados
  I <- model_output[,'I1']+model_output[,'I2']+model_output[,'I3']+model_output[,'I4']
  E <- model_output[,'E1']+model_output[,'E2']+model_output[,'E3']+model_output[,'E4']
  Infected = I + E
  
  #prevalence
  fitted_values <- Infected/N
  
  #cost
  cost<-sum((observed_data-fitted_values)^2)
  
  return(cost)
  
}
  
#initial guess for beta
initial_beta<-0.006
  
# Using optim to estimate beta
fit <- optim(par = initial_beta, fn = cost_function, observed_data = observed_data, hiv_model_output = hiv_model_output, method = "Brent", lower = 0, upper = 1)

# Best fit beta
beta_best_fit <- fit$par
