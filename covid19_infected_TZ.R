## source: https://github.com/bbolker/bbmisc/blob/master/peak_I_simple.rmd
# https://tinu.shinyapps.io/Flatten_the_Curve/
# https://staff.math.su.se/hoehle/blog/2020/03/16/flatteningthecurve.html


install.packages("tidyverse")
library(tidyverse)

install.packages("pillar")
library(pillar)

install.packages("tidyr")
library(tidyr)

install.packages("dplyr")
library(dplyr)

install.packages("magrittr")
library(dplyr)

install.packages("ggplot2")
library(ggplot2); theme_set(theme_bw(base_size=16))

######################################################################
# sources: https://raw.githubusercontent.com/hoehleatsu/hoehleatsu.github.io/master/_source/2020-03-16-flatteningthecurve.Rmd
# Function to compute the derivative of the ODE system
#
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system
#
# Returns:
#  list with one component being a vector of length two containing
#  dS(t)/dt and dI(t)/dt
######################################################################


# SW
# Infected <- c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,7,7,12,14,15,21,35,94,101,161,203,248,355,500,599,814,961,1022,1103,1190,1279,1439,1639,1763,1934,2046,2286,2526,2840,3069,3447,3700,4028,4435,4947,5568,6131,6443,6830,7206,7693,8419,9141,9685,10151,10151)

# FR
# Infected <- c(0,0,2,3,3,3,4,5,5,5,6,6,6,6,6,6,6,11,11,11,11,11,11,11,12,12,12,12,12,12,12,12,12,12,14,18,38,57,100,130,191,204,285,377,653,949,1126,1209,1784,2281,2281,3661,4469,4499,6633,7652,9043,10871,12612,14282,16018,19856,22304,25233,29155,32964,37575,40174,44550,52128,56989,59105,64338,89953,92839,98010,109069,112950,117749,124869,129654,129654)

# US
# Infected <- c(1,1,2,2,5,5,5,5,5,7,8,8,11,11,11,11,11,11,11,11,12,12,13,13,13,13,13,13,13,13,15,15,15,51,51,57,58,60,68,74,98,118,149,217,262,402,518,583,959,1281,1663,2179,2727,3499,4632,6421,7783,13747,19273,25600,33276,43847,53740,65778,83836,101657,121465,140909,161831,188172,213372,243762,275586,308853,337072,366667,396223,429052,461437,496535,526396,526396)

# BE #Note. count starting as of march 1st when cases>1
#Infected <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,8,13,23,50,109,169,200,239,267,314,314,559,689,886,1058,1243,1486,1795,2257,2815,3401,3743,4269,4937,6235,7284,9134,10836,11899,12775,13964,15348,16770,18431,19691,20814,22194,23403,24983,26667,28018,28018)

# TZ
Infected <- c(3,6,6,6,12,12,12,12,13,13,14,14,19,19,20,20,20,20,22,24,24,25,25,32,32,32,49,53,88,94,147,147,170,254,254,284,284,299,299,299,299,299,480,480,480)

Day <- 1:(length(Infected)) # length of the infection period

# Population size: 1e6 (1 milliion) or 1e7 (10 millions)
#N <- 1e6 OR, better, actual countries' population

#N <- 66000000 # pupulation of the UK
#N <- 10230000 # pupulation of the SW
#N <- 66990000 # pupulation of the FR
#N <- 32820000 # pupulation of the US
#N <- 11460000 # pupulation of the BE
N <- 56320000 # pupulation of the TZ

# Plot coutnry's confirmed cases
old <- par(mfrow = c(1, 2))
#plot(Day, Infected, type ="b") # bubbles
plot(Day, Infected, type ="b")
plot(Day, Infected, log = "y") # type b by default
abline(lm(log10(Infected) ~ Day))
title("Confirmed COVID-19 cases in TZ, lin & log scales", outer = TRUE, line = -2)


# Population size: 1e6 (1 milliion) or 1e7 (10 millions)
#N <- 1e6 
# Rate at which person stays in the infectious compartment (disease specific and tracing specific)
# γ=0.2 corresponding here to an average length of the infective period of 5 days.
gamma <- 1/5 
#gamma <- 0.3813767
# Infectious contact rate - beta = R0/N*gamma and when R0 \approx 2.25 then  2.25/N*gamma i.e
# a contact rate of β=0.0000004 means that the contact rate with a given individual is 0.0000004 contacts per day.
beta <- 4.5e-07 
#beta <- 0.6186234
# R0 for the beta and gamma values
# Altogether, this leads to an R0 of 2.25, which roughly corresponds to the R0 of SARS-CoV-2
#R0 <- beta*N/gamma
R0 <- beta/gamma
#[1] 2.25e-06
#[1] 1.62208


## We can now solve the ODE system using the above parameters and an initial number of infectious of, say, 10:

# Load package to numerically solve ODEs
library(deSolve)

# Grid where to evaluate
#max_time <- 150
times <- seq(0, max_time, by=0.01)

# time in days, prediction to up to 90 days (highly unlikely fact will confirm the hypothesis though)
#fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
fit <- data.frame(rk4(y = c(N - 10, 10), times = times, func = sir, parms = c(beta, gamma)))
col <- 1:3 # colour
print(fit)

#sink('fit1_TZ_03052020.txt'); fit; sink()
#write.csv(fit, "fit1_TZ_03052020.csv")

# Plot the ODE system.
## let's try with this sir model and extrapolate to maybe 90 days, via epidemiologic modeling
SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta/N * I * S
    dI <- beta/N * I * S - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}


init <- c(S = N-Infected[1], I = Infected[1], R = 0)
RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((Infected - fit)^2)
}

# general optimazion nethod based on Nelder-Mead
# https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method
Opt <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1)) # optimize with some sensible conditions
Opt$message
## [1] [1] "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH"

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
print(Opt_par)
#beta     gamma 
#0.5592483 0.4407518 


# Grid where to evaluate
max_time <- 300 # here predictions up to 150 days, i.e + (length(Day)-max_time) beyond length(Day)
times <- seq(0, max_time, by=0.01)

# time in days, prediction to up to 90 days (highly unlikely fact will confirm the hypothesis though)
#fit <- data.frame(ode(y = init, times = t, func = SIR, parms = Opt_par))
## We are trying here another fitting method than with RK
fit <- data.frame(ode(y = init, times = times, func = SIR, parms = Opt_par))
col <- 1:3 # colour
print(fit)

#sink('fit0_TZ_03052020.txt'); fit; sink()
#write.csv(fit, "fit0_TZ_03052020.csv")


# Plot infection curve
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col)
legend("topright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 1, col = col, inset = 0.05)
matplot(fit$time, fit[ , 2:4], type = "l", xlab = "Day", ylab = "Number of subjects", lwd = 2, lty = 1, col = col, log = "y")

points(Day, Infected)
#legend("bottomright", c("Susceptibles", "Infecteds", "Recovereds"), lty = 1, lwd = 1, col = col, inset = 0.05)
title("Predicted 2019-nCoV in TZ (worst case), SIR", outer = TRUE, line = -2)


## We can now solve the ODE system using the above parameters and an initial number of infectious of, say, 10:
# Load package to numerically solve ODEs
library(deSolve)

# Solve ODE system using Runge-Kutta numerical method.
# ode_solution <- rk4(y = c(N - 10, 10), times = times, func = sir, parms = c(beta, gamma)) %>%
ode_solution <- ode(y = init, times = times, func = SIR, parms = Opt_par) %>%
  as.data.frame() %>%
  setNames(c("t", "S", "I", "R")) %>%
  mutate(beta = beta, gama = gamma, R0 = N * beta / gamma, s = S / N, i = I / N, type = "without_intervention")

ode_solution_daily <- ode_solution %>%
  filter(t %in% seq(0, max_time, by = 1)) %>%
  mutate(C = if_else(row_number() == 1, 0, lag(S) - S), c = C / N)


## The epidemic curve of new infections per day is shown below:
df_plot <- ode_solution_daily %>% select(t, c) %>% 
  pivot_longer(-t, names_to= "Quantity", values_to = "Proportion") 
ggplot(df_plot, aes(x=t, y=Proportion)) + geom_col() + 
  xlab("Time(d)") + ylab("Daily new cases (%TZ pop.)") + scale_y_continuous(labels=scales::percent)


# Function to compute the final size.
s_inf <- function(R0) {
  
  f_target <- function(x) { x - exp(-R0*(1-x)) }
  
  result <- uniroot(f_target, lower=1e-12, upper=1-1e-12)$root
  
  return(result)
}

# Final proportion of infected.
#1 - s_inf(R0)
#[1] 0.8534382


## We can use the above equation to verify that the larger $R_0$, the larger is the final size of the outbreak:
R0_grid <- c(1.25, 1.5, 1.75,2, 2.25, 2.5, 3)
map_dbl( R0_grid, ~ 1-s_inf(.x)) %>% setNames(R0_grid) %>% scales::percent(accuracy=1)


