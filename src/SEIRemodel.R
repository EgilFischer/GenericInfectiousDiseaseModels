#################################################################################
#                                                                               #
#                     SEIR model with environmental transmission                #
#                     Deterministic version                                   #
#                     Copy right Dr. ir. Egil A.J. Fischer                      #
#                     e.a.j.fischer@uu.nl / egil@egilfischer.nl                                                          #
#                                                                               #
#################################################################################

library(deSolve)
library(ggplot2)
library(reshape2)


# The model
#Ordinary differential Equations
## this model assumes exponentially distributed latent and infectious periods
SEIRe <- function(Time, State, Pars){
  with(as.list(c(State,Pars)),{
    dX <- -1 * beta * X * e/N  + mu * W +mu * Y + mu * Z
    dW <- 1 * beta * X * e/N - mu * W - sigma *W
    dY <- sigma * W  - gamma * Y - mu * Y
    dZ <- gamma * Y - mu * Z
    de <- -decay * e + omega * Y
    return(list(c(dX,dW,dY,dZ,de)))
  })
}

# parameters
meanL <- 0.25;  #mean duration latency period (time) transferred to exponential distribution parameter sigma
meanI <- 1.5; gamma <- 1/meanI #mean duration of infectious period (time) transferred to exponential distribution parameter gamma
pars <- c(N = 20,   #herd/flock size (number)
          beta = 1, #transmission coefficient from environment to individual per excretion  (time^-1)
          sigma = 1/meanL, #transition L-> I rate    (time^-1)
          gamma = 1/meanI, #recovery rate I-> R      (time^-1)
          mu = 0,          #background mortality rate(time^-1)
          omega = 1, #excretion by animal (infectious material per time)
          decay =0.1 #decay of environmental contamination (time^-1)
          )
#initial values
W0 = 1;               #latent infections
Y0 = 0;               #infectious
Z0 = 0;               #recovered
e0 = 0;               #environment
X0 = as.numeric(pars['N'] - W0-Y0-Z0); #initial susceptible
init <- c(X = X0,W = W0, Y = Y0,Z = Z0,e = e0)

#Data storage time
dt = 0.1#timestep for storing data
times <- seq(0, 50, by = dt)



#Solve the ordinary differential equations
ode.out <- ode(init, times, SEIRe,pars) 
ode.out.data <- reshape2::melt(data.frame(ode.out),id.vars =c("time"))
#plot 
ggplot(data = ode.out.data)+geom_path(aes(time, value, colour = variable))
