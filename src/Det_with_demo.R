library(deSolve)

###############################################################################################
#                                                                                             #
#                              Deterministic with demography                                  #
#                                                                                             #
###############################################################################################

#Ordinary differential Equations
sir <- function(Time, State, Pars){
 with(as.list(c(State,Pars)),{
  dX <- -1 * beta * X * Y /sum(State) +mu *Y + mu * Z
  dy <- beta * X * Y / sum(State)  - gamma * Y - mu * Y
  dz <- gamma * Y - mu * Z
  return(list(c(dX,dy,dz)))
 })
}

#Parameters
pars <- c(beta = 0.6, #beta = transmission rate
          gamma = 0.3, #gamma = recovery rate
          mu = 0.05#mu = mortality rate
          )
#initial values
N0 = 30000 #population size
Y0 = 1 #initial fraction infected 
Z0 = 0# initial fraction recovered, 
X0 = N0 - Y0 - Z0
init <- c(X = X0,Y = Y0,Z = Z0)

#Data storage time
dt = 0.5#timestep for storing data
sim_duration  = 150 #length of the simulation
times <- seq(0, sim_duration, by = dt)

#Solve the ordinary differential equations
ode.out <- ode(init, times, sir,pars) 
#plot X, Y, and Z against time
plot(x = ode.out[,1], y = ode.out[,2], type = "l", xlab = "time", ylab = "Count", col = "purple",lwd =2,ylim=c(0,N0))
lines(x = ode.out[,1], y = ode.out[,3], type = "l", xlab = "time", ylab = "Count", col = "red",lwd =2)
lines(x = ode.out[,1], y = ode.out[,4], type = "l", xlab = "time", ylab = "Count", col = "palegreen1",lwd =2)
legend(x = 200, y = 27000., c("X","Y","Z"),lwd =2, col = c("purple", "red", "palegreen1"))
#plot  X against Y
plot(x = ode.out[,2], y = ode.out[,3], type = "l", xlab ="X", ylab ="Y")
s <- seq(length(ode.out[,2])-1)
arrows(x0 = ode.out[s,2], y0 = ode.out[s,3], x1 = ode.out[s + 1,2], y1 = ode.out[s + 1 ,3],length = .075)


