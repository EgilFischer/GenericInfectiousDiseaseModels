###############################################################################################
#                                                                                             #
#                              Stochastic without demography                                  #
#                                                                                             #
###############################################################################################

#####simulation model#######
sir.sim <- function(N0, Y0, Z0, beta, gamma, deltat, duration){
  # start at time = 0
  time = 0 
  
  #set initial population variables
  n = n0 
  Y = Y0
  Z = Z0
  X  = N0 - Y0 - Z0
  
  #Create a data.frame to hold the simulated data
  res = data.frame(Time = c(0),X = c(X),Y = c(Y), Z = c(Z))
  
  #loop as long as the length of the simulation has not been filled
  while(time < duration)
  {
    ##simulate the number of new cases, recoveries and deaths during this time step 
    #Each time check if the variabel is not 0 because rbinom cannot draw from a population of size 0
    ifelse (X * Y > 0, cases <- rbinom(1, X,  1 -  exp(- beta * (Y / n) * deltat)),cases <- 0)
    ifelse (Y > 0 , recoveries <- rbinom(1, Y, 1 -  exp(- gamma * deltat)),recoveries <- 0)
    
    #subtract/add from current values
    X =  X - cases 
    Y = Y + cases - recoveries 
    Z = Z + recoveries 
    
    #advance one time step
    time <- time + deltat
    
    #store every thing
    res = rbind(res, data.frame(Time  = c(time),X = c(X),Y = c(Y), Z= c(Z)))
  }
  return(res)
}

#####function for multiple runs#####
sir.sim.multiple <- function(n0, i0, r0, beta, gamma, deltat, duration, repeats)
{
  res <- data.frame(Time = c(),X =c(), Y =c(), Z=c())
  r <-  1
  while(r <= repeats)
  {
    sim <- sir.sim(n0, i0, r0, beta, gamma, deltat, simlength)
    res <- rbind(res, cbind(sim, sim = rep(r, length(sim$Time))))
    r = r + 1
  }
  return(res)
}

############################################################################################################
#                      
#                             Use of the model
#
#############################################################################################################


####parameters#####
beta = 0.6
gamma = 0.2


####simulation settings####
deltat = 0.05
simlength = 50

#####initial values####
n0 = 200
Y0 = 1
Z0 = 0
X0 = n0 - Y0 - Z0

sim <- sir.sim(n0, Y0,Z0, beta, gamma, deltat, simlength)
plot(sim$Time, sim$X, type = "l", xlab = "Time", ylab = "Number of animals",col = "green",ylim = c(0,n0))
lines(sim$Time, sim$Y, type = "l", xlab = "Time", ylab = "Number of animals",col = "firebrick3")
lines(sim$Time, sim$Z, type = "l", xlab = "Time", ylab = "Number of animals",col = "black")
legend(40,60,c("X","Y","Z"),col = c("green","firebrick3", "black"),lwd = 2,seg.len = 1)

plot(sim$X, sim$Y, type = "l", xlab = "X", ylab = "Y")
s <- seq(length(sim$X)-1)
arrows(x0 = sim$X[s], y0 = sim$Y[s], x1 = sim$X[s + 1], y1 = sim$Y[s + 1],length = .075)


#####################Run the same simulation multiple time#################
##how many times do we want to run the simulation
repeats = 100
#run the simulation multiple times
multi.sim <- sir.sim.multiple(n0, Y0,Z0, beta, gamma, deltat, simlength, repeats)

##plot the end of the epidemic in a histogram
hist(multi.sim$Z[floor(multi.sim$Time) == simlength],breaks = n0, xlim = c(0-1,n0+1), xlab = "Final size", main = "Final size distribution",probability = TRUE, ylab = "frequency")

###plot the different runs
plot(x = multi.sim$Time[multi.sim$sim==1], y = multi.sim$Y[multi.sim$sim==1], xlab = "Time", ylab = "X,Y,Z", type = "l", col = "firebrick3", ylim = c(0,200))
lines(x = multi.sim$Time[multi.sim$sim==1], y = multi.sim$X[multi.sim$sim==1], col = "green")
lines(x = multi.sim$Time[multi.sim$sim==1], y = multi.sim$Z[multi.sim$sim==1], col = "black")
for(r in c(2:repeats)) {
  lines(x = multi.sim$Time[multi.sim$sim==r], y = multi.sim$Y[multi.sim$sim==r], col = "firebrick3")
  lines(x = multi.sim$Time[multi.sim$sim==r], y = multi.sim$X[multi.sim$sim==r], col = "green")
  lines(x = multi.sim$Time[multi.sim$sim==r], y = multi.sim$Z[multi.sim$sim==r], col = "black")
  
}

plot(x = multi.sim$X[multi.sim$sim==1], y = multi.sim$Y[multi.sim$sim==1], xlab = "X", ylab = "Y", type = "l", col = "firebrick3",xlim = c(-1,200), ylim = c(-1,200))
for(r in c(2:repeats)) {
  lines(x = multi.sim$X[multi.sim$sim==r], y = multi.sim$Y[multi.sim$sim==r], col = "firebrick3")
}



