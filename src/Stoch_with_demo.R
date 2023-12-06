###############################################################################################
#                                                                                             #
#                              Stochastic with demography                                  #
#                                                                                             #
###############################################################################################

#####simulation model#######
sir.sim <- function(N0, Y0, Z0, beta, gamma, mu ,deltat, duration){
  # start at time = 0
  time = 0 
 
   #set initial population variables
  N = N0 
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
    ifelse (X * Y > 0, cases <- rbinom(1, X,  1 -  exp(- beta * (Y / N) * deltat)),cases <- 0)#new cases
    ifelse (Y > 0 , recoveries <- rbinom(1, Y, 1 -  exp(- gamma * deltat)),recoveries <- 0) #recoveries
    ifelse (Y > 0 , ideaths <- rbinom(1, Y, 1 -  exp(- mu * deltat)),ideaths <- 0) #deaths in I class
    ifelse (Z > 0 , rdeaths <- rbinom(1, Z, 1 -  exp(- mu * deltat)), rdeaths <- 0)#deaths in R class
    #The simulation assumes only 1 event per capita per time step
    if (ideaths +  recoveries > Y) {
      print("Time step too large, two occurences of an event in the same individual")
      break()
      }
    #subtract/add from current values
    X =  X - cases + ideaths + rdeaths
    Y = Y + cases - recoveries - ideaths
    Z = Z + recoveries - rdeaths
    
    #advance one time step
    time <- time + deltat
    
    #store every thing
    res = rbind(res, data.frame(Time  = c(time),X = c(X),Y = c(Y), Z= c(Z)))
  }
  return(res)
}

#####function for multiple runs#####
sir.sim.multiple <- function(N0, Y0, Z0, beta, gamma, mu ,deltat, duration, repeats)
{
  res <- data.frame(Time = c(),X =c(), Y =c(), Z=c())
  r <-  1
  while(r <= repeats)
  {
    sim <- sir.sim(N0, Y0, Z0, beta, gamma,mu, deltat, simlength)
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
beta = 0.3
gamma = 0.1
mu = 0.05

####simulation settings####
deltat = 0.1
simlength = 100
repeats = 100

#####initial values####
N0 = 100
Y0 = 1
Z0 = 0
X0 = N0 - Y0 - Z0

sim <- sir.sim(N0, Y0, Z0, beta, gamma,mu, deltat, simlength)
plot(sim$Time, sim$X, type = "l", xlab = "Time", ylab = "Number of animals",col = "black", ylim= (c(0,N0)))
lines(sim$Time, sim$Y, type = "l", xlab = "Time", ylab = "Number of animals",col = "red")
lines(sim$Time, sim$Z, type = "l", xlab = "Time", ylab = "Number of animals",col = "blue")
legend(80,100,c("X","Y","Z"),col = c("black","red", "blue"),lwd = 2,seg.len = 1)

plot(sim$X, sim$Y, type = "l", xlab = "X", ylab = "Y")
s <- seq(length(sim$X)-1)
arrows(x0 = sim$X[s], y0 = sim$Y[s], x1 = sim$X[s + 1], y1 = sim$Y[s + 1],length = .075)


#####################Run the same simulation multiple time#################
multi.sim <- sir.sim.multiple(N0, Y0, Z0, beta, gamma,mu, deltat, simlength,repeats)

###plot the different runs
plot(x = multi.sim$Time[multi.sim$sim==1], y = multi.sim$Z[multi.sim$sim==1], xlab = "Time", ylab = "Z", type = "l", col = "firebrick3", ylim = c(0,100))
for(r in c(2:repeats)) {
  lines(x = multi.sim$Time[multi.sim$sim==r], y = multi.sim$Z[multi.sim$sim==r], col = "firebrick3")
}

plot(x = multi.sim$X[multi.sim$sim==1], y = multi.sim$Y[multi.sim$sim==1], xlab = "X", ylab = "Y", type = "l", col = "firebrick3",xlim = c(-1,100), ylim = c(-1,100))
for(r in c(2:repeats)) {
  lines(x = multi.sim$X[multi.sim$sim==r], y = multi.sim$Y[multi.sim$sim==r], col = "firebrick3")
}


##plot the end of the epidemic in a histogram
hist(multi.sim$Z[floor(multi.sim$Time) == simlength], breaks = 50,xlab = "Final size", main = "Final size distribution")

##plot the end of the epidemic (Y == 0)
select<- function(r){
  last = multi.sim$Time[multi.sim$sim==r & multi.sim$Y ==0]
  last = min(last)
  if(last == Inf){
  return(simlength+1)} 
  return(last)

  }
extinctiontime= sapply(FUN = select, c(1:repeats)) 
extinctiontime
hist(extinctiontime)
