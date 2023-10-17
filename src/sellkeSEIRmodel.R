# 
packages <- c("ggplot2","dplyr")


## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)



##parameters##
animals <- 100 #population size
beta <- 0.2 #transmission coefficient
dL <- function(U){return(-log(U)/0.5)} #Duration latency period as function of a random variable U
dI <- function(U){return(-log(U)/0.1)} #Duration infectious period as function of a random variable U
#L0 <- 1 #number of initially latently infected; this is hardcoded later commented out for later implementation
#I0 <- 0 #number of initially infectious; this is hardcoded later commented out for later implementation
runs <- 100 #number of runs

#initialize
output <- c()
state <- data.frame(time = 0, N = animals, S = animals - 1 , L =1, I =0, R=0, D = 0, C=0, Th =0, run = 0)

#loop over number of runs
while(state$run < runs)
{
  #get for every animal the infection threshold Q, the length of the latent period and length of infectious period
  QIRtimes<- data.frame(Q = sort(rexp(animals - 1,1)),
                        E = sapply(FUN = dL, X = runif(animals-1)), 
                        I = sapply(FUN = dI, X = runif(animals-1)))
  
  #set state of the system with initial values and add it to the output
  state <- data.frame(time = 0, N = animals, S = animals - 1 , L =1, I =0, R=0, D = 0, C=0, Th =0, run = state$run + 1)
  output <- rbind(output, state)
  #initialize event list
  events <-data.frame(type = c("LI","IR"), time = c(dL(runif(1)),dI(runif(1))))
  events$time[2] = events$time[1]+events$time[2]
  #set initial values 
  time <- 0 #time
  foi <- 0  #force of infection
  cumInf <- 0 #cumulative force of infection
  
  #counters for output/debugging not necessary for running the simulations
  cLI <- 0 #counter for latent to infectious transitions
  cIR <-0  #counter for infectious to recover transitions
  cSL <-0 #counter for susceptible to latent transitions
  
  #set manual 'handbrake' that prevents endless loops
  handbrake =0
  while(length(events$time) > 0 & state$L + state$I > 0)
  {
    handbrake = handbrake + 1
    if(handbrake > animals*3){break}
    #process the first event in the list
    #update the cumulative infection pressure
    cumInf <- cumInf + foi * (first(events$time) - state$time)
    #set time to current time
    state$time <- first(events$time)
    #determine the next event
    if(first(events$type) == "LI" ){
      cLI <- cLI + 1
      #add one to I type 
      state$I <- state$I + 1
      #subtract one from L type
      state$L <- state$L - 1
    }
    if(first(events$type) == "IR" )
    {
      cIR <- cIR + 1
      #add one to R
      state$R <- state$R + 1
      #subtract one from I
      state$I <- state$I - 1
      
    }
    if(first(events$type) == "SL" )
    {
      cSL <- cSL + 1
      #add one to L
      state$L <- state$L + 1
      #subtract one from S
      state$S <- state$S - 1
      #set transitions
      events<- rbind(events, data.frame(time =c(state$time + QIRtimes[1,]$E), type = c("LI"))) #L -> I
      events<- rbind(events, data.frame(time =c(state$time + QIRtimes[1,]$E + QIRtimes[1,]$I), type = c("IR"))) #I -> R
      
      #remove lowest resistance
      QIRtimes <-QIRtimes[-1,] 
      
    } 
    
    #set force-of-infection
    foi <- beta * state$I /state$N
    
    #order events by time of execution
    events <- events[order(events$time),]
    
    
    #determine next infection event
    if(state$S * state$I > 0){
      infection <- state$time + (first(QIRtimes$Q) - cumInf)/foi
    } else {infection <- 10^10}
    
    #remove first event
    events <- events[-1,]
    
    #if the infection event is previous to other events schedule it
    if(length(events$time) > 0)
    {
      if(infection < first(events$time))
      {
        events <- rbind(data.frame(time = c(infection),type = "SL"),events)
        #order events by time of execution
        events <- events[order(events$time),]
      }
    }
    
    
    #record this moment
    output <- rbind(output, state)
  }
}
##plot the output
ggplot(data = output)+
  geom_path(aes(x = time, y = L, group = run), color = "blue")+
  geom_path(aes(x = time, y = I, group = run), color = "red")+
  geom_path(aes(x = time, y = R, group = run), color = "black")+
  ylab("L,I,R")

output
