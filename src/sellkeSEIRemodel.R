#################################################################################
#                                                                               #
#                     SEIR model with environmental transmission                #
#                     State transitions are stochastic                          #
#                     Environmental contamination is deterministic              #
#                     Copy right Dr. ir. Egil A.J. Fischer                      #
#                     e.a.j.fischer@uu.nl / egil@egilfischer.nl                                                          #
#                                                                               #
#################################################################################

library(dplyr)
library(ggplot2)
library(LambertW)

##parameters##
animals <- 20; # (Number)
beta <- 1 #transmission coefficient from environment to individual per excretion day (time^-1)
decay <- 0.001 #decay of environmental contamination (time^-1)
meanL <- 0.25 #mean duration latency period (time)
meanI <- 1.5 #mean duration of infectious period (time)
varL <- 1; #variance duration latency period (time^2) if 1 exponential distribution
varI <- 1; #variance duration infectious period (time^2) if 1 exponential distribution

### specify distribution of periods functions
dL <- function(n){return(rgamma(n, shape = meanL*sqrt(varL), scale = sqrt(varL/meanL))) } #Duration latency period as function of a random variable U (time)
dI <- function(n){return(rgamma(n, shape = meanI*sqrt(varI), scale = sqrt(varI/meanI)))  } #Duration infectious period as function of a random variable U (time)


##initial states##
L0 <- 1 #initial exposed / latent infections (integer)
I0 <- 0 #initial infectious animals (integer)

##simulation settings
runs <- 3 #number of runs (integer)
chop.env <- 10^-2 #environment less than 10^-10 contaminated  set to 0 (number)

#initialize

## initialize outputs
output <- data.frame(time = 0, N = animals, S = animals - 1 , L =1, I =0, R=0, e = 0, D = 0, C=0, Th =0, run = 0)
state <- output
#
while(state$run < runs)
{
  QIRtimes<- data.frame(Q = sort(rexp(animals - 1,1)),
                        E = dL(animals-1), 
                        I = dI(animals-1))
state <- data.frame(time = 0, N = animals, S = animals - 1 , L =1, I =0, R=0,e =0, D = 0, C=0, Th =0, run = last(state$run) + 1)
first.latency <- dL(runif(1))
events <-data.frame(type = c("LI","IR"), time = c(first.latency,first.latency + dI(runif(1))))
time <- 0
foi <- 0
cumInf <- 0
cLI <- 0
cIR <-0
cSL <-0

handbreak =0
while(length(events$time) > 0 & state$L + state$I + state$e > 0 & first(state$time) < 10^10)
{
  handbreak = handbreak + 1
  if(handbreak > animals){stop}
  #order events by time of execution
  events <- events[order(events$time),]
  #update environmental contamination
  state$e <- state$e * exp(-decay * (first(events$time) - state$time)) + state$I / decay * (1- exp(-decay * (first(events$time) - state$time)) )
  state$e <- chop.env*(floor(state$e/chop.env)) #create a possibility for the environmental contamination to die out
  #process thefirst event in the list
  #update the cummulative infection pressure
  cumInf <- cumInf + 
    (beta/decay^2) * (decay * state$e - state$I + (exp(-decay * (first(events$time) - state$time)))) * ((state$I - decay * state$e) + decay * state$I * (first(events$time) - state$time)) 
  
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
    if(state$I ==0)print(events)
    cIR <- cIR + 1
    #add one to R
    state$R <- state$R + 1
    #subtract one from I
    state$I <- state$I - 1
   
  }
  if(first(events$type) == "SL" )
  {
     #take into account that environmental contamination might have gone extinct.
    
    if(state$e > 0){
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
   
      #order events by time of execution
      events <- events[order(events$time),]
    } 
   
  } 
 

  
  #determine next infection event
  if(state$S > 0){
    #order events by time of execution
    events <- events[order(events$time),]
    if(state$I ==0)
    {
      infection = state$time +log(1 + (decay *  (first(QIRtimes$Q) - cumInf))/(state$e * beta))/decay
    }else{
     invW <- (exp(-1 + (decay * (-decay * (first(QIRtimes$Q) - cumInf) + state$e / beta))/(state$I /beta)) * (decay *state$e - state$I))/state$I
    if(invW < -1/exp(1) ) {  print("invW < -1/E")  
      print(state)
      print((first(QIRtimes$Q) - cumInf))
      }
     infection <- state$time + (1/decay) - (state$e/state$I) + decay * (first(QIRtimes$Q) - cumInf) / (beta  * state$I) +  W(invW)/decay
    }
   
  } else {infection <- 10^10}
  
 
  
  #remove first event
  events <- events[-1,]
  #order events by time of execution
  events <- events[order(events$time),]
  
  #if this event is previous to other events schedule it
  if(length(events$time) > 0)
  {
    if(infection < first(events$time))
    {
      events <- rbind(data.frame(time = c(infection),type = "SL"),events)
      #order events by time of execution
      events <- events[order(events$time),]
    }
  }else {
    if(state$S > 0){
    events <- rbind(data.frame(time = c(infection),type = "SL"),events)
    #order events by time of execution
    events <- events[order(events$time),]
    }
  }
  #order events by time of execution
  events <- events[order(events$time),]
  
  #record this moment
  output <- rbind(output, state)
}
}
##plot the output
ggplot(data = output[output$run>0, ])+
  geom_point(aes(x = time, y = L), color = "blue")+
  geom_path(aes(x = time, y = L), color = "blue")+
  geom_point(aes(x = time, y = I), color = "red")+
  geom_point(aes(x = time, y = R), color = "black")+
  geom_path(aes(x = time, y = R), color = "black")+
  geom_path (aes(x = time, y = e), color = "gray")+
  facet_grid(run~.)
ggplot(data = output[output$run>0, ])+
  geom_path (aes(x = time, y = e), color = "gray")+
  facet_grid(run~.)
#output

