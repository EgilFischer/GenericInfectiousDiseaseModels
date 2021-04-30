#################################################################################
#                                                                               #
#                     SEIR model with environmental transmission                #
#                     State transitions are stochastic                          #
#                     Environmental contamination is deterministic              #
#                     Copy right Dr. ir. Egil A.J. Fischer                      #
#                     e.a.j.fischer@uu.nl / egil@egilfischer.nl                                                          #
#                                                                               #
#################################################################################

packages <- c("dplyr","ggplot2","LambertW")

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
animals <- 10; # (Number)
beta <- 1 #transmission coefficient from environment to individual per excretion day (time^-1)
decay <- 0.1 #decay of environmental contamination (time^-1)
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
R0 <- 0 #initial recovered animals (integer)
e0 <- 0 #initial contamination of the environment (decimal)

##simulation settings
runs <- 10 #number of runs (integer)
chop.env <- 10^-2 #environment less than 10^-10 contaminated  set to 0 (number)
max.length <- 5 #maximum length of the simulation (time)

#initialize
state <- data.frame(run = 0)
output<- NULL

#multiple runs  of the simulation
#select lines 47 until 213 and press run
while(state$run <= runs)
{ 
  #start new run
  print(paste("Now running: ", state$run)); 
  
  #set up a new run
  if(exists("QIRtimes")){rm(QIRtimes)} #remove previous simulation
  if(exists("events")){rm(events)} #remove previous simulation
  seed <- floor(runif(1,min =0,max = 10^5));print(paste("Using Seed:",seed))
  set.seed(seed)
  
  #set time, foi and cumInf to zero
  foi <- 0
  cumInf <- 0
  
  #set counts of events to 0
  cSL <-0
  cLI <- 0
  cIR <-0
  
  #set-up infection thresholds, latency periods, and infectious periods
  QIRtimes<- data.frame(Q = sort(rexp(animals - 1,1)),
                        E = dL(animals-1), 
                        I = dI(animals-1));
  #set initial state 
  state <- data.frame(time = 0, 
                      N = animals, 
                      S = animals - L0 - I0 , 
                      L =L0, 
                      I =I0, 
                      R=R0, 
                      e = e0, 
                      D = 0, 
                      C=0, 
                      Th =0, 
                      run = last(state$run) + 1);
  #set firsst events
  events <-data.frame(type = c(rep(c("LI","IR"),L0),rep("IR",I0)), #name of event (L to I, I to R)
                     time = c(c(mapply(function(l,i)cumsum(c(l,i)),dL(L0),dI(L0))),dI(I0))); #random time of event
  if(!is.infinite(max.length)) events <- rbind(events,data.frame(type = "END",time = max.length)); # if shorter than infinity add end of simulation time
  events <- events[order(events$time),] #sort events starting with the first
  

#run until either still infection, still events, or the next event time is smaller than 10^10
while(length(events$time) > 0 & 
      #state$L + state$I + state$e > 0 &
       first(state$time) < 10^10)
{
  # if(handbreak > animals){stop}
  #order events by time of execution
  events <- events[order(events$time),]
  
  #process thefirst event in the list
  #calculate the lenght of the time step 
  deltat <- (first(events$time) - state$time)
  
  #update the cummulative infection pressure
  cumInf <- cumInf + 
    beta * (state$e*(1-exp(-decay*deltat))/decay) + beta*(state$I*(exp(-decay*deltat)-1+decay*deltat)/decay^2)
  
    #update environmental contamination
  state$e <- state$e * exp(-decay * deltat) + state$I / decay * (1- exp(-decay * deltat) )
  state$e <- chop.env*(floor(state$e/chop.env)) #create a possibility for the environmental contamination to die out
  
  
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
     #take into account that environmental contamination might have gone extinct.
    
    if(state$e > 0){
      cSL <- cSL + 1
      #add one to L
      state$L <- state$L + 1
      #subtract one from S
      state$S <- state$S - 1
      #set transitions
      events<- rbind(events, 
                     data.frame(time =c(state$time + QIRtimes[1,]$E), type = c("LI"))) #L -> I
      events<- rbind(events, 
                     data.frame(time =c(state$time + QIRtimes[1,]$E + QIRtimes[1,]$I), type = c("IR"))) #I -> R
      
      
      #remove lowest resistance
      QIRtimes <-QIRtimes[-1,] 
   
      #order events by time of execution
      events <- events[order(events$time),]
    } 
    
   
  } 
  if(first(events$type) == "END" )
  {
    #add some output to have a nice graph of the environmental contamination
    if(state$I+state$L == 0){
      next.step <- tail(output,1)
      for(dt in seq(deltat/10,deltat-deltat/10,deltat/10))
      {
        next.step$time <- tail(output,1)$time + deltat/10
        next.step$e <- tail(output,1)$e* exp(-decay * deltat/10)
        output<- rbind(output,
                       next.step
                       )
          }
    }
    
    
  }

  
  #determine next infection event
  if(state$S > 0){
    #order events by time of execution
    events <- events[order(events$time),]
    if(state$I ==0) # no infectious indivuals shedding
    {
      #infection time is time + log of decaying infectivity divided b
      infection = state$time +log(1 + (decay *  (first(QIRtimes$Q) - cumInf))/(state$e * beta))/decay
    }else{
     invW <- (exp(-1 + (decay * (-decay * (first(QIRtimes$Q) - cumInf) + state$e / beta))/(state$I /beta)) * (decay *state$e - state$I))/state$I
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
  print("That's all folks");
}


##plot the output
### rearrange data for plotting
plot.data <- reshape2::melt(output,id.vars = c("time","run"))
### plot data of each run
ggplot(data = plot.data[plot.data$variable %in% c("L","I","R","e"),])+
  geom_point(aes(x = time, y = value, colour = variable))+
  geom_path(aes(x = time, y = value, colour = variable))+
  scale_colour_manual(name = "Infection",values =c("blue","red","black","grey"))+
  facet_grid(run~.)+xlab("Time")+ylab("Animals") +
  theme_bw() #theme_bw is the layoout
ggsave("out/SEIRe_run.jpg") #save plot to output directory
### plot data of runs together for each variable
ggplot(data = plot.data[plot.data$variable %in% c("L","I","R","e"),])+
  #geom_point(aes(x = time, y = value, colour = variable, group = run),alpha = 0.2)+
  geom_path(aes(x = time, y = value, colour = variable, group = run),alpha = 0.2)+
  scale_colour_manual(name = "Infection",values =c("blue","red","black","grey"))+
  facet_grid(variable~.)+xlab("Time")+ylab("Animals") +
  theme_bw() #theme_bw is the layoout
ggsave("out/SEIRe_var.jpg") #save plot to output directory

