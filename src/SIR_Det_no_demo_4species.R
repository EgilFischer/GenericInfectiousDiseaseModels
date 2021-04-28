install.packages("deSolve")
install.packages("ggplot2")
install.packages("reshape2")

library(deSolve)
library(ggplot2)
library(reshape2)

###############################################################################################
#                                                                                             #
#                              Deterministic without demography for 4 species                 #
#                             for an non-lethal disease                                       #                                                             #                                                                                             #
###############################################################################################

#Ordinary differential Equations
sir4 <- function(Time, State, Pars){
 with(as.list(c(State,Pars)),{
  dy1 <- (n1 - y1 -z1) * (beta11 * (y1 / n1) +beta21 * (y2 / n2) + beta31 * (y3 /n3) + beta41 * (y4 / n4)) - gamma1 * y1 #change in infections for species 1
  dy2 <- (n2 - y2 -z2) * (beta12 * (y1 / n1) +beta22 * (y2 / n2) + beta32 * (y3 /n3) + beta42 * (y4 / n4))- gamma2 * y2 #change in infections for species 2
  dy3 <- (n3 - y3 -z3) * (beta13 * (y1 / n1) +beta23 * (y2 / n2) + beta33 * (y3 /n3) + beta43 * (y4 / n4))- gamma3 * y3 #change in infections for species 3
   dy4 <- (n4 - y4 -z4) * (beta14 * (y1 / n1) +beta24 * (y2 / n2) + beta34 * (y3 /n3) + beta44 * (y4 / n4))- gamma4 * y4 #change in infections for species 4
   dz1 <- gamma1 * y1
   dz2 <- gamma2 * y2
   dz3 <- gamma3 * y3
   dz4 <- gamma4 * y4
  return(list(c(dy1,dy2,dy3,dy4,dz1,dz2,dz3,dz4 )))
 })
}

####################################################
#          Parameters                              #
####################################################
pars <- c(#beta = transmission rates between species
          beta11 = 1.0, 
          beta12 = 0.1, #transmission constant from species 1 to 2
          beta13 = 0.0, #transmission constant from species 1 to 3
          beta14 = 0.0, #etc.
          beta21 = 0.1, 
          beta22 = 1.0, 
          beta23 = 0.1, 
          beta24 = 0.0, 
          beta31 = 0.0, 
          beta32 = 0.1, 
          beta33 = 1.0, 
          beta34 = 0.1, 
          beta41 = 0.0, 
          beta42 = 0.0, 
          beta43 = 0.1, 
          beta44 = 1.0, 
          #recovery rates
          gamma1 = 1/6 ,# recovery of species 1 = 1/duration of infection
          gamma2 = 1/6 ,
          gamma3 = 1/6 ,
          gamma4 = 1/6 
          )

#population sizes
n1 = 2500
n2 = 2500
n3 = 2500
n4 = 2500

#initial infections
y1_0 = 1
y2_0 = 0
y3_0 = 0
y4_0 = 0
init <- c(y1 = y1_0,
          y2 = y2_0,
          y3 = y3_0,
          y4 = y4_0,
          z1 =0, z2 =0, z3=0, z4=0)

##summary of the parameter settings
#if transmission is only within the same species the basic reproduction numbers are:
pars["beta11"]/pars["gamma1"]
pars["beta22"]/pars["gamma2"]
pars["beta33"]/pars["gamma3"]
pars["beta44"]/pars["gamma4"]

#R0 of the complete system
NGM <- matrix(pars[1:16],nrow = 4) / pars[17:20] #next generation matrix i.e. how many new infections in each species caused by each of the species in the next generations
ev<- eigen(NGM) #calculate the eigenvalues and eigenvectors of the next generation matrix
max(ev$values) #the largest is R0


#############################################
#       simulation                          #
#############################################
#Data storage time
dt = 0.1#timestep for storing data
times <- seq(0, 50, by = dt)


#Solve the ordinary differential equations
ode.out <- ode(init, times, sir4,pars) 
#plot infection against time
plot.data = melt(data.frame(ode.out), id = c("time")) #create plottable data
plot.data$species = substr(plot.data$variable,2,2) #variable indicating the species
plot.data$state = substr(plot.data$variable,1,1) #variable indicating the disease state

ggplot(data = plot.data)+
  geom_path(aes(x = time, y = value, group = variable, colour = species), size =0.8)+
  ylab("Number") + 
  xlab("Time (in days)")+
  facet_grid(state ~.)
  


