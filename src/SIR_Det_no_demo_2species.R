install.packages("deSolve")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("GGally")
install.packages("network")

library(deSolve)
library(ggplot2)
library(reshape2)
library(GGally)
library(network)

###############################################################################################
#                                                                                             #
#                              Deterministic without demography for 4 species                 #
#                             for an non-lethal disease                                       #                                                             #                                                                                             #
###############################################################################################

#Ordinary differential Equations
sir2 <- function(Time, State, Pars){
 with(as.list(c(State,Pars)),{
  dy1 <- (n1 - y1 -z1) * (beta11 * (y1 / n1) +beta21 * (y2 / n2)) - gamma1 * y1 #change in infections for species 1
  dy2 <- (n2 - y2 -z2) * (beta12 * (y1 / n1) +beta22 * (y2 / n2))- gamma2 * y2 #change in infections for species 2
  dz1 <- gamma1 * y1
  dz2 <- gamma2 * y2
  return(list(c(dy1,dy2,dz1,dz2)))
 })
}

####################################################
#          Parameters                              #
####################################################
pars <- c(#beta = transmission rates between species
          beta11 = 0.8, 
          beta12 = 0.3, #transmission constant from species 1 to 2
          beta21 = 0.2, 
          beta22 = 0.5, 
          #recovery rates
          gamma1 = 1 ,# recovery of species 1 = 1/duration of infection
          gamma2 = 2 
          )

#population sizes
n1 = 2500
n2 = 2500

#initial infections
y1_0 = 1
y2_0 = 0
init <- c(y1 = y1_0,
          y2 = y2_0,
          z1 =0, 
          z2 =0)

##summary of the parameter settings
#if transmission is only within the same species the basic reproduction numbers are:
pars["beta11"]/pars["gamma1"]
pars["beta22"]/pars["gamma2"]

#R0 of the complete system
NGM <- matrix(pars[1:4],nrow = 2) / pars[5:6] #next generation matrix i.e. how many new infections in each species caused by each of the species in the next generations
ev<- eigen(NGM) #calculate the eigenvalues and eigenvectors of the next generation matrix
max(Re(ev$values)) #the largest is R0
#directly solvable for 2x2 NGM
0.5*(pars[1]/pars["gamma1"]+pars[4]/pars["gamma2"]+sqrt(diff(pars[c(1,4)]/pars[c("gamma1","gamma2")])^2+4*prod(pars[c(2,3)]/pars[c("gamma1","gamma2")])))

#plot K22  for which R0 > 1
k22 <- function(k11,k21xk12){ifelse(k11<1,1-(k21xk12)/(1-k11),NA)}
rangek21xK12 = c(0.001,0.01,0.1,0.5,.9);
k11seq = seq(0,2,0.001)
ggplot(data = data.frame(k11 = rep(k11seq,length(rangek21xK12)),
                         k22 = (mapply(FUN =k22,rep(k11seq,length(rangek21xK12)),
                                        rep(rangek21xK12,each = length(k11seq)))),
                         k21xK12 = rep(rangek21xK12,each = length(k11seq))))+
  geom_path(aes(x,y),data = data.frame(x = c(0,2),y = c(1,1)))+
  geom_path(aes(y,x),data = data.frame(x = c(0,2),y = c(1,1)))+
  scale_x_continuous(expand = c(0,0),limits = c(0,2))+
  scale_y_continuous(expand = c(0,0),limits = c(0,2))+
  geom_path(aes(k11,k22,colour = factor(k21xK12)))+
  labs(color = "k21 x K12", x = "K11",y = "K22")  +
  theme_bw()

#############################################
#       simulation                          #
#############################################
#Data storage time
dt = 0.1#timestep for storing data
times <- seq(0, 100, by = dt)


#Solve the ordinary differential equations
ode.out <- ode(init, times, sir2,pars) 
#plot infection against time
plot.data = melt(data.frame(ode.out), id = c("time")) #create plottable data
plot.data$species = substr(plot.data$variable,2,2) #variable indicating the species
plot.data$state = substr(plot.data$variable,1,1) #variable indicating the disease state

ggplot(data = plot.data)+
  geom_path(aes(x = time, y = value, group = variable, colour = species), size =0.8)+
  ylab("Number") + 
  xlab("Time (in days)")+
  facet_grid(state ~.,scales = "free_y")
  


