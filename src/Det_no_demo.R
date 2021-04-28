library(deSolve)
library(ggplot2)
library(reshape2)
###############################################################################################
#                                                                                             #
#                              Deterministic without demography                               #
#                                                                                             #
###############################################################################################

#Ordinary differential Equations
sir <- function(Time, State, Pars){
 with(as.list(c(State,Pars)),{
  dx <- -1 * beta * x * y 
  dy <- beta * x * y  - gamma * y
  dz <- gamma * y
  return(list(c(dx,dy,dz)))
 })
}

#Parameters
pars <- c(beta = 1.25, #beta = transmission rate
          gamma = 0.25 #gamma = recovery rate
          )
#Data storage time
dt = 0.1#timestep for storing data
times <- seq(0, 50, by = dt)

#initial values
y0 = 0.01 #initial fraction infected 
z0 = 0# initial fraction recovered, 
x0 = 1 - y0 - z0
init <- c(x = x0,y = y0,z = z0)

#Solve the ordinary differential equations
ode.out <- ode(init, times, sir,pars) 
#plot x, y, and z against time
plot(x = ode.out[,1], y = ode.out[,2], type = "l", xlab = "time", ylab = "proportion", col = "purple",lwd =2)
lines(x = ode.out[,1], y = ode.out[,3], type = "l", xlab = "time", ylab = "proportion", col = "red",lwd =2)
lines(x = ode.out[,1], y = ode.out[,4], type = "l", xlab = "time", ylab = "proportion", col = "palegreen1",lwd =2)
legend(x = 40, y = 0.4, c("x","y","z"),,lwd =2, col = c("purple", "red", "palegreen1"))
#plot  x against y
plot(x = ode.out[,2], y = ode.out[,3], type = "l", xlab ="x", ylab ="y")
s <- seq(length(ode.out[,2])-1)
arrows(x0 = ode.out[s,2], y0 = ode.out[s,3], x1 = ode.out[s + 1,2], y1 = ode.out[s + 1 ,3],length = .075)

#pretty pictures
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#
plot.data <- melt(data.frame(ode.out), "time")
final.size = tail(data.frame(ode.out)$z,1)
peak =max(data.frame(ode.out)$y)
time.peak = data.frame(ode.out)[data.frame(ode.out)$y == peak,]$time
exp.phase = 2.5
#epidemic curve
ggplot()+
  geom_path(aes( x = time, y = value, colour = variable),data = plot.data)+
  geom_path(aes( x = time, y = value, colour = variable), size = 2,data = plot.data[plot.data$time < exp.phase,])+
  geom_path(aes(x = c(30,30)+.4, y=c(0,final.size)), arrow = arrow(length=unit(0.25,"cm"),ends = "last",type = "open"))+
  geom_line(aes(x = c(time.peak,20), y=c(peak,peak)), arrow = arrow(length=unit(0.25,"cm"),ends = "first",type = "open"))+
  annotate("rect", xmin = 0, xmax = exp.phase, ymin = 0, ymax = 1,
           alpha = .2)+
  scale_y_continuous(
     
    sec.axis = sec_axis(~ . , name = "Fraction of population")
  )+
  scale_colour_manual(values = cbbPalette[c(1,7,3)])+
  xlim(c(0,30.5))+
  xlab("Time")+
  ylab("Fraction of population")+
  ggtitle("Epidemic curve")+
  theme(legend.position = "bottom")+
  theme_bw()
#different phse (Re>1,Re<1)
ggplot()+
  geom_path(aes( x = time, y = value, colour = variable),data = plot.data)+
  #geom_path(aes( x = time, y = value, colour = variable), size = 2,data = plot.data[plot.data$time < exp.phase,])+
  #geom_path(aes(x = c(30,30)+.4, y=c(0,final.size)), arrow = arrow(length=unit(0.25,"cm"),ends = "last",type = "open"))+
  #geom_line(aes(x = c(time.peak+5,20), y=c(peak,peak)), arrow = arrow(length=unit(0.25,"cm"),ends = "first",type = "open"))+
  annotate("rect", xmin = 0, xmax = time.peak, ymin = 0, ymax = 1,
           alpha = .2)+
  scale_y_continuous(
    
    sec.axis = sec_axis(~ . , name = "Fraction of population")
  )+
  scale_colour_manual(values = cbbPalette[c(1,7,3)])+
  xlim(c(0,30.5))+
  xlab("Time")+
  ylab("Fraction of population")+
  ggtitle("Epidemic curve")+
  theme(legend.position = "bottom")+
  theme_bw()

#xy plane####
ggplot()+
  geom_path(aes( x = x, y = y), colour = "black",data =data.frame(ode.out))+
  geom_point(aes( x = x[1], y = y[1]), colour = cbbPalette[2],data =data.frame(ode.out))+
  geom_path(aes( x = x+.05, y = y),size =1.2, colour = cbbPalette[3],arrow = arrow(length=unit(0.25,"cm"),ends = "last",type = "open"),data =data.frame(ode.out)[c(1,25),])+
  geom_text(aes (x = 1.05, y = 0.05), label = "Time" , colour = cbbPalette[3], angle = -60)+
  ylim(0,.6)+
  scale_y_continuous(
    
    sec.axis = sec_axis(~ . , name = "proportion infectious")
  )+
  scale_x_continuous(
    
    sec.axis = sec_axis(~ . , name = "proportion susceptible")
  )+
   ggtitle("xy plane")+
  
  theme(legend.position = "none")+
  theme_bw()



#xy plane with x and y peak####
ggplot()+
  geom_path(aes( x = x, y = y), colour = "black",data =data.frame(ode.out))+
  geom_point(aes( x = x[1], y = y[1]), colour = cbbPalette[2],data =data.frame(ode.out))+
  geom_path(aes( x = x+.05, y = y),size =1.2, colour = cbbPalette[3],arrow = arrow(length=unit(0.25,"cm"),ends = "last",type = "open"),data =data.frame(ode.out)[c(1,25),])+
  geom_point(aes( x = c(1), y =c(peak)), colour = cbbPalette[2])+
  geom_text(aes (x = 1.05, y = 0.05), label = "Time" , colour = cbbPalette[3], angle = -60)+
  scale_y_continuous(
    
    sec.axis = sec_axis(~ . , name = "proportion infectious")
  )+
  scale_x_continuous(
    
    sec.axis = sec_axis(~ . , name = "proportion susceptible")
  )+
  ggtitle("xy plane")+
  ylim(0,.6)+
  theme(legend.position = "none")+
  theme_bw()
