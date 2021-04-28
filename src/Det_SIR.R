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
    dx <- -1 * beta * x * y  +mu * y + mu * z
    dy <- beta * x * y / sum(State)  - gamma * y - mu * y
    dz <- gamma * y - mu * z
    return(list(c(dx,dy,dz)))
  })
}

#Parameters
pars <- c(beta = 1.25, #beta = transmission rate
          gamma = 0.25, #gamma = recovery rate
          mu = 0
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

#exponential growth phase
ggplot()+
  geom_path(aes( x = time, y = value, colour = variable),data = plot.data)+
  #geom_path(aes( x = time, y = value, colour = variable), size = 2,data = plot.data[plot.data$time < exp.phase,])+
  #geom_path(aes(x = c(30,30)+.4, y=c(0,final.size)), arrow = arrow(length=unit(0.25,"cm"),ends = "last",type = "open"))+
  #geom_line(aes(x = c(time.peak,20), y=c(peak,peak)), arrow = arrow(length=unit(0.25,"cm"),ends = "first",type = "open"))+
  #annotate("rect", xmin = 0, xmax = exp.phase, ymin = 0, ymax = 1,
   #        alpha = .2)+
  scale_y_continuous(
    
    sec.axis = sec_axis(~ . , name = "Fraction of population")
  )+
  scale_colour_manual(values = cbbPalette[c(1,7,3)])+
  xlim(c(0,1.5))+
  ylim(c(0,0.06))+
  xlab("Time")+
  ylab("Fraction of population")+
 # ggtitle("Epidemic curve")+
  theme(legend.position = "none")+
  #theme_bw()

#xy plane with x and y peak####
ggplot()+
  geom_path(aes( x = x, y = y), colour = "black",data =data.frame(ode.out))+
  geom_point(aes( x = x[1], y = y[1]), colour = cbbPalette[2],data =data.frame(ode.out))+
  geom_path(aes( x = x+.05, y = y),size =1.2, colour = cbbPalette[3],arrow = arrow(length=unit(0.25,"cm"),ends = "last",type = "open"),data =data.frame(ode.out)[c(1,25),])+
  geom_path(aes( x = c(0,pars["gamma"]/pars["beta"],pars["gamma"]/pars["beta"]), y =c(peak,peak,0)), colour = cbbPalette[4])+
  geom_text(aes (x = 0.05, y =peak+.02), label = "y_peak" , colour = cbbPalette[4])+
  geom_text(aes (x = pars["gamma"]/pars["beta"], y = 0), label = "x_peak" , colour = cbbPalette[4])+
  geom_text(aes (x = 1.05, y = 0.05), label = "Time" , colour = cbbPalette[3], angle = -60)+
  xlim(-.2,1)+
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

#xy plane with final size####
ggplot()+
  geom_path(aes( x = x, y = y), colour = "black",data =data.frame(ode.out))+
  geom_point(aes( x = x[1], y = y[1]), colour = cbbPalette[2],data =data.frame(ode.out))+
  geom_point(aes( x = 1-final.size, y = 0), colour = cbbPalette[6],data =data.frame(ode.out))+
  geom_path(aes( x = x+.05, y = y),size =1.2, colour = cbbPalette[3],arrow = arrow(length=unit(0.25,"cm"),ends = "last",type = "open"),data =data.frame(ode.out)[c(1,25),])+
  #geom_path(aes( x = c(0,pars["gamma"]/pars["beta"],pars["gamma"]/pars["beta"]), y =c(peak,peak,0)), colour = cbbPalette[4])+
  geom_text(aes (x = 1-final.size+.23, y =.02), label = "1-x_inf = final size" , colour = cbbPalette[6])+
  #geom_text(aes (x = pars["gamma"]/pars["beta"], y = 0), label = "x_peak" , colour = cbbPalette[4])+
  geom_text(aes (x = 1.05, y = 0.05), label = "Time" , colour = cbbPalette[3], angle = -60)+
  xlim(-.2,1)+
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

#final size graph
R0finalsize <- function(zinf){log(1-zinf)/-zinf}
final.size.table <- data.frame(z = c(seq(from= 0, to = .9,.1),.95,.99,.999))
final.size.table$R0 <- sapply(final.size.table, R0finalsize)
final.size.table$R0[is.nan(final.size.table$R0)]<- 0
final.size.table<- rbind(final.size.table,data.frame(z = c(0),R0=c(0)))
final.size.table<- rbind(final.size.table,data.frame(z = c(0),R0=c(1)))
ggplot(data = final.size.table) +
  geom_line(aes(R0,z), size =1.2)+
  ggtitle("Final size")+
  theme(legend.position = "none")+
  theme_bw()

#expected number of infections if R0 < 1 
R0finalsize <- function(zinf){log(1-zinf)/-zinf}
exp.inf <- data.frame(R0 = c(seq(from= 0, to =.99,.05),1), Expected = c(1/(1-c(seq(from= 0, to = .99,.05))),NA))
ggplot(data = exp.inf) +
  geom_line(aes(R0,Expected), size =1.2)+
  ggtitle("R0 < 1")+
  theme(legend.position = "none")+
  theme_bw()


#with demography
#Parameters
pars <- c(beta = 1., #beta = transmission rate
          gamma = 0.3, #gamma = recovery rate
          mu = 0.02
)
times <- seq(0, 150, by = dt)

#Solve the ordinary differential equations
ode.out <- ode(init, times, sir,pars)
#
plot.data <- melt(data.frame(ode.out), "time")
final.size = tail(data.frame(ode.out)$z,1)
peak =max(data.frame(ode.out)$y)
time.peak = data.frame(ode.out)[data.frame(ode.out)$y == peak,]$time
exp.phase = 2.5
#epidemic curve
ggplot()+
  geom_path(aes( x = time, y = value, colour = variable),data = plot.data)+
  #geom_path(aes(x = c(30,30)+.4, y=c(0,final.size)), arrow = arrow(length=unit(0.25,"cm"),ends = "last",type = "open"))+
  #geom_line(aes(x = c(time.peak,20), y=c(peak,peak)), arrow = arrow(length=unit(0.25,"cm"),ends = "first",type = "open"))+
  #annotate("rect", xmin = 0, xmax = exp.phase, ymin = 0, ymax = 1,
   #        alpha = .2)+
  scale_y_continuous(
    
    sec.axis = sec_axis(~ . , name = "Fraction of population")
  )+
  scale_colour_manual(values = cbbPalette[c(1,7,3)])+
  #xlim(c(0,30.5))+
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
