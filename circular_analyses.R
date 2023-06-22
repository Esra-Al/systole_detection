rm(list=ls())
require(circular)
library(plotrix)
library(ggplot2)
library(extrafont)
main_path='~/Desktop/systole_detection/' # change this to the path of systole_detection folder
setwd(main_path)

ecgdata=read.csv2("filtecg_final.csv")
ecgdata=na.omit(ecgdata)
ecgdata$sysdegree=360*(ecgdata$syslength/ecgdata$RRinterval) #the average length of systole and diastole

#PARTICIPANTS
# Loop through all participants and calculate hits distribution

hit_sys = matrix(data=NA,nrow=41,ncol=1) 
miss_sys = matrix(data=NA,nrow=41,ncol=1)

mean_degree_hit = matrix(data=NA,nrow=41,ncol=1)
mean_degree_miss = matrix(data=NA,nrow=41,ncol=1)

ID_list = unique(ecgdata$subject) 

for (s in ID_list) {
  hit_sys[s]=mean(circular(ecgdata$sysdegree[ecgdata$subject==s & abs(ecgdata$stimuli)==1 & ecgdata$detans==1], type="angles", units="degree", rotation="clock", zero=0))
  miss_sys[s]=mean(circular(ecgdata$sysdegree[ecgdata$subject==s & abs(ecgdata$stimuli)==1 & ecgdata$detans==0], type="angles", units="degree", rotation="clock", zero=0))
  
  mean_degree_hit[s]=mean(circular(ecgdata$stim_degree[ecgdata$subject==s & abs(ecgdata$stimuli)==1 & ecgdata$detans==1], type="angles", units="degree", rotation="clock", zero=0))
  mean_degree_miss[s]=mean(circular(ecgdata$stim_degree[ecgdata$subject==s & abs(ecgdata$stimuli)==1 & ecgdata$detans==0], type="angles", units="degree", rotation="clock", zero=0))
 }


hitsys2= circular(hit_sys[ID_list], type="angles", units="degree", rotation="clock", zero=pi/2)
Hits_secondlevel <- circular(mean_degree_hit[ID_list], type="angles", units="degree", rotation="clock", zero=pi/2)
#svg(filename = "Hits_miss_circular.svg", width=5,height=3)
par(mfrow=c(1,2), mai = c(0.1, 0.2, 0.3, 0.3))
plot(Hits_secondlevel, stack=TRUE, bins = 720, col = "gray47", cex = 0.77, lwd = 2)
draw.arc(x=0,y=0,radius=1, deg1=90,deg2=90-mean(hitsys2)[[1]], col="red")
draw.arc(x=0,y=0,radius=1, deg1=360-mean(hitsys2)[[1]],deg2=90, col="blue")
arrows.circular(mean(Hits_secondlevel), y=rho.circular(Hits_secondlevel),lwd = 3.3,  col = "gray22", length = 0.08)
circ.dens = density(Hits_secondlevel, bw=20)
lines(circ.dens, col="gray", lwd = 2, xpd=TRUE)
mtext('a', side=3, at=-1, font =2, family=fonts()[7], cex=2)

misssys2= circular(miss_sys[ID_list], type="angles", units="degree", rotation="clock", zero=pi/2)
Miss_secondlevel <- circular(mean_degree_miss[ID_list], type="angles", units="degree", rotation="clock", zero=pi/2)
plot(Miss_secondlevel, stack=TRUE, bins = 720, col = "gray47", cex = 0.77, lwd = 2)
draw.arc(x=0,y=0,radius=1, deg1=90,deg2=90-mean(misssys2)[[1]], col="red")
draw.arc(x=0,y=0,radius=1, deg1=360-mean(misssys2)[[1]], deg2=90, col="blue")
arrows.circular(mean(Miss_secondlevel), y=rho.circular(Miss_secondlevel), lwd = 3.3,  col = "gray22", length = 0.08)
circ.dens = density(Miss_secondlevel, bw=20)
lines(circ.dens, col="gray", lwd = 2, xpd=TRUE)
mtext('b', side=3, at=-1, font =2, family=fonts()[7], cex=2)
#dev.off()

rayleigh.test(Hits_secondlevel)
mean(Hits_secondlevel)[[1]]+360

rayleigh.test(Miss_secondlevel)
mean(Miss_secondlevel)[[1]]

