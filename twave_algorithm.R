# The below code for implementation of T-wave end algorithm is written by Esra Al.

require(R.matlab)
rm(list = ls())
main_path='~/Desktop/systole_detection/' # change this to the path of systole_detection folder
setwd(main_path) # set working directory
files <- list.files(pattern = ".mat") #the list of valid blocks of subjects
bigdata=c() #initialize an array to store all trials from each valid block of every subject
# specify ECG sampling frequency (in Hz)
fs=2500 #sampling frequency of ECG
numtrials=120
s=1  #subject 1
b=1  #block 1
#for (s in 1:41){ #loop for subjects
  subdata=c()
  #for (b in 1:8){ #loop for blocks
    dataname=paste('data',s,'_',b, '.mat', sep='')
    data=c() #initialize array to store all valid blocks of each subject
    if (sum(files==dataname)==1){
      #load data containing behavior and stimulus onset latencies
      data=readMat(dataname)
      data=data[[1]]
      # load ECG data that was filtered 1) High-pass 0.5 Hz, and 2) low pass 30Hz with a 4th order of Butterworth filter
      ECGname=paste(main_path, 'filtECGVP', s,'_',b, '.txt', sep='')
      ECG=read.table(ECGname)
      
      #load R peak locations
      kubiosname=paste(main_path,'Kubios_VP', s,'_',b, '_hrv.mat', sep='')
      R_peaks= readMat(kubiosname)
      R_peaks= as.vector(R_peaks$Res[[4]][[2]][[2]]) 
      
      
      data=cbind(matrix(data=NA,nrow=numtrials,ncol=3), data, matrix(data=NA,nrow=numtrials,ncol=7))
      for (ind in 1:numtrials) { 
        # calculate the position of R peak that occured just before the stimulus onset
        pos=max(which(R_peaks < data[ind,8]))
        
        # calculate the position of R peak that occured just before the stimulus onset in data points in the R_peaks vector
        ecgpos1=R_peaks[pos]*fs 
        
        # calculate the duration of RR interval containing the stimulus in data points
        RRint=(R_peaks[pos+1]-R_peaks[pos])*fs #in data points
        
        # Calculate a time interval where the end of t wave would be included (in data points). Cut off 60ms or fs*(150/2500) data points so that t wave max calculation would not detect the next R peak.    
        ecgpos2=ecgpos1+RRint-150   
        
        # create a vector containing all possible data points in the RR interval for visualization purposes. 120 ms or fs*(300/2500) data points are chosen arbitrarily to include the Rpeaks in the plots.
        twave1=ECG[(ecgpos1-300):ecgpos2+300,]
        
        # create a vector (twave) where the maximum of t wave could be calculated. After visually observing a block of each participant, it is decided that tmax should occur at least after 140 ms or (350/2500)*fs datapoints after the R peak. 
        twave=ECG[(ecgpos1+350):ecgpos2,]
        
        #In order to calculate the position of tmax, delete 140 ms or (350/2500)*fs data points from the end of RRint vector 
        #and search it until the 1/3 of the remaining vector. 1/3 ratio was given to compensate for the variation of trial to trial RR interval durations.
        
        tmaxpos=which.max(twave[1:((RRint-350)/3),2]) 
        #calculate a snip of t wave which starts with the tmax
        twave2=twave[tmaxpos:dim(twave)[1],]
        par(mfrow=c(1,2)) #for creating two subplots in a figure
        
        # Determine a point called xm located in the segment after the T peak, which has a minimum value in the first derivative. 
        #The algoritm searces for xm in a 120 ms time window startng from tmax. In case twave2 does not contain 0.12*fs data points, it searches only until the last point of twave2.
        dp=0.12*fs 
        if (dp>dim(twave2)[1]) {
          xm=which(diff(twave2[,2])==min(diff(twave2[,2]))) 
        } else {
          xm=which(diff(twave2[1:dp,2])==min(diff(twave2[1:dp,2]))) 
        }
        xm=xm[1]
        ym=twave2[xm,2]
        
        # determine a point xr which is supposed to happen after the end of twave, I chose it to be 80 ms , if your sampling rate is different than 2500 Hz use: (200/2500)*fs 
        xr=200+xm 
        
        # make a vector starting from xm and goes until xr
        xseq=xm:xr
        yseq=twave2[xm:xr,2]
        
        # write a function find the end of twave: first calculation of the trapeziums areas of all the points located between xseq and yseq,then identification of the point which gives the maximum area and label it as the t-wave end
        
        trapez_area <- function(xm, ym, xseq, yseq, xr) {
          a <- numeric()
          for (i in seq_along(xseq)){
            a[i] <- 0.5 * (ym - yseq[i]) * ((2*xr) - xseq[i] - xm)
          }
          x_tend <- which.max(a)+xm-1
          return(x_tend)
        }
        tend=trapez_area(xm, ym, xseq, yseq, xr)
        
        
        nextR=(R_peaks[pos+1]*fs)-10 # here delete 4 ms from the next R peak (to prevent labeling a stimulus onset occuring during the R peak as both systole and diastole)
        syslen=twave2[tend,1]-R_peaks[pos] #length of systole
        diastole_dp=(nextR-syslen*fs):(nextR-1)
        diastole_ecg=ECG[diastole_dp,]
        diastole_ecg[1,1]
        
        #To check whether twave end detection worked well
        
        jpeg(file = paste('subject ',s,'_block ',b,'_trial ', ind,".jpg"))
        par(mfrow=c(1,2))
        plot(twave1,col='black',xlab='time(ms)', ylab='electrical potential(mV)')
        points(twave[tmaxpos,1],twave[tmaxpos,2],col='magenta',pch='+',cex=2)
        points(twave2[tend,1],twave2[tend,2],col='green',pch='+',cex=2)

        plot(twave2,col='black', xlab='time(ms)', ylab='electrical potential(mV)')
        title(paste('subject ',s,'_block ',b,'_trial ', ind, sep=''),line=-2, outer=TRUE)
        points(twave2[xm,1],twave2[xm,2],col='blue',pch='+',cex=2)
        points(twave2[xr,1],twave2[xr,2],col='blue',pch='+',cex=2)
        points(twave2[tend,1],twave2[tend,2],col='green',pch='+',cex=2)
        points(twave[tmaxpos,1],twave[tmaxpos,2],col='magenta',pch='+',cex=2)
        title(paste('subject ',s,'_block ',b,'_trial ', ind, sep=''),line=-2, outer=TRUE)
        dev.off()
        # 
        # the difference between the onset and previous R peak
        diff2peak=data[ind,8] - R_peaks[pos]
        # relative position of onset in cardiac cycle
        stim_degree=360 * diff2peak/(R_peaks[pos+1] - R_peaks[pos]) # this part is for circular analyses
        data[ind,1]=s
        data[ind,2]=b
        data[ind,3]=ind
        data[ind,9]=stim_degree
        data[ind,10]=(R_peaks[pos+1] - R_peaks[pos]) #RR interval
        data[ind,11]=twave2[tend,1] # end of twave in s
        data[ind,12]=twave2[tend,1]-R_peaks[pos] # systole length
        data[ind,13]=data[ind,8]<twave2[tend,1] # whether the stimulus onset is in systole 
        data[ind,14]=data[ind,8]>diastole_ecg[1,1] # whether in diastole
        data[ind,15]=diff2peak
      }
      #subdata=rbind(subdata,data)
    }
 # }
  #bigdata=rbind(bigdata,subdata)
#}
#bigdata=data.frame(bigdata)
#names(bigdata)=c('subject','block','trial','stimuli','detans','locans','reacttime','stimonset','stim_degree','RRinterval','tend', 'syslength','systole','diastole','diff2peak')
#write.csv2(bigdata,file= "C:/Users/alesra/Desktop/analyses/ecg/cardiac_phase/filtecg_equalsysdys.csv")
#####
