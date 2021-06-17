graphics.off() #
rm(list = ls()) #deleting all variables from work space

## load libraries-------
library(rstudioapi)
library(multcompView)
#load statistics routine adapted from package visStatistics
# https://cran.r-project.org/web/packages/visStatistics/index.html
source("./scripts/Kruskal_Wallis_visualization.R")


## define directories ------
codedirec=c("./scripts")
figdirec=c("./figures") #directory of figures
csvdirec=c("./data") # directory of data

## global graphics parameters ------
cexsize = 1  #variable defining size of cexsize
figtype="pdf" #type of figures to get generated
lwdsize=1
background="white" #figure background
foreground="black"#figure foreground
#set colors and size of axis for all plots----
par(bg=background,fg=foreground,col.lab=foreground,col.axis=foreground,
    cex.axis=cexsize,cex.lab=cexsize, lwd=lwdsize, oma = c(0, 0, 2, 0),mfrow = c(1, 1),
    xaxs="i",family = "Times New Roman")
default_graphical_parameters=par() #makes copy of current settings

## load experimental data ------
experimental_inputfile=paste(csvdirec,"/raw_data_endothelial.csv",sep="") #adapt this line to your directory
endo=read.csv(experimental_inputfile,header=TRUE,sep=",",stringsAsFactors = F)
listm=c(6,23,24,35,37,38) #column numbers in experimental data sheet from FIJI
#6: area
#23: Angle
#24: circularity=4*pi area/perimeter%2
#35: Feret angle
#37: Aspect ratio
#38 Roundness: =4*pi*area/(pi*major axis)^2:

## Create column with words instead of number codes for graphical output ------
endo$FlowWords[which(endo$Flow==1)]="30 dynes/cm^2"
endo$FlowWords[which(endo$Flow==2)]="2 dynes/cm^2"
endo$FlowWords[which(endo$Flow==3)]="80 dynes/cm^2"
endo$FlowWords=as.factor(endo$FlowWords)
names_for_plots=c("area in squared pixels" ,"angle","circularity","angle gamma","aspect ratio","roundness")
#names_for_plots=c("")

##Angle------
## change angle for 2 dynes and 30 dynes/cm^2, correcting for different plate orientation in these flow conditions (change: April 2017)
endo$Angle[which(endo$cellComment=="change")]=endo$Angle[which(endo$cellComment=="change")]-90
endo$FeretAngle[which(endo$cellComment=="change")]=endo$FeretAngle[which(endo$cellComment=="change")]-90
#correct for negative angles----
endo$FeretAngle[which(endo$FeretAngle<0)]=180-endo$FeretAngle[which(endo$FeretAngle<0)]
endo$Angle[which(endo$Angle<0)]=180-endo$Angle[which(endo$Angle<0)]


endo$Angle=projectAngleToFirstQuadrant(endo$Angle)
endo$FeretAngle=projectAngleToFirstQuadrant(endo$FeretAngle)
endo$flow=endo$FlowWords
#make flow subsets #much shorter, splits by all possible Flows
 endoflow=split(endo,f=endo$Flow)
 names(endoflow)

flownamesreduced_=c("2 dynes/cm^2","30 dynes/cm^2","80 dynes/cm^2")

#name change for density plots
names(endo)[2] <- "flow"
names(endo)[3] <- "flowNumeric"
names(endo)[40] <- "flowWords"

# Plotting------
graphics.off()
j=0
for (i in listm){
  #Wilcoxon for all cases
  # setwd(figdirec)
  j=j+1

  openGraph()
  res_krus_3_conditions=vis_Kruskal_Wallis_clusters(endo[,i],endo$flow,alpha=0.05,xlab="",ylab=names_for_plots[j],cex=1,notch=F)
  saveGraph(paste(figdirec, "/Kruskal_3_conditions",variable.names(endo[i]),sep=""),type=figtype)

}
