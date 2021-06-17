#GNU General Public License v3.0
#Copyright (c) 2021 Sabine Schilling
#Feedback highly welcome: sabineschilling@gmx.ch

# clear variables and close graphics----
graphics.off() #
rm(list = ls()) #deleting all variables from work space

## load library-------
library(multcompView)

#load statistics routine adapted from package visStatistics----
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
experimental_inputfile=paste(csvdirec,"/raw_data_endothelial_Frontiers.csv",sep="") #adapt this line to your directory
endo=read.csv(experimental_inputfile,header=TRUE,sep=";",stringsAsFactors = F)
listm=c(2,4,7,8) #column numbers in experimental data:
#2: area
#4: circularity=4*pi area/perimeter^2
#7: aspect ratio=major/minor axis
#8: angle

## Create column with words instead of number codes for graphical output ------
endo$FlowWords=endo$Flow
endo$FlowWords[which(endo$Flow==1)]="30 dynes/cm^2"
endo$FlowWords[which(endo$Flow==2)]="2 dynes/cm^2"
endo$FlowWords[which(endo$Flow==3)]="80 dynes/cm^2"
endo$FlowWords=as.factor(endo$FlowWords)
names_for_plots=c("area in squared pixels" ,"circularity","aspect ratio","angle")
#names_for_plots=c("")

##Angle---------
## change angle for 2 dynes and 30 dynes/cm^2, correcting for different plate orientation in these flow conditions
endo$Angle[which(endo$cellComment=="change")]=endo$Angle[which(endo$cellComment=="change")]-90
#correct for negative angles----
endo$Angle[which(endo$Angle<0)]=180-endo$Angle[which(endo$Angle<0)]
endo$Angle=projectAngleToFirstQuadrant(endo$Angle)

#Sort data by FlowWords----
endosort <- endo[order(endo$FlowWords),]

j=0
for (i in listm){
  #Wilcoxon for all cases
  j=j+1
  openGraph()
  res_krus_3_conditions=vis_Kruskal_Wallis_clusters(endosort[,i],endosort$FlowWords,alpha=0.05,xlab="",ylab=names_for_plots[j],cex=1,notch=F)
  saveGraph(paste(figdirec, "/Kruskal_3_conditions",variable.names(endosort[i]),sep=""),type=figtype)
}

