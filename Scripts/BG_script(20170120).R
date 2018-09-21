#..........................................#
#...........Blue-Green project.............#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#... Contact: eric.harvey@eawag.ch                                                                                           #                                                                       #
#..........................................................................................................................................#


###THIS IS AN OLD SCRIPT THAT WAS SEPARATED INTO ITS USEFUL COMPONENTS INTO DIFFERENT SCRIPTS

#########################################################################
################# DATA STRUCTURE
#########################################################################

##################
#Clear any variables from the R environment 
##################
rm(list=ls())

##################
#Directory paths 
##################
to.data = "./data/"
to.script = "./scripts/"
to.output = "./output/"
to.figs = "./figs/"
to.R = "./output/"


##################
#Load packages
##################
library(tidyverse)

# library(vegan)
# library(plyr)
# library(plotrix)
# library(sciplot)
# library(car)
# library(MASS)

##################
#Load data in R environment
##################
Prot.b0 =  read_csv(paste0(to.data,"BG_protist_data_(20170307).csv")) #Protist
Cyto0 = read_tsv(paste0(to.data,"FULL_Cyto_BlueGreen(20160816).txt")) # Bacteria
Prot.b = Prot.b0
Cyto = Cyto0

##################
#Bacteria data
##################

#.Add experimental day column 
Cyto = Cyto[order(Cyto$Date),] #make sure that data are ordered by dates
Cyto$day = c(rep(0,32),rep(7,524),rep(15,524),rep(21,524),rep(29,524)) #create a new variable "day"

#...Add a bacteria density per mL column 
Cyto$density = ((Cyto$Count*1000)/Cyto$Volume.uL)*1000 #convert from dens/50ul to dens/mL and then multiply by the cytometry dilution factor (1000)

##################
#Protist data
##################

#...Rename dates consistently (they are not)
Prot.b$date = with(Prot.b,revalue(date, c("16-05-02"="20160502", "16-05-23"="20160523","16-05-31"="20160531","5/17/2016"="20160517","5/9/2016"="20160509")))
Prot.b$date = factor(Prot.b$date,levels=c("20160502","20160509","20160517","20160523","20160531"))

#...Create a new variable called day
Prot.b = Prot.b[order(Prot.b$date),] #make sure that data are ordered by dates 
Prot.b$day = c(rep(0,30),rep(7,302),rep(15,302),rep(21,302),rep(29,302))

#...Add T0 to the data 
{  
#1.Extract and Repeat pooled cultures from Day 0 for each patch size
#Protist
Prot.day0 = Prot.b[rep(0:16,times=4),] #so that we have each patch size for both connected and isolated patches
row.names(Prot.day0) = 1:64
Prot.day0$Size[1:16] = 7.5
Prot.day0$Size[17:32] = 13
Prot.day0$Size[33:48] = 22.5
Prot.day0$Size[49:64] = 45

#Bacteria (for bacteria only pooled cultured from blue landscapes, connected and isolated, were measured [see original data] for logistic reasons, we assumed that our initial bacteria meausres are representative of the entire intial conditions for all treatments)
#Repeat for each patch size
Cyto.day0 = Cyto[rep(0:16,times=4),]
row.names(Cyto.day0) = 1:64
Cyto.day0$Size[1:16] = 7.5
Cyto.day0$Size[17:32] = 13
Cyto.day0$Size[33:48] = 22.5
Cyto.day0$Size[49:64] = 45

#2.Remove Day 0 pooled cultures from the original dataset but keep day 0 monoculture info
#Protist
sel.date = c("20160509","20160517","20160523","20160531")
Prot.b = Prot.b[which(Prot.b$date %in% sel.date | Prot.b$Landscape=="monoculture"),]
#Bacteria
Cyto = Cyto[which(Cyto$Date %in% sel.date | Cyto$Landscape=="Monoculture"),]

#3. Merge the new 'Day 0' duplicated for each patch size with dataset
#Protist
Prot.b = rbind(Prot.day0,Prot.b)
str(Prot.b)
row.names(Prot.b) = 1:1286
#data.frame':	1286 obs. of  31 variables:

#Bacteria
Cyto = rbind(Cyto.day0,Cyto)
str(Cyto)
#data.frame':	2176 obs. of  20 variables:
row.names(Cyto) = 1:2176

} 


#...Convert species abundance per ul into abundance per ML
species = c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other")
Prot.b[,species] = Prot.b[,species]*1000

#... add column for total protist abundance per mL per patch
Prot.b$abundance = apply(Prot.b[,species],1,sum)

#... add column for bioarea per ml
Prot.b$bioarea = Prot.b$bioarea_per_ul*1000

#...Calculate total protist abundance per patch (as opposed to abundance per volume [density])
new.COL.names = c("Rot.all","Spi.all","Ble.all","Pca.all","Col.all","Chi.all","Tet.all","Other.all","bioarea.all")
COL.names = c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other","bioarea")
for(i in 1:length(COL.names)){
  Prot.b[,new.COL.names[i]] = 0
  Prot.b[,new.COL.names[i]] = Prot.b[,COL.names[i]]*Prot.b$Size
}

#... add species richness per patch column
species = c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other")
Prot.b$richness = apply(Prot.b[,species]>0,1,sum)

#... add community eveness per patch column 
simpson = function(x){return(sum((x/sum(x))^2))}; 
Prot.b$simpson = apply(Prot.b[,species],1,simpson)

##################
#Specify factors and Drop unused levels
##################
Prot.b$Size = factor(Prot.b$Size)
Prot.b = droplevels(Prot.b)

#...Fix structure and drop unused levels
Cyto$Date = factor(Cyto$Date)
Cyto$Size = factor(Cyto$Size)
Cyto = droplevels(Cyto)


#########################################################################
################# Exploratory FIGURES 
#########################################################################

##################
# Protist each species total abundance per treatment (connected/isolated) and per patch size
##################
{
#Subset data of interest for the figure
x = subset(Prot.b, Prot.b$Treatment %in% c("Connected","Isolated") & Prot.b$Landscape %in% c("Blue"))
x = droplevels(x)
#Species to plot on each graph
ab.all = c("Rot.all","Spi.all","Ble.all","Pca.all","Col.all","Chi.all","Tet.all","Other.all")
#Colors for each species
col.p = rainbow(8)

pdf(paste(graphpath,"Prot_Blue_Treatment_Size.pdf"),width=8,height=8)

for(i in 1:nlevels(x$Treatment)){
  
  for(j in 1:nlevels(x$Size)) { 
  
  max = max(x[x$Treatment==levels(x$Treatment)[i] & x$Size==levels(x$Size)[j],ab.all])
    
    for(k in 1:length(ab.all)){ 
  with(x[x$Treatment==levels(x$Treatment)[i] & x$Size==levels(x$Size)[j],],lineplot.CI(day,eval(parse(text= paste(ab.all[k]))),col=col.p[k],ylim=c(0,max),legend=T,xlab="Time",ylab="Protist total abundance"))
  par(new=T)
            }
  title(main=paste(levels(x$Treatment)[i],levels(x$Size)[j],"mL"))
  col.labels=c("Rot","Spi","Ble","PCA","Col","Chi","Tet","Oth")
  testcol = col.p
  color.legend(4.5,10000,5,max,col.labels,testcol,cex=1,gradient="y")
  #dif between x1 and xr is the width of the rectangle
} }

dev.off()
}

##################
# Protist temporal community patterns (aggregate properties) across patch sizes
##################
{
#colors for each path size
col.p = rainbow(n=4,start=0.1)
#Subset data of interest for the figure
x = subset(Prot.b, Prot.b$Treatment %in% c("Connected","Isolated") & Prot.b$Landscape %in% c("Blue"))
x = droplevels(x)
#Metrics of interest 
metrics = c("abundance","richness","simpson")
#Grouping factors
grouping = c("Size","Treatment")
#Figure
pdf(paste(graphpath,"Prot_Blue_DIV_AB_TIME.pdf"),width=8,height=8)
for(i in 1:length(metrics)){
  for(j in 1:length(grouping)) {
    with(x,lineplot.CI(day,eval(parse(text= paste(metrics[i]))),eval(parse(text= paste(grouping[j]))),col=col.p,legend=T,xlab="Time",ylab=metrics[i]))
  }
}
dev.off()
}

##################
# Bacteria density patterns over time and by patch size
##################
{

#Subset data of interest for the figure
x = subset(Cyto, Cyto$Treatment %in% c("Connected","Isolated") & Cyto$Replicate %in% c("A","B"))
x = droplevels(x)
  
pdf(paste(graphpath,"BACT_AB_SIZE_LANDSCAPE_FLOW_TIME.pdf"),width=8,height=8)
for(i in 1:nlevels(x$Treatment)){
  
  for(j in 1:nlevels(x$Landscape)) {
    
    max = max(x[x$Treatment==levels(x$Treatment)[i] & x$Landscape==levels(x$Landscape)[j],20])
    
    with(x[x$Treatment==levels(x$Treatment)[i] & x$Landscape==levels(x$Landscape)[j],],lineplot.CI(day,density,Size,col="black",ylim=c(0,max),legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
    title(main=paste(levels(Cyto$Treatment)[i],levels(Cyto$Landscape)[j]))

  }
}
  dev.off()
  }

##################
# Bacteria density temporal patterns across treatments (connected vs. isolated) and landscape types
##################
{

#Subset data of interest for the figure
x = subset(Cyto, Cyto$Treatment %in% c("Connected","Isolated") & Cyto$Replicate %in% c("A","B"))
x = droplevels(x)  

pdf(paste(graphpath,"BACT_AB_LANDSCAPE_FLOW_TIME.pdf"),width=8,height=8)
  
with(x[x$Landscape=="Blue",],lineplot.CI(day,density,Treatment,legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
title(main="Blue")
with(x[x$Landscape=="Green",],lineplot.CI(day,density,Treatment,legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
title(main="Green")

with(x[x$Treatment=="Connected",],lineplot.CI(day,density,Landscape,legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
title(main="Connected")
with(x[x$Treatment=="Isolated",],lineplot.CI(day,density,Landscape,legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
title(main="Isolated")
dev.off()
}

##################
#Protist beta-diversity (need to be re-done)
##################
{
sp.mat = Prot.b[Prot.b$day!=0,32:39]
sp.mat[899,1:8] = 0.1
sp.mat[1073,1:8] = 0.1
env = Prot.b[Prot.b$day!=0,c(2,5,8)]

MDS.mod1 = metaMDS(decostand(sp.mat,"hell"), k=3,autotransform=FALSE,distance="bray")
pdf(paste(graphpath,"PROT_BETADIV_FLOW.pdf"),width=8,height=8)
ordiplot(MDS.mod1,type="n",choices=c(1,2),main=NULL,xlab="", ylab="")
ordihull(MDS.mod1,groups=env$Treatment,show.groups="Isolated",col="gray88",label=T,lwd=3,lty=1)
ordihull(MDS.mod1,groups=env$Treatment,show.groups="Connected",col="black",label=T,lwd=3,lty=2)
text(MDS.mod1, display="species", col="black",cex=0.5)
#points(MDS.mod1, display="sites", col=c("red"),cex=0.5)
dev.off()
}


#########################################################################
################# Dendritic network
#########################################################################

#Build the dendritic landscapes from connectivity matrices and calculate network metrics of interest
library(igraph)
source(paste0(to.script,"Network_metric_bg.R"))

#Merged network metrics to data

#.....Protist
net.met = c("degree","dist","diam")
for(i in 1:length(net.met)){
  #Protist
  Prot.b[,net.met[i]] = NA
  Prot.b[which(Prot.b$Landscape == "Blue" & Prot.b$day!=0),net.met[i]] = rep(c(eval(parse(text= paste("RepA",net.met[i],sep="."))),eval(parse(text= paste("RepB",net.met[i],sep="."))),eval(parse(text= paste("RepC",net.met[i],sep="."))),eval(parse(text= paste("RepD",net.met[i],sep=".")))),8)
  #Bacteria
  Cyto[,net.met[i]] = NA
  Cyto[which(Cyto$Landscape == "Green" & Cyto$Treatment == "Isolated" & Cyto$day!=0 & Cyto$Replicate %in% c("A","B")),net.met[i]] = rep(c(eval(parse(text= paste("RepA",net.met[i],sep="."))),eval(parse(text= paste("RepB",net.met[i],sep=".")))),4)
  Cyto[which(Cyto$Landscape == "Green" & Cyto$Treatment == "Connected" & Cyto$day!=0 & Cyto$Replicate %in% c("A","B","C","D")),net.met[i]] = rep(c(eval(parse(text= paste("RepA",net.met[i],sep="."))),eval(parse(text= paste("RepB",net.met[i],sep="."))),eval(parse(text= paste("RepC",net.met[i],sep="."))),eval(parse(text= paste("RepD",net.met[i],sep=".")))),4)
  Cyto[which(Cyto$Landscape == "Blue" & Cyto$day!=0 & Cyto$Replicate %in% c("A","B","C","D")),net.met[i]] = rep(c(eval(parse(text= paste("RepA",net.met[i],sep="."))),eval(parse(text= paste("RepB",net.met[i],sep="."))),eval(parse(text= paste("RepC",net.met[i],sep="."))),eval(parse(text= paste("RepD",net.met[i],sep=".")))),8)
}

##################
#Preliminary figures 
##################
#..Select the data needed
dates = c("20160509","20160517","20160523","20160531" )
sizes = c(7.5,13,22.5,45)
isolated = Prot.b$Label[which(Prot.b$Treatment=="Isolated")] #Labels of isolated patches
connected = Prot.b$Label[which(Prot.b$Treatment=="Connected")] #Labels of connected patches
x=subset(Prot.b,date %in% dates & Label %in% c(connected) & Landscape=="Blue")[,c("date","Label","Size","Replicate","degree","dist","diam","abundance","bioarea","richness","simpson",species)]
x = droplevels((x))

bargraph.CI(Size,richness,group=dist,data=x,legend=T)
bargraph.CI(Size,abundance,group=dist,data=x,legend=T)
bargraph.CI(Size,bioarea,group=dist,data=x,legend=T)

x = subset(Cyto,Landscape=="Blue" & Treatment %in% c("Isolated","Connected") & Date!=20160502 & Replicate %in% c("A","B"))[,c("day","Label","Treatment","Landscape","Replicate","Size","degree","diam","dist","Count","density")] 
x = droplevels(x)

bargraph.CI(Size,density,group=dist,data=x,legend=T)

##################
#Pie chart network plots (need to be re-done)
##################

################
#...Bacteria

#To graph an averaged through time, just remove the loop through Cyto$Date
{
for(i in 1:nlevels(Cyto$Date)) {
  for(j in 1:nlevels(Cyto$Replicate)) {
    for(k in 1:nlevels(Cyto$Treatment)) {
      for(z in 1:nlevels(Cyto$Landscape)) { 
      
      abundance = Cyto$count.ml[which(Cyto$Date == levels(Cyto$Date)[i] & Cyto$Treatment==levels(Cyto$Treatment)[k] & Cyto$Replicate==levels(Cyto$Replicate)[j] & Cyto$Landscape==levels(Cyto$Landscape)[z])]
      patch.size = as.numeric(as.character(Cyto$Size[which(Cyto$Date == levels(Cyto$Date)[i] & Cyto$Treatment==levels(Cyto$Treatment)[k] & Cyto$Replicate==levels(Cyto$Replicate)[j] & Cyto$Landscape==levels(Cyto$Landscape)[z])]))
      magn = 0.6 #amplification magnitude for the patch size vector
      
      #Color palette for vertices
      clrs = colorRampPalette(c('blue',"cyan",'yellow',"red"))
      fine = 36 # this will adjust the resolving power.
      graphCol = clrs(fine)[as.numeric(cut(abundance,breaks = fine))]
      
      #Plot the figure
      
      #Generate figure layout corresponding to replicate landscape
      if(levels(Cyto$Landscape)[z]=="Blue") { 
      coords1 <- layout_on_grid(eval(parse(text= paste("Rep",levels(Cyto$Replicate)[j],".","g",sep=""))),width=6,height=6) #set the rectangular layout 
      } else {coords1 <- layout_on_grid(Green.g,width=6,height=6)}  #set the rectangular layout 
      
      coords1[,2] <- max(coords1[,2])-coords1[,2] 
      
      #Open pdf
      pdf(paste(graphpath,levels(Cyto$Landscape)[z],levels(Cyto$Replicate)[j],levels(Cyto$Treatment)[k],"BACT",levels(Cyto$Date)[i],".pdf",sep="_"),width=8,height=8)
      
      #if(i==4) magn=0.008 else magn=0.001
    
      if(levels(Cyto$Landscape)[z]=="Blue") { 
      plot(eval(parse(text= paste("Rep",levels(Cyto$Replicate)[j],".","g",sep=""))),
           layout=coords1,
           vertex.color=graphCol,
           vertex.size=(patch.size)*magn, 
           vertex.label=NA) } else{plot(Green.g,layout=coords1,vertex.color=graphCol,vertex.size=(patch.size)*magn,vertex.label=NA) }  
            #plot: you can change vertex.label by bacteria count to make sure that the color fits with the count
      
      title(main=paste(levels(Cyto$Landscape)[z],levels(Cyto$Replicate)[j],levels(Cyto$Treatment)[k],"BACT",levels(Cyto$Date)[i],sep="_"))
       #Legend
      library(aqfig)
      vertical.image.legend(col=clrs(fine), zlim=range(abundance))
      
  
      dev.off()
      }
    }
  }
}
}

################
#...Protist
{# Color palett for all pies
V(RepA.g)$pie.color=list(rainbow(7))
V(RepB.g)$pie.color=list(rainbow(7))
V(RepC.g)$pie.color=list(rainbow(7))
V(RepD.g)$pie.color=list(rainbow(7))

for(i in 1:nlevels(Prot.b$date)) {
  for(j in 1:nlevels(Prot.b$Replicate)) {
    for(k in 1:nlevels(Prot.b$Treatment)) {


      rel.ab = Prot.b[which(Prot.b$date == levels(Prot.b$date)[i] & Prot.b$Treatment==levels(Prot.b$Treatment)[k] & Prot.b$Replicate==levels(Prot.b$Replicate)[j]),32:38]
      rel.ab = rel.ab+0.00000000000000000000000000000001
      rel.ab.list <- as.list(as.data.frame(t(rel.ab)))
      abundance = Prot.b$Pcount[which(Prot.b$date == levels(Prot.b$date)[i] & Prot.b$Treatment==levels(Prot.b$Treatment)[k] & Prot.b$Replicate==levels(Prot.b$Replicate)[j])]

      #Plot the figure

      #Generate figure layout corresponding to replicate landscape
      coords1 <- layout_on_grid(eval(parse(text= paste("Rep",levels(Prot.b$Replicate)[j],".","g",sep=""))),width=6,height=6) #set the rectangular layout
      coords1[,2] <- max(coords1[,2])-coords1[,2]

      pdf(paste(graphpath,"Blue",levels(Prot.b$Replicate)[j],levels(Prot.b$Treatment)[k],"Prot",levels(Prot.b$date)[i],".pdf",sep="_"),width=8,height=8)

      if(i==4) magn=0.008 else magn=0.001

      plot(eval(parse(text= paste("Rep",levels(Prot.b$Replicate)[j],".","g",sep=""))),
           layout=coords1,
           vertex.shape="pie",
           vertex.pie=rel.ab.list,
           vertex.size=(abundance)*magn,
           vertex.label=NA)
      #Legend
      col.labels=c("Rot","Spi","Ble","PCA","Col","Chi","Tet")
      testcol = c("#FF0000FF","#FFDB00FF","#49FF00FF","#00FF92FF","#0092FFFF","#4900FFFF","#FF00DBFF")
      color.legend(-1.2,-0.3,-1.4,-1,col.labels,testcol,cex=1,gradient="y")
      #dif between x1 and xr is the width of the rectangle

      dev.off()

    }
  }
}}


#########################################################################
################# Log Response Ratio 
#########################################################################

##################
# Clean the environment
#################
#...Clean the environment 
remove(lala);remove(lala2);remove(lala3);rm(sel.date);rm(sel.rep);rm(sel.treatment);rm(col.p)
remove(Green);remove(Green.m);remove(RepA);remove(RepA.m);remove(RepB);remove(RepB.m);remove(RepC);remove(RepC.m);remove(RepD);remove(RepD.m)
remove(coords1);remove(abundance);remove(fine);remove(graphCol);rm(i);rm(j);rm(k);rm(magn);rm(patch.size);rm(z)
rm(clrs)
detach("package:igraph") #igraph masks several functions in other packages

##################
# Calculate and Plot LRR 
##################

#... Set a minimum value to replace 0, for log response ratios
mini=0.001

##############
#...Protist

#..Replace values 0 by mini value
Prot.b$abundance = ifelse(Prot.b$abundance == 0,mini,Prot.b$abundance)
Prot.b$bioarea = ifelse(Prot.b$bioarea == 0,mini,Prot.b$bioarea)
Prot.b$richness = ifelse(Prot.b$richness == 0,mini,Prot.b$richness)
Prot.b$simpson = ifelse(is.nan(Prot.b$simpson),mini,Prot.b$simpson)
for(i in 1:length(species)){
  Prot.b[,species[i]] = ifelse(Prot.b[,species[i]] == 0,mini,Prot.b[,species[i]])
}

#..Calculate Log response ratio
dates = c("20160509","20160517","20160523","20160531" )
sizes = c(7.5,13,22.5,45)
isolated = Prot.b$Label[which(Prot.b$Treatment=="Isolated")] #Labels of isolated patches
connected = Prot.b$Label[which(Prot.b$Treatment=="Connected")] #Labels of connected patches

#..Select the data needed
x=subset(Prot.b,date %in% dates & Label %in% c(isolated,connected))[,c("date","Label","Size","Replicate","abundance","bioarea","richness","simpson",species)]

#..Calculate the different log response ratios
{ 
logES_ab = logES_ba = logES_rich = logES_simp = logES_species = vector(mode="list",length=4);

for(i in 1:length(dates)){
  x_is=subset(x,date==dates[i] & Label %in% isolated)
  x_con=subset(x,date==dates[i] & Label %in% connected)
  logES_ab[[i]] = cbind(x_is[,c("date","Size","Replicate")],logES=log(x_con$abundance/x_is$abundance))
  logES_ba[[i]] = cbind(x_is[,c("date","Size","Replicate")],logES=log(x_con$bioarea/x_is$bioarea))
  logES_rich[[i]] = cbind(x_is[,c("date","Size","Replicate")],logES=log(x_con$richness/x_is$richness))
  logES_simp[[i]] = cbind(x_is[,c("date","Size","Replicate")],logES=log(x_con$simpson/x_is$simpson))
  logES_species[[i]] = cbind(x_is[,c("date","Size","Replicate")],logRot=log(x_con$Rot/x_is$Rot),logSpi=log(x_con$Spi/x_is$Spi),logBle=log(x_con$Ble/x_is$Ble),logPca=log(x_con$Pca/x_is$Pca),logCol=log(x_con$Col/x_is$Col),logChi=log(x_con$Chi/x_is$Chi),logTet=log(x_con$Tet/x_is$Tet),logOther=log(x_con$Other/x_is$Other))
}
}

#..Put the data together
logES_ab_tot = rbind(logES_ab[[1]],logES_ab[[2]],logES_ab[[3]],logES_ab[[4]])
logES_ba_tot = rbind(logES_ba[[1]],logES_ba[[2]],logES_ba[[3]],logES_ba[[4]])
logES_rich_tot = rbind(logES_rich[[1]],logES_rich[[2]],logES_rich[[3]],logES_rich[[4]])
logES_simp_tot = rbind(logES_simp[[1]],logES_simp[[2]],logES_simp[[3]],logES_simp[[4]])
logES_species_tot = rbind(logES_species[[1]],logES_species[[2]],logES_species[[3]],logES_species[[4]])

#..Define the metrics that will be plotted 
metrics = list(logES_ab_tot,logES_ba_tot,logES_rich_tot,logES_simp_tot)
names(metrics)= c("abundance","bioarea","richness","simpson")

#... Plots
#..Effect size by time and by patch size
{ 
pdf(paste(graphpath,"Effect_size_Protists.pdf"),width=8,height=8)
for(i in 1:length(metrics)){
  d=metrics[[i]]
  #... Plots per patch size for the LRR of each metric
  plot(NA,xlim=c(0.5,4.5),ylim=c(-1,1),main = names(metrics)[i],xlab="PatchSize",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:4,labels=sizes)
  abline(h=0,lty=3)
  for(k in 1:length(sizes)){
    #Calculate CI 95% after Hedges et al., 1999
    x0_con=subset(x[,names(metrics)[i]],x$Size==sizes[k] & x$Label %in% connected)
    x0_is = subset(x[,names(metrics)[i]],x$Size==sizes[k] & x$Label %in% isolated)
    m_con=mean(x0_con)
    m_is = mean(x0_is)
    n_con = length(x0_con)
    n_is = length(x0_is)
    CI= qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
    #Calculate mean LRR
    x0=subset(d,Size==sizes[k])$logES 
    m=mean(x0)
    #add points and arrows
    arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
  #... Plots per date for the LRR of each metric
  plot(NA,xlim=c(0.5,4.5),ylim=c(-1,1),main = names(metrics)[i],xlab="Date",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:4,labels=dates)
  abline(h=0,lty=3)
  for(k in 1:length(dates)){
    #Calculate CI 95% after Hedges et al., 1999
    x0_con=subset(x[,names(metrics)[i]],x$date==dates[k] & x$Label %in% connected)
    x0_is = subset(x[,names(metrics)[i]],x$date==dates[k] & x$Label %in% isolated)
    m_con=mean(x0_con)
    m_is = mean(x0_is)
    n_con = length(x0_con)
    n_is = length(x0_is)
    CI= qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
    #Calculate mean LRR
    x0=subset(d,date==dates[k])$logES
    m=mean(x0)
    #add points and arrows
    arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
  
  #... interactions dates x size
  for(j in 1:length(dates)){
    d=subset(metrics[[i]],date==dates[j])
    plot(NA,xlim=c(0.5,4.5),ylim=range(d$logES),main = paste(names(metrics)[i],dates[j]),xlab="PatchSize",ylab="ln(patch connected / patch isolated)",xaxt="n")
    axis(1,at=1:4,labels=sizes)
    abline(h=0,lty=3)
    for(k in 1:length(sizes)){
      #Calculate CI 95% after Hedges et al., 1999
      x0_con=subset(x[,names(metrics)[i]],x$Size==sizes[k] & x$date==dates[j] & x$Label %in% connected)
      x0_is = subset(x[,names(metrics)[i]],x$Size==sizes[k] & x$date==dates[j] & x$Label %in% isolated)
      m_con=mean(x0_con)
      m_is = mean(x0_is)
      n_con = length(x0_con)
      n_is = length(x0_is)
      CI= qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
      #Calculate mean LRR
      x0=subset(d,Size==sizes[k])$logES
      m=mean(x0)
      #add points and arrows
      arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
      points(k,m,pch=19)
    }
  }
}
dev.off()

}

#..Effect size by species
{ 
pdf(paste(graphpath,"Effect_size_Protists_byspecies.pdf"),width=8,height=8)
ES_species= c("logRot","logSpi","logBle","logPca","logCol","logChi","logTet","logOther")
d=logES_species_tot
plot(NA,xlim=c(0.5,length(species)+0.5),ylim=c(-2,2),main = "Protist species",xlab="Species",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:length(species),labels=species,las=2)
abline(h=0,lty=3)
for(k in 1:length(species)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x[,species[k]],x$Label %in% connected)
  x0_is = subset(x[,species[k]],x$Label %in% isolated)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
  #Calculate mean LRR
  x0=d[,ES_species[k]]
  m=mean(x0)
  #add points and arrows
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}
#Interaction species*dates
ES_species= c("logRot","logSpi","logBle","logPca","logCol","logChi","logTet","logOther")
for(j in 1:length(dates)){ 
  d=subset(logES_species_tot,date==dates[j])
  plot(NA,xlim=c(0.5,length(species)+0.5),ylim=c(-5,5),main =paste("Protist species",dates[j]),xlab="Species",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:length(species),labels=species,las=2)
  abline(h=0,lty=3)
  for(k in 1:length(species)){
    #Calculate CI 95% after Hedges et al., 1999
    x0_con= subset(x[,species[k]],x$Label %in% connected & x$date==dates[j])
    x0_is = subset(x[,species[k]],x$Label %in% isolated & x$date==dates[j])
    m_con=mean(x0_con)
    m_is = mean(x0_is)
    n_con = length(x0_con)
    n_is = length(x0_is)
    CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
    #Calculate mean LRR
    x0=d[,ES_species[k]]
    m=mean(x0)
    #add points and arrows
    arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
   }
#Interaction species*patch size
ES_species= c("logRot","logSpi","logBle","logPca","logCol","logChi","logTet","logOther")
for(j in 1:length(sizes)){ 
  d=subset(logES_species_tot,Size==sizes[j])
  plot(NA,xlim=c(0.5,length(species)+0.5),ylim=c(-5,5),main =paste("Protist species",sizes[j]),xlab="Species",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:length(species),labels=species,las=2)
  abline(h=0,lty=3)
  for(k in 1:length(species)){
    #Calculate CI 95% after Hedges et al., 1999
    x0_con= subset(x[,species[k]],x$Label %in% connected & x$Size==sizes[j])
    x0_is = subset(x[,species[k]],x$Label %in% isolated & x$Size==sizes[j])
    m_con=mean(x0_con)
    m_is = mean(x0_is)
    n_con = length(x0_con)
    n_is = length(x0_is)
    CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
    #Calculate mean LRR
    x0=d[,ES_species[k]]
    m=mean(x0)
    #add points and arrows
    arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
}
dev.off()
}

#############
#..Bacteria

#Select the data needed 
x = subset(Cyto,Treatment %in% c("Isolated","Connected") & Date!=20160502 & Replicate %in% c("A","B"))[,c("day","Label","Treatment","Landscape","Replicate","Size","degree","diam","dist","Count","density")] 
#Add patch size and richness and abundance in the blue landscape
x$Size_b = 10
x$richness_b = 0
x$abundance_b = 0
x$bact_blue = 0
x$Size_b[which(x$Landscape=="Green")] = as.numeric(as.character(Prot.b$Size[which(Prot.b$Replicate %in% c("A","B") & Prot.b$date!=20160502 & Prot.b$Treatment %in% c("Isolated","Connected"))]))
x$richness_b[which(x$Landscape=="Green")] = Prot.b$richness[which(Prot.b$Replicate %in% c("A","B") & Prot.b$date!=20160502 & Prot.b$Treatment %in% c("Isolated","Connected"))]
x$abundance_b[which(x$Landscape=="Green")] = Prot.b$abundance[which(Prot.b$Replicate %in% c("A","B") & Prot.b$date!=20160502 & Prot.b$Treatment %in% c("Isolated","Connected"))]
x$bact_blue[which(x$Landscape=="Green")] = x$density[which(x$Landscape=="Blue")]

#grouping variables
days = c(7,15,21,29)
richness_blue = c(0.001,1,2,3,4,5,6,7,8)
sizes = c(7.5,13,22.5,45)
dist_outlet = c(0,1,2,3,4,5,6)
dist_outlet2 = c(0,1,2,5,6)
degree_cent  = c(1,2,3,4,5,6)
#Separate isolated/connected green/blue
isolated_green = Cyto$Label[which(Cyto$Treatment=="Isolated" & Cyto$Landscape=="Green")] #Labels of isolated patches
isolated_blue = Cyto$Label[which(Cyto$Treatment=="Isolated" & Cyto$Landscape=="Blue")] #Labels of isolated patches
connected_green = Cyto$Label[which(Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")]#Labels of connected patches
connected_blue = Cyto$Label[which(Cyto$Treatment=="Connected" & Cyto$Landscape=="Blue")]#Labels of connected patches


with(x[which(x$Label %in% connected_green & x$day==15),],plot(density ~ abundance_b))
with(x[which(x$Label %in% connected_green & x$day==29),],plot(density ~ bact_blue))
with(x[which(x$Label %in% connected_green & x$day==29),],plot(density ~ as.factor(richness_b)))

with(x[which(x$Label %in% connected_green),],bargraph.CI(Size_b,richness_b,group=dist,legend=T))
with(x[which(x$Label %in% connected_green),],bargraph.CI(Size_b,abundance_b,group=dist,legend=T))
with(x[which(x$Label %in% connected_green),],bargraph.CI(Size_b,bact_blue,group=dist,legend=T))



#Calculate Log Response Ratio
logES_Blue = logES_Green = vector(mode="list",length=4);
{ 
for(i in 1:length(days)){
  x_is_blue=subset(x,day==days[i] & Label %in% isolated_blue)
  x_con_blue=subset(x,day==days[i] & Label %in% connected_blue)
  logES_Blue[[i]] = cbind(x_is_blue[,c("day","Size","Replicate","density")],logES=log(x_con_blue$density/x_is_blue$density))
  x_is_green=subset(x,day==days[i] & Label %in% isolated_green)
  x_con_green=subset(x,day==days[i] & Label %in% connected_green)
  logES_Green[[i]] = cbind(x_con_green[,c("day","Size","Size_b","degree","diam","richness_b","abundance_b","bact_blue","dist","Replicate","density")],logES=log(x_con_green$density/x_is_green$density))
 }
}

#Put LRR back together
logES_Blue_tot = rbind(logES_Blue[[1]],logES_Blue[[2]],logES_Blue[[3]],logES_Blue[[4]])
logES_Green_tot = rbind(logES_Green[[1]],logES_Green[[2]],logES_Green[[3]],logES_Green[[4]])


#Plot
{ 
pdf(paste(graphpath,"Effect_size_Bacteria.pdf"),width=8,height=8)

#BACTERIA IN GREEN LANDSCAPES
  d=logES_Green_tot
  
#... Plot over days for green landscapes
plot(NA,xlim=c(0.5,4.5),ylim=c(-0.5,0.5),main = "Bacteria in Green",xlab="Days",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:4,labels=days)
abline(h=0,lty=3)
for(k in 1:length(days)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x$density,x$day==days[k] & x$Label %in% connected_green)
  x0_is = subset(x$density,x$day==days[k] & x$Label %in% isolated_green)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
  #Calculate mean LRR
  x0=subset(d,day==days[k])$logES
  m=mean(x0)
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}
#... Plot over patch size (in blue connected) for green landscapes
plot(NA,xlim=c(0.5,4.5),ylim=c(-0.5,0.5),main = "Bacteria in Green",xlab="Size in blue connected",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:4,labels=sizes)
abline(h=0,lty=3)
for(k in 1:length(sizes)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x$density,x$Size_b==sizes[k] & x$Label %in% connected_green)
  x0_is = subset(x$density,x$Size_b==sizes[k] & x$Label %in% isolated_green)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))
  #Calculate mean LRR
  x0=subset(d,d$Size_b==sizes[k])$logES
  m=mean(x0)
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}
#... Plot over distance to outlet (in blue connected) for green landscapes
plot(NA,xlim=c(0.5,7.5),ylim=c(-1,1),main = "Bacteria in Green",xlab="Distance to outlet in blue connected",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:7,labels=dist_outlet)
abline(h=0,lty=3)
for(k in 1:length(dist_outlet)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x$density,x$dist==dist_outlet[k] & x$Label %in% connected_green)
  x0_is = subset(x$density,x$dist==dist_outlet[k] & x$Label %in% isolated_green)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))
  #Calculate mean LRR
  x0=subset(d,d$dist==dist_outlet[k])$logES
  m=mean(x0)
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}
#... Plot over degree (in blue connected) for green landscapes
plot(NA,xlim=c(0.5,6.5),ylim=c(-1,1),main = "Bacteria in Green",xlab="Degree centrality in blue connected",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:6,labels=degree_cent)
abline(h=0,lty=3)
for(k in 1:length(degree_cent)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x$density,x$degree==degree_cent[k] & x$Label %in% connected_green)
  x0_is = subset(x$density,x$degree==degree_cent[k] & x$Label %in% isolated_green)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))
  #Calculate mean LRR
  x0=subset(d,d$degree==degree_cent[k])$logES
  m=mean(x0)
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}
#... Plot over protist richness (in blue connected) for green landscapes
plot(NA,xlim=c(0.5,9.5),ylim=c(-1,1),main = "Bacteria in Green",xlab="Protist richness (in blue connected)",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:9,labels=c(0,1,2,3,4,5,6,7,8))
abline(h=0,lty=3)
for(k in 1:length(richness_blue)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x$density,x$richness_b==richness_blue[k] & x$Label %in% connected_green)
  x0_is = subset(x$density,x$richness_b==richness_blue[k] & x$Label %in% isolated_green)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))
  #Calculate mean LRR
  x0=subset(d,d$richness_b==richness_blue[k])$logES
  m=mean(x0)
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}
#Plot over Protist and bacteria density in blue
plot(NA,xlim=range(d$abundance),ylim=range(d$logES),xlab="Protist density in blue connected",ylab="ln(patch connected / patch isolated)",main = "Bacteria in Green")
points(d$abundance_b,d$logES)
plot(NA,xlim=range(d$bact_blue),ylim=range(d$logES),xlab="Bacteria density in blue connected",ylab="ln(patch connected / patch isolated)",main = "Bacteria in Green")
points(d$bact_blue,d$logES)

#Interaction between distanece to outlet and patch size
plot(NA,ylim=c(-1,1.2),xlim=c(0.5,4.5),xaxt="n",xlab="Patch size",ylab="LRR")
axis(1,at=1:4,labels=sizes)
abline(h=0,lty=3)
#library(RColorBrewer)
#col.dist=brewer.pal(7,"Set3")
col.dist = c('#7fc97f','#beaed4','#fdc086','#ffff99','#386cb0','#f0027f','#bf5b17')
#col.dist = rainbow(7)
for(j in 1:length(dist_outlet2)){ 
  for(k in 1:length(sizes)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x$density,x$Size_b==sizes[k] & x$dist == dist_outlet2[j] & x$Label %in% connected_green)
  x0_is = subset(x$density,x$Size_b==sizes[k] & x$dist == dist_outlet2[j] & x$Label %in% isolated_green)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))
  #Calculate mean LRR
  x0=subset(d,Size_b==sizes[k] & dist == dist_outlet2[j])$logES
  m=mean(x0)
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05,col=col.dist[j])
  points(k,m,pch=19,col=col.dist[j])
 }
  }

legend("bottomright",legend=dist_outlet2,col=col.dist,pch=16,title="Distance to outlet",ncol=2,cex=0.8)

#BACTERIA IN BLUE LANDSCAPES
d=logES_Blue_tot

#... Plot over days for blue landscapes
plot(NA,xlim=c(0.5,4.5),ylim=c(-0.5,0.5),main = "Bacteria in Blue",xlab="Days",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:4,labels=days)
abline(h=0,lty=3)
for(k in 1:length(days)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x$density,x$day==days[k] & x$Label %in% connected_blue)
  x0_is = subset(x$density,x$day==days[k] & x$Label %in% isolated_blue)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
  #Calculate mean LRR
  x0=subset(d,day==days[k])$logES
  m=mean(x0)
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}

#... Plot by patch size for blue landscapes
plot(NA,xlim=c(0.5,4.5),ylim=c(-0.5,0.5),main = "Bacteria in Blue",xlab="Patch Size",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:4,labels=sizes)
abline(h=0,lty=3)
for(k in 1:length(sizes)){
  #Calculate CI 95% after Hedges et al., 1999
  x0_con= subset(x$density,x$Size==sizes[k] & x$Label %in% connected_blue)
  x0_is = subset(x$density,x$Size==sizes[k] & x$Label %in% isolated_blue)
  m_con=mean(x0_con)
  m_is = mean(x0_is)
  n_con = length(x0_con)
  n_is = length(x0_is)
  CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
  #Calculate mean LRR
  x0=subset(d,Size==sizes[k])$logES
  m=mean(x0)
  arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}
#... Interactions day x size blue landscapes
for(j in 1:length(days)){
  d=subset(logES_Blue_tot,day==days[j])
  plot(NA,xlim=c(0.5,4.5),ylim=range(d$logES),main = paste("Bacteria in Blue",days[j]),xlab="Patch Size",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:4,labels=sizes)
  abline(h=0,lty=3)
  for(k in 1:length(sizes)){
    #Calculate CI 95% after Hedges et al., 1999
    x0_con= subset(x$density,x$day==days[j] & x$Size == sizes[k] & x$Label %in% connected_blue)
    x0_is = subset(x$density,x$day==days[j] & x$Size == sizes[k] & x$Label %in% isolated_blue)
    m_con=mean(x0_con)
    m_is = mean(x0_is)
    n_con = length(x0_con)
    n_is = length(x0_is)
    CI = qnorm(0.975)*sqrt( (sd(x0_con)^2/(n_con*m_con^2)) + (sd(x0_is)^2/(n_is*m_is^2)))     
    #Calculate mean LRR
    x0=subset(d,Size==sizes[k])$logES
    m=mean(x0)
    arrows(k,m+CI,k,m-CI,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
}
dev.off()

}


####The END#######

#Some other figures that needs to be place somewhere  
x = subset(Cyto,Landscape %in% c("Blue") & Treatment %in% c("Isolated","Connected") & Date!=20160502 & Replicate %in% c("A","B"))[,c("day","Label","Treatment","Landscape","Replicate","Size","Count","density")] 
x = droplevels(x)
lineplot.CI(x$Size,x$density,xlab="Size",ylab="Bacteria density (Blue)")
x = subset(Prot.b,Treatment %in% c("Isolated","Connected") & date!=20160502 & Replicate %in% c("A","B"))[,c("day","Label","Treatment","Landscape","Replicate","Size","abundance","richness","simpson")] 
x = droplevels(x)
lineplot.CI(x$Size,x$richness,xlab="Size",ylab="Protist richness")
lineplot.CI(x$Size,x$abundance,xlab="Size",ylab="Protist density")
lineplot.CI(x$Size,x$simpson,xlab="Size",ylab="Evenness")



library(nlme)
Mod = lme(var.prot.rich ~ as.factor(Size)+day, ~ as.factor(date)|Replicate,data=Prot.b.connected,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
hist(Mod$residuals,breaks=50)
summary(Mod)$tTable

Mod2 = lme(var.bact.green ~ as.factor(Size)+day, ~ as.factor(Date)|Replicate,data=Cyto.g.isolated,method="REML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
hist(Mod2$residuals,breaks=50)
summary(Mod2)$tTable




#########################################################################
################# MAIN ANALYSES (NEED TO RE_DONE)
#########################################################################


##################
# Remove Day 0 
#################

#Protist
sel.date = c("20160509","20160517","20160523","20160531")
Prot.b = Prot.b[Prot.b$date %in% sel.date,] 
Prot.b$date = factor(Prot.b$date)
#Bacteria
Cyto = Cyto[Cyto$Date %in% sel.date,]
Cyto$Date = factor(Cyto$Date)

##################
# Merge network metrics to Prot.b, Cyto.b, Cyto.g data
#################
{ 
  
  #...Extract blue and green landscapes only for bacteria (we don't have protist data for green landscapes)
  selb = c("Blue")
  selg = c("Green")
  Cyto.b = Cyto[Cyto$Landscape %in% selb,]
  Cyto.g = Cyto[Cyto$Landscape %in% selg,]
  Cyto.b= droplevels(Cyto.b)
  Cyto.g= droplevels(Cyto.g)
  
  #...Remove moncultures from protist data
  sel.treat = c("Connected","Isolated")
  Prot.b = Prot.b[which(Prot.b$Treatment %in% sel.treat),]
  Prot.b=droplevels(Prot.b)
  
  #...merge network metrics to data
  #.....Protist
  Prot.b$Degree = rep(c(RepA.degree,RepB.degree,RepC.degree,RepD.degree),8)
  Prot.b$Dist = rep(c(RepA.dist,RepB.dist,RepC.dist,RepD.dist),8)
  Prot.b$Diam = rep(c(RepA.diam,RepB.diam,RepC.diam,RepD.diam),8)
  Prot.b$Drainage = rep(c(RepA.drainage,RepB.drainage,RepC.drainage,RepD.drainage),8)
  #.....Bacteria
  Cyto.b$Degree = rep(c(RepA.degree,RepB.degree),8)
  Cyto.b$Dist = rep(c(RepA.dist,RepB.dist),8)
  Cyto.b$Diam = rep(c(RepA.diam,RepB.diam),8)
  Cyto.b$Drainage = rep(c(RepA.drainage,RepB.drainage),8)
  
  Cyto.g$Size = Prot.b$Size[which(Prot.b$Replicate %in% c("A","B"))]
  Cyto.g$Degree = rep(c(RepA.degree,RepB.degree),8)
  Cyto.g$Dist = rep(c(RepA.dist,RepB.dist),8)
  Cyto.g$Diam = rep(c(RepA.diam,RepB.diam),8)
  Cyto.g$Drainage = rep(c(RepA.drainage,RepB.drainage),8)
  
  #...Fix some variables
  Prot.b$Size = as.numeric(as.character(Prot.b$Size))
  Cyto.b$Size = as.numeric(as.character(Cyto.b$Size))
  Cyto.g$Size = as.numeric(as.character(Cyto.g$Size))
  
}

##################
# Plot effects from Green to blue
##################
{ 
  pdf(paste0(graphpath,"Green_to_Blue.pdf"),width=8,height=8)
  
  plot(Prot.b$Prot.tot.ab[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==7 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")],pch=16,ylab="Protist total abundance (Blue-connected)",xlab="Bacteria density (ind/mL)(Green-connected)",main="Day 7")
  abline(lm(Prot.b$Prot.tot.ab[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==7 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")]), col="red",lwd=2)
  
  plot(Prot.b$Prot.tot.ab[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==29 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")],pch=16,ylab="Protist total abundance (Blue-connected)",xlab="Bacteria density (ind/mL)(Green-connected)",main="Day 29")
  abline(lm(Prot.b$Prot.tot.ab[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==29 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")]),col="red",lwd=2)
  
  plot(Prot.b$Prot.rich[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==7 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")],pch=16,ylab="Protist richness (Blue-connected)",xlab="Bacteria density (ind/mL)(Green-connected)",main='Day 7')
  abline(lm(Prot.b$Prot.rich[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==7 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")]),col="red",lwd=2)
  
  plot(Prot.b$Prot.rich[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==29 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")],pch=16,ylab="Protist richness (Blue-connected)",xlab="Bacteria density (ind/mL)(Green-connected)", main="Day 29")
  abline(lm(Prot.b$Prot.rich[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==29 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")]),col="red",lwd=2)
  
  plot(Cyto$count.ml[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==7 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")],pch=16,ylab="Bacteria density (ind/mL) (Blue-connected)",xlab="Bacteria density (ind/mL)(Green-connected)", main = "Day 7")
  abline(lm(Cyto$count.ml[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==7 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")]),col="red",lwd=2)
  
  plot(Cyto$count.ml[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==29 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")],pch=16,ylab="Bacteria density (ind/mL) (Blue-connected)",xlab="Bacteria density (ind/mL)(Green-connected)", main="Day 29")
  abline(lm(Cyto$count.ml[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] ~ Cyto$count.ml[which(Cyto$day==29 & Cyto$Landscape=="Green" & Cyto$Treatment=="Connected")]),col="red",lwd=2)
  
  dev.off()
}


##################
#Plot  effects from Blue to Green
##################
{ 
  
  pdf(paste0(graphpath,"Blue_to_Green.pdf"),width=8,height=8)
  
  #Effect of patch size in blue-connected on green-connected bacteria
  #day 7
  lineplot.CI(Prot.b$Size[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))],
              Cyto$count.ml[which(Cyto$day==7 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")],
              xlab="Patch size (blue-connected)",ylab="Bacteria density (ind/mL) (Green connected)",
              log="y",main="day 7",col="darkgreen")#ylim=c(5500,65000000)
  
  lineplot.CI(Prot.b$Size[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))],
              Prot.b$Prot.tot.ab[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))],
              col="blue",xlab="Patch size (blue-connected)",ylab="Protist abundance/bacteria density (blue-connected) - LOG",
              ylim=c(5500,65000000),log="y",pch=17,main="Day 7")
  par(new=T)
  lineplot.CI(Prot.b$Size[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))],
              Cyto$count.ml[which(Cyto$day==7 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Blue")],
              col="blue",xlab="",ylab="",ylim=c(5500,65000000),log="y",pch=16)
  
  #day 29
  lineplot.CI(Prot.b$Size[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))],
              Cyto$count.ml[which(Cyto$day==29 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")],
              xlab="Patch size (blue-connected)",ylab="Bacteria density (ind/mL) (Green connected)",
              log="y",col="darkgreen",main="day 29")#ylim=c(5500,65000000)
  
  lineplot.CI(Prot.b$Size[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))],
              Prot.b$Prot.tot.ab[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))],
              col="blue",xlab="Patch size (blue-connected)",ylab="Protist abundance/bacteria density (blue-connected) - LOG",
              ylim=c(1400,70000000),log="y",main ="day 29",pch=17)
  par(new=T)
  lineplot.CI(Prot.b$Size[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))],
              Cyto$count.ml[which(Cyto$day==29 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Blue")],
              col="blue",xlab="",ylab="",ylim=c(1400,70000000),log="y",pch=16)
  
  
  
  #Effect of closeness in blue-connected on green-connected bacteria
  
  #day 7
  plot(Cyto$count.ml[which(Cyto$day==7 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")] 
       ~ Prot.b$Diam[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
       pch=16,col="darkgreen",xlab="Closeness (Blue-connected)",ylab="Bacteria density (ind/mL) (Green connected)",main="Day 7")
  abline(lm(Cyto$count.ml[which(Cyto$day==7 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")] 
            ~ Prot.b$Diam[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))]),
         col="green",lwd=2)
  plot(Prot.b$Prot.tot.ab[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] 
       ~ Prot.b$Diam[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
       pch=17,col="blue",xlab="Closeness (Blue-connected)",ylab="Protist abundance (blue-connected)",main = "Day 7")
  # plot(Cyto$count.ml[which(Cyto$day==7 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Blue")] 
  #      ~ Prot.b$Diam[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
  #      pch=16,col="blue",xlab="Closeness (Blue-connected)",ylab="Protist abundance (blue-connected)",main = "Day 7")
  
  #day 29
  plot(Cyto$count.ml[which(Cyto$day==29 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")] 
       ~ Prot.b$Diam[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
       pch=16,col="darkgreen",xlab="Closeness (Blue-connected)",ylab="Bacteria density (ind/mL) (Green connected)",main="Day 29")
  abline(lm(Cyto$count.ml[which(Cyto$day==29 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")] 
            ~ Prot.b$Diam[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))]),
         col="green",lwd=2)
  plot(Prot.b$Prot.tot.ab[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] 
       ~ Prot.b$Diam[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
       pch=17,col="blue",xlab="Closeness (Blue-connected)",ylab="Protist abundance (blue-connected)",main = "Day 29")
  
  #Effect of distance to outlet in blue-connected on green-connected bacteria
  
  #day 7
  plot(Cyto$count.ml[which(Cyto$day==7 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")] 
       ~ Prot.b$Dist[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
       pch=16,col="darkgreen",xlab="Distance to outlet (Blue-connected)",ylab="Bacteria density (ind/mL) (Green connected)",main="Day 7")
  abline(lm(Cyto$count.ml[which(Cyto$day==7 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")] 
            ~ Prot.b$Dist[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))]),
         col="green",lwd=2)
  plot(Prot.b$Prot.tot.ab[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] 
       ~ Prot.b$Dist[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
       pch=17,col="blue",xlab="Distance to outlet (Blue-connected)",ylab="Protist abundance (blue-connected)",main = "Day 7")
  # plot(Cyto$count.ml[which(Cyto$day==7 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Blue")] 
  #      ~ Prot.b$Dist[which(Prot.b$day==7 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
  #      pch=16,col="blue",xlab="Distance to outlet (Blue-connected)",ylab="Protist abundance (blue-connected)",main = "Day 7")
  
  #day 29
  plot(Cyto$count.ml[which(Cyto$day==29 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")] 
       ~ Prot.b$Dist[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
       pch=16,col="darkgreen",xlab="Distance to outlet (Blue-connected)",ylab="Bacteria density (ind/mL) (Green connected)",main="Day 29")
  abline(lm(Cyto$count.ml[which(Cyto$day==29 & Cyto$Treatment=="Connected" & Cyto$Landscape=="Green")] 
            ~ Prot.b$Dist[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))]),
         col="green",lwd=2)
  plot(Prot.b$Prot.tot.ab[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))] 
       ~ Prot.b$Dist[which(Prot.b$day==29 & Prot.b$Treatment=="Connected" & Prot.b$Replicate %in% c("A","B"))], 
       pch=17,col="blue",xlab="Distance to outlet (Blue-connected)",ylab="Protist abundance (blue-connected)",main = "Day 29")
  
  dev.off()
}





{#..........................................#
#           Blue-green interactions        #
#..........................................#

#Run model
library(nlme)
Mod = lme(log(Count) ~ Treatment*Landscape*as.numeric(Date), ~ Date|Replicate,data=Cyto,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod)$tTable

##################
# Statistics
##################

with(Cyto.b,plot(density(count.ml)))

Mod = lme(log(Count) ~ Treatment*Size+as.numeric(Date), ~ as.factor(Date)|Replicate,data=Cyto.b,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod)$tTable

Mod2 = lme(Prot.tot.ab ~ Treatment*Size+as.numeric(date), ~ as.factor(date)|Replicate,data=Prot.b,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod2)$tTable

Mod3 = lme(Prot.rich ~ Treatment*Size+as.numeric(date), ~ as.factor(date)|Replicate,data=Prot.b,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod3)$tTable

biotic.cor = with(Cyto.b,cor(cbind(Size,Degree,Dist,Diam,Drainage)))
library(corrplot)
corrplot(biotic.cor,method="number",type="full")

#test the effect of bacteria density in green landscape on protist and bacteria in blue landscapes
Prot.b.AB = Prot.b[which(Prot.b$Replicate=="A" | Prot.b$Replicate=="B"),]
Prot.b.AB$BACT.count.ml.g = Cyto.g$count.ml
Prot.b.AB$BACT.count.ml.b = Cyto.b$count.ml

Mod4 = lme(BACT.count.ml.b ~ Treatment*Size*BACT.count.ml.g, ~ as.factor(date)|Replicate,data=Prot.b.AB,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod4)$tTable

Mod5 = lme(Prot.tot.ab ~ Treatment*Size*BACT.count.ml.g, ~ as.factor(date)|Replicate,data=Prot.b.AB,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod5)$tTable

Mod6 = lme(Prot.rich ~ Treatment*Size*BACT.count.ml.g, ~ as.factor(date)|Replicate,data=Prot.b.AB,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod6)$tTable

#Use only last day so taht we can use all four replicates
test = Cyto2[which(Cyto2$Date=="20160531" & Cyto2$Landscape=="Blue"),]
test$Replicate = factor(test$Replicate)
test2 = Prot.b[which(Prot.b$date=="20160531"),]
test2$Count = test$count

Mod6 = lme(Prot.rich ~ Treatment*Size*Count, ~ 1|Replicate,data=test2,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod6)$tTable


#..........................................#
#           Green Landscape only           #
#..........................................#


#...Stat

with(Cyto.g,plot(density(count.ml)))

library(nlme)
Mod = lme(log(count.ml) ~ Treatment*patch.size.blue,random = ~ 1|Replicate,data=Cyto.g,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod)$tTable

mod1.3.1 = with(Cyto.g,update(Mod, ~. - Treatment:patch.size.blue:Date))
summary(mod1.3.1)$tTable

mod1.3.2 = with(Cyto.g,update(mod1.3.1, ~. - patch.size.blue:Date))
summary(mod1.3.2)$tTable

mod1.3.3 = with(Cyto.g,update(mod1.3.2, ~. - Treatment:Date))
summary(mod1.3.3)$tTable

mod1.3.4 = with(Cyto.g,update(mod1.3.3, ~. - Treatment:patch.size.blue))
summary(mod1.3.4)$tTable




#Look at all replicates on the last sampling day
patch.size = Cyto2$Size[which(Cyto2$Date=="20160531" & Cyto2$Landscape=="Blue" & Cyto2$Treatment=="Connected")]
with(Cyto2[Cyto2$Date=="20160531" & Cyto2$Treatment=="Connected" & Cyto2$Landscape=="Green",],lineplot.CI(patch.size,count.ml,Treatment))


}

