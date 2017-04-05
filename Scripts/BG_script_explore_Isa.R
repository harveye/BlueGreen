#..........................................#
#...........Blue-Green project.............#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#... Contact: eric.harvey@eawag.ch                                                                                           #                                                                       #
#..........................................................................................................................................#

#########################################################################
################# DATA STRUCTURE
#########################################################################

##################
#Clear any variables from the R environment 
##################
rm(list=ls())
# search()
# detach(pos=3)

##################
#Directory paths 
##################
datapath = "~/Documents/Recherche/GITKRAKEN/BlueGreen/Data/"
to.script = "~/Documents/Recherche/GITKRAKEN/BlueGreen/Scripts/"
graphpath = "~/Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/4_RESULTS/BlueGreen"

##################
#Load packages
##################
library(vegan)
library(plyr)
library(plotrix)
library(sciplot)
library(car)
library(MASS)

##################
# Protist data
##################
setwd(datapath)
prot0 =  read.delim(paste0(datapath,"BG_protist_data_(20170307).txt")) #Protist
prot = prot0

#... Set a minimum value to replace 0, for log response ratios
mini=0.001

#... add column for abundance per ml
prot$abundance = apply(prot[,c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other")],1,sum)*1000
prot$abundance = ifelse(prot$abundance == 0,mini,prot$abundance)

#... set the minimum in individual species abundance
prot$Rot = ifelse(prot$Rot == 0,mini,prot$Rot*1000); prot$Spi = ifelse(prot$Spi == 0,mini,prot$Spi*1000);
prot$Ble = ifelse(prot$Ble == 0,mini,prot$Ble*1000);prot$Pca = ifelse(prot$Pca == 0,mini,prot$Pca*1000);
prot$Col = ifelse(prot$Col == 0,mini,prot$Col*1000);prot$Chi = ifelse(prot$Chi == 0,mini,prot$Chi*1000);
prot$Tet = ifelse(prot$Tet == 0,mini,prot$Tet*1000);prot$Other = ifelse(prot$Other == 0,mini,prot$Other*1000);
species = c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other")

#... add column for bioarea per ml
prot$bioarea = prot$bioarea_per_ul*1000
prot$bioarea = ifelse(prot$bioarea == 0,mini,prot$bioarea)

#... add column for richness
prot$richness = apply(prot[,species]>0,1,sum)
prot$richness = ifelse(prot$richness == 0,mini,prot$richness)

#... add column for evenness
simpson = function(x){return(sum((x/sum(x))^2))}; 
prot$simpson = apply(prot[,species],1,simpson)
prot$simpson = ifelse(is.nan(prot$simpson),mini,prot$simpson)

dates = c("5/9/2016",  "5/17/2016", "16-05-23",  "16-05-31" )
sizes = c(7.5,13,22.5,45)
isolated = 217:360 #Labels of isolated patches
connected = 361:504 #Labels of connected patches

#... Select the data
x=subset(prot,date %in% dates & Label %in% c(isolated,connected))[,c("date","Label","Size","Replicate","abundance","bioarea","richness","simpson",species)]


#... Calculate the different log response ratios
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

logES_ab_tot = rbind(logES_ab[[1]],logES_ab[[2]],logES_ab[[3]],logES_ab[[4]])
logES_ba_tot = rbind(logES_ba[[1]],logES_ba[[2]],logES_ba[[3]],logES_ba[[4]])
logES_rich_tot = rbind(logES_rich[[1]],logES_rich[[2]],logES_rich[[3]],logES_rich[[4]])
logES_simp_tot = rbind(logES_simp[[1]],logES_simp[[2]],logES_simp[[3]],logES_simp[[4]])
logES_species_tot = rbind(logES_species[[1]],logES_species[[2]],logES_species[[3]],logES_species[[4]])

metrics = list(logES_ab_tot,logES_ba_tot,logES_rich_tot,logES_simp_tot)
names(metrics)= c("abundances","bioarea","richness","simpson")

#... Plots
setwd(graphpath)
pdf("Effect_size_Protists.pdf")
for(i in 1:length(metrics)){
  d=metrics[[i]]
  #... Plots per patch size for the LRR of each metric
  plot(NA,xlim=c(0.5,4.5),ylim=c(-1,1),main = names(metrics)[i],xlab="PatchSize",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:4,labels=sizes)
  abline(h=0,lty=3)
  for(k in 1:length(sizes)){
    x0=subset(d,Size==sizes[k])$logES
    m=mean(x0)
    s=sd(x0)/sqrt(length(x0))*1.96
    arrows(k,m+s,k,m-s,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
  #... Plots per date for the LRR of each metric
  plot(NA,xlim=c(0.5,4.5),ylim=c(-1,1),main = names(metrics)[i],xlab="Date",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:4,labels=dates)
  abline(h=0,lty=3)
  for(k in 1:length(dates)){
    x0=subset(d,date==dates[k])$logES
    m=mean(x0)
    s=sd(x0)/sqrt(length(x0))*1.96
    arrows(k,m+s,k,m-s,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
  
  #... interactions dates x size
  for(j in 1:length(dates)){
    d=subset(metrics[[i]],date==dates[j])
    plot(NA,xlim=c(0.5,4.5),ylim=range(d$logES),main = paste(names(metrics)[i],dates[j]),xlab="PatchSize",ylab="ln(patch connected / patch isolated)",xaxt="n")
    axis(1,at=1:4,labels=sizes)
    abline(h=0,lty=3)
    for(k in 1:length(sizes)){
      x0=subset(d,Size==sizes[k])$logES
      m=mean(x0)
      s=sd(x0)/sqrt(length(x0))*1.96
      arrows(k,m+s,k,m-s,angle=90,code=3,length=0.05)
      points(k,m,pch=19)
    }
  }
}
dev.off()


pdf("Effect_size_Protists_byspecies.pdf")
ES_species= c("logRot","logSpi","logBle","logPca","logCol","logChi","logTet","logOther")
d=logES_species_tot
plot(NA,xlim=c(0.5,length(species)+0.5),ylim=c(-1,1),main = "Protist species",xlab="Species",ylab="ln(patch connected / patch isolated)",xaxt="n")
axis(1,at=1:length(species),labels=species,las=2)
abline(h=0,lty=3)
for(k in 1:length(species)){
  x0=d[,ES_species[k]]
  m=mean(x0)
  s=sd(x0)/sqrt(length(x0))*1.96
  arrows(k,m+s,k,m-s,angle=90,code=3,length=0.05)
  points(k,m,pch=19)
}
dev.off()




##################
# BACTERIA
##################
cyto = read.delim(paste0(datapath,"FULL_Cyto_BlueGreen(20160816).txt")) # Bacteria

#...Arrange structure cytometry data
cyto = cyto[order(cyto$Date),]
cyto$day = c(rep(0,32),rep(7,524),rep(15,524),rep(21,524),rep(29,524)) #create a new variable "day"
cyto$density = ((cyto$Count)*1000/cyto$Volume.uL)*1000 #convert from dens/50ul to dens/mL and then multiply by the cytometry dilution factor (1000)
x = subset(cyto,Treatment %in% c("Isolated","Connected") & Date!=20160502 &Replicate %in% c("A","B"))[,c("day","Label","Treatment","Replicate","Size","Count","density")] 
days = c(7,15,21,29)
sizes = c(7.5,13,22.5,45)

(unique(subset(xg,Treatment == "Connected" & Replicate == "B")$Label))
isolated_green = 1:72
connected_green = 73:144

logES_Blue = logES_Green = vector(mode="list",length=4);

for(i in 1:length(dates)){
  x_is_blue=subset(x,day==dates[i] & Label %in% isolated)
  x_con_blue=subset(x,day==dates[i] & Label %in% connected)
  logES_Blue[[i]] = cbind(x_is_blue[,c("day","Size","Replicate","density")],logES=log(x_con_blue$density/x_is_blue$density))
  x_is_green=subset(x,day==dates[i] & Label %in% isolated_green)
  x_con_green=subset(x,day==dates[i] & Label %in% connected_green)
  logES_Green[[i]] = cbind(x_is_green[,c("day","Size","Replicate","density")],logES=log(x_con_green$density/x_is_green$density))
}

logES_Blue_tot = rbind(logES_Blue[[1]],logES_Blue[[2]],logES_Blue[[3]],logES_Blue[[4]])
logES_Green_tot = rbind(logES_Green[[1]],logES_Green[[2]],logES_Green[[3]],logES_Green[[4]])


#... Plot
setwd(graphpath)
pdf("Effect_size_Bacteria.pdf")

  #... Plot over days for green landscapes
  d=logES_Green_tot
  plot(NA,xlim=c(0.5,4.5),ylim=c(-1,1),main = "Bacteria in Green",xlab="Days",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:4,labels=days)
  abline(h=0,lty=3)
  for(k in 1:length(days)){
    x0=subset(d,day==days[k])$logES
    m=mean(x0)
    s=sd(x0)/sqrt(length(x0))*1.96
    arrows(k,m+s,k,m-s,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
  
  #... Plot over days for blue landscapes
  d=logES_Blue_tot
  plot(NA,xlim=c(0.5,4.5),ylim=c(-1,1),main = "Bacteria in Blue",xlab="Days",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:4,labels=days)
  abline(h=0,lty=3)
  for(k in 1:length(days)){
    x0=subset(d,day==days[k])$logES
    m=mean(x0)
    s=sd(x0)/sqrt(length(x0))*1.96
    arrows(k,m+s,k,m-s,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
  
  #... Plot by patch size for blue landscapes
  plot(NA,xlim=c(0.5,4.5),ylim=c(-1,1),main = "Bacteria in Blue",xlab="Patch Size",ylab="ln(patch connected / patch isolated)",xaxt="n")
  axis(1,at=1:4,labels=sizes)
  abline(h=0,lty=3)
  for(k in 1:length(sizes)){
    x0=subset(d,Size==sizes[k])$logES
    m=mean(x0)
    s=sd(x0)/sqrt(length(x0))*1.96
    arrows(k,m+s,k,m-s,angle=90,code=3,length=0.05)
    points(k,m,pch=19)
  }
  
  #... Interactions day x size blue landscapes
  for(j in 1:length(dates)){
    d=subset(logES_Blue_tot,day==dates[j])
    plot(NA,xlim=c(0.5,4.5),ylim=range(d$logES),main = paste("Bacteria in Blue",dates[j]),xlab="Patch Size",ylab="ln(patch connected / patch isolated)",xaxt="n")
    axis(1,at=1:4,labels=sizes)
    abline(h=0,lty=3)
    for(k in 1:length(sizes)){
      x0=subset(d,Size==sizes[k])$logES
      m=mean(x0)
      s=sd(x0)/sqrt(length(x0))*1.96
      arrows(k,m+s,k,m-s,angle=90,code=3,length=0.05)
      points(k,m,pch=19)
    }
  }
dev.off()

