rm(list = ls())



##################
#Load data in R environment
##################
library(tidyverse)
library(vegan)
library(RColorBrewer)
library(gridExtra)
library(nlme)
library(scales)

##################
#Directory paths 
##################
to.point = "/Users/isabellegounand/Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/7_MANUSCRIPT_BG2/NEW_FIG"
to.data <- paste0(to.point,"/data/")
to.script <- paste0(to.point,"/scripts/")
to.output <- paste0(to.point,"/output/")
to.figs <- paste0(to.point,"/figs/")
to.R <- paste0(to.point,"/output/")

##################
#Load data in R environment
##################
data = readRDS(paste0(to.output,"cl.data.RDS"))
SP <- c("Rot","Spi","Ble","Pca","Col","Chi","Tet") #comm matrix

##################
#Set colour gradient for each species
##################
colvec <- brewer.pal(7,"Set2")
colvec <- c("#ff7f00","#1f78b4","navy","purple3","violetred","orangered2","#e31a1c")

##################
# Ordination (temporal trends in community structure)
##################
{
#Create interaction variable
data$Treat.Size <- interaction(data$Treatment,data$Size) 
data$Treat.Day <- interaction(data$Treatment,data$day)
data$Treat.Size.Day <- interaction(data$Treatment,data$Size,data$day)

#Place vars of interest into vectors for each component of ordination
SP <- c("Rot","Spi","Ble","Pca","Col","Chi","Tet") #comm matrix
ENV <- c("day.cont","centrality","Size","dist.outlet","Treatment","Treat.Size","Treat.Size.Day") #Constrains
CONT <- c("Replicate") #Effect to partial out

#Filter data as wanted (for bacteria need to use only rep A and B)
X <- data %>% 
  filter(day!=0 & Replicate %in% c("A","B","C","D")) %>% 
  mutate(Treatment = as.factor(Treatment),
         day = factor(day,levels=c("0","7","15","21","29")),
         day.cont = as.numeric(as.character(day)),
         Size = as.numeric(as.character(Size)),
         centrality = as.numeric(as.character(centrality)),
         dist.outlet = as.numeric(as.character(dist.outlet)))


#Create matrices for analyses 
C <- decostand(X[,SP],"hell")
E <- X[,ENV]
Z <- as.matrix(X[,CONT])

#Run the ordination model
rda.mod <- rda(C ~ ., as.data.frame(E))
rda.mod
summary(rda.mod)$concont

# Run a permutation ANOVA
anova(rda.mod,by="terms",permu=200)
}



# Ordination plot
colvec.R = c("darkgreen","chartreuse3")

##-- Plot species and variables
##-----------------------------
setwd(to.figs)
pdf("RDA_species.pdf",useDingbats = F)
ordiplot(rda.mod,display="site",type="n",xlim=c(-1.5,1.6),cex.axis=1.4,cex.lab=1.4)
points(rda.mod, display="cn", col="gray20",cex=0.0001,lwd=2)
text(rda.mod, display="sp", col=colvec,cex=1.5,lwd=2)
text(x=-1.3,y=-0.95,"T",col="gray20",cex=1.5)
text(x=-0.55,y=0.98,"C",col="gray20",cex=1.5)
text(x=-0.48,y=1.05,"V",col="gray20",cex=1.5)
text(x=0.48,y=-0.9,"D",col="gray20",cex=1.5)
dev.off()


coord <- scores(rda.mod, display = c("sites"))

library(grDevices)
library(car)
drawHull = function(xy,colo,pch=19,lty0=1,draw=T,ellips = F){
  df = data.frame(X=xy[,1], Y= xy[,2])
  con.hull.pos = chull(df) # find positions of convex hull
  con.hull <- rbind(df[con.hull.pos,],df[con.hull.pos[1],]) # get coordinates for convex hull
  if(draw==T){
    if(ellips == F){
      polygon(con.hull[,1],con.hull[,2],col = alpha(colo,0.1),lty=0)
      lines(con.hull,col=colo,lwd=1,lty=lty0) # add lines for convex hull
    }else{
      x=xy[,1] ; y=xy[,2]
      dataEllipse(x, y,levels=c(0.95), center.pch=19,add = T,plot.points=F,col= colo,center.cex = 0) 
    }
  
  points(df[,1],df[,2], cex=0.8,pch=pch,col=alpha(colo,0.7))
  #points(mean(df[,1]),mean(df[,2]),pch=4,lwd=2,col="black")
  }
  return(c(mean(df[,1]),mean(df[,2])))
}


## Figures decomposing RDA site representation
##............................................
setwd(to.figs)

## Isolated versus connected
##--------------------------
setwd(to.figs)
pdf("RDA_conn_vs_isol_Time_split.pdf",useDingbats = F,width = 10)
{
  coordAi7 = subset(coord[which(X$Treatment == "Isolated"  & X$day == 7),])
  coordAi15 = subset(coord[which(X$Treatment == "Isolated"  & X$day == 15),])
  coordAi21 = subset(coord[which(X$Treatment == "Isolated" & X$day == 21),])
  coordAi29 = subset(coord[which(X$Treatment == "Isolated" & X$day == 29),])
  
  coordAc7 = subset(coord[which(X$Treatment == "Connected"  & X$day == 7),])
  coordAc15 = subset(coord[which(X$Treatment == "Connected"  & X$day == 15),])
  coordAc21 = subset(coord[which(X$Treatment == "Connected" & X$day == 21),])
  coordAc29 = subset(coord[which(X$Treatment == "Connected" & X$day == 29),])
  
  centroidsI = centroidsC =matrix(nrow=4,ncol=2,dimnames = list(NULL,c("x","y")))
  
  
  
  layout(matrix(1:2,ncol=2),widths = c(1.2,1))
  
  par(mar=c(5,5,5,0))
  ordiplot(rda.mod,display="site",type="n",xlim=c(-0.4,0.35),cex.axis=1.4,cex.lab=1.4,main="Isolated")
  
  centroidsI[1,] = drawHull(coordAi7,alpha(colvec.R[2]),lty0=2)
  centroidsI[2,] = drawHull(coordAi15,alpha(colvec.R[2]),draw=FALSE)
  centroidsI[3,] = drawHull(coordAi21,alpha(colvec.R[2]),draw=FALSE)
  centroidsI[4,] = drawHull(coordAi29,alpha(colvec.R[2]),pch=17)
  points(centroidsI[,1],centroidsI[,2],pch=4,cex=1,lwd=2, col = colvec.R[2])
  for(j in 2:nrow(centroidsI)) arrows(x0 =centroidsI[j-1,"x"],x1=centroidsI[j,"x"],y0 =centroidsI[j-1,"y"],y1=centroidsI[j,"y"],length=0.1,lwd=2,col="black")
  legend("topleft",legend=c("day 7","day 29"),pch=c(19,17),lty=c(2,1),bty="n",cex=1)
  
  par(mar=c(5,0,5,1))
  ordiplot(rda.mod,display="site",type="n",xlim=c(-0.4,0.35),cex.axis=1.4,cex.lab=1.4,main="Connected",yaxt="n")
  centroidsC[1,] = drawHull(coordAc7,alpha(colvec.R[1]),lty0=2)
  centroidsC[2,] = drawHull(coordAc15,alpha(colvec.R[1]),draw=FALSE)
  centroidsC[3,] = drawHull(coordAc21,alpha(colvec.R[1]),draw=FALSE)
  centroidsC[4,] = drawHull(coordAc29,alpha(colvec.R[1]),pch=17)
  points(centroidsC[,1],centroidsC[,2],pch=4,cex=1,lwd=2, col = colvec.R[1])
  for(j in 2:nrow(centroidsC)) arrows(x0 =centroidsC[j-1,"x"],x1=centroidsC[j,"x"],y0 =centroidsC[j-1,"y"],y1=centroidsC[j,"y"],length=0.1,lwd=2,col="black")

}
dev.off()

## Decomposition by landscape (replicate)
##---------------------------------------
pdf("RDA_byLand_conn_vs_isol_Time_split.pdf",useDingbats = T,height = 10)
{
land = unique(X$Replicate)
for(i in 1:length(land)){
  
  coordAi7 = subset(coord[which(X$Replicate == land[i] & X$Treatment == "Isolated"  & X$day == 7),])
  coordAi15 = subset(coord[which(X$Replicate == land[i] & X$Treatment == "Isolated"  & X$day == 15),])
  coordAi21 = subset(coord[which(X$Replicate == land[i] & X$Treatment == "Isolated" & X$day == 21),])
  coordAi29 = subset(coord[which(X$Replicate == land[i] & X$Treatment == "Isolated" & X$day == 29),])
  
  coordAc7 = subset(coord[which(X$Replicate == land[i] & X$Treatment == "Connected"  & X$day == 7),])
  coordAc15 = subset(coord[which(X$Replicate == land[i] & X$Treatment == "Connected"  & X$day == 15),])
  coordAc21 = subset(coord[which(X$Replicate == land[i] & X$Treatment == "Connected" & X$day == 21),])
  coordAc29 = subset(coord[which(X$Replicate == land[i] & X$Treatment == "Connected" & X$day == 29),])
  
  centroidsI = centroidsC =matrix(nrow=4,ncol=2,dimnames = list(NULL,c("x","y")))
  
  layout(matrix(1:2,nrow=2),heights = c(1,1))
  
  par(mar=c(0,5,5,1))
  ordiplot(rda.mod,display="site",type="n",xlim=c(-0.4,0.35),main = paste0("Landscape ", land[i], " - Isolated"),cex.axis=1.4,cex.lab=1.4,xaxt="n")
  centroidsI[1,] = drawHull(coordAi7,alpha(colvec.R[2]),lty0=2)
  centroidsI[2,] = drawHull(coordAi15,alpha(colvec.R[2]),draw=FALSE)
  centroidsI[3,] = drawHull(coordAi21,alpha(colvec.R[2]),draw=FALSE)
  centroidsI[4,] = drawHull(coordAi29,alpha(colvec.R[2]),pch=17)
  points(centroidsI[,1],centroidsI[,2],pch=4,cex=1,lwd=2, col = colvec.R[2])
  for(j in 2:nrow(centroidsI)) arrows(x0 =centroidsI[j-1,"x"],x1=centroidsI[j,"x"],y0 =centroidsI[j-1,"y"],y1=centroidsI[j,"y"],length=0.1,lwd=2,col="black")
  legend("topleft",legend=c("day 7","day 29"),pch=c(19,17),lty=c(2,1),bty="n",cex=1)
  
  par(mar=c(5,5,0,1))
  ordiplot(rda.mod,display="site",type="n",xlim=c(-0.4,0.35),cex.axis=1.4,cex.lab=1.4) #,main = paste0("Landscape ", land[i], " - Connected")
  centroidsC[1,] = drawHull(coordAc7,alpha(colvec.R[1]),lty0=2)
  centroidsC[2,] = drawHull(coordAc15,alpha(colvec.R[1]),draw=FALSE)
  centroidsC[3,] = drawHull(coordAc21,alpha(colvec.R[1]),draw=FALSE)
  centroidsC[4,] = drawHull(coordAc29,alpha(colvec.R[1]),pch=17)
  points(centroidsC[,1],centroidsC[,2],pch=4,cex=1,lwd=2, col = colvec.R[1])
  for(j in 2:nrow(centroidsC)) arrows(x0 =centroidsC[j-1,"x"],x1=centroidsC[j,"x"],y0 =centroidsC[j-1,"y"],y1=centroidsC[j,"y"],length=0.1,lwd=2,col="black")
}

}
dev.off()

## Decomposition by Patch size
##----------------------------
pdf("RDA_bySize_conn_vs_isol_Time_split.pdf", useDingbats = T,height = 10)
{
  sizes = unique(X$Size)
  for(i in 1:length(sizes)){
    
    coordAi7 = subset(coord[which(X$Size == sizes[i] & X$Treatment == "Isolated"  & X$day == 7),])
    coordAi15 = subset(coord[which(X$Size == sizes[i] & X$Treatment == "Isolated"  & X$day == 15),])
    coordAi21 = subset(coord[which(X$Size == sizes[i] & X$Treatment == "Isolated" & X$day == 21),])
    coordAi29 = subset(coord[which(X$Size == sizes[i] & X$Treatment == "Isolated" & X$day == 29),])
    
    coordAc7 = subset(coord[which(X$Size == sizes[i] & X$Treatment == "Connected"  & X$day == 7),])
    coordAc15 = subset(coord[which(X$Size == sizes[i] & X$Treatment == "Connected"  & X$day == 15),])
    coordAc21 = subset(coord[which(X$Size == sizes[i] & X$Treatment == "Connected" & X$day == 21),])
    coordAc29 = subset(coord[which(X$Size == sizes[i] & X$Treatment == "Connected" & X$day == 29),])
    
    centroidsI = centroidsC =matrix(nrow=4,ncol=2,dimnames = list(NULL,c("x","y")))
    
    layout(matrix(1:2,nrow=2),heights = c(1,1))
    
    par(mar=c(0,5,5,1))
    ordiplot(rda.mod,display="site",type="n",xlim=c(-0.4,0.4),main = paste0("Patch size ", sizes[i], " - Isolated"),cex.axis=1.4,cex.lab=1.4,xaxt="n")
    centroidsI[1,] = drawHull(coordAi7,alpha(colvec.R[2]),lty0=2)
    centroidsI[2,] = drawHull(coordAi15,alpha(colvec.R[2]),draw=FALSE)
    centroidsI[3,] = drawHull(coordAi21,alpha(colvec.R[2]),draw=FALSE)
    centroidsI[4,] = drawHull(coordAi29,alpha(colvec.R[2]),pch=17)
    points(centroidsI[,1],centroidsI[,2],pch=4,cex=1,lwd=2, col = colvec.R[2])
    for(j in 2:nrow(centroidsI)) arrows(x0 =centroidsI[j-1,"x"],x1=centroidsI[j,"x"],y0 =centroidsI[j-1,"y"],y1=centroidsI[j,"y"],length=0.1,lwd=2,col="black")
    legend("topleft",legend=c("day 7","day 29"),pch=c(19,17),lty=c(2,1),bty="n",cex=1)
    
    
    
    par(mar=c(5,5,0,1))
    ordiplot(rda.mod,display="site",type="n",xlim=c(-0.4,0.4),cex.axis=1.4,cex.lab=1.4) #,main = paste0("Patch size ", sizes[i], " - Connected")
    centroidsC[1,] = drawHull(coordAc7,alpha(colvec.R[1]),lty0=2)
    centroidsC[2,] = drawHull(coordAc15,alpha(colvec.R[1]),draw=FALSE)
    centroidsC[3,] = drawHull(coordAc21,alpha(colvec.R[1]),draw=FALSE)
    centroidsC[4,] = drawHull(coordAc29,alpha(colvec.R[1]),pch=17)
    points(centroidsC[,1],centroidsC[,2],pch=4,cex=1,lwd=2, col = colvec.R[1])
    for(j in 2:nrow(centroidsC)) arrows(x0 =centroidsC[j-1,"x"],x1=centroidsC[j,"x"],y0 =centroidsC[j-1,"y"],y1=centroidsC[j,"y"],length=0.1,lwd=2,col="black")
  }
  
}
dev.off()

#################################
### LRR - days 7 to 29
################################

#-- Get the LRR data for protists for Day 7 to 29 all patch size pooled
#----------------------------------------------------------------------
{
#Create output file
LRR_Res_protist <- data.frame(species=character(0),LRR=numeric(0),CI_upper=numeric(0),CI_lower=numeric(0))

#Loop for each species
for(i in 1:length(SP)) { 
  
  # Extract info for SP[i] and summarise data
  D_Res_protist <- data %>% 
    filter(day !=0) %>% 
    select(SP[i],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
    mutate(Treatment = as.factor(Treatment),
           Replicate = as.factor(Replicate),
           Size = as.factor(Size))  %>% 
    Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment")) #get the mean by treatment level for each species
  
  LRR <- log(D_Res_protist[1,3]) - log(D_Res_protist[2,3]) 
  SE_LRR <- sqrt( (D_Res_protist$se[1]^2/D_Res_protist[1,3]^2) + (D_Res_protist$se[2]^2/D_Res_protist[2,3]^2) ) #From Hedge et al., 1999
  CI_upper <- LRR + 1.96*SE_LRR
  CI_lower <- LRR - 1.96*SE_LRR
  #CI_both <- LRR + c(-1, 1) * qnorm(1 - (1 - 0.95)/2) * SE_LRR #from ARPobservation:::logRespRatio
  species <- SP[i]
  
  bind.out_Res_protist <- data.frame(species,LRR,CI_upper,CI_lower)
  #colnames(bind.out) <- c("time","size","LRR","species","CI_upper","CI_lower")
  
  LRR_Res_protist <- rbind.data.frame(LRR_Res_protist,bind.out_Res_protist)
  
}

#Data calculated later in the document (added here for plotting purpose only - see below)
#LRR_bact_green <- data.frame(species="Bacteria",LRR=0.2313328,CI_upper=0.3173062,CI_lower=0.1453594)


}
LRR_Res_protist



#-- Get the LRR data for BACTERIA IN GREEN for Day 7 to 29 all patches pooled
#-------------------------------------------------------------------
{
  ####Generate data for LRR calculation
  D_Res_green <- data %>% 
    filter(day !=0 & Replicate %in% c("A","B")) %>% 
    Rmisc::summarySE(.,measurevar = "bact.green", groupvars = c("Treatment"))
  
  #Calculate LRR
  species <- "bacteriaG"
  LRR <- log(D_Res_green[1,3]) - log(D_Res_green[2,3]) 
  SE_LRR <- sqrt( (D_Res_green$se[1]^2/D_Res_green[1,3]^2) + (D_Res_green$se[2]^2/D_Res_green[2,3]^2) ) #From Hedge et al., 1999
  CI_upper <- LRR + 1.96*SE_LRR
  CI_lower <- LRR - 1.96*SE_LRR
}
bind.out_Res_green <- data.frame(species,LRR,CI_upper,CI_lower)

#-- Get the LRR data for BACTERIA IN BLUE for Day 7 to 29 all patches pooled
#-------------------------------------------------------------------
{
  ####Generate data for LRR calculation
  D_Res_blue <- data %>% 
    filter(day !=0 & Replicate %in% c("A","B")) %>% 
    Rmisc::summarySE(.,measurevar = "bact.blue", groupvars = c("Treatment"))
  
  #Calculate LRR
  species <- "bacteriaB"
  LRR <- log(D_Res_blue[1,3]) - log(D_Res_blue[2,3]) 
  SE_LRR <- sqrt( (D_Res_blue$se[1]^2/D_Res_blue[1,3]^2) + (D_Res_blue$se[2]^2/D_Res_blue[2,3]^2) ) #From Hedge et al., 1999
  CI_upper <- LRR + 1.96*SE_LRR
  CI_lower <- LRR - 1.96*SE_LRR
}
bind.out_Res_blue <- data.frame(species,LRR,CI_upper,CI_lower)

dat <- rbind(LRR_Res_protist,bind.out_Res_blue,bind.out_Res_green)

##-- Plot the results
##-------------------
setwd(to.figs)
layout(1)
pdf("LRR_species_7to29.pdf",useDingbats = F,family = "ArialMT",width = 9)
{
plot(c(1:8,10),LRR_sp$LRR,type="n",xlim=c(0.5,10.5),xlab = "Species",ylab = "Effect size (log ratio of mean densities)",ylim=c(-0.8,0.32), yaxt="n",xaxt="n",cex.lab=1.4,main ="Day 7 to 29")
axis(2,at = c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4),labels = c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4),cex.axis=1.2,las=2)
axis(1,at = c(1:8,10),labels = levels(LRR_sp$species),cex.axis=1.2)
abline(h=0,lty=3); abline(v=9,lwd=1)
segments(x0=c(1:8,10),x1=c(1:8,10),y0=LRR_sp$CI_upper,y1=LRR_sp$CI_lower,col=c(colvec,"grey60","darkgreen"),lwd=3)
points(c(1:8,10),LRR_sp$LRR,col=c(colvec,"grey60","darkgreen"),pch=20,cex=2.5)
}
dev.off()




###############################################################
#-- BY PATCH SIZE
#--------------------------------------------------------------
{
  #Create output file
LRR_Res_Vol_protist <- data.frame(size=numeric(0),species=character(0),LRR=numeric(0),CI_upper=numeric(0),CI_lower=numeric(0))

#Loop for each species

for(i in 1:length(SP)) { 
  
  # Extract info for SP[i] and summarise data
  D_Res_Vol_protist <- data %>% 
    filter(day !=0) %>% 
    select(SP[i],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
    mutate(Treatment = as.factor(Treatment),
           Replicate = as.factor(Replicate),
           Size = as.factor(Size))  %>% 
    Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment","Size"))
  
  LRR.7 <- log(D_Res_Vol_protist[1,4]) - log(D_Res_Vol_protist[5,4]) 
  LRR.13 <- log(D_Res_Vol_protist[2,4]) - log(D_Res_Vol_protist[6,4]) 
  LRR.22 <- log(D_Res_Vol_protist[3,4]) - log(D_Res_Vol_protist[7,4]) 
  LRR.45 <- log(D_Res_Vol_protist[4,4]) - log(D_Res_Vol_protist[8,4]) 
  LRR <- rbind(LRR.7,LRR.13,LRR.22,LRR.45)
  
  SE_LRR.7 <- sqrt( (D_Res_Vol_protist$se[1]^2/D_Res_Vol_protist[1,4]^2) + (D_Res_Vol_protist$se[5]^2/D_Res_Vol_protist[5,4]^2) )
  SE_LRR.13 <- sqrt( (D_Res_Vol_protist$se[2]^2/D_Res_Vol_protist[2,4]^2) + (D_Res_Vol_protist$se[6]^2/D_Res_Vol_protist[6,4]^2) )
  SE_LRR.22 <- sqrt( (D_Res_Vol_protist$se[3]^2/D_Res_Vol_protist[3,4]^2) + (D_Res_Vol_protist$se[7]^2/D_Res_Vol_protist[7,4]^2) )
  SE_LRR.45 <- sqrt( (D_Res_Vol_protist$se[4]^2/D_Res_Vol_protist[4,4]^2) + (D_Res_Vol_protist$se[8]^2/D_Res_Vol_protist[8,4]^2) )
  SE_LRR <- rbind(SE_LRR.7,SE_LRR.13,SE_LRR.22,SE_LRR.45)
  
  CI_upper.7 <- LRR.7 + 1.96*SE_LRR.7
  CI_upper.13 <- LRR.13 + 1.96*SE_LRR.13
  CI_upper.22 <- LRR.22 + 1.96*SE_LRR.22
  CI_upper.45 <- LRR.45 + 1.96*SE_LRR.45
  CI_upper <- rbind(CI_upper.7,CI_upper.13,CI_upper.22,CI_upper.45)
  
  CI_lower.7 <- LRR.7 - 1.96*SE_LRR.7
  CI_lower.13 <- LRR.13 - 1.96*SE_LRR.13
  CI_lower.22 <- LRR.22 - 1.96*SE_LRR.22
  CI_lower.45 <- LRR.45 - 1.96*SE_LRR.45
  CI_lower <- rbind(CI_lower.7,CI_lower.13,CI_lower.22,CI_lower.45)
  
  species <- rep(SP[i],4)
  size <- unique(D_Res_Vol_protist$Size)
  
  bind.out_Res_Vol_protist <- data.frame(size,species,LRR,CI_upper,CI_lower)
  #colnames(bind.out) <- c("time","size","LRR","species","CI_upper","CI_lower")
  
  LRR_Res_Vol_protist <- rbind.data.frame(LRR_Res_Vol_protist,bind.out_Res_Vol_protist)
  
}
}
LRR_Res_Vol_protist


#-- Get the LRR data for BACTERIA for Day 7 to 29 by patch size
#---------------------------------------------------------------

{
D_Res_Vol_green <- data %>% 
  filter(day !=0 & Replicate %in% c("A","B")) %>% 
  Rmisc::summarySE(.,measurevar = "bact.green", groupvars = c("Treatment","Size"))


#Calculate LRR
species <- "bacteria"
LRR.7 <- log(D_Res_Vol_green[1,4]) - log(D_Res_Vol_green[5,4]) 
LRR.13 <- log(D_Res_Vol_green[2,4]) - log(D_Res_Vol_green[6,4])
LRR.22 <- log(D_Res_Vol_green[3,4]) - log(D_Res_Vol_green[7,4])
LRR.45 <- log(D_Res_Vol_green[4,4]) - log(D_Res_Vol_green[8,4])
LRR <- rbind(LRR.7,LRR.13,LRR.22,LRR.45)

SE_LRR.7 <- sqrt( (D_Res_Vol_green$se[1]^2/D_Res_Vol_green[1,4]^2) + (D_Res_Vol_green$se[5]^2/D_Res_Vol_green[5,4]^2) ) #From Hedge et al., 1999
SE_LRR.13 <- sqrt( (D_Res_Vol_green$se[2]^2/D_Res_Vol_green[2,4]^2) + (D_Res_Vol_green$se[6]^2/D_Res_Vol_green[6,4]^2) )
SE_LRR.22 <- sqrt( (D_Res_Vol_green$se[3]^2/D_Res_Vol_green[3,4]^2) + (D_Res_Vol_green$se[7]^2/D_Res_Vol_green[7,4]^2) )
SE_LRR.45 <- sqrt( (D_Res_Vol_green$se[4]^2/D_Res_Vol_green[4,4]^2) + (D_Res_Vol_green$se[8]^2/D_Res_Vol_green[8,4]^2) )
SE_LRR <- rbind(SE_LRR.7,SE_LRR.13,SE_LRR.22,SE_LRR.45)

CI.7 <- 1.96*SE_LRR.7
CI.13 <- 1.96*SE_LRR.13
CI.22 <- 1.96*SE_LRR.22
CI.45 <- 1.96*SE_LRR.45
CI <- rbind(CI.7,CI.13,CI.22,CI.45)

size <- unique(D_Res_Vol_green$Size)
}
bind.out_Res_Vol_green <- data.frame(size,LRR,CI)


#-- Get the LRR data for BACTERIA BLUE for Day 7 to 29 by patch size
#-------------------------------------------------------------------

{
  D_Res_Vol_blue <- data %>% 
    filter(day !=0 & Replicate %in% c("A","B")) %>% 
    Rmisc::summarySE(.,measurevar = "bact.blue", groupvars = c("Treatment","Size"))
  
  
  #Calculate LRR
  species <- "bacteria"
  LRR.7 <- log(D_Res_Vol_blue[1,4]) - log(D_Res_Vol_blue[5,4]) 
  LRR.13 <- log(D_Res_Vol_blue[2,4]) - log(D_Res_Vol_blue[6,4])
  LRR.22 <- log(D_Res_Vol_blue[3,4]) - log(D_Res_Vol_blue[7,4])
  LRR.45 <- log(D_Res_Vol_blue[4,4]) - log(D_Res_Vol_blue[8,4])
  LRR <- rbind(LRR.7,LRR.13,LRR.22,LRR.45)
  
  SE_LRR.7 <- sqrt( (D_Res_Vol_blue$se[1]^2/D_Res_Vol_blue[1,4]^2) + (D_Res_Vol_blue$se[5]^2/D_Res_Vol_blue[5,4]^2) ) #From Hedge et al., 1999
  SE_LRR.13 <- sqrt( (D_Res_Vol_blue$se[2]^2/D_Res_Vol_blue[2,4]^2) + (D_Res_Vol_blue$se[6]^2/D_Res_Vol_blue[6,4]^2) )
  SE_LRR.22 <- sqrt( (D_Res_Vol_blue$se[3]^2/D_Res_Vol_blue[3,4]^2) + (D_Res_Vol_blue$se[7]^2/D_Res_Vol_blue[7,4]^2) )
  SE_LRR.45 <- sqrt( (D_Res_Vol_blue$se[4]^2/D_Res_Vol_blue[4,4]^2) + (D_Res_Vol_blue$se[8]^2/D_Res_Vol_blue[8,4]^2) )
  SE_LRR <- rbind(SE_LRR.7,SE_LRR.13,SE_LRR.22,SE_LRR.45)
  
  CI.7 <- 1.96*SE_LRR.7
  CI.13 <- 1.96*SE_LRR.13
  CI.22 <- 1.96*SE_LRR.22
  CI.45 <- 1.96*SE_LRR.45
  CI <- rbind(CI.7,CI.13,CI.22,CI.45)
  
  size <- unique(D_Res_Vol_blue$Size)
}
bind.out_Res_Vol_blue <- data.frame(size,LRR,CI)



#--- Pool protist and bacteria LRR data per patch size
dat=rbind(LRR_Res_Vol_protist,
          data.frame("size" = bind.out_Res_Vol_blue$size,"species"=rep("BacteriaB",4),"LRR"=bind.out_Res_Vol_blue$LRR,"CI_upper"=bind.out_Res_Vol_blue$LRR+bind.out_Res_Vol_blue$CI,"CI_lower"=bind.out_Res_Vol_blue$LRR-bind.out_Res_Vol_blue$CI),
          data.frame("size" = bind.out_Res_Vol_green$size,"species"=rep("BacteriaG",4),"LRR"=bind.out_Res_Vol_green$LRR,"CI_upper"=bind.out_Res_Vol_green$LRR+bind.out_Res_Vol_green$CI,"CI_lower"=bind.out_Res_Vol_green$LRR-bind.out_Res_Vol_green$CI))

dat

##-- Plot the results with detailed on patch size
##-----------------------------------------------
setwd(to.figs)

#--- positions on the x-axis
pos7.5 = c(1,7,13,19,25,31,37,43,50)
labelPos = c(2,8,14,20,26,32,38,44,51)+0.5

#--- patch sizes and size of point to display patch size
sizes = levels(dat$size)
cexVol= c(0.8,1.3,1.8,2.5)

pdf("LRR_species_Sizes_7to29.pdf",useDingbats = F,family = "ArialMT",width = 10)
{
plot(NA,type="n",xlim=c(0,54),col=c(colvec,"black"),xlab = "Species",ylab = "Effect size (log ratio of mean densities)",ylim=c(-1.4,0.8), yaxt="n",xaxt="n",cex.lab=1.4,main ="Day 7 to 29")
axis(2,at = c(-1.2,-0.8,-0.4,0,0.4,0.8),labels = c(-1.2,-0.8,-0.4,0,0.4,0.8),cex.axis=1.2,las=2)
axis(1,at = labelPos,labels = levels(dat$species),cex.axis=1.2)
abline(h=0,lty=3); abline(v=48,lwd=1)

for(i in 1:length(levels(LRR_size$size))){
  subdat = subset(LRR_size,size==sizes[i])
  segments(x0=c(pos7.5+i-1),c(pos7.5+i-1),y0=subdat$CI_upper,y1=subdat$CI_lower,col=c(colvec,"grey60","darkgreen"),lwd=2)
  points(c(pos7.5+i-1),subdat$LRR,col=c(colvec,"grey60","darkgreen"),pch=20,cex=cexVol[i])
}
legend("bottomleft",lty=rep(1,4),pch=20,pt.cex=cexVol,lwd=2,legend = sizes,title="Patch sizes",bty = "n",cex=0.8)
}
dev.off()



#################################
### LRR - days 29
################################

#-- Get the LRR data for protists for day 29 all patch size pooled
#-----------------------------------------------------------------
{
#Create output file
LRR_Res_protist_29 <- data.frame(species=character(0),LRR=numeric(0),CI_upper=numeric(0),CI_lower=numeric(0))

#Loop for each species

for(i in 1:length(SP)) { 
  
  # Extract info for SP[i] and summarise data
  D_Res_protist <- data %>% 
    filter(day ==29) %>% 
    select(SP[i],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
    mutate(Treatment = as.factor(Treatment),
           Replicate = as.factor(Replicate),
           Size = as.factor(Size))  %>% 
    Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment")) #get the mean by treatment level for each species
  
  LRR <- log(D_Res_protist[1,3]) - log(D_Res_protist[2,3]) 
  SE_LRR <- sqrt( (D_Res_protist$se[1]^2/D_Res_protist[1,3]^2) + (D_Res_protist$se[2]^2/D_Res_protist[2,3]^2) ) #From Hedge et al., 1999
  CI_upper <- LRR + 1.96*SE_LRR
  CI_lower <- LRR - 1.96*SE_LRR
  #CI_both <- LRR + c(-1, 1) * qnorm(1 - (1 - 0.95)/2) * SE_LRR #from ARPobservation:::logRespRatio
  species <- SP[i]
  
  bind.out_Res_protist <- data.frame(species,LRR,CI_upper,CI_lower)
  #colnames(bind.out) <- c("time","size","LRR","species","CI_upper","CI_lower")
  
  LRR_Res_protist_29 <- rbind.data.frame(LRR_Res_protist_29,bind.out_Res_protist)
  
}
}
LRR_Res_protist_29

#-- Get the LRR data for BACTERIA for Day 29 all patches pooled
#--------------------------------------------------------------

{
  ####Generate data for LRR calculation
D_Res_green <- data %>% 
  filter(day ==29 & Replicate %in% c("A","B")) %>% 
  Rmisc::summarySE(.,measurevar = "bact.green", groupvars = c("Treatment"))

#Calculate LRR
species <- "bacteriaG"
LRR <- log(D_Res_green[1,3]) - log(D_Res_green[2,3]) 
SE_LRR <- sqrt( (D_Res_green$se[1]^2/D_Res_green[1,3]^2) + (D_Res_green$se[2]^2/D_Res_green[2,3]^2) ) #From Hedge et al., 1999
CI_upper <- LRR + 1.96*SE_LRR
CI_lower <- LRR - 1.96*SE_LRR

bind.out_Res_green_29 <- data.frame(species,LRR,CI_upper,CI_lower)
}
bind.out_Res_green_29

#-- Get the LRR data for BACTERIA IN BLUE for Day 7 to 29 all patches pooled
#-------------------------------------------------------------------
{
  ####Generate data for LRR calculation
  D_Res_blue <- data %>% 
    filter(day ==29 & Replicate %in% c("A","B")) %>% 
    Rmisc::summarySE(.,measurevar = "bact.blue", groupvars = c("Treatment"))
  
  #Calculate LRR
  species <- "bacteriaB"
  LRR <- log(D_Res_blue[1,3]) - log(D_Res_blue[2,3]) 
  SE_LRR <- sqrt( (D_Res_blue$se[1]^2/D_Res_blue[1,3]^2) + (D_Res_blue$se[2]^2/D_Res_blue[2,3]^2) ) #From Hedge et al., 1999
  CI_upper <- LRR + 1.96*SE_LRR
  CI_lower <- LRR - 1.96*SE_LRR
}
bind.out_Res_blue_29 <- data.frame(species,LRR,CI_upper,CI_lower)

dat <- rbind(LRR_Res_protist_29,bind.out_Res_blue_29,bind.out_Res_green_29)

dat

##-- Plot the results
##-------------------
max(dat$CI_upper)
min(dat$CI_lower)

setwd(to.figs)
layout(1)
pdf("LRR_species_day29.pdf",useDingbats = F,family = "ArialMT",width = 9)
{
  plot(c(1:8,10),dat$LRR,type="n",xlim=c(0.5,10.5),xlab = "Species",ylab = "Effect size (log ratio of mean densities)",ylim=c(-1.7,0.9), yaxt="n",xaxt="n",cex.lab=1.4,main ="day 29")
  axis(2,at = c(-1.6,-1.2,-0.8,-0.4,0,0.4,0.8),labels = c(-1.6,-1.2,-0.8,-0.4,0,0.4,0.8),cex.axis=1.2,las=2)
  axis(1,at = c(1:8,10),labels = levels(dat$species),cex.axis=1.2)
  abline(h=0,lty=3); abline(v=9,lwd=1)
  segments(x0=c(1:8,10),x1=c(1:8,10),y0=dat$CI_upper,y1=dat$CI_lower,col=c(colvec,"grey60","darkgreen"),lwd=3)
  points(c(1:8,10),dat$LRR,col=c(colvec,"grey60","darkgreen"),pch=20,cex=2.5)
}
dev.off()



###############################################################
#-- BY PATCH SIZE
#--------------------------------------------------------------


#-- Get the LRR data for protists for day 29 by patch size
#----------------------------------------------------------
{
#Create output file
LRR_Res_Vol_protist_29 <- data.frame(size=numeric(0),species=character(0),LRR=numeric(0),CI_upper=numeric(0),CI_lower=numeric(0))

#Loop for each species

for(i in 1:length(SP)) { 
  
  # Extract info for SP[i] and summarise data
  D_Res_Vol_protist <- data %>% 
    filter(day ==29) %>% 
    select(SP[i],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
    mutate(Treatment = as.factor(Treatment),
           Replicate = as.factor(Replicate),
           Size = as.factor(Size))  %>% 
    Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment","Size"))
  
  LRR.7 <- log(D_Res_Vol_protist[1,4]) - log(D_Res_Vol_protist[5,4]) 
  LRR.13 <- log(D_Res_Vol_protist[2,4]) - log(D_Res_Vol_protist[6,4]) 
  LRR.22 <- log(D_Res_Vol_protist[3,4]) - log(D_Res_Vol_protist[7,4]) 
  LRR.45 <- log(D_Res_Vol_protist[4,4]) - log(D_Res_Vol_protist[8,4]) 
  LRR <- rbind(LRR.7,LRR.13,LRR.22,LRR.45)
  
  SE_LRR.7 <- sqrt( (D_Res_Vol_protist$se[1]^2/D_Res_Vol_protist[1,4]^2) + (D_Res_Vol_protist$se[5]^2/D_Res_Vol_protist[5,4]^2) )
  SE_LRR.13 <- sqrt( (D_Res_Vol_protist$se[2]^2/D_Res_Vol_protist[2,4]^2) + (D_Res_Vol_protist$se[6]^2/D_Res_Vol_protist[6,4]^2) )
  SE_LRR.22 <- sqrt( (D_Res_Vol_protist$se[3]^2/D_Res_Vol_protist[3,4]^2) + (D_Res_Vol_protist$se[7]^2/D_Res_Vol_protist[7,4]^2) )
  SE_LRR.45 <- sqrt( (D_Res_Vol_protist$se[4]^2/D_Res_Vol_protist[4,4]^2) + (D_Res_Vol_protist$se[8]^2/D_Res_Vol_protist[8,4]^2) )
  SE_LRR <- rbind(SE_LRR.7,SE_LRR.13,SE_LRR.22,SE_LRR.45)
  
  CI_upper.7 <- LRR.7 + 1.96*SE_LRR.7
  CI_upper.13 <- LRR.13 + 1.96*SE_LRR.13
  CI_upper.22 <- LRR.22 + 1.96*SE_LRR.22
  CI_upper.45 <- LRR.45 + 1.96*SE_LRR.45
  CI_upper <- rbind(CI_upper.7,CI_upper.13,CI_upper.22,CI_upper.45)
  
  CI_lower.7 <- LRR.7 - 1.96*SE_LRR.7
  CI_lower.13 <- LRR.13 - 1.96*SE_LRR.13
  CI_lower.22 <- LRR.22 - 1.96*SE_LRR.22
  CI_lower.45 <- LRR.45 - 1.96*SE_LRR.45
  CI_lower <- rbind(CI_lower.7,CI_lower.13,CI_lower.22,CI_lower.45)
  
  species <- rep(SP[i],4)
  size <- unique(D_Res_Vol_protist$Size)
  
  bind.out_Res_Vol_protist <- data.frame(size,species,LRR,CI_upper,CI_lower)
  #colnames(bind.out) <- c("time","size","LRR","species","CI_upper","CI_lower")
  
  LRR_Res_Vol_protist_29 <- rbind.data.frame(LRR_Res_Vol_protist_29,bind.out_Res_Vol_protist)
  
}
}
LRR_Res_Vol_protist_29

##-- the to a define value the extinct species for visualization purpose
LRR_Res_Vol_protist_29$LRR = ifelse(is.infinite(LRR_Res_Vol_protist_29$LRR),-4.5,LRR_Res_Vol_protist_29$LRR)
#LRR_Res_Vol_protist_29$CI_upper = ifelse(is.nan(LRR_Res_Vol_protist_29$CI_upper),NA,LRR_Res_Vol_protist_29$CI_upper)
LRR_Res_Vol_protist_29

#-- Get the LRR data for BACTERIA for Day 29 by patch size
#---------------------------------------------------------
{
D_Res_Vol_green <- data %>% 
  filter(day ==29 & Replicate %in% c("A","B")) %>% 
  Rmisc::summarySE(.,measurevar = "bact.green", groupvars = c("Treatment","Size"))

#Calculate LRR
species <- "bacteria"
LRR.7 <- log(D_Res_Vol_green[1,4]) - log(D_Res_Vol_green[5,4]) 
LRR.13 <- log(D_Res_Vol_green[2,4]) - log(D_Res_Vol_green[6,4])
LRR.22 <- log(D_Res_Vol_green[3,4]) - log(D_Res_Vol_green[7,4])
LRR.45 <- log(D_Res_Vol_green[4,4]) - log(D_Res_Vol_green[8,4])
LRR <- rbind(LRR.7,LRR.13,LRR.22,LRR.45)

SE_LRR.7 <- sqrt( (D_Res_Vol_green$se[1]^2/D_Res_Vol_green[1,4]^2) + (D_Res_Vol_green$se[5]^2/D_Res_Vol_green[5,4]^2) ) #From Hedge et al., 1999
SE_LRR.13 <- sqrt( (D_Res_Vol_green$se[2]^2/D_Res_Vol_green[2,4]^2) + (D_Res_Vol_green$se[6]^2/D_Res_Vol_green[6,4]^2) )
SE_LRR.22 <- sqrt( (D_Res_Vol_green$se[3]^2/D_Res_Vol_green[3,4]^2) + (D_Res_Vol_green$se[7]^2/D_Res_Vol_green[7,4]^2) )
SE_LRR.45 <- sqrt( (D_Res_Vol_green$se[4]^2/D_Res_Vol_green[4,4]^2) + (D_Res_Vol_green$se[8]^2/D_Res_Vol_green[8,4]^2) )
SE_LRR <- rbind(SE_LRR.7,SE_LRR.13,SE_LRR.22,SE_LRR.45)

CI.7 <- 1.96*SE_LRR.7
CI.13 <- 1.96*SE_LRR.13
CI.22 <- 1.96*SE_LRR.22
CI.45 <- 1.96*SE_LRR.45
CI <- rbind(CI.7,CI.13,CI.22,CI.45)

size <- unique(D_Res_Vol_green$Size)
}
bind.out_Res_Vol_green_29 <- data.frame(size,LRR,CI)
bind.out_Res_Vol_green_29

#-- Get the LRR data for BACTERIA BLUE for Day 7 to 29 by patch size
#-------------------------------------------------------------------

{
  D_Res_Vol_blue <- data %>% 
    filter(day ==29 & Replicate %in% c("A","B")) %>% 
    Rmisc::summarySE(.,measurevar = "bact.blue", groupvars = c("Treatment","Size"))
  
  
  #Calculate LRR
  species <- "bacteria"
  LRR.7 <- log(D_Res_Vol_blue[1,4]) - log(D_Res_Vol_blue[5,4]) 
  LRR.13 <- log(D_Res_Vol_blue[2,4]) - log(D_Res_Vol_blue[6,4])
  LRR.22 <- log(D_Res_Vol_blue[3,4]) - log(D_Res_Vol_blue[7,4])
  LRR.45 <- log(D_Res_Vol_blue[4,4]) - log(D_Res_Vol_blue[8,4])
  LRR <- rbind(LRR.7,LRR.13,LRR.22,LRR.45)
  
  SE_LRR.7 <- sqrt( (D_Res_Vol_blue$se[1]^2/D_Res_Vol_blue[1,4]^2) + (D_Res_Vol_blue$se[5]^2/D_Res_Vol_blue[5,4]^2) ) #From Hedge et al., 1999
  SE_LRR.13 <- sqrt( (D_Res_Vol_blue$se[2]^2/D_Res_Vol_blue[2,4]^2) + (D_Res_Vol_blue$se[6]^2/D_Res_Vol_blue[6,4]^2) )
  SE_LRR.22 <- sqrt( (D_Res_Vol_blue$se[3]^2/D_Res_Vol_blue[3,4]^2) + (D_Res_Vol_blue$se[7]^2/D_Res_Vol_blue[7,4]^2) )
  SE_LRR.45 <- sqrt( (D_Res_Vol_blue$se[4]^2/D_Res_Vol_blue[4,4]^2) + (D_Res_Vol_blue$se[8]^2/D_Res_Vol_blue[8,4]^2) )
  SE_LRR <- rbind(SE_LRR.7,SE_LRR.13,SE_LRR.22,SE_LRR.45)
  
  CI.7 <- 1.96*SE_LRR.7
  CI.13 <- 1.96*SE_LRR.13
  CI.22 <- 1.96*SE_LRR.22
  CI.45 <- 1.96*SE_LRR.45
  CI <- rbind(CI.7,CI.13,CI.22,CI.45)
  
  size <- unique(D_Res_Vol_blue$Size)
}
bind.out_Res_Vol_blue_29 <- data.frame(size,LRR,CI)
bind.out_Res_Vol_blue_29



#--- Pool protist and bacteria LRR data per patch size
dat=rbind(LRR_Res_Vol_protist_29,
          data.frame("size" = bind.out_Res_Vol_blue_29$size,"species"=rep("BacteriaB",4),"LRR"=bind.out_Res_Vol_blue_29$LRR,"CI_upper"=bind.out_Res_Vol_blue_29$LRR+bind.out_Res_Vol_blue_29$CI,"CI_lower"=bind.out_Res_Vol_blue_29$LRR-bind.out_Res_Vol_blue_29$CI),
          data.frame("size" = bind.out_Res_Vol_green_29$size,"species"=rep("BacteriaG",4),"LRR"=bind.out_Res_Vol_green_29$LRR,"CI_upper"=bind.out_Res_Vol_green_29$LRR+bind.out_Res_Vol_green_29$CI,"CI_lower"=bind.out_Res_Vol_green_29$LRR-bind.out_Res_Vol_green_29$CI))

dat


##-- Plot the results
##-------------------
max(dat$CI_upper,na.rm=T)
min(dat$CI_lower,na.rm=T)

#--- positions on the x-axis
pos7.5 = c(1,7,13,19,25,31,37,43,50)
labelPos = c(2,8,14,20,26,32,38,44,51)+0.5

#--- patch sizes and size of point to display patch size
sizes = levels(dat$size)
cexVol= c(0.8,1.3,1.8,2.5)

setwd(to.figs)
pdf("LRR_species_Sizes_day29.pdf",useDingbats = F,family = "ArialMT",width = 10)
{
plot(NA,type="n",xlim=c(0,54),col=c(colvec,"black"),xlab = "Species",ylab = "Effect size (log ratio of means)",ylim=c(-4.5,3.3), yaxt="n",xaxt="n",cex.lab=1.4,main ="Day 29")
axis(2,at = c(-4.5,-3,-2,-1,0,1,2,3),labels = c("Extinction",-3,-2,-1,0,1,2,3),cex.axis=1.2,las=2)
axis(1,at = labelPos,labels = levels(dat$species),cex.axis=1.2)
abline(h=0,lty=3); abline(v=48,lwd=1)

for(i in 1:length(levels(dat$size))){
  positions = c(pos7.5+i-1)
  subdat = subset(dat,size==sizes[i])
  segments(x0=positions,positions,y0=subdat$CI_upper,y1=subdat$CI_lower,col=c(colvec,"grey60","darkgreen"),lwd=2)
  points(positions,subdat$LRR,col=c(colvec,"grey60","darkgreen"),pch=20,cex=cexVol[i])
  for(k in 1:length(positions)){if((subdat$LRR)[k] == -4.5)points(positions[k],-4.5,pch = 3,cex=1.4)}
}
legend("bottomleft",lty=rep(1,4),pch=20,pt.cex=cexVol,lwd=2,legend = sizes,title="Patch sizes",bty = "n",cex=0.8)
}
dev.off()





# 
# 
# 
# ## Figure Connectance
# 
# D_connectance <- data %>%
#   filter(day !=0) %>% 
#   mutate(Size = as.numeric(as.character(Size))) %>% 
#   Rmisc::summarySE(.,measurevar = "degree", groupvars = c("Size"))
# 
# ggplot(D_connectance,mapping=aes(x=Size,y=degree)) +
#   ylab("Connectance (degree)") + xlab("Patch volume (mL)") +
#   geom_errorbar(aes(ymin=degree-se, ymax=degree+se), width=0) +
#   geom_point() + theme_classic()
# 
# 
# # Statistical model
# 
# #Extract data 
# 
# X.bact <- data %>% 
#   filter(day!=0 & Replicate %in% c("A","B")) %>% 
#   mutate(Treatment = as.factor(Treatment),
#          day = factor(day,levels=c("7","15","21","29")),
#          day.cont = log(as.numeric(as.character(day))),
#          Size = as.numeric(as.character(Size)),
#          centrality = as.numeric(as.character(centrality)),
#          dist.outlet = as.numeric(as.character(dist.outlet)))
# 
# #Run model (no log transformation based on histogramm of residuals)
# #Also: here Time was not significant so it was removed from the fixed variables
# Mod4 <- nlme:::lme(bact.green ~ Treatment*Size + degree + bact.blue + prot.ab + prot.rich, 
#                    random = ~ as.factor(day)|Replicate, data=X.bact,
#                    method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
# summary(Mod4)$tTable
# 
# 
# ## Effects of pulse on Blue
# # Filter data (bacteria were only measured for replicates A and B)
# 
# X.bact <- data %>% 
#   filter(day!=0 & Replicate %in% c("A","B")) %>% 
#   mutate(Treatment = as.factor(Treatment),
#          day = factor(day,levels=c("7","15","21","29")),
#          day.cont = log(as.numeric(as.character(day))),
#          Size = as.numeric(as.character(Size)),
#          centrality = as.numeric(as.character(centrality)),
#          dist.outlet = as.numeric(as.character(dist.outlet)))
# 
# X.prot <- data %>% 
#   filter(day!=0 & Replicate %in% c("A","B","C","D")) %>% 
#   mutate(Treatment = as.factor(Treatment),
#          day = factor(day,levels=c("7","15","21","29")),
#          day.cont = log(as.numeric(as.character(day))),
#          Size = as.numeric(as.character(Size)),
#          centrality = as.numeric(as.character(centrality)),
#          dist.outlet = as.numeric(as.character(dist.outlet)))
# 
# #Run model (no log transformation based on histogramm of residuals)
# Mod1 <- nlme:::lme(bact.blue ~ Treatment*Size*day.cont + bact.green, 
#                    random = ~ 1|Replicate, data=X.bact,
#                    method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
# summary(Mod1)$tTable
# 
# 
# # Generate data frame with each predictor to extract predictions
# newdat <-  with(X.bact,
#                 expand.grid(Treatment=unique(Treatment),
#                             Size = c(min(Size),max(Size)),
#                             bact.green=c(min(bact.green),max(bact.green)),
#                             day.cont = c(min(day.cont),max(day.cont))))
# 
# X.1 <- X.bact %>% 
#   filter(day==29)
# 
# ggplot(mapping=aes(x=as.factor(Size),y=log(bact.blue)),data=X.1) +
#   geom_boxplot(col="blue",alpha=0.5) -> p1
# 
# ggplot(mapping=aes(x=day.cont,y=bact.blue),data=X.bact) +
#   geom_point(col="blue",alpha=0.5) +
#   geom_line(data=newdat,aes(y=predict(Mod1,newdata=newdat,level=0)),col="red") -> p2
# 
# grid.arrange(p1,p2,ncol=2)
# 
# # Effects of pulse on blue density
# #Run model (no log transformation based on histogramm of residuals)
# Mod2 <- nlme:::lme(prot.ab ~ Treatment*Size*day.cont, 
#                    random = ~ 1|Replicate, data=X.prot,
#                    method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
# summary(Mod2)$tTable
# 
# # Plot results
# 
# ggplot(mapping=aes(x=as.factor(Size),y=prot.ab),data=X.prot) +
#   geom_boxplot() + facet_wrap(~ day)   
# 
# 
# 
# ggplot(mapping=aes(x=day,y=prot.ab),data=X.prot) +
#   geom_boxplot(col="blue",alpha=0.5) 
# 
# 
# 
# ## Effects of pulse on richness in blue
# 
# #Run model (no log transformation based on histogramm of residuals)
# Mod3 <- nlme:::lme(prot.rich ~ Treatment*Size*day.cont, 
#                    random = ~ 1|Replicate, data=X.prot,
#                    method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
# summary(Mod3)$tTable
# 
# # Plot results
# 
# X1.1 <- X.prot %>%
#   filter(day==29)
# ggplot(mapping=aes(x=as.factor(Size),y=prot.rich),data=X1.1) +
#   geom_boxplot(col="blue",alpha=0.5) -> p1
# 
# ggplot(mapping=aes(x=day,y=prot.rich),data=X.prot) +
#   geom_boxplot(col="blue",alpha=0.5)  -> p2
# 
# grid.arrange(p1,p2,ncol=2)
# 


####-----------------------------------------------------
#### Plot LRR in landscapes
####-----------------------------------------------------


lands = LETTERS[1:4]

#... Patch coordinates
#.....................

# Download connectivity matrices and patch coordinate
pTOL = cbind(rep(1:6,6),rep(6:1,each=6))
dimnames(pTOL) = list(1:36,c("x","y"))

setwd(to.script)
source("BG2_functions_20190301.R")
setwd(to.output)
connect.list = vector(mode = "list",length = 4); names(connect.list) = unique(data$Replicate) 
for(i in 1:4)connect.list[[i]] = read.table(paste("connectivity",i,".txt",sep=""))
patchsizes = read.csv2("patchsize.csv",row.names = 1)
minmaxPie = c(0.05,0.4)
minmaxPie = c(0.12,0.45)
legendPie=c(7.5,13,22.5,45)
spnames = SP

setwd(to.figs)

pdf("LRR_spatial.pdf")

for(i in 1:length(lands)){
  
  ## get the LRR data
  data_connected = subset(data,Treatment == "Connected" & Replicate == lands[i] & day == 29)
  data_isolated = subset(data,Treatment == "Isolated" & Replicate == lands[i] & day == 29)
  LRR_totdensit = log(data_connected$prot.ab)-log(data_isolated$prot.ab)
  LRR_totrich = log(data_connected$prot.rich)-log(data_isolated$prot.rich)
  if(i %in% c(1,2))LRR_bactdens = log(data_connected$bact.green)-log(data_isolated$bact.green)
  
  LRR_totdensit = ifelse(is.infinite(LRR_totdensit),NA,LRR_totdensit)
  LRR_totrich = ifelse(is.infinite(LRR_totrich),NA,LRR_totrich)
  
  rangeab = range(unlist(LRR_totdensit),na.rm=T)
  rangeab = range(unlist(LRR_totrich),na.rm=T)
  if(i %in% c(1,2))rangeab = range(unlist(LRR_bactdens),na.rm=T)
  
  
  colour.gradient=rainbow(100,start=0.1,end=0.9)
  colour.gradient= c("#f9f8e4","#f9f6c2","#fce99a","#f9d88f","#f9bb69","#f95341","#db0e18","#93052e","#7c0303","#560202","black")
  colour.gradient= c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30")
  

  
  
  legendCol = c(-4.5,-4,-3,-2,-1,0,1,2,3,4,4.5);

  #colour.pal = rainbow(p["ns"],start=0.1,end=0.8)
  
  plot.landscape.twoVariables(connect.list[[i]],LRR_totdensit,legendCol,colour.gradient,patchsizes[,lands[i]], legendPie,minmaxPie,tit=paste0("Landscape ",lands[i]," - Density"))
  plot.landscape.twoVariables(connect.list[[i]],LRR_totrich,legendCol,colour.gradient,patchsizes[,lands[i]], legendPie,minmaxPie,tit=paste0("Landscape ",lands[i]," - Richness"))
  
  if(i %in% c(1,2))plot.landscape.twoVariables(connect.list[[i]],LRR_bactdens,legendCol,colour.gradient,patchsizes[,lands[i]], legendPie,minmaxPie,tit=paste0("Landscape ",lands[i]," - Bacteria"))
  
  
}

dev.off()


col.morNEG = c("#fcfbfd","#efedf5","#dadaeb","#bcbddc","#9e9ac8","#807dba","#6a51a3","#54278f","#3f007d")
col.morPOS=c("#fff5eb","#fee6ce","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#7f2704")
colour.gradient= rev(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))

colour.gradient =c(col.morNEG[9:2],"white",col.morPOS[c(2,3,5,7,9)])

setwd(to.figs)

maxi = mini = 0

pdf("LRR_spatial_species_day29.pdf")

SP2 = c(SP,"bact.blue","bact.green")
head(data)
for(j in 1:length(SP2)){
  for(i in 1:length(lands)){
  
    ## get the LRR data
    data_connected = subset(data,Treatment == "Connected" & Replicate == lands[i] & day == 29)[,SP2[j]]
    if(j < 9){
      data_isolated = subset(data,Treatment == "Isolated" & Replicate == lands[i] & day == 29)[,SP2[j]]
    }else{
      data_isolatedA = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[1] & day == 29)[,SP2[j]])
      data_isolatedB = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[2] & day == 29)[,SP2[j]])
      data_isolated = apply(rbind(data_isolatedA,data_isolatedB),2,mean)
      
    }
    
    LRR_spdensit = unlist(log(data_connected))-unlist(log(data_isolated))
    
    LRR_spdensit = ifelse(is.infinite(LRR_spdensit),NA,LRR_spdensit)

    rangeab = range(unlist(LRR_spdensit),na.rm=T)
    maxi = ifelse(maxi<rangeab[2],rangeab[2],maxi)
    mini = ifelse(mini>rangeab[1],rangeab[1],mini)

    legendCol = c(-4.8, -4.2, -3.6, -3.0, -2.4, -1.8, -1.2, -0.6,  0.0, 0.6, 1.2, 1.8, 2.4, 3.0);
    
    
    plot.landscape.twoVariables(connect.list[[i]],LRR_spdensit,legendCol,colour.gradient,patchsizes[,lands[i]], legendPie,minmaxPie,tit=paste0("Landscape ",lands[i]," - Density of ",SP2[j]))

  }
}

dev.off()

mini
maxi
seq(-4.8,0,by=0.6)
seq(0,3,by=0.6)





colour.gradient =c(col.morNEG[9:2],"white",col.morPOS[c(3,5,7,9)])

setwd(to.figs)

maxi = mini = 0

pdf("LRR_spatial_species_days7to29.pdf")

SP2 = c(SP,"bact.blue","bact.green")
head(data)
for(j in 1:length(SP2)){
  for(i in 1:length(lands)){
    
    ## get the LRR data
    data_connected7 = unlist(subset(data,Treatment == "Connected" & Replicate == lands[i] & day == 7)[,SP2[j]])
    data_connected15 = unlist(subset(data,Treatment == "Connected" & Replicate == lands[i] & day == 15)[,SP2[j]])
    data_connected21 = unlist(subset(data,Treatment == "Connected" & Replicate == lands[i] & day == 21)[,SP2[j]])
    data_connected29 = unlist(subset(data,Treatment == "Connected" & Replicate == lands[i] & day == 29)[,SP2[j]])
    
    data_connected = apply(rbind(data_connected7,data_connected15,data_connected21,data_connected29),2,mean,na.rm=T)
    
    if(j < 9){
      data_isolated7 = subset(data,Treatment == "Isolated" & Replicate == lands[i] & day == 7)[,SP2[j]]
      data_isolated15 = subset(data,Treatment == "Isolated" & Replicate == lands[i] & day == 15)[,SP2[j]]
      data_isolated21 = subset(data,Treatment == "Isolated" & Replicate == lands[i] & day == 21)[,SP2[j]]
      data_isolated29 = subset(data,Treatment == "Isolated" & Replicate == lands[i] & day == 29)[,SP2[j]]
      
      data_isolated = apply(rbind(data_isolated7,data_isolated15,data_isolated21,data_isolated29),2,mean,na.rm=T)
    }else{
      data_isolatedA7 = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[1] & day == 7)[,SP2[j]])
      data_isolatedA15 = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[1] & day == 15)[,SP2[j]])
      data_isolatedA21 = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[1] & day == 21)[,SP2[j]])
      data_isolatedA29 = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[1] & day == 29)[,SP2[j]])
      data_isolatedB7 = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[2] & day == 7)[,SP2[j]])
      data_isolatedB15 = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[2] & day == 15)[,SP2[j]])
      data_isolatedB21 = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[2] & day == 21)[,SP2[j]])
      data_isolatedB29 = unlist(subset(data,Treatment == "Isolated" & Replicate == lands[2] & day == 29)[,SP2[j]])
      data_isolated = apply(rbind(data_isolatedA7,data_isolatedA15,data_isolatedA21,data_isolatedA29,
                                  data_isolatedB7,data_isolatedB15,data_isolatedB21,data_isolatedB29),2,mean,na.rm=T)
      
    }
    
    LRR_spdensit = log(data_connected)-log(data_isolated)
    
    LRR_spdensit = ifelse(is.infinite(LRR_spdensit),NA,LRR_spdensit)
    
    rangeab = range(unlist(LRR_spdensit),na.rm=T)
    maxi = ifelse(maxi<rangeab[2],rangeab[2],maxi)
    mini = ifelse(mini>rangeab[1],rangeab[1],mini)
    
    legendCol = c(-4.8, -4.2, -3.6, -3.0, -2.4, -1.8, -1.2, -0.6,  0.0, 0.6, 1.2, 1.8, 2.4);
    
    tit = paste0("Landscape ",lands[i]," - Density of ",SP2[j], " - mean effects over days 7 to 29")
    if(j %in% c(8,9) & i %in% c(3,4)) tit = paste0("Landscape ",lands[i]," - Density of ",SP2[j], " - only on Day 29")
    
    plot.landscape.twoVariables(connect.list[[i]],LRR_spdensit,legendCol,colour.gradient,patchsizes[,lands[i]], legendPie,minmaxPie,tit=tit)
    
  }
}

dev.off()

mini
maxi
seq(-4.8,0,by=0.6)

##########################################
##-- Barplots of species composition
##########################################

{
  # Summary for first species
  D_Res_Vol_protist <- data %>% 
    filter(day !=0) %>% 
    select(SP[1],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
    mutate(Treatment = as.factor(Treatment),
           Replicate = as.factor(Replicate),
           Size = as.factor(Size))  %>% 
    Rmisc::summarySE(.,measurevar = SP[1], groupvars = c("Treatment","Size"))
  

  Comm_Vol_protist = D_Res_Vol_protist[,c("Treatment","Size","N","Rot")]
  
  #Loop for each species
  
  for(i in 2:length(SP)) { 
    
    # Extract info for SP[i] and summarise data
    D_Res_Vol_protist <- data %>% 
      filter(day !=0) %>% 
      select(SP[i],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
      mutate(Treatment = as.factor(Treatment),
             Replicate = as.factor(Replicate),
             Size = as.factor(Size))  %>% 
      Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment","Size"))
    
    Comm_Vol_protist = cbind(Comm_Vol_protist,D_Res_Vol_protist[,SP[i]])
 
  }
  
  colnames(Comm_Vol_protist) = c("Treatment","Size","N",SP)
}
Comm_Vol_protist


{
  # Summary for first species
  D_Res_Vol_protistt0 <- data %>% 
    filter(day ==0) %>% 
    select(SP[1],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
    mutate(Treatment = as.factor(Treatment),
           Replicate = as.factor(Replicate),
           Size = as.factor(Size))  %>% 
    Rmisc::summarySE(.,measurevar = SP[1], groupvars = c("Treatment","Size"))
  
  
  Comm_Vol_protistt0 = D_Res_Vol_protistt0[,c("Treatment","Size","N","Rot")]
  
  #Loop for each species
  
  for(i in 2:length(SP)) { 
    
    # Extract info for SP[i] and summarise data
    D_Res_Vol_protistt0 <- data %>% 
      filter(day !=0) %>% 
      select(SP[i],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
      mutate(Treatment = as.factor(Treatment),
             Replicate = as.factor(Replicate),
             Size = as.factor(Size))  %>% 
      Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment","Size"))
    
    Comm_Vol_protistt0 = cbind(Comm_Vol_protistt0,D_Res_Vol_protistt0[,SP[i]])
    
  }
  
  colnames(Comm_Vol_protistt0) = c("Treatment","Size","N",SP)
}
Comm_Vol_protistt0

t0 = subset(data,day==0)
dim(t0)

Connt0 = apply(subset(Comm_Vol_protistt0,Treatment == "Connected")[,4:10],2,mean)
propConnt0 = Connt0 / sum(Connt0,na.rm=T)

Isolt0 = apply(subset(Comm_Vol_protistt0,Treatment == "Isolated")[,4:10],2,mean)
propIsolt0 = Isolt0 / sum(Isolt0,na.rm=T)


Allt0 = apply(Comm_Vol_protistt0[,4:10],2,mean)
propAllt0 = Allt0 / sum(Allt0,na.rm=T)

size7.5 = subset(Comm_Vol_protist,Size == 7.5)
size13 = subset(Comm_Vol_protist,Size == 13)
size22.5 = subset(Comm_Vol_protist,Size == 22.5)
size45 = subset(Comm_Vol_protist,Size == 45)

prop7.5i = ((size7.5[2,4:10] / apply(size7.5[2,4:10],1,sum,na.rm=T)))
prop13i = ((size13[2,4:10] / apply(size13[2,4:10],1,sum,na.rm=T)))
prop22.5i = ((size22.5[2,4:10] / apply(size22.5[2,4:10],1,sum,na.rm=T)))
prop45i= ((size45[2,4:10] / apply(size45[2,4:10],1,sum,na.rm=T)))

prop7.5c = ((size7.5[1,4:10] / apply(size7.5[1,4:10],1,sum,na.rm=T)))
prop13c = ((size13[1,4:10] / apply(size13[1,4:10],1,sum,na.rm=T)))
prop22.5c = ((size22.5[1,4:10] / apply(size22.5[1,4:10],1,sum,na.rm=T)))
prop45c= ((size45[1,4:10] / apply(size45[1,4:10],1,sum,na.rm=T)))

colvec2 = makeTransparent(colvec,alpha=0.85)



setwd(to.figs)

pdf("Barplot_bySize_Day7to29.pdf")

layout(matrix(1:9,ncol =9))


par(mar=c(1,0.1,3,1))
barplot(as.matrix(propAllt0),col=colvec,axes=F,main="t0")

par(mar=c(1,1,3,0.1))
barplot(t(as.matrix(prop7.5i)),names.arg=NA,density=15,axes=F,col="black")
barplot(t(as.matrix(prop7.5i)),col=colvec2,main="7.5",axes=F,add=T)
par(mar=c(1,0.1,3,1))
barplot(t(as.matrix(prop7.5c)),col=colvec,axes=F)


par(mar=c(1,1,3,0.1))
barplot(t(as.matrix(prop13i)),names.arg=NA,density=15,axes=F,col="black")
barplot(t(as.matrix(prop13i)),col=colvec2,main="13",axes=F,add=T)
par(mar=c(1,0.1,3,1))
barplot(t(as.matrix(prop13c)),col=colvec,axes=F)

par(mar=c(1,1,3,0.1))
barplot(t(as.matrix(prop22.5i)),names.arg=NA,density=15,axes=F,col="black")
barplot(t(as.matrix(prop22.5i)),col=colvec2,main="22.5",axes=F,add=T)
par(mar=c(1,0.1,3,1))
barplot(t(as.matrix(prop22.5c)),col=colvec,axes=F)

par(mar=c(1,1,3,0.1))
barplot(t(as.matrix(prop45i)),names.arg=NA,density=15,axes=F,col="black")
barplot(t(as.matrix(prop45i)),col=colvec2,main="45",axes=F,add=T)
par(mar=c(1,0.1,3,1))
barplot(t(as.matrix(prop45c)),col=colvec,axes=F)

dev.off()


######
##-- Day 29
######


colvec2 = makeTransparent(colvec,alpha=0.85)

{
  # Summary for first species
  D_Res_Vol_protist <- data %>% 
    filter(day ==29) %>% 
    select(SP[1],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
    mutate(Treatment = as.factor(Treatment),
           Replicate = as.factor(Replicate),
           Size = as.factor(Size))  %>% 
    Rmisc::summarySE(.,measurevar = SP[1], groupvars = c("Treatment","Size"))
  
  
  Comm_Vol_protist = D_Res_Vol_protist[,c("Treatment","Size","N","Rot")]
  
  #Loop for each species
  
  for(i in 2:length(SP)) { 
    
    # Extract info for SP[i] and summarise data
    D_Res_Vol_protist <- data %>% 
      filter(day ==29) %>% 
      select(SP[i],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
      mutate(Treatment = as.factor(Treatment),
             Replicate = as.factor(Replicate),
             Size = as.factor(Size))  %>% 
      Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment","Size"))
    
    Comm_Vol_protist = cbind(Comm_Vol_protist,D_Res_Vol_protist[,SP[i]])
    
  }
  
  colnames(Comm_Vol_protist) = c("Treatment","Size","N",SP)
}
Comm_Vol_protist

day29i = subset(Comm_Vol_protist,Treatment == "Isolated")
day29c = subset(Comm_Vol_protist,Treatment == "Connected")

day29propi = (day29i[,4:10] / apply(day29i[,4:10],1,sum,na.rm=T))
day29propc = (day29c[,4:10] / apply(day29c[,4:10],1,sum,na.rm=T))

setwd(to.figs)

pdf("Barplot_bySize_Day29.pdf")
{
layout(matrix(1:3,ncol =3),widths = c(0.95,2,2))

par(mar=c(3,3,3,3))
barplot(as.matrix(propAllt0),col=colvec,axes=F,main="t0")

par(mar=c(3,1,3,1))
barplot(t(as.matrix(day29propi)),density=15,axes=F,col="black",names.arg =sizes,cex.names = 1.5)
barplot(t(as.matrix(day29propi)),col=colvec2,axes=F,add=T,main="Isolated - day 29",names.arg =rep("",4),cex.names = 1.5)
#barplot(t(as.matrix(day29propi)),col=colvec,axes=F,main="Isolated")
barplot(t(as.matrix(day29propc)),col=colvec,axes=F,main="Connected - day 29",names.arg =sizes,cex.names = 1.5)


size7.5 = subset(Comm_Vol_protist,Size == 7.5)
size13 = subset(Comm_Vol_protist,Size == 13)
size22.5 = subset(Comm_Vol_protist,Size == 22.5)
size45 = subset(Comm_Vol_protist,Size == 45)

prop7.5i = ((size7.5[2,4:10] / apply(size7.5[2,4:10],1,sum,na.rm=T)))
prop13i = ((size13[2,4:10] / apply(size13[2,4:10],1,sum,na.rm=T)))
prop22.5i = ((size22.5[2,4:10] / apply(size22.5[2,4:10],1,sum,na.rm=T)))
prop45i= ((size45[2,4:10] / apply(size45[2,4:10],1,sum,na.rm=T)))

prop7.5c = ((size7.5[1,4:10] / apply(size7.5[1,4:10],1,sum,na.rm=T)))
prop13c = ((size13[1,4:10] / apply(size13[1,4:10],1,sum,na.rm=T)))
prop22.5c = ((size22.5[1,4:10] / apply(size22.5[1,4:10],1,sum,na.rm=T)))
prop45c= ((size45[1,4:10] / apply(size45[1,4:10],1,sum,na.rm=T)))



layout(matrix(1:8,ncol =8))

par(mar=c(1,1,3,0.1))
barplot(t(as.matrix(prop7.5i)),names.arg=NA,density=15,axes=F,col="black")
barplot(t(as.matrix(prop7.5i)),col=colvec2,main="7.5",axes=F,add=T)
par(mar=c(1,0.1,3,1))
barplot(t(as.matrix(prop7.5c)),col=colvec,axes=F)

par(mar=c(1,1,3,0.1))
barplot(t(as.matrix(prop13i)),names.arg=NA,density=15,axes=F,col="black")
barplot(t(as.matrix(prop13i)),col=colvec2,main="13",axes=F,add=T)
par(mar=c(1,0.1,3,1))
barplot(t(as.matrix(prop13c)),col=colvec,axes=F)

par(mar=c(1,1,3,0.1))
barplot(t(as.matrix(prop22.5i)),names.arg=NA,density=15,axes=F,col="black")
barplot(t(as.matrix(prop22.5i)),col=colvec2,main="22.5",axes=F,add=T)
par(mar=c(1,0.1,3,1))
barplot(t(as.matrix(prop22.5c)),col=colvec,axes=F)

par(mar=c(1,1,3,0.1))
barplot(t(as.matrix(prop45i)),names.arg=NA,density=15,axes=F,col="black")
barplot(t(as.matrix(prop45i)),col=colvec2,main="45",axes=F,add=T)
par(mar=c(1,0.1,3,1))
barplot(t(as.matrix(prop45c)),col=colvec,axes=F)
}
dev.off()


######
##-- Over time
######




{
  # Summary for first species
  D_Protists <- data %>% 
    Rmisc::summarySE(.,measurevar = SP[1], groupvars = c("Treatment","day"))
  
  
  Comm_protist = D_Protists[,c("Treatment","day","N","Rot")]
  
  #Loop for each species
  
  for(i in 2:length(SP)) { 
    
    # Extract info for SP[i] and summarise data
    D_Protists <- data %>% 
      Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment","day"))
    
    Comm_protist = cbind(Comm_protist,D_Protists[,SP[i]])
    
  }
  
  colnames(Comm_protist) = c("Treatment","day","N",SP)
}
Comm_protist

samp_days = c(7,15,21,29)
paste0("t",samp_days)
propCommc = (subset(Comm_protist,Treatment == "Connected" & day %in% samp_days)[,4:10] / apply(subset(Comm_protist,Treatment == "Connected" & day %in% samp_days)[,4:10],1,sum,na.rm=T))
propCommi = (subset(Comm_protist,Treatment == "Isolated" & day %in% samp_days)[,4:10] / apply(subset(Comm_protist,Treatment == "Isolated" & day %in% samp_days)[,4:10],1,sum,na.rm=T))

Allt0 = apply(Comm_Vol_protistt0[,4:10],2,mean)
propAllt0 = Allt0 / sum(Allt0,na.rm=T)

setwd(to.figs)

pdf("Barplot_overTime.pdf")

layout(matrix(1:3,ncol =3),widths = c(0.95,2,2))

par(mar=c(3,3,3,3))
barplot(as.matrix(propAllt0),col=colvec,axes=F,names.arg ="t0",cex.names = 1.5)

par(mar=c(3,1,3,1))
barplot(t(as.matrix(propCommi)),density=15,axes=F,col="black",names.arg =paste0("t",samp_days),cex.names = 1.5)
barplot(t(as.matrix(propCommi)),col=colvec2,axes=F,add=T,main="Isolated ",names.arg =rep("",4))
barplot(t(as.matrix(propCommc)),col=colvec,axes=F,main="Connected",names.arg =paste0("t",samp_days),cex.names = 1.5)

dev.off()
