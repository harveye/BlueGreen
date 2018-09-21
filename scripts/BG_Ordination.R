#..........................................#
#...........Blue-Green project.............#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#... Contact: eric.harvey@eawag.ch                                                                                           #                                                                       #
#..........................................................................................................................................#

##################
#Clear any variables from the R environment 
##################
rm(list=ls())

##################
#Load data in R environment
##################
library(tidyverse)
library(vegan)

##################
#Directory paths 
##################
to.data <- "./data/"
to.script <- "./scripts/"
to.output <- "./output/"
to.figs <- "./figs/"
to.R <- "./output/"

##################
#Load data in R environment
##################
data = readRDS(paste0(to.output,"cl.data.RDS"))

head(data)

#########################################################################
################# Analysis
#########################################################################

##################
# Preliminary figures
##################

ggplot(data=data,mapping=aes(x=Treatment,y=prot.rich)) +
  geom_boxplot(mapping=aes(color=Treatment)) + 
  facet_wrap(~ day,nrow=2)

ggplot(data=data,mapping=aes(x=bact.blue,y=prot.ab)) +
  geom_point(mapping=aes(color=Treatment)) + 
  facet_wrap(~ day,nrow=2)

ggplot(data=data,mapping=aes(x=bact.blue,y=prot.rich)) +
  geom_point(mapping=aes(color=Treatment)) + 
  facet_wrap(~ day,nrow=2)

ggplot(data=data,mapping=aes(x=bact.green,y=prot.ab)) +
  geom_point(mapping=aes(color=Treatment)) + 
  facet_wrap(~ day,nrow=2)

ggplot(data=data,mapping=aes(x=as.factor(day),y=prot.ab)) +
  geom_boxplot(mapping=aes(color=Treatment))  
  


##################
# Ordination (temporal trends in community structure)
##################

#Place vars of interest into vectors for each component of ordination
SP <- c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other") #comm matrix
ENV <- c("centrality","dist.outlet","Treatment") #Constrains
CONT <- c("bact.green") #Effect to partial out

#Filter data as wanted (for bacteria need to use only rep A and B)
X <- data %>% 
  filter(day!=0 & Replicate %in% c("A","B","C","D")) %>% 
  mutate(Treatment = as.numeric(as.factor(Treatment)))

#Create matrices for analyses 
C <- decostand(X[,SP],"hell")
E <- as.matrix(X[,ENV])
Z <- as.matrix(X[,CONT])


########## REDUNDANCY ANALYSIS (RDA)

rda.mod <- rda(C ~ ., as.data.frame(E))
#rda.mod <- rda(C ~ X$Treatment + X$centrality +  X$dist.outlet + Condition(Z))
rda.mod
anova(rda.mod)
anova(rda.mod,by="terms",permu=200)

#define color vectors for sampling days
colvec <- c("blue","green","red","black")
ordiplot(rda.mod,display="site",type="n")
text(rda.mod, display="cn", col="blue",cex=0.5)
points(rda.mod, display="site", cex=0.5,pch=c(16, 17)[X$Treatment],col=colvec[as.factor(X$day)])
#points(rda.mod, display="site", cex=0.8,pch=c(16, 2, 10)[as.factor(X$day)],col=c("red","blue")[as.factor(X$Treatment)])
text(rda.mod, display="sp", col="blue",cex=1)
legend("bottomleft",legend=c("day7","day15","day21","day29"),col=c("blue","green","red","black"),pch=16,bty="n",cex=0.5)
legend("bottomright",legend=c("connected","isolated"),pch=c(1,2),bty="n",cex=0.5)
text(5,2,labels = "Hellinger transformed",cex=0.5)
# ordihull(rda.mod,groups=X$day,show.groups=c("7","15","21","29"),col=colvec,label=F,lwd=1,lty=1)
# ordihull(rda.mod,groups=X$Treatment,show.groups=c("1","2"),col=c("blue","green"),label=F,lwd=1,lty=1)


############# Tb-PCA
C <- decostand(X[,SP],"hell")

pca.mod <- rda(C)
pca.mod

colvec <- c("blue","green","red","black")
ordiplot(pca.mod,display="site",type="n")
points(pca.mod, display="site", cex=0.8,pch=c(16, 2, 10)[as.factor(X$day)],col=c("red","blue")[as.factor(X$Treatment)])

points(pca.mod, display="site", cex=0.5,pch=c(16, 17)[X$Treatment],col=colvec[as.factor(X$day)])


############# NMDS

#Generate distance-based matrix 
# C <- decostand(X[,SP],"hell")
# C[which(rowSums(C)==0)[1],1:8] = 0.000001#does not allow rowSums of 0
# C[which(rowSums(C)==0)[1],1:8] = 0.000002#does not allow rowSums of 0
# C <- vegdist(C,"bray")
# nmds.mod <- metaMDS(C,k=3)#engine="isoMDS"



