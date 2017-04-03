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
datapath = "~/Documents/Research/Eawag/Projects/8.Blue-Green/BlueGreen/Data/"
to.script = "~/Documents/Research/Eawag/Projects/8.Blue-Green/BlueGreen/Scripts/"
graphpath = "~/Documents/Research/Eawag/Projects/8.Blue-Green/4.Results/"


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
#Load data in R environment
##################
# Prot.b =  read.delim("BG_protist_data_(20161124).txt") #Green were removed
Prot.b =  read.delim(paste0(datapath,"BG_protist_data_(20170307).txt")) #Protist
Cyto2 = read.delim(paste0(datapath,"FULL_Cyto_BlueGreen(20160816).txt")) # Bacteria

##################
#Fix bacteria data structure for further preliminary exploration and analyses
##################

#...Arrange structure cytometry data
Cyto2 = Cyto2[order(Cyto2$Date),]
Cyto2$day = c(rep(0,32),rep(7,524),rep(15,524),rep(21,524),rep(29,524)) #create a new variable "day"
Cyto2$count.ml = ((Cyto2$Count*1000)/Cyto2$Volume.uL)*1000 #convert from dens/50ul to dens/mL and then multiply by the cytometry dilution factor (1000)
sel.rep = c("A","B") #A,B for the whole duration - ABCD only on the last day
Cyto = Cyto2[Cyto2$Replicate %in% sel.rep,] #Extract only replicate A and B 
sel.treatment = c("Isolated","Connected")
Cyto = Cyto[Cyto$Treatment %in% sel.treatment,] #Drop "monocultures"


##################
#Fix Protist data structure for further preliminary exploration and analyses
##################
Prot.b$date = with(Prot.b,revalue(date, c("16-05-02"="20160502", "16-05-23"="20160523","16-05-31"="20160531","5/17/2016"="20160517","5/9/2016"="20160509")))
Prot.b$date = factor(Prot.b$date,levels=c("20160502","20160509","20160517","20160523","20160531"))
Prot.b = Prot.b[order(Prot.b$date),]
Prot.b$day = c(rep(0,30),rep(7,302),rep(15,302),rep(21,302),rep(29,302))


##################
#Aproximate initial protist total abundances per microcosm size from initial pooled cultures that were used to fill each microcosm at day 0
#(Lines 1:30 from dataset)
################### 
{  #1.Extract and Repeat pooled cultures from Day 0 for each patch size  
#Protist
Prot.day0 = Prot.b[rep(0:30,times=4), ]
row.names(Prot.day0) = 1:120
Prot.day0$Size[1:16] = 7.5
Prot.day0$Size[31:46] = 13
Prot.day0$Size[61:76] = 22.5
Prot.day0$Size[91:106] = 45

#Bacteria (for bacteria only pooled cultured from blue landscapes, connected and isolated, were measured [see original data] for logistic reasons, we assumed that our initial bacteria meausres are representative of the entire intial conditions for all treatments)
#Repeat for each patch size
Cyto.day0 = Cyto[rep(0:8,times=4),]
row.names(Cyto.day0) = 1:32
Cyto.day0$Size[1:8] = 7.5
Cyto.day0$Size[9:16] = 13
Cyto.day0$Size[17:24] = 22.5
Cyto.day0$Size[25:32] = 45
#Repeat for each landscape type (Blue and Green)
Cyto.day0.green = Cyto.day0[rep(0:32,times=2),]
row.names(Cyto.day0.green) = 1:64
Cyto.day0.green$Landscape[33:64] = "Green"
Cyto.day0.green$Size[Cyto.day0.green$Landscape=="Green"]=10

  #2.Remove Day 0 from the original dataset
#Protist
sel.date = c("20160509","20160517","20160523","20160531")
Prot.b = Prot.b[Prot.b$date %in% sel.date,] 
#Bacteria
Cyto = Cyto[Cyto$Date %in% sel.date,] 

  #3. Merge the new 'Day 0' repeated for each patch size with dataset
#Protist
Prot.b = rbind(Prot.day0,Prot.b)
str(Prot.b)
#'data.frame':	1328 obs. of  31 variables:

#Bacteria
Cyto = rbind(Cyto.day0.green,Cyto)
str(Cyto)
#'data.frame':	1216 obs. of  20 variables:
row.names(Cyto) = 1:1216

  #4. Remove monocultures from the Protist dataset
sel.treatment = c("Isolated","Connected")#,"monoculture"
Prot.b = Prot.b[Prot.b$Treatment %in% sel.treatment,] #Drop "monocultures"
str(Prot.b)
#'data.frame':	1216 obs. of  31 variables:
row.names(Prot.b) = 1:1216 #Same number of rows with Cyto but very different data structure - WARNING: do not merge! (Protist data contains only Blue treatments for each replicate ABCD, Cyto data contains only replicate AB for both Blue and Green treatments )

}

##################
#Calculate protist TOTAL abundance per microcosm (as opposed to density per volume)
##################
colnames(Prot.b)[23:30] = c("Rot.ul","Spi.ul","Ble.ul","Pca.ul","Col.ul","Chi.ul","Tet.ul","Other.ul")
Prot.b$Rot.all = Prot.b$Rot.ul*Prot.b$Size*1000
Prot.b$Spi.all = Prot.b$Spi.ul*Prot.b$Size*1000
Prot.b$Ble.all = Prot.b$Ble.ul*Prot.b$Size*1000
Prot.b$Pca.all = Prot.b$Pca.ul*Prot.b$Size*1000
Prot.b$Col.all = Prot.b$Col.ul*Prot.b$Size*1000
Prot.b$Chi.all = Prot.b$Chi.ul*Prot.b$Size*1000
Prot.b$Tet.all = Prot.b$Tet.ul*Prot.b$Size*1000
Prot.b$Other.all = Prot.b$Other.ul*Prot.b$Size*1000
Prot.b$bioarea.all = Prot.b$bioarea_per_ul*Prot.b$Size*1000

Prot.b$Tet.all[Prot.b$day==7 & Prot.b$Size==7.5 & Prot.b$Treatment=="Connected"]
# [1]  830.23255 1196.51163 1304.65116  165.69767    0.00000  231.97675  568.60465  436.04651  559.88372


##################
#Specify factors and Drop unused levels
##################
Prot.b$Label = factor(Prot.b$Label)
Prot.b$Size = factor(Prot.b$Size)
Prot.b = droplevels(Prot.b)

#...Fix structure and drop unused levels
Cyto$Date = factor(Cyto$Date)
Cyto$Label = factor(Cyto$Label)
Cyto$Size = factor(Cyto$Size)
Cyto = droplevels(Cyto)


#...Calculate some Protist community aggregate properties
Prot.b$Prot.tot.ab = rowSums(Prot.b[,32:39])
Prot.b$Prot.rich = specnumber(Prot.b[,32:39])
Prot.b$Prot.div = diversity(Prot.b[,32:39],"simpson")


#########################################################################
################# PRELIMINARY FIGURES 
#########################################################################

##################
# Protist species total abundance per treatment (connected/isolated) and per patch size
##################
{col.p = rainbow(8)

pdf(paste(graphpath,"Prot_Blue_Treatment_Size.pdf"),width=8,height=8)


for(i in 1:nlevels(Prot.b$Treatment)){
  
  for(j in 1:nlevels(Prot.b$Size)) { 
  
  max = max(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],32:38])

  with(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],],lineplot.CI(day,Rot.all,col=col.p[1],ylim=c(0,max),legend=T,xlab="Time",ylab="Protist total abundance"))
  par(new=T)
  with(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],],lineplot.CI(day,Chi.all,col=col.p[6],yaxt="n",ylim=c(0,max),legend=T,xlab="Time",ylab="Protist total abundance"))
  par(new=T)
  with(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],],lineplot.CI(day,Tet.all,col=col.p[7],ylim=c(0,max),yaxt="n",legend=T,xlab="Time",ylab="Protist total abundance"))
  par(new=T)
  with(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],],lineplot.CI(day,Col.all,col=col.p[5],ylim=c(0,max),yaxt="n",legend=T,xlab="Time",ylab="Protist total abundance"))
  par(new=T)
  with(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],],lineplot.CI(day,Pca.all,col=col.p[4],ylim=c(0,max),yaxt="n",legend=T,xlab="Time",ylab="Protist total abundance"))
  par(new=T)
  with(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],],lineplot.CI(day,Ble.all,col=col.p[3],ylim=c(0,max),yaxt="n",legend=T,xlab="Time",ylab="Protist total abundance"))
  par(new=T)
  with(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],],lineplot.CI(day,Spi.all,col=col.p[2],ylim=c(0,max),yaxt="n",legend=T,xlab="Time",ylab="Protist total abundance"))
  par(new=T)
  with(Prot.b[Prot.b$Treatment==levels(Prot.b$Treatment)[i] & Prot.b$Size==levels(Prot.b$Size)[j],],lineplot.CI(day,Other.all,col=col.p[2],ylim=c(0,max),yaxt="n",legend=T,xlab="Time",ylab="Protist total abundance"))
  
  
  title(main=paste(levels(Prot.b$Treatment)[i],levels(Prot.b$Size)[j],"mL"))
  col.labels=c("Rot","Spi","Ble","PCA","Col","Chi","Tet","Oth")
  testcol = col.p
  color.legend(4.5,10000,5,max,col.labels,testcol,cex=1,gradient="y")
  #dif between x1 and xr is the width of the rectangle
} 
  }

dev.off()}

##################
# Protist temporal community patterns (aggregate properties) across patch sizes
##################
{col.p = rainbow(n=4,start=0.1)

pdf(paste(graphpath,"Prot_Blue_DIV_AB_TIME.pdf"),width=8,height=8)

with(Prot.b,lineplot.CI(day,Prot.tot.ab,Size,col=col.p,legend=T,xlab="Time",ylab="Protist total abundance"))

with(Prot.b,lineplot.CI(day,Prot.tot.ab,Treatment,legend=T,xlab="Time",ylab="Protist total abundance"))

with(Prot.b,lineplot.CI(day,Prot.rich,Size,col=col.p,legend=T,xlab="Time",ylab="Protist species richness"))

with(Prot.b,lineplot.CI(day,Prot.rich,Treatment,col=col.p,legend=T,xlab="Time",ylab="Protist species richness"))

with(Prot.b,lineplot.CI(day,Prot.div,Size,col=col.p,legend=T,xlab="Time",ylab="Protist diversity (Simpson)"))

with(Prot.b,lineplot.CI(day,Prot.div,Treatment,col=col.p,legend=T,xlab="Time",ylab="Protist diversity (Simpson)"))

with(Prot.b,lineplot.CI(day,bioarea.all,Treatment,col=col.p,legend=T,xlab="Time",ylab="Total bioarea"))
dev.off()
}

##################
# Bacteria abundance patterns over time and by patch size*landscape*connected/isolated
##################
{

pdf(paste(graphpath,"BACT_AB_SIZE_LANDSCAPE_FLOW_TIME.pdf"),width=8,height=8)
  

for(i in 1:nlevels(Cyto$Treatment)){
  
  for(j in 1:nlevels(Cyto$Landscape)) {
    
    max = max(Cyto[Cyto$Treatment==levels(Cyto$Treatment)[i] & Cyto$Landscape==levels(Cyto$Landscape)[j],20])

    
    
    with(Cyto[Cyto$Treatment==levels(Cyto$Treatment)[i] & Cyto$Landscape==levels(Cyto$Landscape)[j],],lineplot.CI(day,count.ml,Size,col="black",ylim=c(0,max),legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
    
    
    title(main=paste(levels(Cyto$Treatment)[i],levels(Cyto$Landscape)[j]))

  }
}
  dev.off()
  }

##################
# Bacteria temporal patterns across treatment (connected vs. isolated) and landscape
##################
{
pdf(paste(graphpath,"BACT_AB_LANDSCAPE_FLOW_TIME.pdf"),width=8,height=8)
  
with(Cyto[Cyto$Landscape=="Blue",],lineplot.CI(day,count.ml,Treatment,legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
title(main="Blue")
with(Cyto[Cyto$Landscape=="Green",],lineplot.CI(day,count.ml,Treatment,legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
title(main="Green")

with(Cyto[Cyto$Treatment=="Connected",],lineplot.CI(day,count.ml,Landscape,legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
title(main="Connected")
with(Cyto[Cyto$Treatment=="Isolated",],lineplot.CI(day,count.ml,Landscape,legend=T,xlab="Time",ylab="Bacteria density (/mL)"))
title(main="Isolated")
dev.off()
}

##################
#Protist beta-diversity
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


##################
# Network figures
##################
#Build the dendritic landscapes from connectivity matrices and calculate network metrics of interest
library(igraph)
source(paste0(to.script,"Network_metric_bg.R"))


################
#...Bacteria
#The first date is based on pooled cultures therefore the graph they produce cannot be interpreted
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


##################
# Clean the environment before going further
#################
#...Clean the environment 
remove(lala);remove(lala2);remove(lala3);rm(sel.date);rm(sel.rep);rm(sel.treatment);rm(col.p)
remove(Green);remove(Green.m);remove(RepA);remove(RepA.m);remove(RepB);remove(RepB.m);remove(RepC);remove(RepC.m);remove(RepD);remove(RepD.m)
remove(coords1);remove(abundance);remove(fine);remove(graphCol);rm(i);rm(j);rm(k);rm(magn);rm(patch.size);rm(z)
rm(clrs)
detach("package:igraph") #igraph masks several functions in other packages

#Remove Day 0 for further analyses 
#Protist
sel.date = c("20160509","20160517","20160523","20160531")
Prot.b = Prot.b[Prot.b$date %in% sel.date,] 
Prot.b$date = factor(Prot.b$date)
#Bacteria
Cyto = Cyto[Cyto$Date %in% sel.date,]
Cyto$Date = factor(Cyto$Date)


##################
# Log Response Ratio
##################
{ 


#...Extract blue landscape only for bacteria
selb = c("Blue")
Cyto.b = Cyto[Cyto$Landscape %in% selb,]

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

#...Fix some variables
Prot.b$Size = as.numeric(as.character(Prot.b$Size))
Cyto.b$Size = as.numeric(as.character(Cyto.b$Size))

#...Separate connected and isolated landscapes
#.......Connected
#..........Protist
Prot.b.connected = Prot.b[which(Prot.b$Treatment=="Connected"),]
Prot.b.connected$indiv_per_ul[497] = 1
Prot.b.connected$Prot.rich[497] = 1
Prot.b.connected$indiv_per_ul[467] = 1
Prot.b.connected$Prot.rich[467] = 1
#..........Bacteria
Cyto.b.connected = Cyto.b[which(Cyto.b$Treatment=="Connected"),]
#.......Isolated
#..........Protist
Prot.b.isolated = Prot.b[which(Prot.b$Treatment=="Isolated"),]
Prot.b.isolated$indiv_per_ul[467] = 1
Prot.b.isolated$Prot.rich[467] = 1
Prot.b.isolated$indiv_per_ul[497] = 1
Prot.b.isolated$Prot.rich[497] = 1
#..........Bacteria
Cyto.b.isolated = Cyto.b[which(Cyto.b$Treatment=="Isolated"),]

#...Calculate log Response Ratio
var.prot.rich = log(Prot.b.connected$Prot.rich/Prot.b.isolated$Prot.rich)
var.prot.ab = log(Prot.b.connected$indiv_per_ul/Prot.b.isolated$indiv_per_ul)
var.bact.ab = log(Cyto.b.connected$count.ml/Cyto.b.isolated$count.ml)

#Calculate 95% CI for each factor levels of interest
err.rich.size = c(qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==7.5)])/sqrt(340),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==13)])/sqrt(136),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==22.5)])/sqrt(56),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==45)])/sqrt(44))
err.rich.day = c(qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$day==7)])/sqrt(144),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$day==15)])/sqrt(144),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$day==21)])/sqrt(144),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$day==29)])/sqrt(144))

err.ab.size = c(qnorm(0.975)*sd(var.prot.ab[which(Prot.b.isolated$Size==7.5)])/sqrt(340),qnorm(0.975)*sd(var.prot.ab[which(Prot.b.isolated$Size==13)])/sqrt(136),qnorm(0.975)*sd(var.prot.ab[which(Prot.b.isolated$Size==22.5)])/sqrt(56),qnorm(0.975)*sd(var.prot.ab[which(Prot.b.isolated$Size==45)])/sqrt(44))
err.ab.day = c(qnorm(0.975)*sd(var.prot.ab[which(Prot.b.isolated$day==7)])/sqrt(144),qnorm(0.975)*sd(var.prot.ab[which(Prot.b.isolated$day==15)])/sqrt(144),qnorm(0.975)*sd(var.prot.ab[which(Prot.b.isolated$day==21)])/sqrt(144),qnorm(0.975)*sd(var.prot.ab[which(Prot.b.isolated$day==29)])/sqrt(144))

err.bact.size = c(qnorm(0.975)*sd(var.bact.ab[which(Cyto.b.isolated$Size==7.5)])/sqrt(340),qnorm(0.975)*sd(var.bact.ab[which(Cyto.b.isolated$Size==13)])/sqrt(136),qnorm(0.975)*sd(var.bact.ab[which(Cyto.b.isolated$Size==22.5)])/sqrt(56),qnorm(0.975)*sd(var.bact.ab[which(Cyto.b.isolated$Size==45)])/sqrt(44))
err.bact.day = c(qnorm(0.975)*sd(var.bact.ab[which(Cyto.b.isolated$day==7)])/sqrt(144),qnorm(0.975)*sd(var.bact.ab[which(Cyto.b.isolated$day==15)])/sqrt(144),qnorm(0.975)*sd(var.bact.ab[which(Cyto.b.isolated$day==21)])/sqrt(144),qnorm(0.975)*sd(var.bact.ab[which(Cyto.b.isolated$day==29)])/sqrt(144))

#Plot Figures
 
library(Hmisc)

pdf(paste(graphpath,"LRR.pdf"),width=8,height=8)

#...Patch size
#Species richness
plot(Prot.b.isolated$Size,var.prot.rich,type="n",ylab="Effect size",ylim=c(-0.5,0.5),xaxt="n", main="Species richness")
axis(side = 1, at =c(7.5,13,22.5,45))
abline(h=0,lwd=1,col="gray")
x = as.numeric(levels(as.factor(Prot.b.isolated$Size)))
y = tapply(var.prot.rich,Prot.b.isolated$Size,mean)
ymax = y + err.rich.size
yminus = y - err.rich.size
errbar(x,y,ymax,yminus,add=T,errbar.col="blue",col="blue",lty=1,pch=16)

#Density per ul 
plot(Prot.b.isolated$Size,var.prot.ab,type="n",ylab="Effect size",ylim=c(-0.5,0.5),xaxt="n",main="Protist density")
axis(side = 1, at =c(7.5,13,22.5,45))
abline(h=0,lwd=1,col="gray")
x = as.numeric(levels(as.factor(Prot.b.isolated$Size)))
y = tapply(var.prot.ab,Prot.b.isolated$Size,mean)
ymax = y + err.ab.size
yminus = y - err.ab.size
errbar(x,y,ymax,yminus,add=T,errbar.col="blue",col="blue",lty=1,pch=16)

#Bacteria density 
plot(Cyto.b.isolated$Size,var.bact.ab,type="n",ylab="Effect size",ylim=c(-0.5,0.5),xaxt="n",main="Bacteria density")
axis(side = 1, at =c(7.5,13,22.5,45))
abline(h=0,lwd=1,col="gray")
x = as.numeric(levels(as.factor(Cyto.b.isolated$Size)))
y = tapply(var.bact.ab,Cyto.b.isolated$Size,mean)
ymax = y + err.bact.size
yminus = y - err.bact.size
errbar(x,y,ymax,yminus,add=T,errbar.col="blue",col="blue",lty=1,pch=16)

#...Time
#Species richness
plot(Prot.b.isolated$day,var.prot.rich,type="n",ylab="Effect size",ylim=c(-0.5,0.5),xaxt="n", main="Species richness")
axis(side = 1, at =c(7,15,21,29))
abline(h=0,lwd=1,col="gray")
x = as.numeric(levels(as.factor(Prot.b.isolated$day)))
y = tapply(var.prot.rich,Prot.b.isolated$day,mean)
ymax = y + err.rich.day
yminus = y - err.rich.day
errbar(x,y,ymax,yminus,add=T,errbar.col="blue",col="blue",lty=1,pch=16)

#Density per ul 
plot(Prot.b.isolated$day,var.prot.ab,type="n",ylab="Effect size",ylim=c(-0.5,0.5),xaxt="n",main="Protist density")
axis(side = 1, at =c(7,15,21,29))
abline(h=0,lwd=1,col="gray")
x = as.numeric(levels(as.factor(Prot.b.isolated$day)))
y = tapply(var.prot.ab,Prot.b.isolated$day,mean)
ymax = y + err.ab.day
yminus = y - err.ab.day
errbar(x,y,ymax,yminus,add=T,errbar.col="blue",col="blue",lty=1,pch=16)

#Bacteria density 
plot(Cyto.b.isolated$day,var.bact.ab,type="n",ylab="Effect size",ylim=c(-0.5,0.5),xaxt="n",main="Bacteria density")
axis(side = 1, at =c(7,15,21,29))
abline(h=0,lwd=1,col="gray")
x = as.numeric(levels(as.factor(Cyto.b.isolated$day)))
y = tapply(var.bact.ab,Cyto.b.isolated$day,mean)
ymax = y + err.bact.day
yminus = y - err.bact.day
errbar(x,y,ymax,yminus,add=T,errbar.col="blue",col="blue",lty=1,pch=16)


#Patch size at day 7
err.rich.size.day7 = c(qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==7.5 & Prot.b.isolated$day==7)])/sqrt(85),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==13 & Prot.b.isolated$day==7)])/sqrt(34),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==22.5 & Prot.b.isolated$day==7)])/sqrt(14),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==45 & Prot.b.isolated$day==7)])/sqrt(11))
plot(Prot.b.isolated$Size[which(Prot.b.isolated$day==7)],var.prot.rich[which(Prot.b.isolated$day==7)],type="n",ylab="Effect size",ylim=c(-0.5,0.5),xaxt="n", main="Species richness at day 7")
axis(side = 1, at =c(7.5,13,22.5,45))
abline(h=0,lwd=1,col="gray")
x = as.numeric(levels(as.factor(Prot.b.isolated$Size)))
y = tapply(var.prot.rich[which(Prot.b.isolated$day==7)],Prot.b.isolated$Size[which(Prot.b.isolated$day==7)],mean)
ymax = y + err.rich.size
yminus = y - err.rich.size
errbar(x,y,ymax,yminus,add=T,errbar.col="blue",col="blue",lty=1,pch=16)

#Patch size at day 29
err.rich.size.day29 = c(qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==7.5 & Prot.b.isolated$day==29)])/sqrt(85),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==13 & Prot.b.isolated$day==29)])/sqrt(34),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==22.5 & Prot.b.isolated$day==29)])/sqrt(14),qnorm(0.975)*sd(var.prot.rich[which(Prot.b.isolated$Size==45 & Prot.b.isolated$day==29)])/sqrt(11))
plot(Prot.b.isolated$Size[which(Prot.b.isolated$day==29)],var.prot.rich[which(Prot.b.isolated$day==29)],type="n",ylab="Effect size",ylim=c(-0.5,0.5),xaxt="n", main="Species richness at day 29")
axis(side = 1, at =c(7.5,13,22.5,45))
abline(h=0,lwd=1,col="gray")
x = as.numeric(levels(as.factor(Prot.b.isolated$Size)))
y = tapply(var.prot.rich[which(Prot.b.isolated$day==29)],Prot.b.isolated$Size[which(Prot.b.isolated$day==29)],mean)
ymax = y + err.rich.size
yminus = y - err.rich.size
errbar(x,y,ymax,yminus,add=T,errbar.col="blue",col="blue",lty=1,pch=16)


dev.off()
detach("package:Hmisc") 

}

#########################################################################
################# MAIN ANALYSES
#########################################################################


#..........................................#
#           Blue-green interactions        #
#..........................................#

#..Identify probability distribution
{ 
with(Cyto,plot(density(count.ml)))
with(Cyto,plot(density(log(count.ml))))
summary(Cyto$count.ml)
with(Cyto,qqp(count.ml,"norm"))
with(Cyto,qqp(count.ml,"lnorm"))
# qqp requires estimates of the parameters of the negative binomial, Poisson
# and gamma distributions. You can generate estimates using the fitdistr
# function. Save the output and extract the estimates of each parameter as I
# have shown below.
nbinom <- with(Cyto,fitdistr(count.ml, "negative binomial"))
with(Cyto,qqp(count.ml, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]]))
poisson <- with(Cyto,fitdistr(count.ml, "Poisson"))
with(Cyto,qqp(count.ml, "pois",poisson$estimate))
gamma = with(Cyto,fitdistr(count.ml, "gamma"))
with(Cyto,qqp(Count, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]]))
}

#Run model
library(nlme)
Mod = lme(log(Count) ~ Treatment*Landscape*as.numeric(Date), ~ Date|Replicate,data=Cyto,method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod)$tTable

#..........................................#
#           Blue Landscape only            #
#..........................................#


##################
# Extract blue landscapes only from Bacteria data
##################

#############
#...Bacteria
selb = c("Blue")
Cyto.b = Cyto[Cyto$Landscape %in% selb,]
#.Fix structure and drop unused levels
Cyto.b$Date = factor(Cyto.b$Date)
Cyto.b$Label = factor(Cyto.b$Label)
Cyto.b$ID_TUBE = factor(Cyto.b$ID_TUBE)
Cyto.b$Treatment = factor(Cyto.b$Treatment)
Cyto.b$Landscape = factor(Cyto.b$Landscape)
Cyto.b$Size = factor(Cyto.b$Size)
Cyto.b$Replicate = factor(Cyto.b$Replicate)

##################
# Merge network metrics to data
##################

#############
#...Bacteria
Cyto.b$Degree = rep(c(RepA.degree,RepB.degree),8)
Cyto.b$Dist = rep(c(RepA.dist,RepB.dist),8)
Cyto.b$Diam = rep(c(RepA.diam,RepB.diam),8)
Cyto.b$Drainage = rep(c(RepA.drainage,RepB.drainage),8)

#############
#...Protist
Prot.b$Degree = rep(c(RepA.degree,RepB.degree,RepC.degree,RepD.degree),8)
Prot.b$Dist = rep(c(RepA.dist,RepB.dist,RepC.dist,RepD.dist),8)
Prot.b$Diam = rep(c(RepA.diam,RepB.diam,RepC.diam,RepD.diam),8)
Prot.b$Drainage = rep(c(RepA.drainage,RepB.drainage,RepC.drainage,RepD.drainage),8)

##################
# Plot effects from green to blue
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

#...Extract green landscapes
selg = c("Green")
Cyto.g = Cyto[Cyto$Landscape %in% selg,]

#Extract network metrics of the connected blue landscapes

patch.size.blue = as.numeric(as.character(Cyto$Size[which(Cyto$Date !="20160502" & Cyto$Landscape=="Blue")]))

#...Fix structure and drop unused levels
Cyto.g$Date = as.numeric(Cyto.g$Date)
Cyto.g$Label = factor(Cyto.g$Label)
Cyto.g$ID_TUBE = factor(Cyto.g$ID_TUBE)
Cyto.g$Treatment = factor(Cyto.g$Treatment)
Cyto.g$Landscape = factor(Cyto.g$Landscape)
Cyto.g$Size = factor(Cyto.g$Size)
Cyto.g$Replicate = factor(Cyto.g$Replicate)
str(Cyto.g)
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

#Look at all replicates on the last sampling day
patch.size = Cyto2$Size[which(Cyto2$Date=="20160531" & Cyto2$Landscape=="Blue" & Cyto2$Treatment=="Connected")]
with(Cyto2[Cyto2$Date=="20160531" & Cyto2$Treatment=="Connected" & Cyto2$Landscape=="Green",],lineplot.CI(patch.size,count.ml,Treatment))




