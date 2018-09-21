#..........................................#
#...........Blue-Green project.............#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#... Date latest modifications: 2017-01-23                                                                                               #                                                                       #
#..........................................................................................................................................#

#########################################################################
################# This script generate all the network variables for the analyses of the Blue-Green data
#########################################################################

##################
#Clear any variables from the R environment 
##################
rm(list=ls())

##################
#Load data in R environment
##################
library(igraph)

##################
#Directory paths 
##################
to.data <- "./data/"
to.script <- "./scripts/"
to.output <- "./output/"
to.figs <- "./figs/"
to.R <- "./output/"

Replicate = as.factor(c(rep("A",36),rep("B",36),rep("C",36),rep("D",36)))


##################
#Extract each connectivity matrix 
# from Am.Nat. paper (Carrara et al. 2014) 
##################

#The row names from the original text file needs to be changed 
RepA = read.delim(paste0(to.data,"rep1.txt"),header=F)
rownames(RepA) = c(1,7,13,19,25,31,2,8,14,20,26,32,3,9,15,21,27,33,4,10,16,22,28,34,5,11,17,23,29,35,6,12,18,24,30,36) #Rename rows from the adjacency matrix to fit with row numbers in the experimental landscapes
colnames(RepA) = rownames(RepA)
RepA = RepA[order(as.numeric(as.character(rownames(RepA)))),order(as.numeric(as.character(colnames(RepA))))]

RepB = read.delim(paste0(to.data,"rep2.txt"),header=F)
rownames(RepB) = c(1,7,13,19,25,31,2,8,14,20,26,32,3,9,15,21,27,33,4,10,16,22,28,34,5,11,17,23,29,35,6,12,18,24,30,36)
colnames(RepB) = rownames(RepB)
RepB = RepB[order(as.numeric(as.character(rownames(RepB)))),order(as.numeric(as.character(colnames(RepB))))]

RepC = read.delim(paste0(to.data,"rep3.txt"),header=F)
rownames(RepC) = c(1,7,13,19,25,31,2,8,14,20,26,32,3,9,15,21,27,33,4,10,16,22,28,34,5,11,17,23,29,35,6,12,18,24,30,36)
colnames(RepC) = rownames(RepC)
RepC = RepC[order(as.numeric(as.character(rownames(RepC)))),order(as.numeric(as.character(colnames(RepC))))]

RepD = read.delim(paste0(to.data,"rep4.txt"),header=F)
rownames(RepD) = c(1,7,13,19,25,31,2,8,14,20,26,32,3,9,15,21,27,33,4,10,16,22,28,34,5,11,17,23,29,35,6,12,18,24,30,36)
colnames(RepD) = rownames(RepD)
RepD = RepD[order(as.numeric(as.character(rownames(RepD)))),order(as.numeric(as.character(colnames(RepD))))]

Green = read.delim(paste0(to.data,"Green.txt"),header=F)
colnames(Green) = rownames(Green)

##################
#Convert each data frame into an adjacency matrix
##################

RepA.m = as.matrix(RepA)
RepB.m = as.matrix(RepB)
RepC.m = as.matrix(RepC)
RepD.m = as.matrix(RepD)
Green.m = as.matrix(Green)

##################
#Convert each matrix into an igraph object 
##################
RepA.g = graph_from_adjacency_matrix(RepA.m, mode = "undirected", weighted = NULL, diag = F,
                                     add.colnames = NULL, add.rownames = NA)

RepB.g = graph_from_adjacency_matrix(RepB.m, mode = "undirected", weighted = NULL, diag = F,
                                     add.colnames = NULL, add.rownames = NA)

RepC.g = graph_from_adjacency_matrix(RepC.m, mode = "undirected", weighted = NULL, diag = F,
                                     add.colnames = NULL, add.rownames = NA)

RepD.g = graph_from_adjacency_matrix(RepD.m, mode = "undirected", weighted = NULL, diag = F,
                                     add.colnames = NULL, add.rownames = NA)

Green.g = graph_from_adjacency_matrix(Green.m, mode = "undirected", weighted = NULL, diag = F,
                                      add.colnames = NULL, add.rownames = NA)

# #. Visualize each replicate 
#Rep1
png(paste(to.figs,"RepA.jpg"))
coords1 <- layout_on_grid(RepA.g,width=6,height=6) #set the rectangular layout
coords1[,2] <- max(coords1[,2])-coords1[,2] #Inverse the Y position
plot.igraph(RepA.g,layout=coords1) #plot
dev.off()

#Rep2
png(paste(to.figs,"RepB.png"))
coords2 <- layout_on_grid(RepB.g,width=6,height=6) #set the rectangular layout
coords2[,2] <- max(coords2[,2])-coords2[,2] #Inverse the Y position
plot.igraph(RepB.g,layout=coords2) #plot
dev.off()

#Rep3
png(paste(to.figs,"RepC.png"))
coords3 <- layout_on_grid(RepC.g,width=6,height=6) #set the rectangular layout
coords3[,2] <- max(coords3[,2])-coords3[,2] #Inverse the Y position
plot.igraph(RepC.g,layout=coords3) #plot
dev.off()

#Rep4
png(paste(to.figs,"RepD.png"))
coords4 <- layout_on_grid(RepD.g,width=6,height=6) #set the rectangular layout
coords4[,2] <- max(coords4[,2])-coords4[,2] #Inverse the Y position
plot.igraph(RepD.g,layout=coords4) #plot
dev.off()

#Green
png(paste(to.figs,"Green_landscape.png"))
coords6 <- layout_on_grid(Green.g,width=6,height=6) #set the rectangular layout
coords6[,2] <- max(coords6[,2])-coords6[,2] #Inverse the Y position
plot.igraph(Green.g,layout=coords6) #plot
dev.off()

##################
#Calculate Network Metrics of interest
##################

#..Degree of each local community

RepA.degree = degree(RepA.g)
RepB.degree = degree(RepB.g)
RepC.degree = degree(RepC.g)
RepD.degree = degree(RepD.g)
Green.degree = degree(Green.g)

degree = c(RepA.degree,RepB.degree,RepC.degree,RepD.degree)

#...Distance to the outlet

RepA.dist = distances(RepA.g)[36,]
RepB.dist = distances(RepB.g)[36,]
RepC.dist = distances(RepC.g)[36,]
RepD.dist = distances(RepD.g)[36,]

dist.outlet = c(RepA.dist,RepB.dist,RepC.dist,RepD.dist)

#...Closeness centrality 

RepA.diam = closeness(RepA.g)
RepB.diam = closeness(RepB.g)
RepC.diam = closeness(RepC.g)
RepD.diam = closeness(RepD.g)

centrality = c(RepA.diam,RepB.diam,RepC.diam,RepD.diam)

#...Drainage volume
RepA.drainage = 
  c(0,0,0,0,0,0,
    0,22.5,0,0,22.5,0,
    0,65.5,0,0,15,43,
    0,0,95.5,0,0,84,
    0,15,43,189,0,114,
    0,0,0,0,241.5,430.5)

RepB.drainage = 
  c(0,0,0,0,0,0,
    0,15,0,22.5,0,7.5,
    0,7.5,28,50.5,0,20.5,
    0,22.5,50.5,121.5,104.5,0,
    0,0,0,33.5,198,134.5,
    0,7.5,20.5,0,0,398)

RepC.drainage = 
  c(0,0,0,0,0,0,
    0,22.5,0,15,0,7.5,
    0,0,58,0,28,28,
    0,7.5,0,71,89.5,0,
    0,0,63.5,185,0,119.5,
    0,15,0,0,237.5,387)

RepD.drainage = 
  c(0,0,0,0,0,0,
    0,30,0,0,22.5,0,
    0,0,50.5,0,50.5,0,
    0,22.5,78.5,0,78.5,0,
    0,0,35.5,209.5,366.5,0,
    0,7.5,28,89.5,0,419)

drainage_vol = c(RepA.drainage,RepB.drainage,RepC.drainage,RepD.drainage)


#Merge data in one dataframe

net_metric = data.frame(ID=1:36,Replicate,degree,centrality,dist.outlet,drainage_vol)
str(net_metric)
head(net_metric)

#Save data
saveRDS(net_metric,paste0(to.output,"net_metric.RDS"))
