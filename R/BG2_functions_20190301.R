#... Blue-Green Project                                              
#... co-authors: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt
#... Author of the script: Isabelle Gounand
#... date: 06-07-2017
#... Script description: functions for the analysis of the C output: plot the landscapes (of 6x6 patches)

#install.packages("plotrix")
library(plotrix)
#install.packages("mapplots")
library(mapplots)


#... Plot the vertices (connectivity between patches)
#.....................
  # provide the connectivity matrix of size patch number ^2

plot.vertices = function(connectivity.mat)
{
  co = connectivity.mat
  co[upper.tri(co)]=0
  for(i in 1:nrow(co))for(j in 1:ncol(co))if(co[i,j]==1)segments(pTOL[i,"x"],pTOL[i,"y"],pTOL[j,"x"],pTOL[j,"y"])
}



#... Scale an array
#..................
  # x0 the array to scale
  # minmax the boundaries (array of 2 number) between which the array has to be scaled (uniform transformation)

scale.func=function(x0,minmax){
  ran.x=max(x0,na.rm=T)-min(x0,na.rm=T)
  ran.s=max(minmax,na.rm=T)-min(minmax,na.rm=T)
  x1=x0 * ran.s / ran.x
  return(x1+(min(minmax,na.rm=T)-min(x1,na.rm=T)))
}


#... Plot the patches with patch size and colour gradients
#.........................................................
  # value1: array of patch values to scale to the colour gradient
  # minmaxV1: min and max values to scale to the colour extrema (array of 2 numbers)
  # value2: array of patch values to scale to the pie size
  # minmaxV2: min and max values to scale to the pie size extrema (array of 2 numbers)
  # colour.grad: colour gradient
  # minmaxRay: min and max of pie radius (array of 2 numbers)

plot.patch.twoVariables = function(values1,minmaxV1,values2,minmaxV2,colour.grad,minmaxRay)
{
  values.scaled.col = scale.func(x0 = c(minmaxV1,values1),minmax = c(1:length(colour.grad)))[-c(1,2)]
  values.scaled.ray = scale.func(x0 = sqrt(sqrt(c(minmaxV2,values2))),minmax = minmaxRay)[-c(1,2)]
  for(i in 1:length(values1)){
    diffcol = abs(c(1:length(colour.grad))-values.scaled.col[i])
    colo = colour.grad[which(diffcol == min(diffcol))]
    draw.circle(x=pTOL[i,"x"],y=pTOL[i,"y"],radius=values.scaled.ray[i],nv=100,border=NULL,col=colo,lty=1,lwd=1)
  }
}



#... Plot the landscape with colour gradient (species richness)
#..............................................................
  # L.connectivity: matrix of landscape connectivity
  # values1: array of patch values to scale to the colour gradient
  # legendCol: array of the values to put in colour legend (including min and max)
  # colour.gradient: colour gradient
  # values2 : array of patch values to scale to the pie radius scale
  # legendPie: array of 6 values to put in pie size legend (including min and max)
  # minmaxPie: min and max values for the size scale of pies
  # tit: title of the plot
  # Test values
  #L.connectivity = connect.list[[1]]
  #values1 = runif(36,0.2,5.9); values2 = runif(36,0,89); values1 = rich ; values2 = ab
  

plot.landscape.twoVariables = function(L.connectivity,values1,legendCol,colour.gradient,values2, legendPie,minmaxPie,tit)
{
  par(mar=c(0,0,0,0))
  plot(NA,xlim=c(0,9),ylim=c(-1,8),xlab ="",ylab="", xaxt="n",yaxt="n",frame.plot = F)
  text(3.5,7.5,labels = tit)
  plot.vertices(L.connectivity)
  plot.patch.twoVariables(values1,range(legendCol),values2,range(legendPie),colour.gradient,minmaxPie)

  gradient.rect(xleft=7.1,ybottom = 1-2.5/length(colour.gradient), xright = 7.5, ytop=6+2.5/length(colour.gradient),col = colour.gradient, gradient = "y")
  yticksCol = scale.func(x0 = legendCol,minmax = c(1,6))
  for(i in 1:length(yticksCol)){
    segments(7.5,yticksCol[i],7.7,yticksCol[i])
    text(7.8,yticksCol[i],labels=legendCol[i],adj=0,cex=0.8)
  }
  text(7.5,6.5,labels = "LRR")
  
  centerLegendPie = 1:6
  radiusLegendPie = scale.func(x0 = sqrt(sqrt(legendPie)),minmax = minmaxPie)
  for(i in 1:length(legendPie)){
    draw.circle(x=centerLegendPie[i],y=-0.2,radius=radiusLegendPie[i],nv=100,border=NULL,col="grey",lty=1,lwd=1)
    text(x=centerLegendPie[i],y=-0.8,labels=legendPie[i],cex=0.8)
  }
}

makeTransparent = function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}

