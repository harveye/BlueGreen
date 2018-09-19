#...................................#
#... Project Blue-Green (EMERGe 3) ...#
#...................................#

#.........................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt                   #
#... Leading Author: Eric Harvey                                                                          #
#... Author of the script: Isabelle Gounand                                                               #
#... Date: 18 November 2016                                                                               #
#... the script identify the species in mixtures and calculate the densities of the respective specie     #
#.........................................................................................................#


#... Clear variables
rm(list = ls())


#... Load libraries
#..................
#install.packages("devtools",dependencies = T)
library(devtools)
#install_github("pennekampster/bemovi", ref="master")
library(bemovi)
#install.packages("party",dependencies = T)
library(e1071)
#install.packages("moments",dependencies = T)
library(moments)
#install.packages("FactoMineR",dependencies = T)
library(FactoMineR)


#... Directory and file names
#............................
video.description.folder = "0_video_description/"
video.description.file = "video_description.txt"
merged.data.folder = "5_merged_data/"
to.data.folders = c("T0_20160502_Blue_Lanscapes/","T1_20160509_Blue_Landscapes/","T2_20160517_Blue_Landscapes/","T3_20160523_Blue_Landscapes/","T4_20160531_Blue_Landscapes/")

## Folders for Mac computer
general.folder="~/Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/2_DATA/"
to.data.mono = "~/Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/2_DATA/Monoculture_test_20160502/"
to.script ="~/Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/3_ANALYSIS"
to.plots = "~/Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/4_RESULTS"

## Folders for Linux computer
general.folder="~/Documents/FreeFileSync_Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/2_DATA/"
to.data.mono ="~/Documents/FreeFileSync_Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/2_DATA/Monoculture_test_20160502/"
to.script = "/home/isabelle/Documents/FreeFileSync_Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/3_ANALYSIS/"
to.plots = "/home/isabelle/Documents/FreeFileSync_Documents/Recherche/PROJECTS/EW10_BLUE-GREEN/4_RESULTS/"


#... Load home-made functions
#............................
setwd(to.script)
source("BG_functions_20161117.R")


#... Parameters (those used to obtain the traits)
#..............
today = Sys.Date() 
fps = 25 #... Video frame rate (in frames per second)
nsv = 5 #... Number of seconds per videos
measured_volume = 34.4 #... Volume in µL, for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
pixel_to_scale = 4.05 #... Size of a pixel in µm, for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
filter_min_net_disp = 20 #... minimum net displacement (µm)
filter_min_duration = 0.5 #... duration (s)
filter_detection_freq = 0.1 #... detection frequency
filter_median_step_length = 3 #... median step length [µm]


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
#...........................................................#
#... IDENTIFICATION MODELS: DETERMINATION AND VALIDATION ...#
#...........................................................#

{
    ##########################################################################################################################################
    #..................................................#
    #... Create the training dataset (Monocultures) ...#
    #..................................................#
    {
      #... Select folders where to pick the monocultures (sampling date folders; we chose to pool all the monocultures at all dates)
      #..................
      to.data.folders.sel = to.data.folders
      
      #... List to record the monoculture data
      monoculture.list = vector(mode = "list", length = length(to.data.folders.sel))
      
      
      #... Get the individual-based data of the monocultures and add momments of the trait distributions
      #.................................
      for(i in 1:length(to.data.folders.sel))
      {
        to.data = paste(general.folder,to.data.folders.sel[i],sep="")
        
        #... Load the merged data
        load(paste0(to.data, merged.data.folder, "Master.RData"))
        
        #... Filter data
        trajectory.data.filtered = filter_data(trajectory.data, filter_min_net_disp, filter_min_duration, filter_detection_freq, filter_median_step_length)
        
        #... Summarize trajectory data to individual-based data
        morph_mvt0 = summarize_trajectories(trajectory.data.filtered, write = T, to.data, calculate.median=F, merged.data.folder)
      
        #... Add the moments of morphology trait distributions (k, sk, q25, q75)  
        morph_mvt1 = add.momentsMorphology.to.indData(fdat = trajectory.data.filtered, mm = morph_mvt0)
        
        #... Add the moments of movement trait distributions (k, sk, q25, q75) 
        morph_mvt = add.momentsMovement.to.indData(fdat = trajectory.data.filtered, mm = morph_mvt1)
        
        monoculture.list[[i]] = subset(morph_mvt, Composition != "mixture")
      }
      
      #... Pool all the monocultures in a single dataset
      #.................................................
      morph_mvt_monos = monoculture.list[[1]] ; if(length(to.data.folders.sel)>1)for(i in 2:length(to.data.folders))morph_mvt_monos = rbind(morph_mvt_monos,monoculture.list[[i]])
      
      
      # ... Differenciate Ble from Chi in monocultures and merge with other monocultures in a single training dataset
      #..............................................................................................................
      
      #... plot to choose discriminating criteria
      ble = subset(morph_mvt_monos,Composition == "Ble")
      plot(ble$mean_minor,ble$mean_major)
      abline(h=45,v=35,lty=3)
      
      morph_mvt_noBle = subset(morph_mvt_monos,Composition != "Ble")
      morph_mvt_Ble = subset(morph_mvt_monos,Composition == "Ble" & mean_major >45)
      morph_mvt_Chi1 = subset(morph_mvt_monos,Composition == "Ble" & mean_major < 45)
      morph_mvt_Chi1$Composition = "Chi"  # small species are identified as Chi in the "Composition" column
      
      #... Create training dataset with no NAs in data
      #...............................................
      training_data_complete = rbind(morph_mvt_noBle,morph_mvt_Ble,morph_mvt_Chi1)
      training_data = training_data_complete[complete.cases(training_data_complete), ]
    }
    
    
    ##########################################################################################################################################
    #...........................................................#
    #... Preliminary analysis: Species trait distribution    ...#
    #...........................................................#
    
    #....................................#
    #... Plot the trait distributions ...#
    #....................................#
    
    #... Positions of traits in the dataset
    traits = c(2:16,18:28,38:75) # with all new traits
    traits.names = colnames(training_data)[traits]
    
    {
      #... Species names, total and the ones selected for the graphs
      speciesNames = c("Tet","Chi","Col","Pca","Rot","Ble","Spi")
      speciesNames2 = c("Tet","Chi")
      
      
      #... Species data
      species.list = vector(mode = "list",length(speciesNames2))
      for(i in 1:length(speciesNames2))species.list[[i]] = subset(training_data,Composition == speciesNames2[i])
      names(species.list) = speciesNames2
      
      #... Species colours
      col.list=c("yellow","chartreuse3", "cyan","purple","darkred","brown2","blue")[which(speciesNames %in% speciesNames2)]
      col.list2=col.list ; for(i in 1:length(col.list))col.list2[i] = makeTransparent(col.list[i],alpha=0.3)
      
      
      #... change working directory
      setwd(to.plots)
      
      #... bandwidth for each trait (degree of smoothing )
      bws=c(1,0.5,5,1,2,1,1,0.3,0.5,0.2,
            0.05,0.01,0.05,0.04,0.09,10,10,7,10,3,
            0.5,0.1,5,0.4,1,0.5,10,10,10,10,
            0.1,0.2,2,0.02,0.08,0.01,0.1,0.1,0.5,0.5,
            0.01,1,5,1,0.5,0.05,0.3,0.1,1,1,
            0.1,0.05,0.3,0.1,10,10,0.5,10,25,25,
            0.1,20,0.1,0.1) # for Tet Chi combination
    
      
      #... define a list for each trait with the values for the different species
      trait.species.list = vector(mode="list",length = length(traits))
      
      #... Plot series of trait distributions
      #......................................
      
      pdf(paste("BG_trait_distributions_NewTraitsChosen_",paste0(speciesNames2,collapse = "_"),today,".pdf"),width=12)
      ct=1 #... counter
      for(i in traits)
      {
        density.list = vector(mode="list",length = length(species.list))
        traitvalue.list = vector(mode="list",length = length(species.list))
        xlims = range(species.list[[1]][,i],na.rm=T)
        ylims = numeric(2)
        
        #... get the density distribution data
        for(j in 1:length(species.list))
        { 
          trait = species.list[[j]][,i]
          xy= density(trait[which(!is.na(trait))],from=min(trait,na.rm=T),to=max(trait,na.rm=T),bw=bws[ct])
          if(min(xy$x) < xlims[1]) xlims[1] = min(xy$x)
          if(max(xy$x) > xlims[2]) xlims[2] = max(xy$x)
          if(max(xy$y) > ylims[2]) ylims[2] = max(xy$y)    
          density.list[[j]] = xy
          traitvalue.list[[j]] = trait
        }  
        trait.species.list[[ct]] = traitvalue.list
    
        #... plot them
        plot(NA,type="n", xlim = xlims, ylim = ylims, xlab = "trait value", ylab = "density", main = colnames(training_data)[i])
        for(j in 1:length(species.list))
        {
          xx=density.list[[j]]$x
          yy=density.list[[j]]$y
          polygon(c(xx[1],xx,xx[length(xx)],xx[length(xx)],rev(xx),xx[1]),c(yy[1],yy,yy[length(yy)],0,numeric(length(yy)),0),col=col.list2[j],lty=0)
          lines(xx,yy,col=col.list[j],lwd = 2)
          abline(v = median(traitvalue.list[[j]],na.rm=T),lty = 3,col = col.list[j])
        }
        legend("topright",lty=rep(1,length(species.list)),col=col.list,legend=speciesNames2,bty="n",cex=0.8)
        ct = ct + 1
      }
      dev.off()
    }
    
    
    #.....................#
    #... PCA on traits ...#
    #.....................#
    
    {
      speciesNames2 = c("Tet","Chi","Col","Pca","Rot","Ble","Spi")
      col.list=c("yellow","chartreuse3", "cyan","purple","darkred","brown2","blue")[which(speciesNames %in% speciesNames2)]
      col.list2=col.list ; for(i in 1:length(col.list))col.list2[i] = makeTransparent(col.list[i],alpha=0.3)
      
      
      #... Select data
      selected.traits = c("mean_grey","sd_grey","mean_area","sd_area","mean_perimeter","sd_perimeter","mean_major","sd_major","mean_minor","sd_minor","mean_ar",
                          "sd_ar","mean_turning","sd_turning","duration","max_net","net_disp","net_speed","gross_disp" ,"gross_speed","max_step","min_step",
                          "sd_step","sd_gross_speed","max_gross_speed","min_gross_speed","k_area","k_perimeter","k_major","k_minor","k_ar","sk_area","sk_perimeter",
                          "sk_major","sk_minor","sk_ar", "q25_area","q25_perimeter","q25_major","q25_minor","q25_ar","q75_area","q75_perimeter","q75_major","q75_minor","q75_ar")
      
      selected.data =  training_data[,c("Composition",selected.traits)]
      
      #... PCA
      acp.traits = PCA(selected.data,ncp=length(selected.traits),quali.sup=1)
      
      #... PCA ouputs 
      v = acp.traits$var$coord #... Coordinates of variables (end of arrows)
      inertie = acp.traits$eig$`percentage of variance`#... Percent of variance explained on each dimentsion
      ind = acp.traits$ind$coord #... Coordinates of individuals (with names in row.names) -coord, cos2, contrib, dist
      qs = acp.traits$quali.sup$coord #... Qualitative variables, coordinates of all - coord, cos2, v.test, dist, eta2
      dims = c(1,2) #... dimensions to be plotted
      
      pdf("PCA_species-traits_newtraits.pdf",width=12)
      layout(matrix(1:2,nrow =1))
    
      #... 1st Plot: 2-dimension graph with individuals (colours are species)
      #................................................
      plot(ind[,paste("Dim",dims[1],sep=".")]*1.2,ind[,paste("Dim",dims[2],sep=".")]*1.2,type="n",xlab=paste("Axe 1 (",round(inertie[dims[1]],1),"%)",sep=""),ylab=paste("Axe 2 (",round(inertie[dims[2]],1),"%)",sep=""))
      abline(h=0,v=0,lty=2)
      for(i in 1:nrow(ind)){
        col.e = col.list[which(speciesNames == (selected.data$Composition)[i])]
        points(ind[i,dims[1]],ind[i,dims[2]],pch=21,col = col.e, bg = makeTransparent(col.e,alpha=0.3),cex=1)
      }
      
      #... Trace the ellipses of confidence
          #... Calculate the ellipse
          concat.species = cbind.data.frame(selected.data$Composition,acp.traits$ind$coord)
          ellipses.species = coord.ellipse(concat.species,bary=TRUE)
          
          #... colours for centroids
          col.list3 = c("gray30","white","black","black","white","black","black","white")
          
          #... Trace the ellipses and their centroids
          for(i in 1: length(speciesNames)){
            el =ellipses.species$res[ which((ellipses.species$res)[,1] == speciesNames[i]),]
            col.e = col.list[i]
            polygon(el[,2],el[,3],col=makeTransparent("white",alpha=0.3),lty=0)
            polygon(el[,2],el[,3],col=makeTransparent(col.e,alpha=0.3),lty=0)
            lines(el[,2],el[,3],col=col.e)
            points(qs[speciesNames[i],1],qs[speciesNames[i],2],pch=3,col=col.list3[i])
            text(qs[speciesNames[i],1],qs[speciesNames[i],2]*1.1,adj = -0.2,cex=0.8,label=speciesNames[i],col=col.list3[i])
          }
    
      legend("topright",legend = speciesNames, fill = col.list, bty="n",border="white",cex=1)
     
       
      #... 2nd Plot: Traits 
      #....................
      graph.var(acp.traits,draw=c("var"),cex=1)
      
      dev.off()
    }  
    
    
    ##########################################################################################################################################
    #........................................................................#
    #... First set of models, to identify large species (not Tet and Chi) ...#
    #........................................................................#
    
    #... The models
    #..............
    {
    
        #... Select traits (we keep all, it's better)
        #.................
        expl.vars = traits.names
        
        #... Complete model (with all the individuals included)
        #..................
        svm0 = svm(as.formula(paste("as.factor(Composition) ~ ", paste(expl.vars, collapse="+"))),
                   data=training_data, probability=T,na.action=na.pass)
        
        
                   
        #... Partial models specific to some species identification (select typical individuals for training datasets)
        #..........................................................
        
        #... The different complete monocultures
        #.......................................
        morph_mvt_Tet0 = subset(training_data,Composition == "Tet")  
        morph_mvt_Chi0 = subset(training_data,Composition == "Chi")
        morph_mvt_Col0 = subset(training_data,Composition == "Col")
        morph_mvt_Pca0 = subset(training_data,Composition == "Pca")
        morph_mvt_Rot0 = subset(training_data,Composition == "Rot")
        morph_mvt_Ble0 = subset(training_data,Composition == "Ble")
        morph_mvt_Spi0 = subset(training_data,Composition == "Spi")
        
        
        #... Tet specific model
        #......................
        morph_mvt_Chi = subset(morph_mvt_Chi0, mean_minor < 20 & mean_major<32&net_speed < 100 & max_gross_speed < 450 & min_gross_speed <70)
        training_data_Tet = rbind(morph_mvt_Tet0,morph_mvt_Chi,morph_mvt_Col0,morph_mvt_Pca0,morph_mvt_Rot0,morph_mvt_Ble0,morph_mvt_Spi0)
        svmTet = svm(as.formula(paste("as.factor(Composition) ~ ", paste(expl.vars, collapse="+"))),
                     data=training_data_Tet, probability=T,na.action=na.pass)
        
        
        #... Chi specific model (1)
        #......................
        morph_mvt_Tet = subset(morph_mvt_Tet0, mean_turning > -0.3 & mean_turning < 0.2 & sd_turning < 0.6 & net_speed >50)
        morph_mvt_Chi = subset(morph_mvt_Chi0, mean_minor < 14)
        morph_mvt_Spi =  subset(morph_mvt_Spi0, mean_grey >50 & sd_perimeter > 50& mean_ar > 3.7&mean_major > 120& gross_speed <450)
        training_data_Chi = rbind(morph_mvt_Tet,morph_mvt_Chi,morph_mvt_Col0,morph_mvt_Pca0,morph_mvt_Rot0,morph_mvt_Ble0,morph_mvt_Spi)
        svmChi = svm(as.formula(paste("as.factor(Composition) ~ ", paste(expl.vars, collapse="+"))),
                     data=training_data_Chi, probability=T,na.action=na.pass)
        
        
        #... Chi specific model (2)
        #......................
        morph_mvt_Tet = subset(morph_mvt_Tet0, mean_turning > -0.3 & mean_turning < 0.2 & sd_turning < 0.7 & net_speed >50)
        morph_mvt_Spi =  subset(morph_mvt_Spi0, mean_grey >50 & sd_perimeter > 50& mean_ar > 3.7&mean_major > 120& gross_speed <450)
        training_data_Chi2 = rbind(morph_mvt_Tet,morph_mvt_Chi0,morph_mvt_Col0,morph_mvt_Pca0,morph_mvt_Rot0,morph_mvt_Ble0,morph_mvt_Spi)
        svmChi2 = svm(as.formula(paste("as.factor(Composition) ~ ", paste(expl.vars, collapse="+"))),
                      data=training_data_Chi2, probability=T,na.action=na.pass)
        
        
        #... Pca (-Rot) specific training dataset
        #..................................
        morph_mvt_Pca = subset(morph_mvt_Pca0, mean_ar > 1.9&q75_major>100)
        morph_mvt_Rot = subset(morph_mvt_Rot0, mean_major < 100& sk_step_speed>0.5&q75_major<120 &sk_perimeter>1)
        morph_mvt_Spi =  subset(morph_mvt_Spi0, mean_grey >50 & sd_perimeter > 50& mean_ar > 3.7&mean_major > 120& gross_speed <450)
        training_data_Pca = rbind(morph_mvt_Tet0,morph_mvt_Chi0,morph_mvt_Col0,morph_mvt_Pca,morph_mvt_Rot,morph_mvt_Ble0,morph_mvt_Spi)
        svmPca = svm(as.formula(paste("as.factor(Composition) ~ ", paste(expl.vars, collapse="+"))),
                     data=training_data_Pca, probability=T,na.action=na.pass)
        
    }
    
    #... MODEL COMBINATION VALIDATION
    #................................
    #... Principle of the analysis: we run each model and for each individual we keep the prediction which is predicted with the highest confidence (proba per species)
    #... we apply the model on the complete set of monocultures and the error rates is the difference between known composition and prediction
    
    {
        #... Species and data on which we test the model accuracy
        speciesNames = c("Tet","Chi","Col","Pca","Rot","Ble","Spi")
        data.to.predict = training_data
        
        
        #... Record the probabilities of attribution of each species for each individual, for each model
        x0_prob = attr(predict(svm0,OOB=F,data.to.predict, probability=T),"prob")
        xTet_prob = attr(predict(svmTet,OOB=F,data.to.predict, probability=T),"prob")
        xChi_prob = attr(predict(svmChi,OOB=F,data.to.predict, probability=T),"prob")
        xChi2_prob = attr(predict(svmChi2,OOB=F,data.to.predict, probability=T),"prob")
        xPca_prob = attr(predict(svmPca,OOB=F,data.to.predict, probability=T),"prob")
        
        #... Record the prediction of the complete model to error rates with our combination of models
        x0_resp = predict(svm0,OOB=F,data.to.predict, type="response")
        
        #... Objects to record outputs of model predictions
        model.probs = list(x0_prob,xTet_prob,xChi_prob,xPca_prob,xChi2_prob ) # list of the model proba outputs
        nmodels = length(model.probs) # the number of models used
        x_best = vector(mode="character",length = nrow(x0_prob)) # array to record the species attribution based on the better proba on all models for each individual
        model_congruence = best_model = vector(mode="numeric",length = nrow(x0_prob)) # array to record if models attribute the same species (1) or not (0) for each individual
        model_predictions = matrix(NA,ncol = nmodels,nrow = nrow(x0_prob),dimnames = list(NULL,paste("M",1:nmodels,sep=""))) # matrix to record the predictions for each individual per model
        model_best_probs = matrix(NA,ncol = nmodels,nrow = nrow(x0_prob),dimnames = list(NULL,paste("pM",1:nmodels,sep=""))) # matrix to record the best proba for each individual per model
        
        
        #... Record predictions
        #......................
        
        #... Per individual
        for(i in 1:nrow(x0_prob)){
          
            #... Record the prediction and best proba for each model
            for(j in 1:nmodels){
              mp = model.probs[[j]][i,]
              model_best_probs[i,j] = max(mp)
              model_predictions[i,j] = names(mp)[which(mp == max(mp))]
            }
            
            #... Determine the better model (the one with the highest best proba)
            whichmaxp = which(model_best_probs[i,] == max(model_best_probs[i,]))
            
            #... Determine the better prediction (species predicted by the best model)
            if(length(whichmaxp) == 1){
              best_model[i] = whichmaxp
              x_best[i] = model_predictions[i,whichmaxp]
            }else{
              best_model[i] = NA
              x_best[i] = unique(model_predictions[i,which(model_best_probs[i,] == max(model_best_probs[i,]))])
            }
            
            #... Determine if the model predictions are congruent
            model_congruence[i] = ifelse(length(unique(model_predictions[i,])) == 1,1,0)
        }
    
        #... Add the outputs to the to-predict dataset
        predictions.step1 = cbind(data.to.predict,model_predictions,model_best_probs,x_best,best_model,model_congruence)
        
        
        
        #... Analyse predictions
        #.......................
        
        #... Define and format the error tables
        error.table.bestModels = matrix (NA,nrow = (length(speciesNames)+3), ncol = length(speciesNames), dimnames = list(c("n",paste("pred",speciesNames,sep="-"),"nError","pError"),speciesNames))
        error.table0 = error.table.bestModels
        
        #... Select the known and predicted species
        x = x_best
        y = data.to.predict$Composition
        
        #... Fill the error tables
        for(i in 1:length(speciesNames))
        {
          #... Attribution and error rates for the complete model alone
          pos0 = which(y == speciesNames[i])
          error.table0[1,i] = length(pos0) 
          for(j in 1:length(speciesNames))error.table0[j+1,i] = sum(x0_resp[pos0] == speciesNames[j])
          error.table0["nError",i] = sum(x0_resp[pos0] != speciesNames[i])
          error.table0["pError",i] = sum(x0_resp[pos0] != speciesNames[i])/length(pos0) 
          
          #... Attribution and error rates for the combination of models
          pos = which(y == speciesNames[i])
          error.table.bestModels[1,i] = length(pos)
          for(j in 1:length(speciesNames))error.table.bestModels[j+1,i] = sum(x[pos] == speciesNames[j])
          error.table.bestModels["nError",i] = sum(x[pos] != speciesNames[i])
          error.table.bestModels["pError",i] = sum(x[pos] != speciesNames[i])/length(pos)
        }
    
        #... Display the error tables
        error.table0
        error.table.bestModels
    }
   
     
    ##########################################################################################################################################
    #............................................................#
    #... Second set of models, to identify better Tet and Chi ...#
    #............................................................#
    #... Trait space overlap between Tet and Chi leads to cross misclassifications that could be reduced a bit if treated with separate models
    #... We select the individuals which have been classified as Tet or Chi (predicted Tet and Chi are rarely another species than one another)
    #... Then we applied other models for which all other species are grouped in a class "Other" in the training dataset.
    
    #... The models
    #..............
    
    {
        #... Selected traits
        #...................
        expl.vars = traits.names
        
        
        #... New complete model (with all the individuals included / all other species than Tet and Chi are identified as "Other" in the Composition column)
        #......................
        training_data_pTC_complete = training_data
        training_data_pTC_complete$Composition = ifelse(training_data_pTC_complete$Composition=="Chi","Chi",ifelse(training_data_pTC_complete$Composition=="Tet","Tet","Other"))
        svm_pTC = svm(as.formula(paste("as.factor(Composition) ~ ", paste(expl.vars, collapse="+"))),
                      data=training_data_pTC_complete, probability=T,na.action=na.pass)
        
        #... Partial model (1)
        #....................
        morph_mvt_Tet = subset(training_data_pTC_complete,Composition == "Tet"  )
        morph_mvt_Chi = subset(training_data_pTC_complete,Composition == "Chi" & mean_minor < 18 & mean_major<35&net_speed < 270 & gross_speed<300 & min_step<3 )
        morph_mvt_Oth = subset(training_data_pTC_complete, Composition == "Other")
        training_data_p1 = rbind(morph_mvt_Tet,morph_mvt_Chi,morph_mvt_Oth)
        svm_p1 = svm(as.formula(paste("as.factor(Composition) ~ ", paste(expl.vars, collapse="+"))),
                      data=training_data_p1, probability=T,na.action=na.pass)
        
        #... Partial model (2)
        #....................
        morph_mvt_Tet = subset(training_data_pTC_complete,Composition == "Tet" & mean_turning > -0.3& mean_minor < 25 & mean_turning < 0.2 & sd_turning < 0.7 & net_speed >45)
        morph_mvt_Chi = subset(training_data_pTC_complete,Composition == "Chi" & mean_minor> 13)
        morph_mvt_Oth = subset(training_data_pTC_complete, Composition == "Other")
        training_data_p2 = rbind(morph_mvt_Tet,morph_mvt_Chi,morph_mvt_Oth)
        svm_p2 = svm(as.formula(paste("as.factor(Composition) ~ ", paste(expl.vars, collapse="+"))),
                     data=training_data_p2, probability=T,na.action=na.pass)
    }
    
    
    #... MODEL COMBINATION VALIDATION
    #................................
    {
        #... select the individuals that have been identified as Tet or Chi
        #..................................................................
        speciesNames2 = c("Tet","Chi","Other")
        data.to.predict2 =  subset(predictions.step1, x_best %in% c("Tet","Chi"))
        
        #... Record the probabilities of attribution of each species for each individual, for each model
        x_pTC_prob = attr(predict(svm_pTC,OOB=F,data.to.predict2, probability=T),"prob")
        x_p1_prob = attr(predict(svm_p1,OOB=F,data.to.predict2, probability=T),"prob")
        x_p2_prob = attr(predict(svm_p2,OOB=F,data.to.predict2, probability=T),"prob")
        
        
        #... Objects to record outputs of model predictions
        model.probs = list(x_p1_prob,x_pTC_prob,x_p2_prob) # list of the model proba outputs
        nmodels = length(model.probs) # the number of models used
        x_best = vector(mode="character",length = nrow(x_pTC_prob))  # array to record the species attribution based on the better proba on all models for each individual
        model_congruence = best_model = vector(mode="numeric",length = nrow(x_pTC_prob)) # array to record if models attribute the same species (1) or not (0) for each individual
        model_predictions = matrix(NA,ncol = nmodels,nrow = nrow(x_pTC_prob),dimnames = list(NULL,paste("M",1:nmodels,sep=""))) # matrix to record the predictions for each individual per model
        model_best_probs = matrix(NA,ncol = nmodels,nrow = nrow(x_pTC_prob),dimnames = list(NULL,paste("pM",1:nmodels,sep=""))) # matrix to record the best proba for each individual per model
        
        
        
        #... Record predictions
        #......................
        
        #... Per individual
        for(i in 1:nrow(x_pTC_prob)){
          
          #... Record the prediction and best proba for each model
          for(j in 1:nmodels){
            mp = model.probs[[j]][i,]
            model_best_probs[i,j] = max(mp)
            model_predictions[i,j] = names(mp)[which(mp == max(mp))]
          }
          
          #... Determine the better model (the one with the highest best proba)
          whichmaxp = which(model_best_probs[i,] == max(model_best_probs[i,]))
          
          #... Determine the better prediction (species predicted by the best model)
          if(length(whichmaxp) == 1){
            best_model[i] = whichmaxp
            x_best[i] = model_predictions[i,whichmaxp]
          }else{
            best_model[i] = NA
            x_best[i] = unique(model_predictions[i,which(model_best_probs[i,] == max(model_best_probs[i,]))])
          }
          
          #... Determine if the model predictions are congruent
          model_congruence[i] = ifelse(length(unique(model_predictions[i,])) == 1,1,0)
        }
        
        #... Add the outputs to the to-predict dataset
        predictions.step2 = cbind(data.to.predict2[,1:75],model_predictions,model_best_probs,x_best,best_model,model_congruence)
        
        
        #... Analyse predictions
        #.......................
        
        #... Define and format the error table
        error.table.bestModels2 = matrix (NA,nrow = (length(speciesNames2)+3), ncol = length(speciesNames2), dimnames = list(c("n",paste("pred",speciesNames2,sep="-"),"nError","pError"),speciesNames2))
        
        #... Select the known and predicted species
        x = x_best
        y = data.to.predict2$Composition
        
        #... Fill the error table
        for(i in 1:length(speciesNames2))
        {
          pos = which(y == speciesNames2[i])
          error.table.bestModels2[1,i] = length(pos)
          for(j in 1:length(speciesNames2))error.table.bestModels2[j+1,i] = sum(x[pos] == speciesNames2[j])
          error.table.bestModels2["nError",i] = sum(x[pos] != speciesNames2[i])
          error.table.bestModels2["pError",i] = sum(x[pos] != speciesNames2[i])/length(pos)
        }
        
        #... Display the error table
        error.table.bestModels2
    
    }
    
    
    ##########################################################################################################################################
    #........................#
    #... FINAL VALIDATION ...#
    #........................#
    {
        #... Combine the predicted data
        #..............................
        Large.species.predicted = subset(predictions.step1[,c(1:75,86)], x_best != c("Tet") & x_best != c("Chi"))
        Tet.Chi.predicted = predictions.step2[,c(1:75,82)]
        final_validation.data = rbind(Large.species.predicted,Tet.Chi.predicted)
        
        
        #... Analyse predictions
        #.......................
        
        #... Define and format the error table
        speciesNames = c("Tet","Chi","Col","Pca","Rot","Ble","Spi","Other")
        error.table.final = matrix (NA,nrow = (length(speciesNames)+3), ncol = length(speciesNames), dimnames = list(c("n",paste("pred",speciesNames,sep="-"),"nError","pError"),speciesNames))
        
        #... Select the known and predicted species
        x = final_validation.data$x_best
        y = final_validation.data$Composition
        
        #... Fill the error table
        for(i in 1:length(speciesNames))
        {
          pos = which(y == speciesNames[i])
          error.table.final[1,i] = length(pos)
          for(j in 1:length(speciesNames))error.table.final[j+1,i] = sum(x[pos] == speciesNames[j])
          error.table.final["nError",i] = sum(x[pos] != speciesNames[i])
          error.table.final["pError",i] = sum(x[pos] != speciesNames[i])/length(pos)
        }
        
        #... Display the error table
        error.table.final
    }
    
}

##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
#........................................................................................................................................................
#................................................................#
#... IDENTIFY SPECIES AND CALCULATE DENSITIES IN VIDEO SERIES ...#
#................................................................#
{
  #... for each sampling date
  measure.date = 0
  for(m in 1:length(to.data.folders))
  {
    
    #... LOAD AND PREPARE TRAIT DATA (ADD THE NEW TRAIT DISTRIBUTION MOMENTS)
    #...............................
    {
      to.data = paste(general.folder,to.data.folders[m],sep="")
      
      #... Load the merged data
      load(paste0(to.data, merged.data.folder, "Master.RData"))
      
      #... Filter data
      trajectory.data.filtered = filter_data(trajectory.data, filter_min_net_disp, filter_min_duration, filter_detection_freq, filter_median_step_length)
      
      #... Summarize trajectory data to individual-based data
      morph_mvt0 = summarize_trajectories(trajectory.data.filtered, write = T, to.data, calculate.median=F, merged.data.folder)
      
      #... Add the moments of morphology trait distributions (k, sk, q25, q75)  
      morph_mvt1 = add.momentsMorphology.to.indData(fdat = trajectory.data.filtered, mm = morph_mvt0)
      
      #... Add the moments of movement trait distributions (k, sk, q25, q75) 
      morph_mvt = add.momentsMovement.to.indData(fdat = trajectory.data.filtered, mm = morph_mvt1)
    }
    
    
    #... FIRST STEP OF PREDICTIONS (LARGE SPECIES)
    #.............................   
    {    
        #... Species names
        speciesNames = c("Tet","Chi","Col","Pca","Rot","Ble","Spi")
        
        #... Data to predict
        data.to.predict = morph_mvt[complete.cases(morph_mvt), ]
        
        
        #... Record the probabilities of attribution of each species for each individual, for each model
        x0_prob = attr(predict(svm0,OOB=F,data.to.predict, probability=T),"prob")
        xTet_prob = attr(predict(svmTet,OOB=F,data.to.predict, probability=T),"prob")
        xChi_prob = attr(predict(svmChi,OOB=F,data.to.predict, probability=T),"prob")
        xChi2_prob = attr(predict(svmChi2,OOB=F,data.to.predict, probability=T),"prob")
        xPca_prob = attr(predict(svmPca,OOB=F,data.to.predict, probability=T),"prob")
        
        #... Objects to record outputs of model predictions
        model.probs = list(x0_prob,xTet_prob,xChi_prob,xPca_prob,xChi2_prob ) # list of the model proba outputs
        nmodels = length(model.probs) # the number of models used
        x_best = vector(mode="character",length = nrow(x0_prob)) # array to record the species attribution based on the better proba on all models for each individual
        model_predictions = matrix(NA,ncol = nmodels,nrow = nrow(x0_prob),dimnames = list(NULL,paste("M",1:nmodels,sep=""))) # matrix to record the predictions for each individual per model
        model_best_probs = matrix(NA,ncol = nmodels,nrow = nrow(x0_prob),dimnames = list(NULL,paste("pM",1:nmodels,sep=""))) # matrix to record the best proba for each individual per model
        
        #... Record predictions
        #......................
        
        #... Per individual
        for(i in 1:nrow(x0_prob)){
          
          #... Record the prediction and best proba for each model
          for(j in 1:nmodels){
            mp = model.probs[[j]][i,]
            model_best_probs[i,j] = max(mp)
            model_predictions[i,j] = names(mp)[which(mp == max(mp))]
          }
          
          #... Determine the better model (the one with the highest best proba)
          whichmaxp = which(model_best_probs[i,] == max(model_best_probs[i,]))
          
          #... Determine the better prediction (species predicted by the best model)
          if(length(whichmaxp) == 1){
            best_model[i] = whichmaxp
            x_best[i] = model_predictions[i,whichmaxp]
          }else{
            best_model[i] = NA
            x_best[i] = unique(model_predictions[i,which(model_best_probs[i,] == max(model_best_probs[i,]))])
          }
        }
        
        #... Add the outputs to the to-predict dataset
        predictions.step1 = cbind(data.to.predict,x_best)
    }
 
  
         
    #... SECOND STEP OF PREDICTIONS (TET AND CHI)
    #..............................   
    {
      
      #... select the individuals that have been identified as Tet or Chi
      data.to.predict2 =  subset(predictions.step1, x_best %in% c("Tet","Chi"))
      
      
      #... Record the probabilities of attribution of each species for each individual, for each model
      x_pTC_prob = attr(predict(svm_pTC,OOB=F,data.to.predict2, probability=T),"prob")
      x_p1_prob = attr(predict(svm_p1,OOB=F,data.to.predict2, probability=T),"prob")
      x_p2_prob = attr(predict(svm_p2,OOB=F,data.to.predict2, probability=T),"prob")
      
      
      #... Objects to record outputs of model predictions
      model.probs = list(x_pTC_prob,x_p1_prob,x_p2_prob) # list of the model proba outputs
      nmodels = length(model.probs) # the number of models used
      x_best = vector(mode="character",length = nrow(x_pTC_prob))  # array to record the species attribution based on the better proba on all models for each individual
      model_predictions = matrix(NA,ncol = nmodels,nrow = nrow(x_pTC_prob),dimnames = list(NULL,paste("M",1:nmodels,sep=""))) # matrix to record the predictions for each individual per model
      model_best_probs = matrix(NA,ncol = nmodels,nrow = nrow(x_pTC_prob),dimnames = list(NULL,paste("pM",1:nmodels,sep=""))) # matrix to record the best proba for each individual per model
      
      
      
      #... Record predictions
      #......................
     
      #... Per individual
      for(i in 1:nrow(x_pTC_prob)){
        
        #... Record the prediction and best proba for each model
        for(j in 1:nmodels){
          mp = model.probs[[j]][i,]
          model_best_probs[i,j] = max(mp)
          model_predictions[i,j] = names(mp)[which(mp == max(mp))]
        }
        
        #... Determine the better model (the one with the highest best proba)
        whichmaxp = which(model_best_probs[i,] == max(model_best_probs[i,]))
        
        #... Determine the better prediction (species predicted by the best model)
        if(length(whichmaxp) == 1){
          best_model[i] = whichmaxp
          x_best[i] = model_predictions[i,whichmaxp]
        }else{
          best_model[i] = NA
          x_best[i] = unique(model_predictions[i,which(model_best_probs[i,] == max(model_best_probs[i,]))])
        }
      }
      
      #... Add the outputs to the to-predict dataset
      predictions.step2 = cbind(data.to.predict2[,1:75],x_best)
      
    }
 
       
    #... COMBINE THE PREDICTED DATA
    #..............................   
    Large.species.predicted = subset(predictions.step1, x_best != c("Tet") & x_best != c("Chi"))
    Tet.Chi.predicted = predictions.step2
    final.predicted.data = rbind(Large.species.predicted,Tet.Chi.predicted)
    final.predicted.data$predicted_species = final.predicted.data$x_best
    
    
    
    #... CALCULATE THE DENSITY OF EACH SPECIES
    #.........................................   
    #... Species names
    species.names = as.character(unique(final.predicted.data$x_best))
    
    #... Population level data
    pop.data = summarize_populations2(trajectory.data.filtered, morph_mvt, write=T, to.data, merged.data.folder, video.description.folder, video.description.file,total_frame = fps * nsv)
    
    #... Add the population density for each species
    output = species.density(pop.data,final.predicted.data,species.names,total_frames = fps*nsv, mv = measured_volume)
    

    #... RECORD THE OUTPUT BY VIDEO SERIE
    #....................................
    setwd(to.data)
    write.table(output,paste("BG_protist_blue_",measure.date,".txt",sep=""),sep="\t")
    
    measure.date = measure.date + 1
  }
  
  
  #... COMPILE THE SAMPLING DATE DATA INTO A SINGLE FILE
  #.....................................................
  bg.data = output[1,]
  measure.date = 0
  for(i in  1:length(to.data.folders))
  {
    to.data = paste(general.folder,to.data.folders[i],sep="")
    setwd(to.data)
    sd.data = read.table(paste("BG_protist_blue_",measure.date,".txt",sep=""))
    bg.data = rbind(bg.data,sd.data)
    measure.date = measure.date + 1
  }
  bg.data = bg.data[-1,]

  
  #... RECORD THE DATA
  #...................
  setwd(general.folder)
  write.table(bg.data,paste("BG_protist_data_",today,".txt"),col.names=T,row.names=F,sep="\t")
  
}

head(bg.data)
