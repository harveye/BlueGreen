#...................................#
#... Project Blue-Green (EMERGe 3) ...#
#...................................#

#.........................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt                   #
#... Leading Author: Eric Harvey                                                                          #
#... Author of the script: Isabelle Gounand                                                               #
#... Date: 17 November 2016                                                                               #
#... functions used in the BG_dat_processing_SVM_final.Rscript                                            #
#.........................................................................................................#



#... Function to add other moments of morphology trait distribution per individual (based on trajectory data)
#.................................................................................
#... fdat are the filtered data per individual per frame (the one from which we want to extract new traits)
#... mm is the morph_mvt data with traits per individuals (the one to which we want to add new traits)

#fdat = as.data.frame(trajectory.data.filtered) ; mm = morph_mvt0; 

library(data.table)
add.momentsMorphology.to.indData = function(fdat,mm){
  mm2 = as.data.table(mm)
  data <- as.data.table(fdat)
  mm2[,id_:=id]
  data[,id_:=id]
  new.traits =    data[, list(
    k_area = kurtosis(Area),
    k_perimeter = kurtosis(Perimeter),
    k_major = kurtosis(Major),
    k_minor = kurtosis(Minor),
    k_ar = kurtosis(AR),
    sk_area = skewness(Area),
    sk_perimeter = skewness(Perimeter),
    sk_major = skewness(Major),
    sk_minor = skewness(Minor),
    sk_ar = skewness(AR), 
    q25_area = quantile(Area,probs = 0.25),
    q25_perimeter = quantile(Perimeter,probs = 0.25),
    q25_major = quantile(Major,probs = 0.25),
    q25_minor = quantile(Minor,probs = 0.25),
    q25_ar = quantile(AR,probs = 0.25), 
    q75_area = quantile(Area,probs = 0.75),
    q75_perimeter = quantile(Perimeter,probs = 0.75),
    q75_major = quantile(Major,probs = 0.75),
    q75_minor = quantile(Minor,probs = 0.75),
    q75_ar = quantile(AR,probs = 0.75)),
    by=id_]
  setkey(new.traits, id_)
  setkey(mm2, id_)
  mm3 <- merge(mm2,new.traits,by=c("id_"), all=T)
  # mm3$id <- mm3$id_
  # mm3$id_ <- NULL
  return(as.data.frame(mm3)[,-1])
}


#... Function to add other moments of morphology trait distribution per individual (based on trajectory data)
#.................................................................................
#... fdat are the filtered data per individual per frame (the one from which we want to extract new traits)
#... mm is the morph_mvt data with traits per individuals (the one to which we want to add new traits)

library(circular)
add.momentsMovement.to.indData = function(fdat,mm){
  mm2 = as.data.table(mm)
  data <- as.data.table(fdat)
  mm2[,id_:=id]
  data[,id_:=id]
  new.traits.mvt =    data[, list(k_net_disp = round(kurtosis(sqrt(net_disp), na.rm=T), digits=2),
                                  k_gross_disp = round(kurtosis(gross_disp, na.rm=T), digits=2),
                                  k_step_length = round(kurtosis(step_length, na.rm=T), digits=2),
                                  k_step_speed = round(kurtosis(step_speed, na.rm=T), digits=2),
                                  sk_net_disp = round(skewness(sqrt(net_disp), na.rm=T), digits=2),
                                  sk_gross_disp = round(skewness(gross_disp, na.rm=T), digits=2),
                                  sk_step_length = round(skewness(step_length, na.rm=T), digits=2),
                                  sk_step_speed = round(skewness(step_speed, na.rm=T), digits=2),
                                  q25_net_disp = round(quantile(sqrt(net_disp),probs = 0.25, na.rm=T), digits=2),
                                  q25_gross_disp = round(quantile(gross_disp,probs = 0.25, na.rm=T), digits=2),
                                  q25_step_length = round(quantile(step_length,probs = 0.25, na.rm=T), digits=2),
                                  q25_step_speed = round(quantile(step_speed,probs = 0.25, na.rm=T), digits=2),
                                  q75_net_disp = round(quantile(sqrt(net_disp),probs = 0.75, na.rm=T), digits=2),
                                  q75_gross_disp = round(quantile(gross_disp,probs = 0.75, na.rm=T), digits=2),
                                  q75_step_length = round(quantile(step_length,probs = 0.75, na.rm=T), digits=2),
                                  q75_step_speed = round(quantile(step_speed,probs = 0.75, na.rm=T), digits=2)),
                           by=id_]
  
  turning <- data[!is.na(rel_angle), list(q25_turning= round(as.numeric(quantile.circular(as.circular(rel_angle, control.circular=list(type='angles', units="radians", template='none', modulo='asis',zero=0, rotation='counter')),probs = 0.25)),2), 
                                          q75_turning= round(as.numeric(quantile.circular(as.circular(rel_angle, control.circular=list(type='angles', units="radians", template='none', modulo='asis',zero=0, rotation='counter')),probs = 0.75)),2)), by=id_]
  
  
  setkey(new.traits.mvt, id_)
  setkey(turning, id_)
  setkey(mm2, id_)
  all.mvt.newTraits <- merge(new.traits.mvt,turning,by=c("id_"), all=T)
  mm3 <- merge(mm2,all.mvt.newTraits,by=c("id_"), all=T)
  # mm3$id <- mm3$id_
  # mm3$id_ <- NULL
  return(as.data.frame(mm3)[,-1])
}


#... Make a colour transparent with a transparency set by alpha
#..............................................................
#... From Ricardo Oliveros-Ramos (http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color)
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





summarize_populations2 <- function(traj.data, sum.data, write=FALSE, to.data, merged.data.folder, video.description.folder, video.description.file,total_frames = 125){
  
  # checks whether frames per second are specified
  if(!exists("fps") ) stop("frames per second not specified (fps)")
  # checks whether the sample volume is specified
  if(!exists("measured_volume") ) stop("measured volume not specified (measured_volume)")
  
  #... results object from video description file
  pop_output <- read.table(paste(to.data, video.description.folder, video.description.file, sep = ""), sep = "\t", header = TRUE)
  
  #... now add the population densities
  pop_count_table <- tapply(traj.data$Major,list(as.factor(traj.data$file),as.factor(traj.data$frame)),length)
  
  help_cnt_rep <- which(is.element(pop_output$file,dimnames(as.matrix(rowMeans(pop_count_table)))[[1]]))
  
  pop_output$indiv_per_frame <- 0
  pop_output$indiv_per_frame[help_cnt_rep] <- as.numeric(apply(pop_count_table,1,sum,na.rm=T))/total_frames
  
  pop_output$indiv_per_µl <- 0
  pop_output$indiv_per_µl[help_cnt_rep] <- as.numeric(apply(pop_count_table,1,sum,na.rm=T))/total_frames /measured_volume
  
  #... add the mean of total area per frame (bioarea by frame)
  pop_count_table2 <- tapply(traj.data$Area,list(as.factor(traj.data$file),as.factor(traj.data$frame)),sum)
  help_cnt_rep <- which(is.element(pop_output$file,dimnames(as.matrix(rowMeans(pop_count_table2)))[[1]]))
  pop_output$bioarea_per_frame <- 0
  pop_output$bioarea_per_frame[help_cnt_rep] <- as.numeric(apply(pop_count_table2,1,sum,na.rm=T))/total_frames 
  
  #... add the mean of total area per volume (bioarea by µl)
  pop_output$bioarea_per_µl <- 0
  pop_output$bioarea_per_µl[help_cnt_rep] <- as.numeric(apply(pop_count_table2,1,sum,na.rm=T))/total_frames/measured_volume
  
  # first get file from id
  sum.data$file <- sub("-.*$", "", sum.data$id )
  
  # get morphology
  pop_output$major_mean <- NA
  pop_output$major_sd <- NA
  pop_output$minor_mean <- NA
  pop_output$minor_sd <- NA
  
  pop_output$major_mean[help_cnt_rep] <- as.numeric(tapply(sum.data$mean_major,sum.data$file,mean,na.rm=T))
  pop_output$major_sd[help_cnt_rep] <- as.numeric(tapply(sum.data$mean_major,sum.data$file,sd,na.rm=T))
  pop_output$minor_mean[help_cnt_rep] <- as.numeric(tapply(sum.data$mean_minor,sum.data$file,mean,na.rm=T))
  pop_output$minor_sd[help_cnt_rep] <- as.numeric(tapply(sum.data$mean_minor,sum.data$file,sd,na.rm=T))
  
  # calculate gross speed
  sum.data$gross_speed <- sum.data$gross_disp/sum.data$duration 
  
  # get movement
  pop_output$gross_speed_mean <- NA
  pop_output$gross_speed_sd <- NA
  pop_output$net_speed_mean <- NA
  pop_output$net_speed_sd <- NA
  pop_output$sd_turning_mean <- NA
  
  # mean gross speed is not returned yet --> calc here
  sum.data$gross_speed <- sum.data$gross_disp/sum.data$duration
  
  pop_output$gross_speed_mean[help_cnt_rep] <- as.numeric(tapply(sum.data$gross_speed,sum.data$file,mean,na.rm=T))
  pop_output$gross_speed_sd[help_cnt_rep] <- as.numeric(tapply(sum.data$gross_speed,sum.data$file,sd,na.rm=T))
  
  pop_output$net_speed_mean[help_cnt_rep] <- as.numeric(tapply(sum.data$net_speed,sum.data$file,mean,na.rm=T))
  pop_output$net_speed_sd[help_cnt_rep] <- as.numeric(tapply(sum.data$net_speed,sum.data$file,sd,na.rm=T))
  
  pop_output$sd_turning_mean[help_cnt_rep] <- as.numeric(tapply(sum.data$sd_turning,sum.data$file,mean,na.rm=T))
  
  
  #output population summary data
  if (write==TRUE){save(pop_output, file = paste0(to.data, merged.data.folder,"Population_Data.RData"))}
  return(as.data.frame(pop_output))
}







#... Function to add a column for the density of each species from a morph_mvt object having a predicted_species column
#......................................................................................................................
#... sample_output is the video description file (line = a sample)
#... indiv_predicted is the data.frame with predicted species of individuals, their sample, and the number of frames on which each appears (N_frames)
#... species_names gives the names of the species (as specified in the predicted_species column)
#... total_frames is the total number of frames on each videos (frames per second x number of seconds)

species.density = function(sample_output,indiv_predicted,species_names,total_frames,mv = measured_volume)
{
  #... get the names of the samples where there are individuals
  samples = unique(indiv_predicted$file)
  
  #... create the matrix of species densities
  sp.dens = matrix(0,nrow(sample_output),length(species_names))
  colnames(sp.dens) = species_names
  
  for(i in 1:length(samples))
  {
    #... select the data for each sample
    indiv = subset(indiv_predicted,file == samples[i])
    
    #... get the species present in there
    spec = unique(indiv$predicted_species)
    for(j in 1:length(spec))
    {
      #... select the data for one species in the sample
      all.indiv.sp = subset(indiv,predicted_species == spec[j])
      
      #... calculate its density from the total number of frames on which each individuals of the species is present
      dens = sum(all.indiv.sp$N_frames)/total_frames/mv
      sp.dens[which(sample_output$file == as.character(samples[i])) ,which(species_names == spec[j])] = dens
    }
  }
  return(cbind(sample_output,sp.dens)) 
}
