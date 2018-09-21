#..........................................#
#...........Blue-Green project.............#
#..........................................#

#..........................................................................................................................................#
#... Collaborators: Eric Harvey, Isabelle Gounand, Emanuel Fronhofer, Florian Altermatt                                                    #
#... Author of the script: Eric Harvey                                                                                                     #
#... Contact: eric.harvey@eawag.ch                                                                                           #                                                                       #
#..........................................................................................................................................#


# This script clean the experimental data to generate the data for the analysis. 

#Experimental design: 
# Blue - isolated (A,B,C,D) or connected (A,B,C,D) --> 2 T x 4 rep x 36 = 288 microcosms (rows) per sampling day
# Green - isolated (A,B) or connected (A,B,C,D) --> ( (1 T x 2 rep) + (1T x 4 rep) ) x 36 = 216 microcosms (rows) per sampling day
# Protist - Blue only --> 288 + 14 monocultures = 302 microcosms per sampling day
# Bacteria - Blue and Green + 20 monocultures = 288 + 216 + 20 = 524 rows (Bacteria were measured only for replicates A and B) per sampling day

##################
#Clear any variables from the R environment 
##################
rm(list=ls())

##################
#Directory paths 
##################
to.data <- "./data/"
to.script <- "./scripts/"
to.output <- "./output/"
to.figs <- "./figs/"
to.R <- "./output/"

##################
#Load packages
##################
library(tidyverse)
library(vegan)

##################
#Load data in R environment
##################
Prot.b <-  read_csv(paste0(to.data,"BG_protist_data_(20170307).csv")) #Protist
Cyto <- read_tsv(paste0(to.data,"FULL_Cyto_BlueGreen(20160623).txt")) # Bacteria


#########################################################################
################# DATA STRUCTURE
#########################################################################

#Create a day variable (experimental day)
Prot.b$day <- ifelse(Prot.b$date=="16-05-02",0, ifelse(Prot.b$date=="5/9/2016",7, ifelse(Prot.b$date=="5/17/2016",15, 
        ifelse(Prot.b$date=="16-05-23",21, ifelse(Prot.b$date=="16-05-31",29, NA)))))
Cyto$day <- ifelse(Cyto$Date=="20160502",0, ifelse(Cyto$Date=="20160509",7, ifelse(Cyto$Date=="20160517",15,
            ifelse(Cyto$Date=="20160523",21, ifelse(Cyto$Date=="20160531",29, NA)))))


#clean protist data

species = c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other") #place species names in a vector


Prot.b[,species] = Prot.b[,species]*1000 #Convert species abundance per ul into abundance per ML

cl.prot <- Prot.b %>% 
  mutate(prot.rich = specnumber(Prot.b[,species]),
         prot.ab = rowSums(Prot.b[,species])) %>% 
  filter(Replicate %in% c("A","B","C","D")) %>% #remove monoculture from data
  select(day,Treatment,Replicate,Size,Rot,Spi,Ble,Pca,Col,Chi,Tet,Other,prot.rich,prot.ab) %>% #select only protist density info
  arrange(day,Treatment,Replicate) #order rows

#Clean cyto data to match protist data 

cl.cyto <- Cyto %>% 
  filter(Replicate %in% c("A","B","C","D") & Landscape=="Blue") %>% #Instead of landscape blue and green we will have bacteria density in blue and green as different columns
  mutate(bact.blue = Cyto$Count[Cyto$Landscape=="Blue"],
         bact.green = Cyto$Count[Cyto$Landscape=="Green"]) %>%  
  mutate(bact.blue = (bact.blue/Volume.uL*1000) * 1000,
         bact.green = (bact.green/Volume.uL*1000) * 1000) %>%  #convert from dens/50ul to dens/mL and then multiply by the cytometry dilution factor (1000)
  select(day,Treatment,Replicate,bact.blue,bact.green)  %>%        
  arrange(day,Treatment,Replicate)

#Merge Protist and Bacteria data together 
merged.data = bind_cols(cl.prot,cl.cyto[,4:5])

#Add network metrics to Protist and Bacteria data
net_metric <- read_rds(paste0(to.output,"net_metric.RDS"))
net_metric <- as.tibble(net_metric) #Transform into Tibble 
net_metric <- net_metric[rep(seq_len(nrow(net_metric)),times=2*4),] #Repeat metric networks for each sampling day
net_metric$day <- c(rep(7,288),rep(15,288),rep(21,288),rep(29,288))  #Create a day var (excluding 0)
dummy_d0 <- data.frame(day=rep(0,16),centrality=rep(NA,16)) #Generate a dummy day 0 dataset (see below)

net_metric <- bind_rows(dummy_d0,net_metric) #now bind the two datasets 

#Merge all clean data together
cl.data <- bind_cols(merged.data,net_metric[,c(2,5:7)])

#Save data
saveRDS(cl.data,paste0(to.output,"cl.data.RDS"))


