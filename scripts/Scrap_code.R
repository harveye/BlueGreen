### Plot effect of protist density in blue ecosysten on bacteria density in green ecosystem
test <- data %>%
  filter(day!=0 & Replicate %in% c("A","B") & Size == 13) %>% 
  select(bact.green,Rot,Spi,Ble,Pca,Col,Chi,Tet) 
  
  test2 <- reshape2:::melt(test,id.var=1)
  
  ggplot(data=test2,mapping=aes(y=bact.green,x=decostand(value,'log'),colour=variable)) + 
  geom_point()


######### Verify calculation of log ratio of means

##Generate data for LRR calculation

t <- c("7","15","21","29") 
SP <- c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other")
S <- c("7.5","13.0","22.5","45.0")

**95CI = LRR +- 1.96 x SE(LRR)**

LRR_treatment <- data.frame(species=character(0),LRR=numeric(0),CI_upper=numeric(0),CI_lower=numeric(0))

for(i in 1:length(SP)) { 
  
#
  D <- data %>% 
    filter(day !=0) %>% 
    select(SP[i],day, Treatment, Replicate, Size, centrality, degree, dist.outlet) %>% 
    mutate(Treatment = as.factor(Treatment),
           Replicate = as.factor(Replicate),
           Size = as.factor(Size))  %>% 
    Rmisc::summarySE(.,measurevar = SP[i], groupvars = c("Treatment"))
  
  LRR <- log(D[1,3]) - log(D[2,3])
  SE_LRR <- sqrt( (D$se[1]^2/D[1,3]^2) + (D$se[2]^2/D[2,3]^2) ) #From Hedge et al., 19999
  CI_upper <- LRR + 1.96*SE_LRR
  CI_lower <- LRR - 1.96*SE_LRR
  #CI_both <- LRR + c(-1, 1) * qnorm(1 - (1 - 0.95)/2) * SE_LRR #from ARPobservation:::logRespRatio
  species <- SP[i]
  
  bind.out <- data.frame(species,LRR,CI_upper,CI_lower)
  #colnames(bind.out) <- c("time","size","LRR","species","CI_upper","CI_lower")
  
  LRR_treatment <- rbind.data.frame(LRR_treatment,bind.out)
  
  }



LRR_treatment




###TEST --> Manuel log ratio of the means (using SummarySE) and 
#ARPobservation::logRespRatio() functions lead to extremely similar results. 
#Bootstraped results are different by a small margin, but with same interpretation.
#The issue with manual --> CI can be hard to calculate manually 
#logRespRatio() offers a great alternative with no information on CI calculation = blackbox
#Bootstrap is intuitive and transparent both in LRR calculation and CI
#Issue with Bootstrap: results vary (marginally) depending on how zeros are delt with (zeros are not an issue when dealing with ratio of means)


#Bootstrap
dat.test <- data.frame(treatment = D$Chi[D$Treatment=="Connected"], control = D$Chi[D$Treatment=="Isolated"] )
#dat.test <- dat.test[-which(dat.test$control == 0 | dat.test$treatment == 0),]
for(i in 1:nrow(dat.test)){
  if(dat.test[i,1] == 0) dat.test[i,1] = summary(dat.test[,1])[2]
  if(dat.test[i,2] == 0) dat.test[i,2] = summary(dat.test[,2])[2]
 }

library(boot)
set.seed(2)
yourFun <- function(x, i) {
  xSub <- x[i, ] #resample x
  LnRR <- log(xSub[, 1]/xSub[ ,2])
  return(mean(LnRR))
}

b <- boot(dat.test[-453,], yourFun, R=999)
plot(b)
ci <-boot.ci(b,type="norm")
b$t0 # meanLRR
ci$normal[2] #lower ci
ci$normal[3] #upper ci


#from logRespRatio

#with zeros
lala = with(D,
            ARPobservation::logRespRatio(observations=Chi,phase=Treatment,base_level="Isolated")
)

$lRR
[1] -0.01376953

$V_lRR
[1] 0.003749426

$CI
[1] -0.1337831  0.1062441

lala$lRR
lala$CI[1]






####Create vectors with variables of interest 
t <- c("7","15","21","29") 
SP <- c("Rot","Spi","Ble","Pca","Col","Chi","Tet","Other")
S <- c("7.5","13.0","22.5","45.0")

####Replace zeros by the same very small value (to avoid INF and -INF in log ratio)
min = 0.001
D[,SP] <- lapply(D[,SP],function(x) ifelse(x==0,min,x))

####Calculate LRR 

output <- data.frame(time=character(0),size=character(0),LRR=numeric(0),SP=character(0))#[1:(144*length(SP)),]

for(i in 1:length(t)){
  for(j in 1:length(SP)) { 
    
    LRR <- log( D[which(D$Treatment=="Connected" & D$day==t[i]),SP[j]] / D[which(D$Treatment=="Isolated" & D$day==t[i]),SP[j]] )
    
    species <- rep(SP[j],144)
    
    time <- rep(t[i],144)
    
    size <- D$Size[which(D$Treatment=="Connected" & D$day==t[i])]
    
    lala <- cbind(time,size,LRR,species)
    colnames(lala) <- c("time","size","LRR","species")
    
    output <- rbind.data.frame(output,lala)
    
    
  }
}

library(metafor)
?escalc()





test <- data.frame(treatment = rnorm(10,10,5), control=rnorm(10,10,5))
test

LRRs_ind <- mean(log(test$treatment/test$control))

LRR_mean_ratio <- log( mean(test$treatment)/mean(test$control) )

LRR_sum <- log(sum(test$treatment/test$control))

LRR <- log(mean(test$treatment/test$control))

LRR_mean_ratio
LRR_ind
LRR_sum


LRR_bact_green <- data.frame(species="bacteria",LRR=0.2313328,CI_upper=0.3173062,CI_lower=0.1453594)

LRR_bact_green <- cbind("Green.Bact",0.2313328,0.3173062,0.1453594)
LRR_Res_protist 

rbind(LRR_Res_protist,LRR_bact_green)





####Generate data for LRR calculation
D_Res_blue_bact <- data %>% 
  filter(day !=0 & Replicate %in% c("A","B")) %>% 
  Rmisc::summarySE(.,measurevar = "bact.blue", groupvars = c("Treatment"))

#Calculate LRR
species <- "bact.blue"
LRR <- log(D_Res_blue_bact[1,3]) - log(D_Res_blue_bact[2,3]) 
SE_LRR <- sqrt( (D_Res_blue_bact$se[1]^2/D_Res_blue_bact[1,3]^2) + (D_Res_blue_bact$se[2]^2/D_Res_blue_bact[2,3]^2) ) #From Hedge et al., 1999
CI_upper <- LRR + 1.96*SE_LRR
CI_lower <- LRR - 1.96*SE_LRR

bind.out_Res_blue_bact <- data.frame(species,LRR,CI_upper,CI_lower)

# Plot 
ggplot(bind.out_Res_blue_bact,mapping=aes(x=species,y=LRR)) +
  ylab("Effect size (log ratio of means)") + xlab("Bacteria in blue ecosystem") +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), width=0) +
  geom_point() +
  geom_hline(yintercept=0,lty=2) +
  annotate("text", x=0.7, y=0.7,label="log(R+/R-) > 0",size=2.5,colour="darkgreen") + 
  annotate("text", x=0.7, y=-1.5,label="log(R+/R-) < 0",size=2.5,colour="darkred") + theme_classic()




