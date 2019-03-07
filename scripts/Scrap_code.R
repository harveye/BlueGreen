### Plot effect of protist density in blue ecosysten on bacteria density in green ecosystem
test <- data %>%
  filter(day!=0 & Replicate %in% c("A","B") & Size == 13) %>% 
  select(bact.green,Rot,Spi,Ble,Pca,Col,Chi,Tet) 
  
  test2 <- reshape2:::melt(test,id.var=1)
  
  ggplot(data=test2,mapping=aes(y=bact.green,x=decostand(value,'log'),colour=variable)) + 
  geom_point()


######### Verify different calculations of log ratio of means

###TEST --> Manuel log ratio of the means (using SummarySE) and 
#ARPobservation::logRespRatio() functions lead to extremely similar results. 
#Bootstraped results are different by a small margin, but with same interpretation.
#The issue with manual --> CI can be hard to calculate manually 
#logRespRatio() offers a great alternative with no information on CI calculation = blackbox
#Bootstrap is intuitive and transparent both in LRR calculation and CI
#Issue with Bootstrap: results vary (marginally) depending on how zeros are delt with (zeros are not an issue when dealing with ratio of means)


## Bootstrap methods
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


#from logRespRatio function

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




##### Some tests on different ways to calculate ratios
test <- data.frame(treatment = sample(c(1:1000000),size=1000), control=sample(c(1:1000000),size=1000))
test

#Bootstrap 
library(boot)
set.seed(2)
yourFun <- function(x, i) {
  xSub <- x[i, ] #resample x
  LnRR <- log(xSub[, 1]/xSub[ ,2])
  return(mean(LnRR))
}

b <- boot(test, yourFun, R=999)
plot(b)
ci <-boot.ci(b,type="norm")
b$t0 # meanLRR
ci$normal[2] #lower ci
ci$normal[3] #upper ci

b$t0 #Boot

(LRRs_ind <- mean(log(test$treatment/test$control))) #mean of each log ratio

(LRR_mean_ratio <- log( mean(test$treatment)/mean(test$control) ) ) #log ratio OF the mean

(LRR_sum <- log(sum(test$treatment/test$control))) #log of the sum of ratio

(LRR <- log(mean(test$treatment/test$control))) #log of the mean ratio


####Generate data for LRR calculation (for bacteria in blue ecosystem)
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



display.brewer.all(n=7,colorblindFriendly = F)
c(brewer.pal(7,"Dark2"),"black")

library(colorspace)

pal <- choose_palette()

display.brewer.all(n=7)

display.brewer.all(n=8,type="qual", colorblindFriendly=T)

display.brewer.pal(n=8,"Set2")

