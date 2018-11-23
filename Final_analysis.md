---
title: "Analysis from Harvey et al., XXX"
author: "Eric"
date: '2018-09-21'
output: 
  html_document:
    highlight: tango
    theme: readable
    toc: true
    toc_float:
      collapsed: true
      smooth_scroll: false
    keep_md: true
---



# Analytical pipeline

The main objective of this project is to identify spatial feedbacks between two landscapes connected by the pulse exchange of resource (see Methods section in the main article). More specifically we are interested in effects of the main treatment (resource pulse) that are dependant on the position in the spatial network. Thus although we are interested in the main effect of our resource pulse treatment, our focus is on the interaction term between patch size (indicates position in the dendritic network) and resource pulse treatment. 

The analytical pipeline to achieve this objective can be synthesize as follow: 

- Effects from green ecosystem to blue ecosystem
  - On community composition
  - On aggregate community metrics (species richness, total density)
  
- Effects from blue ecosystem to green ecosystem
  - On aggregate bacteria metrics (total bacteria density)

Visual explorations from the `Overview.Rmd` already strongly suggests no important effects on aggregate metrics, but significant effects on community composition. For the manuscript, we will thus start by exploring the effect of resource pulse from the green ecosystem to community composition in the blue ecosystem, and then how resource pulse from the blue ecosystem affects bacteria density in the green ecosystem. Those results constitute the core findings presented in the manuscript. At the end we will use statistical modelling to formally test effects on aggregate metrics in the blue ecosystem. 


# Effect from Green to Blue ecosystem

## Ordination



We will use an RDA analysis with those three main components:

* Community matrix (C - Hellinger transformed)

```
##         Rot       Spi       Ble       Pca        Col       Chi       Tet
## 1 0.5860927 0.1684091 0.2783666 0.2835491 0.12630679 0.6277132 0.2454963
## 2 0.4466228 0.2398793 0.3229581 0.2508649 0.10696917 0.6991657 0.2747287
## 3 0.4643459 0.1297158 0.2269056 0.3227714 0.00000000 0.7324779 0.2745273
## 4 0.5782057 0.2695530 0.1790308 0.2037629 0.04272032 0.7118896 0.1040965
## 5 0.6048435 0.1613525 0.2151684 0.1673305 0.11986730 0.7207389 0.0000000
## 6 0.4776683 0.2891105 0.3146273 0.1957177 0.11406729 0.7239454 0.1176608
```
* Environmental matrix (E)

```
## # A tibble: 6 x 7
##   day.cont centrality Size  dist.outlet Treatment Treat.Size Treat.Size.D…
##      <dbl>      <dbl> <ord>       <dbl> <fct>     <fct>      <fct>        
## 1        7    0.00465 7.5             6 Connected Connected… Connected.7.…
## 2        7    0.00465 7.5             6 Connected Connected… Connected.7.…
## 3        7    0.00465 7.5             6 Connected Connected… Connected.7.…
## 4        7    0.00405 7.5             5 Connected Connected… Connected.7.…
## 5        7    0.00405 7.5             5 Connected Connected… Connected.7.…
## 6        7    0.00405 7.5             5 Connected Connected… Connected.7.…
```

So we ran the model:

```r
rda.mod <- rda(C ~ ., as.data.frame(E))
rda.mod
```

```
## Call: rda(formula = C ~ day.cont + centrality + Size + dist.outlet
## + Treatment + Treat.Size + Treat.Size.Day, data =
## as.data.frame(E))
## 
##               Inertia Proportion Rank
## Total         0.24320    1.00000     
## Constrained   0.07771    0.31953    7
## Unconstrained 0.16549    0.68047    7
## Inertia is variance 
## Some constraints were aliased because they were collinear (redundant)
## 
## Eigenvalues for constrained axes:
##    RDA1    RDA2    RDA3    RDA4    RDA5    RDA6    RDA7 
## 0.05551 0.01210 0.00319 0.00284 0.00202 0.00152 0.00053 
## 
## Eigenvalues for unconstrained axes:
##     PC1     PC2     PC3     PC4     PC5     PC6     PC7 
## 0.04076 0.03462 0.02853 0.02347 0.01737 0.01483 0.00591
```

The model show that our environmental matrix explains roughly 31% of the variance, which is not amazing but not so bad.

We then ran a permutation ANOVA on the terms: 

```r
anova(rda.mod,by="terms",permu=200)
```

```
## Permutation test for rda under reduced model
## Terms added sequentially (first to last)
## Permutation: free
## Number of permutations: 999
## 
## Model: rda(formula = C ~ day.cont + centrality + Size + dist.outlet + Treatment + Treat.Size + Treat.Size.Day, data = as.data.frame(E))
##                  Df Variance        F Pr(>F)    
## day.cont          1 0.036435 246.1496  0.001 ***
## centrality        1 0.009496  64.1510  0.001 ***
## Size              3 0.002291   5.1587  0.001 ***
## dist.outlet       1 0.000480   3.2415  0.014 *  
## Treatment         1 0.007661  51.7566  0.001 ***
## Treat.Size        3 0.000993   2.2362  0.011 *  
## Treat.Size.Day   23 0.020352   5.9779  0.001 ***
## Residual       1118 0.165488                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Here we can see that three variables are really standing out:

* Change in time
* Position in the landscape (centrality)
* Treatment (connected vs. isolated)

So it turns out that both the position in the landscape and resource pulse influence community position. Note that the interaction terms are also significant here indicating that the influence of the resource pulse on community composition vary as a function of the position in the network. The importance of these interaction terms in the model, however, remains fairly low compared to the main effects. 

Let's explore this model visually

![](Final_analysis_files/figure-html/ordination figs-1.png)<!-- -->

First looking at the predictors (arrows) one can see that higher patch volume and high centrality values co-varies, which make sense since larger patches are generally more connected. In direct opposition (180 degrees on the figure) higher distance to outlet values are found which also make sense since smaller patches are mainly found upstream and thus further from the outlet. Orthogonal to this gradient (90 degree) is the effect of time. The species scores are well distributed on the plot as a function of the different constraints. Each point on the figure is a microcosm with the two colors identifying them as being part of a connected or an isolated blue ecosystem. From the model we know that the effect of the treatment is significant and relatively important. Here we can see that there is indeed some level of difference between sites in connected and isolated ecosystems.  

Now in terms of interpretation, we can see that PCA and COL are more strongly associated with sites in connected ecosystem while Spi, Tet and Chi are associated with sites in isolated ecosystems. There are some interactions with other preditors though. For instance the figure suggest that PCA also became more abundant later in the experiment, while SPI was more abundant earlier suggesting some temporal community dynamics. Most interestingly, as mentionned, TET and CHI are associated with sites in isolated ecosystem, BUT also sites that tend to be more upstream in these ecosystems (smaller patch size more distant from the outlet). On the opposite ROT seems to be more associated to downstream sites (high centrality, high patch volume) independant of treatment or temporal aspect.

## Log-ratio 

Another useful way of representing those results is by using the log response ratio here defined has: 

**LRR = log(Abundance~connected~/Abundance~isolated~)**

with confidence interval: 

**95CI = LRR +- 1.96 x SE(LRR)**

In complementarity with the RDA model, let's look at the log response ratio of the different protist species to the main treatment. 



![](Final_analysis_files/figure-html/LRR figure 1-1.png)<!-- -->

This figure conveys very similar information to the RDA figure: PCA and COL tend to be more abundant in microcosms from connected ecosystems while CHI, TET and BLE tend to be more abundant in microcosms from isolated ecosystems. GREAT! 

Now let's see if we can still detect the interaction with position in the network using patch size (which according to the RDA is the same as using distance to outlet or centrality): 


![](Final_analysis_files/figure-html/LRR figure 2-1.png)<!-- -->

Again, the same information is conveyed. We do see the interaction with PCA and COL being more abundant in the upstream (smaller volume) of the connected ecosystem, and CHI and TET being more abundant in the upstream (smaller volume) patches of the isolated ecosystem. For CHI this is also true for intermediate patch sizes.

## Preliminary conclusion
The results presented above suggest that: 

* Resource pulse from the green ecosystem does affect the structre of the community in the blue ecosystem by selecting for and againts some species
* This effect depend on the position in the blue spatial network and tend to be stronger in smaller, more isolated patches (based on patch volume, centrality and distance to the outlet)

# Effect from Blue to Green ecosystem

Here the main question of interest is whether we can see imprint of what was going on in the blue ecosystem. Using the log ratio approach seems to be one of the best way to visualize those effects: 

## Log-ratio



![](Final_analysis_files/figure-html/LRR Green figure 1-1.png)<!-- -->


So there is a tendancy for bacteria density in the green ecosystem to be higher when connected to the blue ecosystem versus isolated controls. The strongest effect is observed on day 21. Let's see if we can detect any interactions, even weak, with the position in the blue ecosystem:

![](Final_analysis_files/figure-html/LRR Green figure 2-1.png)<!-- -->

Interestingly, the positive effect of being connected to the blue ecosystem is especially strong for the green patches that were connected with a small volume upstream site in the blue ecosystem. Here this figure is only showing day 21 which is the day where the effect of the treatment was strongest according to the previous figure, but the effect is also significant (LRR > 0) although smaller when averaging across all time periods. 

# Conclusion

Main conclusions:

* The effect of resource pulse on the blue ecosystem depends on the position in the spatial network
* This network mediated effect also feed back on the green ecosystem 
* Upstream (smaller volume, higher distance to outlet and lower centrality) sites are key to this spatial feedback as they are the most responsive to resource pulse and probably most dependant on those resource pulse. 

# Addendum: aggregate community metrics

Visual exploration from the `Overview.Rmd` already strongly suggests no effects on aggregate metrics. We will here use mixed effects modeling to directly test effects on bacteria and protist densities, and protist species richness. We will then see if those results can inform or not our main results presented above.  

* Effect of resource pulse from Green ecosystem on Bacteria density in Blue ecosystem


Let's print out the t Table

```r
# Run model
Mod1 <- nlme:::lme(log(bact.blue) ~ Treatment*Size*day.cont + bact.green, 
                   random = ~ 1|Replicate, data=X.bact,
                   method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod1)$tTable
```

```
##                                         Value    Std.Error  DF    t-value
## (Intercept)                      1.639611e+01 2.996754e-01 566 54.7129148
## TreatmentIsolated                2.042969e-02 4.158309e-01 566  0.0491298
## Size                             2.215569e-02 1.796103e-02 566  1.2335426
## day.cont                         4.950361e-01 1.043469e-01 566  4.7441390
## bact.green                      -1.850781e-09 1.304054e-09 566 -1.4192522
## TreatmentIsolated:Size          -9.301013e-03 2.534468e-02 566 -0.3669809
## TreatmentIsolated:day.cont       6.331961e-02 1.478129e-01 566  0.4283766
## Size:day.cont                   -1.212724e-02 6.372805e-03 566 -1.9029667
## TreatmentIsolated:Size:day.cont  3.153082e-03 8.998941e-03 566  0.3503837
##                                       p-value
## (Intercept)                     3.686880e-228
## TreatmentIsolated                9.608332e-01
## Size                             2.178854e-01
## day.cont                         2.654378e-06
## bact.green                       1.563760e-01
## TreatmentIsolated:Size           7.137704e-01
## TreatmentIsolated:day.cont       6.685399e-01
## Size:day.cont                    5.755201e-02
## TreatmentIsolated:Size:day.cont  7.261811e-01
```
_please note that because Size is ordered, the intercept here represents the MEAN factor level and not the baseline level (Size=7.5) (L: linear, Q: Quadratic, C: Cubic)_

We can see that three effects come as statistically significant: 
- Size (only on day 29): decline of bacteria density with increasing patch size
- Day: average increases of bacteria density over experimental time

![](Final_analysis_files/figure-html/nlme_blue_bact_figure-1.png)<!-- -->

* Effect of resource pulse from Green ecosystem on protist density in Blue ecosystem

```r
# Run model
Mod2 <- nlme:::lme(prot.ab ~ Treatment*Size*day.cont, 
                   random = ~ 1|Replicate, data=X1,
                   method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod2)$tTable
```

```
##                                       Value  Std.Error   DF     t-value
## (Intercept)                     3485.647088 110.006558 1141  31.6858118
## TreatmentIsolated                223.725766 147.270324 1141   1.5191504
## Size                             -30.301551   6.256787 1141  -4.8429892
## day.cont                        -988.410654  36.975119 1141 -26.7317779
## TreatmentIsolated:Size             7.071487   8.848422 1141   0.7991806
## TreatmentIsolated:day.cont       -26.899427  52.290714 1141  -0.5144207
## Size:day.cont                      7.912443   2.221571 1141   3.5616433
## TreatmentIsolated:Size:day.cont   -2.756024   3.141776 1141  -0.8772188
##                                       p-value
## (Intercept)                     1.383199e-158
## TreatmentIsolated                1.290017e-01
## Size                             1.454989e-06
## day.cont                        1.237515e-122
## TreatmentIsolated:Size           4.243520e-01
## TreatmentIsolated:day.cont       6.070575e-01
## Size:day.cont                    3.837329e-04
## TreatmentIsolated:Size:day.cont  3.805525e-01
```

We can see that three effects come as statistically significant: 
- Size: Decline of protist abundance with increasing patch size
  - The effect change with time
- Day: Strong decline of protist density over experimental time

![](Final_analysis_files/figure-html/nlme_prot_dens_figure-1.png)<!-- -->![](Final_analysis_files/figure-html/nlme_prot_dens_figure-2.png)<!-- -->

* Effect of resource pulse from Green ecosystem on protist richness in Blue ecosystem


```r
# Run model
Mod3 <- nlme:::lme(prot.rich ~ Treatment*Size*day.cont, 
                   random = ~ 1|Replicate, data=X1,
                   method="ML",control=lmeControl(optimMethod="BFGS",maxIter=100,opt="optim"))
summary(Mod3)$tTable
```

```
##                                         Value   Std.Error   DF
## (Intercept)                     11.6990005776 0.435453433 1141
## TreatmentIsolated               -1.0877085152 0.615815015 1141
## Size                             0.0427391823 0.026162899 1141
## day.cont                        -1.9314745005 0.154612501 1141
## TreatmentIsolated:Size           0.0088842752 0.036999927 1141
## TreatmentIsolated:day.cont       0.3817514883 0.218655096 1141
## Size:day.cont                   -0.0274476112 0.009289561 1141
## TreatmentIsolated:Size:day.cont  0.0002421182 0.013137423 1141
##                                      t-value       p-value
## (Intercept)                      26.86624952 1.351208e-123
## TreatmentIsolated                -1.76629100  7.761440e-02
## Size                              1.63357975  1.026230e-01
## day.cont                        -12.49235662  1.170148e-33
## TreatmentIsolated:Size            0.24011602  8.102834e-01
## TreatmentIsolated:day.cont        1.74590712  8.109623e-02
## Size:day.cont                    -2.95467257  3.194251e-03
## TreatmentIsolated:Size:day.cont   0.01842966  9.852993e-01
```

We can see that two effects come as statistically significant: 
- Size: Decline of protist richness over patch size but only significant at day 29
- Day: average decline in protist richness over experimental time

![](Final_analysis_files/figure-html/nlme_prot_rich_figure-1.png)<!-- -->

## Conclusion on aggregate metrics

We did not detect any effects of our main treatments in the Blue ecosystem on bacteria and protist densities or on protist richness. We did detect effects of patch size and changes over time. Those effects of patch size on aggegate metrics were thus independant of being connected or not and represent constraints of the dendritic spatial networks on richnes and density. Those patterns found here, especially on species richness, have already been analysed comprehensibly and published elsewhere (please see Harvey et al., Proc.B., in press)  
