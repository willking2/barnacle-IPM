### Will King
### Ch 3 - population model for B. glandula, IPM
### data from field survey of barnacles after 1 year (Jun 2017 - May 2018)

# ~ Data Prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- load data and packages ----

library(lme4)
library(car)
library(lmerTest)
library(MuMIn)
library(merTools)
library(scales)
library(nlme)
library(doBy)
library(boot)
library(fitdistrplus)

#source('MatrixImage.R') # from Ellner, Childs, and Rees 2016, IPM book
# load probability of reproduction model that you made in ch 2. you made a separate version of it to use in ch 3, as follows:
source('ch3_probability-reproduction_20190102.R')

photos <- read.csv('ch3_field-photos.csv', header = T)
rec <- read.csv('ch3_recruits_sizedist.csv', header = T)
totrec <- read.csv('ch3_recruits_total.csv', header = T)




# ---- define domEig ----

# defining function, domEig, to compute just the dominant eigenvalue and eigenvector of the IPM matrix by iteration
# for a size-quality model, matrix is huge; using eigen() will crash R
# based on IPM book pg. 161, code from IPM book code, chapter 6, domEig.R

domEig=function(A,tol=1e-8) {
  qmax=10*tol; lam=1; x=rep(1,nrow(A));   
  while(qmax>tol) {
    x1=A%*%x;
    qmax=sum(abs(x1-lam*x));  
    lam=sum(x1); 
    x=x1/lam; 
  } 
  # Having found w (within tol), get lambda 
  x1 = A%*%x; lam=sum(x1); x=x1/lam;   
  return(list(lambda=lam,w=x/sum(x)))
}  	

# domEig <- function(A, tol = 1e-8){
#   qmax <- 10*tol; lam <- 1;
#   x <- rep(1, nrow(A))/nrow(A);
#   while(qmax > tol){
#     x1 <- A%*%x;
#     qmax <- sum(abs(x1-lam*x));
#     lam <- sum(x1);
#     x <- x1/lam;
#   }
#   return(list(lambda = lam, w = x/sum(x)))
# }



# ---- determine recruits size distribution ----

summary(rec)
rec$z <- rec$operc_length_mm
rec$z <- round(rec$z, 1)

h <- hist(rec$z
          , freq = F
          , breaks = seq( min(rec$z), max(rec$z), 0.1)
)

# rec_sizedist <- h$density * 0.1 # size distribution of recruits (probability)
# rec_totalnew <- mean(totrec$number_of_recruits) # total number of recruits
# howmany <- round(rec_totalnew * rec_sizedist, 0) # how many of each size recruited
# 
# rec_bins <- h$breaks[h$breaks != max(h$breaks)] # size bins of new recruits, excluding final bin
# rec_z1 <- data.frame(z1 = rep(NA, sum(howmany))) # empty data frame for z1 recruits
# rec_z1$z1 <- rep(rec_bins, howmany)
# 
# hist(rec_z1$z1
#      , freq = F
#      , breaks = seq( min(rec_z1$z1), max(rec_z1$z1), 0.1)
# ) # checked to see that the calculated distribution matches empirical data 11/14/2018

### fit candidate distributions

mean(rec$z)
sd(rec$z)
median(rec$z) # mean > median, slightly right skewed

recruit.x <- seq(0, 3, length = 100)


## fit Gaussian distribution

normalfit <- fitdist(data = rec$z
                     , distr = 'norm'
                     , method = 'mle'
)
recruit.Gaussian <- dnorm(recruit.x
                          , m = mean(rec$z) # same as normalfit$estimate[1]
                          , sd = sd(rec$z) # same as normalfit$estimate[2]
)
lines(recruit.x
      , recruit.Gaussian
)

## fit gamma distribution

gammafit <- fitdist(data = rec$z
                    , distr = 'gamma'
                    , method = 'mle'
)
recruit.Gamma <- dgamma(recruit.x
                        , shape = gammafit$estimate[1] # k
                        , rate = gammafit$estimate[2]
)
lines(recruit.x
      , recruit.Gamma
      , col = 'blue'
)

## comparing fits
plot(normalfit)
plot(gammafit) # looks maybe slightly better
c(normalfit$aic, gammafit$aic) # favors normal
c(normalfit$bic, gammafit$bic) # favors normal
# based on AIC and BIC, go w/ Gaussian distribution


# ---- adjust photos dataset ----

### remove unwanted rows ###
photos2 <- photos[photos$died_yn != '', ] # without survival data
photos3 <- photos2[photos2$exclude != 'exclude', ] # rows marked exclude


### calculate percent touch ###
photos3$touch_pct <- round(photos3$touch_raw/360, digits = 2)


### add column for survival ###
# mortality: Y = died, N = lived
# survival: 0 = died, 1 = lived
photos3$Surv <- ''
photos3$Surv[photos3$died_yn == 'N'] <- 1
photos3$Surv[photos3$died_yn == 'Y'] <- 0

# ---- make dataset ----

### associated information
dat <- as.data.frame(photos3[ , c('site'
                                  , 'subsite'
                                  , 'density_trt'
                                  , 'individual_id'
                                  , 'touch_pct'
                                  , 'touch_spp'
                                  , 'day_elapsed'
)
]
)

### initial size, survival, final size
dat$z <- photos3$operc_length_initial_mm
dat$Surv <- photos3$Surv
dat$z1 <- photos3$operc_length_final_mm

### reproduction
# s(z)*Pb(z)*b(z)*Pr*C0(z')
# see surval * from ch 2 * from Strathmann * unknown * ch 3 data but not dependent on z
# see below 


# ---- adjust data: subsite for size, color by touch, surv to factor ----

## remove data for Balanus touching mostly Chthamalus 
#dat <- dat[dat$touch_spp != 'C', ]


## response surface
plot(touch_pct ~ z
     , data = dat
     , las = 1
     , xlab = 'Body size (mm), year t'
     , ylab = 'Crowding'
     , xlim = c(0, 7)
     , type = 'n'
)

rect(xleft = -0.27
     , ybottom = 0.8
     , xright = 7.28
     , ytop = 1.25
     , col = alpha('black', 0.25)
     , border = NA
)
rect(xleft = 5.5
     , ybottom = -0.25
     , xright = 7.28
     , ytop = 0.8
     , col = alpha('black', 0.25)
     , border = NA
)


# abline(h = 0.8, lty = 2)
# abline(v = 5.5, lty = 2)

## color points by survived or no
# lived
points(touch_pct ~ z
       , data = dat[dat$Surv == 1, ]
       , col = alpha('aquamarine4', 0.55)
       , pch = 19
)
# died
points(touch_pct ~ z
       , data = dat[dat$Surv == 0, ]
       , col = alpha('salmon', 0.55)
       , pch = 17
)



# limit data to initial operc size =< 5.5 and touch =< 0.8 
dat <- dat[dat$z <= 5.5 & dat$touch_pct <= 0.8, ]


# make column for colors, based on touch
colPal <- colorRampPalette(c('gray80','black'))
# darker means more touch (lightest gray is touch = 0, black is touch = 1)
dat$colors <- colPal(100)[as.numeric(cut(dat$touch_pct,breaks = 100))]

# keep a numerical version of survival for plotting purposes
dat$Surv.numeric <- dat$Surv
# change survival from character to factor
dat$Surv <- as.factor(dat$Surv)




# ---- calculate growth rate ----
dat$growthrate <- (dat$z1 - dat$z)/dat$day_elapsed

# ---- explore data ----

### general size structure

### adults (z1) and recruits together
hist(c(rec$z, dat$z1), freq = F)

### adults only
hist(dat$z)
hist(dat$z1, freq = F)

### growth
plot(z1 ~ z
     , data = dat
     , xlim = c(0, 8)
     , ylim = c(0, 8)
)
abline(a = 0, b = 1, lty = 2)
abline(lm(z1 ~ z, data = dat))

# color points by touch
points(z1 ~ z
       , data = dat
       , pch = 19
       , col = dat$colors
)


# spp competition?
points(z1 ~ z
       , data = dat[dat$touch_spp == '', ]
       , col = 'black'
)
points(z1 ~ z
       , data = dat[dat$touch_spp == 'B', ]
       , col = 'green'
)
points(z1 ~ z
       , data = dat[dat$touch_spp == 'C', ]
       , col = 'orange'
)

# growth ~ size
plot((z1-z)~ z
     , data = dat
) # no evidence of hump shape, just linear. ok to NOT ln transform size data... see Ellner book pg 20, box 2.2. 


# growth ~ touch
plot((z1-z)~ touch_pct
     , data = dat
)
# points((z1-z)~ touch_pct
#        , data = dat[dat$touch_spp == '', ]
#        , col = 'black'
# )
# points((z1-z)~ touch_pct
#        , data = dat[dat$touch_spp == 'B', ]
#        , col = 'green'
# )
# points((z1-z)~ touch_pct
#        , data = dat[dat$touch_spp == 'C', ]
#        , col = 'orange'
# )




### survival
plot(as.numeric(Surv) ~ z
     , data = dat
)
# color points by touch
points(as.numeric(Surv) ~ z
       , data = dat
       , pch = 19
       , col = alpha(dat$colors, 0.9)
)



plot(Surv ~ touch_pct
     , data = dat
)

# ---- field observed extremes in crowding distributions (supplement) ----
# low crowding dist observed in field: FH2, 5/26/2017
hist(photos$touch_raw[photos$subsite == 'FH2' & photos$density_trt == 'H']/360
     , freq = F
     , xlim = c(0, 1)
     , main = 'FH2'
     , xlab = 'crowding'
)

#high crowding dist observed in field: CP1, 6/12/2017
hist(photos$touch_raw[photos$subsite == 'CP1' & photos$density_trt == 'H']/360
     , freq = F
     , xlim = c(0, 1)
     , main = 'CP1'
     , xlab = 'crowding'
)

# ---- determine field adult and recruit touch distribution for iteration ----

h_t <- hist(dat$touch_pct
            , freq = F
)


### fit candidate distributions

touchy.x <- seq(0, 0.8, length = 100)

## fit Beta distribution

betafit_t <- fitdist(data = dat$touch_pct
                     , distr = 'beta'
                     , method = 'mle'
                     , lower = c(0,0)
                     , start = list(shape1 = 1, shape2 = 1)
)
touchy.beta <- dbeta(touchy.x
                     , shape1 = betafit_t$estimate[1]
                     , shape2 = betafit_t$estimate[2]
)
lines(touchy.x
      , touchy.beta
      , col = 'green'
)

## fit Gaussian distribution

normalfit_t <- fitdist(data = dat$touch_pct
                       , distr = 'norm'
                       , method = 'mle'
)
touchy.Gaussian <- dnorm(touchy.x
                         , m = mean(dat$touch_pct) # same as normalfit_t$estimate[1]
                         , sd = sd(dat$touch_pct) # same as normalfit_t$estimate[2]
)
lines(touchy.x
      , touchy.Gaussian
)

## fit gamma distribution

gammafit_t <- fitdist(data = dat$touch_pct
                      , distr = 'gamma'
                      , method = 'mle'
                      , lower = c(0,0)
                      , start = list(scale = 1, shape = 1)
)
touchy.Gamma <- dgamma(touchy.x
                       , shape = gammafit_t$estimate[1] # k
                       , rate = gammafit_t$estimate[2]
)
lines(touchy.x
      , touchy.Gamma
      , col = 'blue'
)

## fit weibull distribution

weibullfit_t <- fitdist(data = dat$touch_pct
                        , distr = 'weibull'
                        , method = 'mle'
                        , lower = c(0,0)
                        , start = list(scale = 1, shape = 1)
)
touchy.Weibull <- dweibull(touchy.x
                           , shape = weibullfit_t$estimate[1] # k
                           , scale = weibullfit_t$estimate[2]
)
lines(touchy.x
      , touchy.Weibull
      , col = 'red'
)

## comparing fits
plot(betafit_t)
plot(gammafit) # looks maybe slightly better
c(normalfit_t$aic, gammafit_t$aic, weibullfit_t$aic, betafit_t$aic) # favors beta
c(normalfit_t$bic, gammafit_t$bic, weibullfit_t$bic, betafit_t$bic) # favors beta
# based on AIC and BIC, go w/ beta distribution

dev.off()


# ---- determine field initial recruit touch distribution ----

# consider recruits as RT3 HLM only

h_t.r <- hist(dat$touch_pct[dat$subsite == 'RT3']
              , freq = F
)

### fit candidate distributions

touchy.x.r <- seq(0, 0.8, length = 100)

## fit Beta distribution

betafit_t.r <- fitdist(data = dat$touch_pct[dat$subsite == 'RT3']
                     , distr = 'beta'
                     , method = 'mle'
                     , lower = c(0,0)
                     , start = list(shape1 = 1, shape2 = 1)
)
touchy.beta.r <- dbeta(touchy.x.r
                     , shape1 = betafit_t.r$estimate[1]
                     , shape2 = betafit_t.r$estimate[2]
)
lines(touchy.x.r
      , touchy.beta.r
      , col = 'green'
)

# ---- set theoretical touch distributions for IPMsizetouch treatments ----



crowdy.x <- seq(0, 1, length = 100)

# # field observed touch distribution
# hist(dat$touch_pct
#      , freq = F
#      , xlim = c(0, 1)
# )
# lines(crowdy.x
#       , dbeta(crowdy.x
#               , shape1 = betafit_t$estimate[1] # 1
#               , shape2 = betafit_t$estimate[2] # 3.61
#       )
#       , col = 'green'
# )

# mostly low crowd
plot(crowdy.x
      , dbeta(crowdy.x
              , shape1 = 2
              , shape2 = 50
      )
      , col = 'blue'
     , type = 'l'
     , ylab = 'density'
     , xlab = 'crowding'
)

# mostly high crowd
lines(crowdy.x
      , dbeta(crowdy.x
              , shape1 = 50
              , shape2 = 2
      )
      , col = 'red'
)

# # mostly medium crowd
# lines(crowdy.x
#       , dbeta(crowdy.x
#               , shape1 = 300
#               , shape2 = 300
#       )
#       , col = 'purple'
# )

# uniform
lines(crowdy.x
      , dbeta(crowdy.x
              , shape1 = 1
              , shape2 = 1
      )
      , col = 'orange'
)

# # field observed adult touch distribution
# lines(crowdy.x
#       , dbeta(crowdy.x
#               , shape1 = betafit_t$estimate[1] # 1
#               , shape2 = betafit_t$estimate[2] # 3.61
#       )
#       , col = 'green'
# )

### quantifying amount of area under the curve for crowding scenarios

crowdy.uniform <- function(x)dbeta(x, shape1 = 1, shape2 = 1)
integrate(function(x)crowdy.uniform(x), 0.2, 0.3)
integrate(function(x)crowdy.uniform(x), 0.1, 0.2)

crowdy.low <- function(x)dbeta(x, shape1 = 2, shape2 = 50)
integrate(function(x)crowdy.low(x), 0, 0.1)
integrate(function(x)crowdy.low(x), 0, 0.2)

crowdy.high <- function(x)dbeta(x, shape1 = 50, shape2 = 2)
integrate(function(x) crowdy.high(x), 0.8, 1)
integrate(function(x) crowdy.high(x), 0.9, 1)


# ~ Regression analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- analyze: growth, z1 vs. z ----

### determine random structure

# no random intercept or slope
growth.1 <- lm(z1 ~ z * touch_pct
               , data = dat
)

# random intercept for subsite only
growth.2 <- lmer(z1 ~ z * touch_pct
                 + (1 | subsite)
                 , data = dat
)

# random intercept and slope for initial size
growth.3 <- lmer(z1 ~ z * touch_pct
                 + (1 + z | subsite)
                 , data = dat
)

# random intercept and slope for touch
growth.4 <- lmer(z1 ~ z * touch_pct
                 + (1 + touch_pct | subsite)
                 , data = dat
)

# random intercept for subsite nested in site
growth.5 <- lmer(z1 ~ z * touch_pct
                 + (1 | site/subsite)
                 , data = dat
)

# random nested interecept and slope for size 
growth.6 <- lmer(z1 ~ z * touch_pct
                 + (1 + z | site/subsite)
                 , data = dat
)

# random nested interecept and slope for touch
growth.7 <- lmer(z1 ~ z * touch_pct
                 + (1 + touch_pct | site/subsite)
                 , data = dat
)

# compare random structure models
AIC(growth.1
    , growth.2
    , growth.3
    , growth.4
    , growth.5
    , growth.6
    , growth.7
) # prefer growth.2
BIC(growth.1
    , growth.2
    , growth.3
    , growth.4
    , growth.5
    , growth.6
    , growth.7
) # prefer growth.2


### determine fixed structure
summary(growth.2) # R^2
Anova(growth.2)
anova(growth.2)

### rsquared using MuMIn::r.squaredGLMM
MuMIn::r.squaredGLMM(growth.2)

### model validation

# resid vs. fitted ### note non-constant variance ###
plot(x = fitted(growth.2),
     y = resid(growth.2),
     xlab = 'fitted values',
     ylab = 'residuals')
lines(lowess(fitted(growth.2), resid(growth.2)), col="blue", lwd=2)


# resid vs. explanatory
plot(x = dat$z[is.na(dat$z1) == F],
     y = resid(growth.2, type = 'deviance')
)

plot(x = dat$touch_pct[is.na(dat$z1) == F],
     y = resid(growth.2, type = 'deviance')
)

qqnorm(resid(growth.2))
qqline(resid(growth.2))

# 
# ### [archive] try to model the non-constant variance ###
# 
# # hetergeneity of variance; see Zurr pg. 74 and 90
# # variance seems to be greater at smaller body sizes; model fixed variance structure where variance is inversely proportional to body size
# 
# dat.growth <- dat[is.na(dat$z1) == F, ]
# z.growth <- dat.growth$z
# 
# growth.1_lm <- gls(z1 ~ z * touch_pct
#                    , data = dat.growth
# )
# vfix <- varFixed(~z.growth)
# 
# growth.1_gls <- gls(z1 ~ z * touch_pct
#                     , weights = vfix
#                     , data = dat.growth
#                     
# )
# 
# AIC(growth.1_lm, growth.1_gls) # prefers gls
# BIC(growth.1_lm, growth.1_gls) # prefers gls
# 
# ## determine fixed structure of gls
# Anova(growth.1_gls) # effect of size, and size*touch
# 
# ## model validation
# 
# # resid vs. fitted
# plot(x = fitted(growth.1_gls),
#      y = resid(growth.1_gls
#                , type = 'normalized'),
#      xlab = 'fitted values',
#      ylab = 'residuals')
# lines(lowess(fitted(growth.1_gls), resid(growth.1_gls)), col="blue", lwd=2)
# 
# 
# # resid vs. explanatory
# plot(x = dat$z[is.na(dat$z1) == F],
#      y = resid(growth.1_gls, type = 'normalized')
# )
# 
# plot(x = dat$touch_pct[is.na(dat$z1) == F],
#      y = resid(growth.1_gls, type = 'normalized')
# )
# 
# qqnorm(resid(growth.1_gls))
# qqline(resid(growth.1_gls))

# ---- plot: growth, z1 vs. z ----



### calculate predictions and CIs

critval <- 1.96

## touch = 0
# make blank area for new data
newgrowth.00 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0
  , z1 = 0
)
# add predictions
growthpreds.00 <- predict(growth.1, newgrowth.00, se.fit = T)
# CIs
growthCI.upper.00 <- growthpreds.00$fit + critval*growthpreds.00$se.fit
growthCI.lower.00 <- growthpreds.00$fit - critval*growthpreds.00$se.fit
growthfit.00 <- growthpreds.00$fit


## touch = 0.4
# make blank area for new data
newgrowth.04 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0.4
  , z1 = 0
)
# add predictions
growthpreds.04 <- predict(growth.1, newgrowth.04, se.fit = T)
# CIs
growthCI.upper.04 <- growthpreds.04$fit + critval*growthpreds.04$se.fit
growthCI.lower.04 <- growthpreds.04$fit - critval*growthpreds.04$se.fit
growthfit.04 <- growthpreds.04$fit

## touch = 0.8
# make blank area for new data
newgrowth.08 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0.8
  , z1 = 0
)
# add predictions
growthpreds.08 <- predict(growth.1, newgrowth.08, se.fit = T)
# CIs
growthCI.upper.08 <- growthpreds.08$fit + critval*growthpreds.08$se.fit
growthCI.lower.08 <- growthpreds.08$fit - critval*growthpreds.08$se.fit
growthfit.08 <- growthpreds.08$fit

### actually plot


pdf('plots/ch3_growth.absolute.pdf', width = 5, height = 5)
par(cex = 1.2
    , mar = c(4, 4, 0.5, 0.5)
)

## base plot
plot(z1 ~ z
     , data = dat
     , xlim = c(0, 6)
     , ylim = c(0, 8)
     , type = 'n'
     , axes = F
     , xlab = ''
     , ylab = ''
)

# add line of zero growth (1:1)
lines(seq(0, 6, 1), seq(0, 6, 1), col = 'gray')

## color points by touch
points(z1 ~ z
       , data = dat
       , pch = 19
       , col = dat$colors
)

## prediction lines and CIs
# touch = 0
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0
        y = c(growthCI.lower.00, rev(growthCI.upper.00)),
        col = alpha('gray', 0.5),
        border = NA)
lines(growthfit.00 ~ seq(0, 6, 0.1)
      , lwd = 2
      , lty = 3
)
# touch = 0.4
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0
        y = c(growthCI.lower.04, rev(growthCI.upper.04)),
        col = alpha('gray', 0.5),
        border = NA)
lines(growthfit.04 ~ seq(0, 6, 0.1)
      , lwd = 2
      , lty = 5
)
# touch = 0.8
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0
        y = c(growthCI.lower.08, rev(growthCI.upper.08)),
        col = alpha('gray', 0.5),
        border = NA)
lines(growthfit.08 ~ seq(0, 6, 0.1)
      , lwd = 2
      , lty = 1
)



legend(x = 0
       , y = 8
       , legend = c(0, 0.4, 0.8)
       , lty = c(3, 5, 1)
       , lwd = c(2, 2, 2)
       , y.intersp = 0.75
       , box.lty= 0
       , bg = 'transparent'
       , title = 'Crowding'
)

# axes and labels
axis(1 
     , at = seq(from = 0, to = 6, by = 1)
     , pos = 0
     , tck = 0.02
)
axis(2
     , las = 1
     , at = seq(from = 0, to = 8, by = 1)
     , pos = 0
     , tck = 0.02
)
axis(3
     , at = seq(from = 0, to = 6, by = 1)
     , pos = 8
     , tck = 0.02
     , labels = F
)
axis(4
     , seq(from = 0, to = 8, by = 1)
     , tck = 0.02
     , pos = 6
     , labels = F
)


mtext('Body size (operculum length, mm), year t'
      , side = 1
      , line = 2
      , cex = 1.2
)
mtext('Body size (operculum length, mm), year t+1'
      , side = 2
      , line = 2.5
      , cex = 1.2
)


dev.off()

# ---- analyze: growth rate (supplement) ----
gr.1 <- lm(growthrate ~ z * touch_pct
           , data = dat
)

summary(gr.1)
anova(gr.1)
Anova(gr.1) # same pattern as z1 vs. z: interaciton effect btwn z and crowd

# ---- plot: growth rate ----

### calculate predictions and CIs

critval <- 1.96

## touch = 0
# make blank area for new data
gr.newgrowth.00 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0
  , growthrate = 0
)
# add predictions
gr.growthpreds.00 <- predict(gr.1, gr.newgrowth.00, se.fit = T)
# CIs
gr.growthCI.upper.00 <- gr.growthpreds.00$fit + critval*gr.growthpreds.00$se.fit
gr.growthCI.lower.00 <- gr.growthpreds.00$fit - critval*gr.growthpreds.00$se.fit
gr.growthfit.00 <- gr.growthpreds.00$fit


## touch = 0.4
# make blank area for new data
gr.newgrowth.04 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0.4
  , growthrate = 0
)
# add predictions
gr.growthpreds.04 <- predict(gr.1, gr.newgrowth.04, se.fit = T)
# CIs
gr.growthCI.upper.04 <- gr.growthpreds.04$fit + critval*gr.growthpreds.04$se.fit
gr.growthCI.lower.04 <- gr.growthpreds.04$fit - critval*gr.growthpreds.04$se.fit
gr.growthfit.04 <- gr.growthpreds.04$fit

## touch = 0.8
# make blank area for new data
gr.newgrowth.08 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0.8
  , growthrate = 0
)
# add predictions
gr.growthpreds.08 <- predict(gr.1, gr.newgrowth.08, se.fit = T)
# CIs
gr.growthCI.upper.08 <- gr.growthpreds.08$fit + critval*gr.growthpreds.08$se.fit
gr.growthCI.lower.08 <- gr.growthpreds.08$fit - critval*gr.growthpreds.08$se.fit
gr.growthfit.08 <- gr.growthpreds.08$fit


### actually plot


pdf('plots/ch3_growth.rate.pdf', width = 5, height = 5)
par(cex = 1.2
    , mar = c(4, 4, 0.5, 0.5)
)

## base plot
plot(growthrate ~ z
     , data = dat
     # , xlim = c(0, 6)
     # , ylim = c(0, 8)
     # , type = 'n'
     #, axes = F
     , xlab = 'Body size (mm), year t'
     , ylab = 'Growth rate (mm per day)'
)

## color points by touch
points(growthrate ~ z
       , data = dat
       , pch = 19
       , col = dat$colors
)

## prediction lines and CIs
# touch = 0
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0
        y = c(gr.growthCI.lower.00, rev(gr.growthCI.upper.00)),
        col = alpha('gray', 0.5),
        border = NA)
lines(gr.growthfit.00 ~ seq(0, 6, 0.1)
      , lwd = 2
      , lty = 3
)
# touch = 0.4
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0
        y = c(gr.growthCI.lower.04, rev(gr.growthCI.upper.04)),
        col = alpha('gray', 0.5),
        border = NA)
lines(gr.growthfit.04 ~ seq(0, 6, 0.1)
      , lwd = 2
      , lty = 5
)
# touch = 0.8
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0
        y = c(gr.growthCI.lower.08, rev(gr.growthCI.upper.08)),
        col = alpha('gray', 0.5),
        border = NA)
lines(gr.growthfit.08 ~ seq(0, 6, 0.1)
      , lwd = 2
      , lty = 1
)

legend(x = 4
       , y = 0.01
       , legend = c(0, 0.4, 0.8)
       , lty = c(3, 5, 1)
       , lwd = c(2, 2, 2)
       , y.intersp = 0.75
       , box.lty= 0
       , bg = 'transparent'
       , title = 'Crowding'
)

dev.off()

# ---- analyze: survival ----

### determine random structure

# no random intercept or slope
survival.1 <- glm(Surv ~ z * touch_pct
                  , family = binomial
                  , data = dat
)

# random intercept for subsite only
survival.2 <- glmer(Surv ~ z * touch_pct
                    + (1 | subsite)
                    , family = binomial
                    , data = dat
)

# random intercept and slope for size
survival.3 <- glmer(Surv ~ z * touch_pct
                    + (1 + z | subsite)
                    , family = binomial
                    , data = dat
)

# random intercept and slope for touch
survival.4 <- glmer(Surv ~ z * touch_pct
                    + (1 + touch_pct | subsite)
                    , family = binomial
                    , data = dat
)

# random intercept for subsite nested in site
survival.5 <- glmer(Surv ~ z * touch_pct
                    + (1 | site/subsite)
                    , family = binomial
                    , data = dat
)

# random nested intercept and slope for size
survival.6 <- glmer(Surv ~ z * touch_pct
                    + (1 + z | site/subsite)
                    , family = binomial
                    , data = dat
)

# random nested intercept and slope for touch
survival.7 <- glmer(Surv ~ z * touch_pct
                    + (1 + touch_pct | site/subsite)
                    , family = binomial
                    , data = dat
)

# compare random structure models
AIC(survival.1
    , survival.2
    , survival.3
    , survival.4
    , survival.5
    , survival.6
    , survival.7
) # prefer survival.1
BIC(survival.1
    , survival.2
    , survival.3
    , survival.4
    , survival.5
    , survival.6
    , survival.7
) # prefer survival.1


### determine fixed structure
summary(survival.1) # everything NS
Anova(survival.1) # only initial size is significant


### rsquared using MuMIn::r.squaredGLMM
MuMIn::r.squaredGLMM(survival.1)

### model validation

# resid vs. fitted
plot(x = fitted(survival.1),
     y = resid(survival.1),
     xlab = 'fitted values',
     ylab = 'residuals')
lines(lowess(fitted(survival.1), resid(survival.1)), col="blue", lwd=2)


# resid vs. explanatory
plot(x = dat$z,
     y = resid(survival.1, type = 'deviance')
)

plot(x = dat$touch_pct,
     y = resid(survival.1, type = 'deviance')
)

qqnorm(resid(survival.1))
qqline(resid(survival.1))


### trim down fixed effects for parsimony, get it ready for IPM

mod.Surv <- glm(Surv ~ z
                , family = binomial
                , data = dat
)


### compare body size effect to null
survival.null <- glm(Surv ~ 1
                     , family = binomial
                     , data = dat
)

AIC(mod.Surv, survival.null) # AIC prefers mod.Surv (which includes body size effect)
BIC(mod.Surv, survival.null) # BIC prefers mod.Surv.go w/ mod.Surv

# ---- plot: survival ----

# make the model again, here, for convenience
mod.Surv <- glm(Surv ~ z
                , family = binomial
                , data = dat
)

### calculate predictions and CIs

# make blank area for new data
newsurv <- expand.grid(
  z = seq(0, 6, 0.1)
  , z1 = 0
)
# add predictions
survpreds <- predict(mod.Surv, newsurv, type = 'link', se.fit = T)
# CIs
critval <- 1.96
survCI.upper <- mod.Surv$family$linkinv(survpreds$fit + critval*survpreds$se.fit)
survCI.lower <- mod.Surv$family$linkinv(survpreds$fit - critval*survpreds$se.fit)
survfit <- mod.Surv$family$linkinv(survpreds$fit)



### actually plot


pdf('plots/ch3_survival.pdf', width = 5, height = 5)
par(cex = 1.2
    , mar = c(4, 4, 0.5, 0.5)
)



## base plot
plot(Surv.numeric ~ z
     , data = dat
     , xlim = c(0, 6)
     , ylim = c(0, 1)
     , type = 'n'
     , axes = F
     , xlab = ''
     , ylab = ''
)

## color points by touch
points(Surv.numeric ~ z
       , data = dat
       , pch = 19
       , col = dat$colors
)

## prediction lines and CIs
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1)))
        , y = c(survCI.lower, rev(survCI.upper))
        , col = alpha('gray', 0.5)
        , border = NA
)
lines(survfit ~ seq(0, 6, 0.1)
      , lwd = 2
)

# axes and labels
axis(1 
     , at = seq(from = 0, to = 6, by = 1)
     , pos = 0
     , tck = 0.02
)
axis(2
     , las = 1
     , pos = 0
     , tck = 0.02
)
axis(3
     , at = seq(from = 0, to = 6, by = 1)
     , pos = 1
     , tck = 0.02
     , labels = F
)
axis(4
     , tck = 0.02
     , pos = 6
     , labels = F
)


mtext('Operculum length (mm), year t'
      , side = 1
      , line = 2
      , cex = 1.2
)
mtext('Probability of survival'
      , side = 2
      , line = 2.5
      , cex = 1.2
)


dev.off()

# ~ IPM simple ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM simple: general parameters ----

### ceiling 
# U1 refers to upper limit of size z, beyond which vital rates are the same as U1
U1 = 6

### number of 'bins' for iteration matrix
nBigMatrix <- 500

### Probability of recruitment, Pr
# a) IPM book pg. 23. Can estimate by dividing number of recruits by total larval production. A problem with this is that I didn't survey every adult barnacle. But I'm assuming closed population and other assumptions anyway, so just go ahead and do it.
# b) can also adjust as needed to achieve desired lambda? Ken Sebens style

# a) 
# estimate total larval production
total_larval_est <- sum( 0.147 * (10 * dat$z[dat$Surv == 1])^2.74 )
# total number of recruits
total_rec <- mean(totrec$number_of_recruits) # mean total number of recruits for RT3

# calculate and set Pr
p.r.est <- total_rec / total_larval_est

# ---- IPM simple: prepare regression models ----

### survival

# final model from analysis only has body size anyway, so was prepared in section "analyze: survival", above.
# it is named mod.Surv
# repeated here for ease of use
mod.Surv <- glm(Surv ~ z
                , family = binomial
                , data = dat
)



### growth

# code repeated here for ease of use. see "analyze: growth"

## model w/ body size only (lumping all crowding); keep it simple to start with
mod.Grow_simple <- lmer(z1 ~ z
                        + (1 | subsite)
                        , data = dat
)

mod.Grow_simple
## fixed effects only:
# intercept = 2.2127
# slope for z = 0.732

# # calculate average intercept (averaged across sites
# 
# mean(coef(mod.Grow_simple)$subsite[1,1] 
#      , coef(mod.Grow_simple)$subsite[2,1]
#      , coef(mod.Grow_simple)$subsite[3,1]
#      , coef(mod.Grow_simple)$subsite[4,1]
#      , coef(mod.Grow_simple)$subsite[5,1]
# )

### reproduction
# s(z)*Pb(z)*b(z)*Pr*C0(z')

## probability of reproduction Pb(z)
# use ch 2 field regressions, see ch2>expt2.4.1_field-reproduction_bodyweights_20181106.R
reprodyn.1_sizeonly
# note that this model's coefficients are for VOLUME, not operc length. need to convert in IPM before use:
# # predict volume (v) from operc length (z), using regression from ch3_operctovol_20190102.R
# dat$v <- dat$z * 131.18 - 152.52


## fecundity b(z)
# offspring: # embryos = 0.147 * (operculum length in mm^2.74). Based on Strathmann et al. 1981 Oecologia (see Schubart thesis pg. 20 and email w/ Strathmann)
# operculum length multiplied by 10 to make it mm
# 0.147 * (10 * z)^2.74

## probability recruitment Pr
# unknown constant, set in "IPM simple: general parameters"

## recruit size distribution
# intercpet only; recruit size not dependent on parent size
# recruit size has Gaussian distribution w/ constant mean
mod.Rcsz <- lm(z ~ 1 
               , data = rec
)


# ---- IPM simple: prepare parameter estimates ----

# IPM book pg. 23-24 (regressions), 48 (ceiling)

m.par_simple <- c(
  surv = coef(mod.Surv)
  , grow.int = summary(mod.Grow_simple)$coefficients[1]
  , grow.z = summary(mod.Grow_simple)$coefficients[2]
  , grow.sd = summary(mod.Grow_simple)$sigma
  , rcsz = coef(mod.Rcsz)
  , rcsz.sd = summary(mod.Rcsz)$sigma
  , p.r = p.r.est
  , repr = coef(reprodyn.1_sizeonly)
  #, U1 = U1
)

names(m.par_simple) <- c(
  'surv.int'
  , 'surv.z'
  , 'grow.int'
  , 'grow.z'
  , 'grow.sd'
  , 'rcsz.int'
  , 'rcsz.sd'
  , 'p.r'
  , 'repr.int' # model is for volume, NOT operc length
  , 'repr.v' # v here b/c it's volume, NOT operc length
  #, 'U1'
)

# ---- IPM simple: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_simple <- function(z, m.par_simple){
  
  ## linear predictor
  linear.p <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * z # for cases below ceiling
  # linear.p.ceiling <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * U1 # for cases > U1
  
  ## back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  
  
  ## output
  return( unname(p) )
  
  # # code for manually doing ceiling
  # p <- ifelse(z <= U1
  #             , 1/(1+exp(-linear.p)) # below ceiling; use regression coeff of z
  #             , 1/(1+exp(-linear.p.ceiling)) # above ceiling; use vital rates of U1
  # )
}

### growth: G(z',z)

G_z1z_simple <- function(z1, z, m.par_simple){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  mu <- m.par_simple['grow.int'] + m.par_simple['grow.z'] * z # below ceiling
  #mu.ceiling <- m.par_simple['grow.int'] + m.par_simple['grow.z'] * U1 # above ceiling
  sig <- m.par_simple['grow.sd']
  
  ## calculate growth PDF
  
  p.den.grow <- dnorm(z1, mean = mu, sd = sig)
  
  ## output
  return(p.den.grow)
  
  
  
  
  ## graveyard
  
  # p.den.grow <- dnorm(z1, mean = mu, sd = sig)
  
  # #code for prevent shrinking. makes model output 'wonky'. not using.
  # p.den.grow <- ifelse(z1 > z # condition
  #                      , dnorm(z1, mean = mu, sd = sig) # if true. pdf of new size z1
  #                      , 0 # if false. forbid barnacles from shrinking
  # )
  
  # ##code for doing ceiling manually within function. However, do not use, b/c it messes up
  ##the numerical implementation of the kernel. Do the ceiling when forming mk_k.
  # # make blank vector
  # p.den.grow <- c()
  # 
  # loop through z1
  #   for( j in 1:length(z1) ){
  #     p.den.grow[j] <- ifelse(z <= U1
  #                             , dnorm(z1[j], mean = mu, sd = sig) # below ceiling
  #                             , dnorm(z1[j], mean = mu.ceiling, sd = sig) # above ceiling
  #     )
  #   }
  # }
  
}

### prob reproduce: p_bz

p_bz_simple <- function(z, m.par_simple){
  # convert z to v
  v <- z * 131.18 - 152.52
  # linear predictor
  linear.p <- m.par_simple['repr.int'] + m.par_simple['repr.v'] * v
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  # output
  return( unname(p) )
}

### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  p.bz <- 0.147 * (10 * z)^2.74
  return(p.bz)
}

### recruitment: C_0(z') = R(z')

C_0z1 <- function(z1, m.par_simple){
  mu <- m.par_simple['rcsz.int']
  sig <- m.par_simple['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}


### define P component of kernel (survival and growth): P(z',z) = s(z)*G(z',z)

P_z1z <- function(z1, z, m.par_simple){
  return( s_z_simple(z, m.par_simple) * G_z1z_simple(z1, z, m.par_simple) 
  )
}

### define F component of kernel (reproduction): F(z', z) = s(z)*Pb(z)*b(z)*Pr*Co(z')

F_z1z <- function(z1, z, m.par_simple){
  return( s_z_simple(z, m.par_simple) * p_bz_simple(z, m.par_simple) * b_z(z) * m.par_simple['p.r'] * C_0z1(z1, m.par_simple) 
  )
}

# # ---- IPM simple: explore kernel components ----
# 
# # Growth: see MonocarpGrowthEviction.R and IPM book pg. 45-48
# 
# # use a modified growth function here that implements the ceiling manually, for exploration purposes. However, this function is not used for numerical implementation.
# 
# ### create growth ceiling demo function
# G_z1z_simple.ceiling <- function(z1, z, m.par_simple){
#   
#   ## based on growth regression coefficients, define mu and sigma of growth PDF
#   mu <- m.par_simple['grow.int'] + m.par_simple['grow.z'] * z # below ceiling
#   mu.ceiling <- m.par_simple['grow.int'] + m.par_simple['grow.z'] * U1 # above ceiling
#   sig <- m.par_simple['grow.sd']
#   
#   ## calculate growth PDF
#   
#   # make blank vector
#   p.den.grow <- c()
#   
#   #loop through z1
#   for( j in 1:length(z1) ){
#     p.den.grow[j] <- ifelse(z <= U1
#                             , dnorm(z1[j], mean = mu, sd = sig) # below ceiling
#                             , dnorm(z1[j], mean = mu.ceiling, sd = sig) # above ceiling
#     )
#   }
#   
#   ## output
#   return(p.den.grow)
# }
# 
# 
# ### explore growth
# 
# zvals <- seq(0, 10, length = 1000)
# g1 <- G_z1z_simple.ceiling(zvals, min(zvals), m.par_simple) # small individuals
# g2 <- G_z1z_simple.ceiling(zvals, mean(zvals), m.par_simple) # medium individuals
# g3 <- G_z1z_simple.ceiling(zvals, max(zvals), m.par_simple) # large individuals
# 
# plot(x = zvals
#      , y = g1
#      , type = 'l'
#      , lty = 1
#      , lwd = 1
#      , xlab = "Final size z'"
#      , ylab = 'Probability density'
#      , ylim = c(0, 1)
# )
# lines(x = zvals
#       , y = g2
#       , lwd = 2
# )
# lines(x = zvals
#       , y = g3
#       , lwd = 3
# )
# 
# ### create survival ceiling demo function
# s_z_simple.ceiling <- function(z, m.par_simple){
#   
#   ## linear predictor
#   linear.p <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * z # for cases below ceiling
#   linear.p.ceiling <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * U1 # for cases > U1
#   
#   ## back transform from logit
#   # code for manually doing ceiling
#   p <- ifelse(z <= U1
#               , 1/(1+exp(-linear.p)) # below ceiling; use regression coeff of z
#               , 1/(1+exp(-linear.p.ceiling)) # above ceiling; use vital rates of U1
#   )
#   
#   
#   ## output
#   return( unname(p) )
#   
#   
# }
# 
# 
# 
# ### explore fraction wrongfully evicted
# 
# WrongPlace <- function(z, U) {
#   fac1 <- s_z_simple.ceiling(z, m.par_simple)
#   fac2 <- integrate(function(x) G_z1z_simple.ceiling(x, z, m.par_simple), U, Inf)$value
#   return(fac1 * fac2)
# }
# 
# zvals <- Wvals <- seq(0, 9, length = 900)
# 
# for (j in seq_along(zvals)) Wvals[j] <- WrongPlace(zvals[j], 9)
# 
# plot(zvals
#      , Wvals
#      , type = "l"
#      , lty = 1
#      , lwd = 2
#      , col = "black"
#      , xlab = "Initial size z"
#      , ylab = "Fraction wrongfully evicted"
#      , ylim = c(0, 1.1 * max(Wvals))
# )
# 
# ---- IPM simple: numerical implementation ----

### use the mid-point rule

# m = number of mesh points = nBigMatrix
# L = 0.1  = lower bound; minimum observed recruit size in field is 0.2 mm operc length
# U = 10 = upper bound; max observed size in field is 7.637 mm operc length, but I'll do ceiling

mk_K <- function(m, m.par, L, U) {
  
  # compute mesh points
  h <- (U - L)/m
  meshpts <- L + ((1:m) - 1/2) * h # equivalent to L + (1:m)*h - h/2
  
  # compute iteration matrix 
  P <- h * (outer(meshpts, pmin(meshpts, U1), P_z1z, m.par = m.par)) # w/ceiling. IPM book pg. 48
  F <- h * (outer(meshpts, pmin(meshpts, U1), F_z1z, m.par = m.par)) # w/ceiling. IPM book pg. 48
  
  # F <- matrix(0, nrow = m, ncol = m)
  # F[1:m, ] <- F_z1(meshpts, m.par)
  
  K <- P + F
  
  #matrix.image(K) # look at kernel
  
  return( list(K = K, meshpts = meshpts, P = P, F = F) )
  
}
#dev.off()

### create IPM kernel

# IPM_simple <- mk_K(nBigMatrix, m.par_simple, 0.1, 10)
IPM_simple <- mk_K(100, m.par_simple, 0.1, 10)


### dominant eigenvalue, lambda, population growth rate
Re(eigen(IPM_simple$K)$values[1]) 

### eigenvector, w, stable size distribution
meshpts_simple <- IPM_simple$meshpts
w_simple <- Re(eigen(IPM_simple$K)$vectors[,1])
stable.z.dist_simple <- w_simple/sum(w_simple) # scale eigenvector to discrete probability density

plot(stable.z.dist_simple/diff(meshpts_simple)[1]  ~ meshpts_simple
     , type = 'l'
     , xlab = 'Operculum length, mm'
     , ylab = 'Probability'
     , ylim = c(0, 1)
) 

#from Monocarp Calculations.R



# ---- IPM simple: perturbation analysis: set up ----

# see IPM book section 4.4, starting on pg. 96
# following their naming convention

### set up

IPM.sys <- IPM_simple
K <- IPM.sys$K
meshpts <- IPM.sys$meshpts
h <- diff(meshpts[1:2])

IPM.eig.sys <- eigen(K)
lambda <- Re(IPM.eig.sys$values[1]) # dominant right eigenvalue
w.z <- Re(IPM.eig.sys$vectors[,1]) # dominant right eigenvector
v.z1 <- Re(eigen(t(K))$vectors[,1]) # dominant left eigenvector

# ---- IPM simple: perturbation analysis: kernel-level ----

## calculate sensitivity for K
K.sens <- outer(v.z1, w.z, "*") / sum(v.z1 * w.z * h)
# plot
image(x = meshpts
      , y = meshpts
      , z =  t(K.sens)
      , col = grey( seq(0.6, 1, length = 100) )
      , xlab = "Operculum length, mm (t), z"
      , ylab = "Operculum length, mm (t+1), z\'"
)
contour(x = meshpts
        , y = meshpts
        , z = t(K.sens)
        , add = TRUE
)

## calculate elasticity for K
K.elas <- K.sens * (K / h) / lambda
# check elasticity: integral of elasticity should 'sum to 1', see pg. 97
sum(K.elas) * h^2
# plot
image(x = meshpts
      , y = meshpts
      , z =  t(K.elas)
      , col = grey( seq(0.6, 1, length = 100) )
      , xlab = "Operculum length, mm (t), z"
      , ylab = "Operculum length, mm (t+1), z\'"
)
contour(x = meshpts
        , y = meshpts
        , z = t(K.elas)
        , add = TRUE
)

## total elasticity (all transitions to anywhere)
K.elas.total <- colSums(K.elas)
plot(K.elas.total ~ meshpts
     , xlab = 'Operculum length, mm'
     , ylab = 'Total elasticity'
     , type = 'l'
     )


## calculate elasticity for P
Pvals <- IPM.sys$P / h
P.elas <- Pvals * K.sens / lambda
# plot
image(x = meshpts
      , y = meshpts
      , z =  t(P.elas)
      , col = grey( seq(0.6, 1, length = 100) )
      , xlab = "Operculum length, mm (t), z"
      , ylab = "Operculum length, mm (t+1), z\'"
)
contour(x = meshpts
        , y = meshpts
        , z = t(P.elas)
        , add = TRUE
)

## calculate elasticity for F
Fvals <- IPM.sys$F / h
F.elas <- Fvals * K.sens / lambda
# plot
image(x = meshpts
      , y = meshpts
      , z =  t(F.elas)
      , col = grey( seq(0.6, 1, length = 100) )
      , xlab = "Operculum length, mm (t), z"
      , ylab = "Operculum length, mm (t+1), z\'"
)
contour(x = meshpts
        , y = meshpts
        , z = t(F.elas)
        , add = TRUE
)

## calculate relative contributions of P and F
sum(P.elas) * h^2 # 0.68
sum(F.elas) * h^2 # 0.32

# ---- IPM simple: perturbation analysis: vital rate functions ---- 

# see IPM book pg. 99 and ch 3 notebook 1/8/2019 entry

### local perturbation analysis for s(z)

## calculate sensitivity
dK_by_ds_z1z <- outer(meshpts
                      , meshpts
                      , function(z1, z, m.par_simple){
                        G_z1z_simple(z1, z, m.par_simple) 
                        + p_bz_simple(z, m.par_simple)*b_z(z)*p.r.est*C_0z1(z1, m.par_simple)
                      }
                      , m.par_simple
)
s.sens.z <- apply(K.sens * dK_by_ds_z1z, 2, sum) * h

## calculate elasticity
s.elas.z <- s.sens.z * s_z_simple(meshpts, m.par_simple) / lambda

## plot
plot(s.sens.z ~ meshpts
     , type = 'l'
     , lty = 1
     , xlab = 'operculum length, mm'
     , ylab = 'sensitivity or elasticity of s(z)'
)
lines(s.elas.z ~ meshpts
      , type = 'l'
      , lty = 2
)
legend('topleft'
       , legend = c('sensitivity', 'elasticity')
       , lty = c(1, 2)
       , bty = 'n'
)

### local perturbation analysis for pb(z)

## calculate sensitivity
dK_by_dpb_z1z <- outer(meshpts
                       , meshpts
                       , function(z1, z, m.par_simple){
                         s_z_simple(z, m.par_simple)*b_z(z)*p.r.est*C_0z1(z1, m.par_simple)
                       }
                       , m.par_simple
)
pb.sens.z <- apply(K.sens * dK_by_dpb_z1z, 2, sum) * h

## calculate elasticity
pb.elas.z <- pb.sens.z * p_bz_simple(meshpts, m.par_simple) / lambda

## plot
plot(pb.sens.z ~ meshpts
     , type = 'l'
     , lty = 1
     , xlab = 'operculum length, mm'
     , ylab = 'sensitivity or elasticity of Pb(z)'
)
lines(pb.elas.z ~ meshpts
      , type = 'l'
      , lty = 2
)
legend('topleft'
       , legend = c('sensitivity', 'elasticity')
       , lty = c(1, 2)
       , bty = 'n'
)

### local perturbation analysis for b(z)

## calculate sensitivity
dK_by_db_z1z <- outer(meshpts
                      , meshpts
                      , function(z1, z, m.par_simple){
                        s_z_simple(z, m.par_simple)*p_bz_simple(z, m.par_simple)*p.r.est*C_0z1(z1, m.par_simple)
                      }
                      , m.par_simple
)
b.sens.z <- apply(K.sens * dK_by_db_z1z, 2, sum) * h

## calculate elasticity
b.elas.z <- b.sens.z * b_z(meshpts) / lambda

## plot
plot(b.elas.z ~ meshpts # note plotting elasticity first here, b/c much larger scale
     , type = 'l'
     , lty = 2
     , xlab = 'operculum length, mm'
     , ylab = 'sensitivity or elasticity of b(z)'
)
lines(b.sens.z ~ meshpts
      , type = 'l'
      , lty = 1
)
legend('topleft'
       , legend = c('sensitivity', 'elasticity')
       , lty = c(1, 2)
       , bty = 'n'
)

### local perturbation analysis for G(z',z)

## calculate sensitivity
dK_by_dg_z1z <- outer(meshpts
                      , meshpts
                      , function(z1, z, m.par_simple){
                        s_z_simple(z, m.par_simple)
                      }
                      , m.par_simple
)
g.sens.z1z <- K.sens * dK_by_dg_z1z

## calculate elasticity
g.elas.z1z <- g.sens.z1z * outer(meshpts, meshpts, G_z1z_simple, m.par_simple) / lambda

## plot

# sensitivity
image(x = meshpts
      , y = meshpts
      , z =  t(g.sens.z1z)
      , col = grey( seq(0.6, 1, length = 100) )
      , xlab = "Operculum length, mm (t), z"
      , ylab = "Operculum length, mm (t+1), z\'"
)
contour(x = meshpts
        , y = meshpts
        , z = t(g.sens.z1z)
        , add = TRUE
)

# elasticity
image(x = meshpts
      , y = meshpts
      , z =  t(g.elas.z1z)
      , col = grey( seq(0.6, 1, length = 100) )
      , xlab = "Operculum length, mm (t), z"
      , ylab = "Operculum length, mm (t+1), z\'"
)
contour(x = meshpts
        , y = meshpts
        , z = t(g.elas.z1z)
        , add = TRUE
)


### local perturbation analysis for C0(z')

## calculate sensitivity
dK_by_dc0_z1z <- outer(meshpts
                       , meshpts
                       , function(z1, z, m.par_simple){
                         s_z_simple(z, m.par_simple)*p_bz_simple(z, m.par_simple)*b_z(z)*p.r.est
                       }
                       , m.par_simple
)
c0.sens.z1 <- K.sens * dK_by_dc0_z1z

## calculate elasticity
c0.elas.z1 <- c0.sens.z1 * C_0z1(meshpts, m.par_simple) / lambda

## plot
# sensitivity
image(x = meshpts
      , y = meshpts
      , z =  t(c0.sens.z1)
      , col = grey( seq(0.6, 1, length = 100) )
      , xlab = "Operculum length, mm (t), z"
      , ylab = "Operculum length, mm (t+1), z\'"
)
contour(x = meshpts
        , y = meshpts
        , z = t(c0.sens.z1)
        , add = TRUE
)

# elasticity
image(x = meshpts
      , y = meshpts
      , z =  t(c0.elas.z1)
      , col = grey( seq(0.6, 1, length = 100) )
      , xlab = "Operculum length, mm (t), z"
      , ylab = "Operculum length, mm (t+1), z\'"
)
contour(x = meshpts
        , y = meshpts
        , z = t(c0.elas.z1)
        , add = TRUE
)

### perturbation at all sizes 

# see IPM book pg. 102

## survival
# elasticity
sum(s.elas.z * h) # 0.365... would've expected it to be 1? see pg. 102...
# sensitivity
sum(s.sens.z * h) # 0.980

# pb(z)
# elasticity
sum(pb.elas.z * h) # 0.365
# sensitivity
sum(pb.sens.z * h) # 0.594

# b(z)
# elasticity
sum(b.elas.z * h) # 0.365
# sensitivity
sum(b.sens.z * h) # 6.955e-05

## growth
# elasticity
sum(g.elas.z1z * h^2) # 0.690
# sensitivity
sum(g.sens.z1z * h^2) # 12.721

## c0(z',z)
# elasticity
sum(c0.elas.z1 * h^2) # 0.365
# sensitivity
sum(c0.sens.z1 * h^2) # 60.583




# ---- IPM simple: perturbation analysis: parameters ----

# see IPM book pg 103 and Ungulate Pert Lambda Calc.R

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ perturbation at all sizes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### survival intercept 
ds_by_dbeta0_z1z <- outer(meshpts
                          , meshpts
                          , function(z1, z, m.par_simple) {
                            nu <- m.par_simple["surv.int"] + m.par_simple["surv.z"] * z 
                            exp(nu)/(1+exp(nu))^2
                          }
                          , m.par_simple
)

s.int.sens <- sum(K.sens * dK_by_ds_z1z * ds_by_dbeta0_z1z) * h^2
s.int.elas <- s.int.sens * m.par_simple["surv.int"] / lambda


### survival size slope
ds_by_dbetaz_z1z <- outer(meshpts
                          , meshpts
                          , function(z1, z, m.par_simple) {
                            nu <- m.par_simple["surv.int"] + m.par_simple["surv.z"] * z 
                            z * exp(nu)/(1+exp(nu))^2
                          }
                          , m.par_simple
)

s.slp.sens <- sum(K.sens * dK_by_ds_z1z * ds_by_dbetaz_z1z) * h^2
s.slp.elas <- s.slp.sens * m.par_simple["surv.z"] / lambda

### growth intercept
dg_by_dbeta0_z1z <- outer(meshpts
                          , meshpts
                          , function(z1, z, m.par) {
                            G_z1z_simple(z1, z, m.par_simple) * (z1 - m.par_simple["grow.int"] - m.par_simple["grow.z"] * z) / m.par_simple["grow.sd"]^2
                          }
                          , m.par_simple)

g.int.sens <- sum(K.sens * dK_by_dg_z1z * dg_by_dbeta0_z1z) * h^2
g.int.elas <- g.int.sens * m.par_simple["grow.int"] / lambda

### growth size slope
dg_by_dbetaz_z1z <- outer(meshpts
                          , meshpts
                          , function(z1, z, m.par) {
                            G_z1z_simple(z1, z, m.par_simple) * z * (z1 - m.par_simple["grow.int"] - m.par_simple["grow.z"] * z) / m.par_simple["grow.sd"]^2
                          }
                          , m.par_simple
)

g.slp.sens <- sum(K.sens * dK_by_dg_z1z * dg_by_dbetaz_z1z) * h^2
g.slp.elas <- g.slp.sens * m.par_simple["grow.z"] / lambda

### growth standard deviation
dg_by_dsigma_z1z <- outer(meshpts
                          , meshpts
                          , function(z1, z, m.par_simple) {
                            nu <- m.par_simple["grow.int"] + m.par_simple["grow.z"] * z
                            G_z1z_simple(z1, z, m.par_simple) * ((z1 - nu)^2 - m.par_simple["grow.sd"]^2) / m.par_simple["grow.sd"]^3
                          }
                          , m.par_simple)

g.sig.sens <- sum(K.sens * dK_by_dg_z1z * dg_by_dsigma_z1z) * h^2
g.sig.elas <- g.sig.sens * m.par_simple["grow.sd"] / lambda

### offspring size intercept
dc_by_dbeta0_z1z <- outer(meshpts
                          , meshpts
                          , function(z1, z, m.par_simple) {
                            C_0z1(z1, m.par_simple) * (z1 - m.par_simple["rcsz.int"]) / m.par_simple["rcsz.sd"]^2
                          }
                          , m.par_simple
)

c.int.sens <- sum(K.sens * dK_by_dc0_z1z * dc_by_dbeta0_z1z) * h^2
c.int.elas <- c.int.sens * m.par_simple["rcsz.int"] / lambda

### offsping size standard deviation
dc_by_dsigma_z1z <- outer(meshpts
                          , meshpts
                          , function(z1, z, m.par_simple) {
                            nu <- m.par_simple["rcsz.int"]
                            C_0z1(z1, m.par_simple) * ((z1 - nu)^2 - m.par_simple["rcsz.sd"]^2) / m.par_simple["rcsz.sd"]^3
                          }
                          , m.par_simple)

c.sig.sens <- sum(K.sens * dK_by_dc0_z1z * dc_by_dsigma_z1z) * h^2
c.sig.elas <- c.sig.sens * m.par_simple["rcsz.sd"] / lambda


### probability of reproduction intercept 
dpb_by_dbeta0_z1z <- outer(meshpts
                           , meshpts
                           , function(z1, z, m.par_simple) {
                             # convert z to v
                             v <- z * 131.18 - 152.52
                             nu <- m.par_simple["repr.int"] + m.par_simple["repr.v"] * v 
                             exp(nu)/(1+exp(nu))^2
                           }
                           , m.par_simple
)

pb.int.sens <- sum(K.sens * dK_by_dpb_z1z * dpb_by_dbeta0_z1z) * h^2
pb.int.elas <- pb.int.sens * m.par_simple["repr.int"] / lambda


### probability of reproduction size slope
dpb_by_dbetaz_z1z <- outer(meshpts
                           , meshpts
                           , function(z1, z, m.par_simple) {
                             # convert z to v
                             v <- z * 131.18 - 152.52
                             nu <- m.par_simple["repr.int"] + m.par_simple["repr.v"] * v 
                             v * exp(nu)/(1+exp(nu))^2
                           }
                           , m.par_simple
)

pb.slp.sens <- sum(K.sens * dK_by_dpb_z1z * dpb_by_dbeta0_z1z) * h^2
pb.slp.elas <- s.slp.sens * m.par_simple["repr.v"] / lambda

# ~~~~~~~~~~~~~~~~~ perturbation at local sizes of expected size functions ~~~~~~~~~~~~~~~~~~~

# see IPM book pg. 105 and Ungulate PErt Lambda Calc.R

### expected growth function
mu.g.sens.z <- apply(K.sens * dK_by_dg_z1z * dg_by_dbeta0_z1z, 2, sum) * h
mu.g <- m.par_simple["grow.int"] + m.par_simple["grow.z"] * meshpts
mu.g.elas.z <- mu.g.sens.z * mu.g / lambda

## plot
plot(mu.g.elas.z ~ meshpts
     , type = 'l'
     , lty = 2
     , xlab = 'operculum length, mm'
     , ylab = 'sensitivity or elasticity of expected growth function'
)
lines(mu.g.sens.z ~ meshpts
      , type = 'l'
      , lty = 1
)
legend('topright'
       , legend = c('sensitivity', 'elasticity')
       , lty = c(1, 2)
       , bty = 'n'
)

### expected size of offspring function
mu.c.sens.z <- apply(K.sens * dK_by_dc0_z1z * dc_by_dbeta0_z1z, 2, sum) * h
mu.c <- m.par_simple["rcsz.int"]
mu.c.elas.z <- mu.c.sens.z * mu.c / lambda

## plot
plot(mu.c.elas.z ~ meshpts
     , type = 'l'
     , lty = 2
     , xlab = 'operculum length, mm'
     , ylab = 'sensitivity or elasticity of expected offspring size function'
)
lines(mu.c.sens.z ~ meshpts
      , type = 'l'
      , lty = 1
)
legend('topright'
       , legend = c('sensitivity', 'elasticity')
       , lty = c(1, 2)
       , bty = 'n'
)

# ---- *slow, ~ 7 min* IPM simple: lambda CI  ----

# from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30

# NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.

### function to compute lambda from a bootstrapped data set in format required by library(boot)
boot.lam <- function(dataset, sample.index) {

  ### extract the data used to make this fit
  boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.

  ### fit the functions

  ## survival
  mod.Surv <- glm(Surv ~ z
                  , family = binomial
                  , data = boot.data
  )

  ## growth
  boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]

  mod.Grow <- lmer(z1 ~ z
                   + (1 | subsite)
                   , data = boot.data.growth
  )

  ## recruit size distribution
  mod.Rcsz <- lm(z ~ 1
                 , data = rec
  )

  ### Store the estimated parameters

  m.par_simple <- c(
    surv = coef(mod.Surv)
    , grow.int = summary(mod.Grow)$coefficients[1]
    , grow.z = summary(mod.Grow)$coefficients[2]
    , grow.sd = summary(mod.Grow)$sigma
    , rcsz = coef(mod.Rcsz)
    , rcsz.sd = summary(mod.Rcsz)$sigma
    , p.r = p.r.est
    , repr = coef(reprodyn.1_sizeonly)
    #, U1 = U1
  )

  names(m.par_simple) <- c(
    'surv.int'
    , 'surv.z'
    , 'grow.int'
    , 'grow.z'
    , 'grow.sd'
    , 'rcsz.int'
    , 'rcsz.sd'
    , 'p.r'
    , 'repr.int' # model is for volume, NOt operc length
    , 'repr.v' # v here b/c it's volume, NOT operc length
    #, 'U1'
  )

  ### implement IPM and calculate lambda
  # also see IPM: general parameters

  IPM.est <- mk_K(nBigMatrix, m.par_simple, 0.1, 10)

  lam.boot <- Re(eigen(IPM.est$K, only.values = TRUE)$values[1])
  cat(lam.boot, "\n")

  return(lam.boot)
}
### do the bootstrap (code takes ~ 7 min to run)
#starttime <- Sys.time()
boot.out <- boot(data = dat, statistic = boot.lam, simple = TRUE,
                 R = 1000)
#endtime <- Sys.time()

boot.ci(boot.out, type = c('norm', 'basic', 'perc'))



# ~ IPM size only (= simple) ~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size only: general parameters ----

### ceiling 
# Uz1 refers to upper limit of size z, beyond which vital rates are the same as Uz1
Uz1 = 6
# Ut1 refers to upper limit of touch t, beyond withich vital rates are the same as Ut1
Ut1 = 0.8

### number of 'bins' for iteration matrix
nBigMatrix_z <- 500
nBigMatrix_t <- 50

### Probability of recruitment, Pr
# a) IPM book pg. 23. Can estimate by dividing number of recruits by total larval production. A problem with this is that I didn't survey every adult barnacle. But I'm assuming closed population and other assumptions anyway, so just go ahead and do it.
# b) can also adjust as needed to achieve desired lambda? Ken Sebens style

# a) 
# estimate total larval production
total_larval_est <- sum( 0.147 * (10 * dat$z[dat$Surv == 1])^2.74 )
# total number of recruits
total_rec <- mean(totrec$number_of_recruits) # mean total number of recruits for RT3

# calculate and set Pr
p.r.est <- total_rec / total_larval_est

# ---- IPM size only: prepare regression models ----

### survival

# final model from analysis only has body size anyway, so was prepared in section "analyze: survival", above.
# it is named mod.Surv
# repeated here for ease of use
# this is the exact same model as used in IPM_simple
mod.Surv <- glm(Surv ~ z
                , family = binomial
                , data = dat
)



### growth

## model w/ body size * touch
dat.growth <- dat[is.na(dat$z1) == F, ]

mod.Grow_st <- lm(z1 ~ z * touch_pct
                  , data = dat.growth
                  
)

# # code for modeling non-constant variance repeated here for ease of use. see "analyze: growth"

# z.growth <- dat.growth$z
# vfix <- varFixed(~(1/z.growth))



### reproduction
# s(z)*Pb(z)*b(z)*Pr*C0(z')

## probability of reproduction Pb(z)
# use ch 2 field regressions, see ch2>expt2.4.1_field-reproduction_bodyweights_20181106.R, but here using the version that you created for use in ch 3 (see 'load data and packages' section, above)
# here using the model that includes size + touch (no intxn effect)
reprodyn.1b
# note that this model's coefficients are for VOLUME, not operc length. need to convert in IPM before use:
# # predict volume (v) from operc length (z), using regression from ch3_operctovol_20190102.R
# dat$v <- dat$z * 131.18 - 152.52


## fecundity b(z)
# offspring: # embryos = 0.147 * (operculum length in mm^2.74). Based on Strathmann et al. 1981 Oecologia (see Schubart thesis pg. 20 and email w/ Strathmann)
# operculum length multiplied by 10 to make it mm
# 0.147 * (10 * z)^2.74

## probability recruitment Pr
# unknown constant, set in "IPM_st: general parameters"

## recruit size distribution
# intercpet only; recruit size not dependent on parent size
# recruit size has Gaussian distribution w/ constant mean
mod.Rcsz <- lm(z ~ 1 
               , data = rec
)

## recruit touch distribution
# don't need a regression; just use beta distribution as determined above


# ---- IPM size only: prepare parameter estimates ----

# IPM book pg. 23-24 (regressions), 48 (ceiling)

m.par_st <- c(
  surv = coef(mod.Surv)
  , grow = coef(mod.Grow_st)
  , grow.sd = summary(mod.Grow_st)$sigma
  , rcsz = coef(mod.Rcsz)
  , rcsz.sd = summary(mod.Rcsz)$sigma
  , rcst = coef(betafit_t)
  , p.r = p.r.est
  , repr = coef(reprodyn.1b)
  #, U1 = U1
)

names(m.par_st) <- c(
  'surv.int'
  , 'surv.z'
  , 'grow.int'
  , 'grow.z'
  , 'grow.t'
  , 'grow.zt'
  , 'grow.sd'
  , 'rcsz.int'
  , 'rcsz.sd'
  , 'rcst.s1'
  , 'rcst.s2'
  , 'p.r'
  , 'repr.int' # model is for volume, NOT operc length
  , 'repr.t' # model is for volume, NOT operc length
  , 'repr.v' # v here b/c it's volume, NOT operc length
  #, 'U1'
)


# ---- IPM size only (regression coefficients for touch = 0): define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + 0 * t + 0 * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + 0 * t + 0 * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = m.par_st['rcst.s1']
                       , shape2 = m.par_st['rcst.s2']
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + 0 * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
                      , shape1 = 1
                      , shape2 = 1
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size only: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size only: perturbation analysis ----

### Compute elasticity matrix 
repro.val <- matrix(v, mz, mt)
stable.state <- matrix(w, mz, mt) 
stable.state.vector <- apply(stable.state, 1, sum)
v.dot.w <- sum(hz*ht*stable.state*repro.val)
sens <- outer(repro.val,stable.state)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state.vector/sum(stable.state.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size only: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets. 
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
#   
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
#   
#   ### fit the functions
#   
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
#   
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lm(z1 ~ z * touch_pct
#                     , data = boot.data.growth
#                     
#   )
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1 
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow = coef(mod.Grow_st)
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
#   
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
#   
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
#   
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points 
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
#   
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
#   
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
#   
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
#         
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
#   
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
#   
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ~ IPM size touch: prep ~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch: general parameters ----

### ceiling 
# Uz1 refers to upper limit of size z, beyond which vital rates are the same as Uz1
Uz1 = 6

### number of 'bins' for iteration matrix
nBigMatrix_z <- 500
nBigMatrix_t <- 50

### Probability of recruitment, Pr
# a) IPM book pg. 23. Can estimate by dividing number of recruits by total larval production. A problem with this is that I didn't survey every adult barnacle. But I'm assuming closed population and other assumptions anyway, so just go ahead and do it.
# b) can also adjust as needed to achieve desired lambda? Ken Sebens style

# a) 
# estimate total larval production
total_larval_est <- sum( 0.147 * (10 * dat$z[dat$Surv == 1])^2.74 )
# total number of recruits
total_rec <- mean(totrec$number_of_recruits) # mean total number of recruits for RT3

# calculate and set Pr
p.r.est <- total_rec / total_larval_est

# ---- IPM size touch: prepare regression models ----

### survival

# final model from analysis only has body size anyway, so was prepared in section "analyze: survival", above.
# it is named mod.Surv
# repeated here for ease of use
# this is the exact same model as used in IPM_simple
mod.Surv <- glm(Surv ~ z
                , family = binomial
                , data = dat
)



### growth

## model w/ body size * touch
dat.growth <- dat[is.na(dat$z1) == F, ]

mod.Grow_st <- lmer(z1 ~ z * touch_pct
                    + (1 | subsite)
                    , data = dat.growth
)


# # code for modeling non-constant variance repeated here for ease of use. see "analyze: growth"

# z.growth <- dat.growth$z
# vfix <- varFixed(~(1/z.growth))



### reproduction
# s(z)*Pb(z)*b(z)*Pr*C0(z')

## probability of reproduction Pb(z)
# use ch 2 field regressions, see ch2>expt2.4.1_field-reproduction_bodyweights_20181106.R, but here using the version that you created for use in ch 3 (see 'load data and packages' section, above)
# here using the model that includes size + touch (no intxn effect)
reprodyn.1b
# note that this model's coefficients are for VOLUME, not operc length. need to convert in IPM before use:
# # predict volume (v) from operc length (z), using regression from ch3_operctovol_20190102.R
# dat$v <- dat$z * 131.18 - 152.52


## fecundity b(z)
# offspring: # embryos = 0.147 * (operculum length in mm^2.74). Based on Strathmann et al. 1981 Oecologia (see Schubart thesis pg. 20 and email w/ Strathmann)
# operculum length multiplied by 10 to make it mm
# 0.147 * (10 * z)^2.74

## probability recruitment Pr
# unknown constant, set in "IPM_st: general parameters"

## recruit size distribution
# intercpet only; recruit size not dependent on parent size
# recruit size has Gaussian distribution w/ constant mean
mod.Rcsz <- lm(z ~ 1 
               , data = rec
)

## recruit touch distribution
# don't need a regression; just use beta distribution as determined above


# ---- IPM size touch: prepare parameter estimates ----

# IPM book pg. 23-24 (regressions), 48 (ceiling)

m.par_st <- c(
  surv = coef(mod.Surv)
  , grow.int = summary(mod.Grow_st)$coefficients[1]
  , grow.z = summary(mod.Grow_st)$coefficients[2]
  , grow.t = summary(mod.Grow_st)$coefficients[3]
  , grow.zt = summary(mod.Grow_st)$coefficients[4]
  , grow.sd = summary(mod.Grow_st)$sigma
  , rcsz = coef(mod.Rcsz)
  , rcsz.sd = summary(mod.Rcsz)$sigma
  , rcst = coef(betafit_t)
  , p.r = p.r.est
  , repr = coef(reprodyn.1b)
  #, U1 = U1
)

names(m.par_st) <- c(
  'surv.int'
  , 'surv.z'
  , 'grow.int'
  , 'grow.z'
  , 'grow.t'
  , 'grow.zt'
  , 'grow.sd'
  , 'rcsz.int'
  , 'rcsz.sd'
  , 'rcst.s1'
  , 'rcst.s2'
  , 'p.r'
  , 'repr.int' # model is for volume, NOT operc length
  , 'repr.t' # model is for volume, NOT operc length
  , 'repr.v' # v here b/c it's volume, NOT operc length
  #, 'U1'
)


# ~ aFrF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----

# ---- IPM size touch - aFrF: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = m.par_st['rcst.s1']
                       , shape2 = m.par_st['rcst.s2']
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
                      , shape1 = betafit_t.r$estimate[1]
                      , shape2 = betafit_t.r$estimate[2]
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aFrF: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# # ---- IPM size touch - aFrF: perturbation analysis ----
# 
# ### Compute elasticity matrix 
# repro.val_aFrF <- matrix(v, mz, mt)
# stable.state_aFrF <- matrix(w, mz, mt) 
# stable.state_aFrF.vector <- apply(stable.state_aFrF, 1, sum)
# 
# v.dot.w <- sum(hz*ht*stable.state_aFrF*repro.val_aFrF)
# sens <- outer(repro.val_aFrF,stable.state_aFrF)/v.dot.w
# elas <- sens*Kvals/lam.stable
# 
# ### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
# total.elas <- hz*ht*apply(elas,c(3,4),sum) 
# 
# ### Checks
# cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
# cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 
# 
# ### Plots
# 
# ## stable size distribution
# 
# plot(yz
#      , (stable.state_aFrF.vector/sum(stable.state_aFrF.vector))*10
#      , xlab = "Operculum length, mm"
#      , ylab = "Probability"
#      , type = "l"
#      , ylim = c(0, 1)
# )
# 
# # # code following IPM book (I modified this to make it a probability)
# # plot(yz
# #      , stable.state_aFrF.vector
# #      , xlab = "Size x"
# #      , ylab = "Frequency"
# #      , type = "l"
# # )
# 
# # # stable size dist with small ones excluded
# # plot(yz
# #      , (stable.state_aFrF.vector/sum(stable.state_aFrF.vector))*10
# #      , xlab = "Operculum length, mm"
# #      , ylab = "Probability"
# #      , type = "l"
# #      , ylim = c(0, 0.1)
# #      , xlim = c(2.5, 10)
# # )
# 
# 
# ## stable size and crowding distribution
# image(yz
#       , yt
#       , stable.state_aFrF
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , stable.state_aFrF
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## reproductive value
# image(yz
#       , yt
#       , repro.val_aFrF
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , repro.val_aFrF
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## total elasticity
# # image(yt
# #       , yz
# #       , t(total.elas)
# #       , col = grey(seq(0.5, 1, length=100))
# #       , xlab = "Crowding"
# #       , ylab = "Operculum length, mm"
# # )
# # contour(yt
# #         , yz
# #         , t(total.elas)
# #         , add = TRUE
# #         , nlevels = 6
# #         , labcex = 0.8
# # )
# 
# image(yz
#       , yt
#       , total.elas
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , total.elas
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

# # ---- *slow, ~ 34 hrs* IPM size touch - aFrF: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
# 
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
# 
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
# 
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
# 
#   ### Store the estimated parameters
# 
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
# 
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aUrU ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aUrU: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 1
                       , shape2 = 1
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
                      , shape1 = 1
                      , shape2 = 1
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aUrU: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aUrU: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aUrU <- matrix(v, mz, mt)
stable.state_aUrU <- matrix(w, mz, mt) 
stable.state_aUrU.vector <- apply(stable.state_aUrU, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aUrU*repro.val_aUrU)
sens <- outer(repro.val_aUrU,stable.state_aUrU)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aUrU.vector/sum(stable.state_aUrU.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aUrU.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aUrU
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aUrU
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aUrU
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aUrU
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aUrU: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
#   
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aLrU ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aLrU: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 2
                       , shape2 = 50
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
                      , shape1 = 1
                      , shape2 = 1
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aLrU: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aLrU: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aLrU <- matrix(v, mz, mt)
stable.state_aLrU <- matrix(w, mz, mt) 
stable.state_aLrU.vector <- apply(stable.state_aLrU, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aLrU*repro.val_aLrU)
sens <- outer(repro.val_aLrU,stable.state_aLrU)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aLrU.vector/sum(stable.state_aLrU.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aLrU.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aLrU
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aLrU
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aLrU
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aLrU
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aLrU: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
# 
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aHrU ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aHrU: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 50
                       , shape2 = 2
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
                      , shape1 = 1
                      , shape2 = 1
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aHrU: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aHrU: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aHrU <- matrix(v, mz, mt)
stable.state_aHrU <- matrix(w, mz, mt) 
stable.state_aHrU.vector <- apply(stable.state_aHrU, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aHrU*repro.val_aHrU)
sens <- outer(repro.val_aHrU,stable.state_aHrU)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aHrU.vector/sum(stable.state_aHrU.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aHrU.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aHrU
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aHrU
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aHrU
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aHrU
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aHrU: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
# 
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aUrL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aUrL: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 1
                       , shape2 = 1
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with low crowd
                      , shape1 = 2
                      , shape2 = 50
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aUrL: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aUrL: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aUrL <- matrix(v, mz, mt)
stable.state_aUrL <- matrix(w, mz, mt) 
stable.state_aUrL.vector <- apply(stable.state_aUrL, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aUrL*repro.val_aUrL)
sens <- outer(repro.val_aUrL,stable.state_aUrL)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aUrL.vector/sum(stable.state_aUrL.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aUrL.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aUrL
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aUrL
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aUrL
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aUrL
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aUrL: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
# 
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aUrH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aUrH: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 1
                       , shape2 = 1
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with low crowd
                      , shape1 = 50
                      , shape2 = 2
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aUrH: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aUrH: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aUrH <- matrix(v, mz, mt)
stable.state_aUrH <- matrix(w, mz, mt) 
stable.state_aUrH.vector <- apply(stable.state_aUrH, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aUrH*repro.val_aUrH)
sens <- outer(repro.val_aUrH,stable.state_aUrH)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aUrH.vector/sum(stable.state_aUrH.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aUrH.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aUrH
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aUrH
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aUrH
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aUrH
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aUrH: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
# 
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aLrL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aLrL: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 2
                       , shape2 = 50
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with low crowd
                      , shape1 = 2
                      , shape2 = 50
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aLrL: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aLrL: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aLrL <- matrix(v, mz, mt)
stable.state_aLrL <- matrix(w, mz, mt) 
stable.state_aLrL.vector <- apply(stable.state_aLrL, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aLrL*repro.val_aLrL)
sens <- outer(repro.val_aLrL,stable.state_aLrL)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aLrL.vector/sum(stable.state_aLrL.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aLrL.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aLrL
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aLrL
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aLrL
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aLrL
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aLrL: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
# 
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aLrH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aLrH: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 2
                       , shape2 = 50
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with low crowd
                      , shape1 = 50
                      , shape2 = 2
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aLrH: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aLrH: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aLrH <- matrix(v, mz, mt)
stable.state_aLrH <- matrix(w, mz, mt) 
stable.state_aLrH.vector <- apply(stable.state_aLrH, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aLrH*repro.val_aLrH)
sens <- outer(repro.val_aLrH,stable.state_aLrH)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas_aLrH <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aLrH.vector/sum(stable.state_aLrH.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aLrH.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aLrH
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aLrH
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aLrH
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aLrH
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas_aLrH
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas_aLrH
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aLrH: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
# 
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aHrL ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aHrL: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 50
                       , shape2 = 2
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
                      , shape1 = 2
                      , shape2 = 50
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aHrL: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aHrL: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aHrL <- matrix(v, mz, mt)
stable.state_aHrL <- matrix(w, mz, mt) 
stable.state_aHrL.vector <- apply(stable.state_aHrL, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aHrL*repro.val_aHrL)
sens <- outer(repro.val_aHrL,stable.state_aHrL)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas_aHrL <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aHrL.vector/sum(stable.state_aHrL.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aHrL.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aHrL
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aHrL
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aHrL
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aHrL
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas_aHrL
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas_aHrL
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aHrL: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
#   
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~ aHrH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ---- IPM size touch - aHrH: define kernel components ----

# IPM book pg. 24-25

### survival: s(z)

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

### growth: G(z', z, t)

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
  
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}
### touch dynamics: T(t', t)

T_t1t_st <- function(t1, t, m.par_st){
  
  p.den.touch <- dbeta(t1
                       , shape1 = 50
                       , shape2 = 2
  )
  
  ## output
  return(p.den.touch)
}

### prob reproduce: p_bz

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}


### fecundity: b_z

b_z <- function(z){
  # from Strathmann
  
  p.bz <- ifelse(z <= Uz1
                 , 0.147 * (10 * z)^2.74
                 , 0.147 * (10 * Uz1)^2.74
  )
  
  return(p.bz)
}

### recruit size dist: C_0(z')

C_0z1 <- function(z1, m.par_st){
  mu <- m.par_st['rcsz.int']
  sig <- m.par_st['rcsz.sd']
  p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
  return(p.den.rcsz)
}

### recruit touch dist: C_0(t')

C_0t1 <- function(t1, m.par_st){
  
  p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
                      , shape1 = 50
                      , shape2 = 2
  )
  
  # shp1 <- m.par_st['rcst.s1']
  # shp2 <- m.par_st['rcst.s2']
  # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
  
  
  return(p.den.rcst)
}


### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)

P_z1z <- function(z1, z, t1, t, m.par_st){
  return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
  )
}

### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')

F_z1z <- function(z1, z, t1, t, m.par_st){
  return(
    unname(
      s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
    )
  )
}


### define K kernel
k_st <- function(z1, t1,  z, t, m.par_st){
  P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
}

# ---- *slow, ~ 2 min* IPM size touch - aHrH: numerical implementation ----

# following IPM book pg. 162-164 and SizeQualityExample.R

### compute meshpoints
mz <- 100 # number of size mesh points 
mt <- 50  # number of touch mesh points
Lz <- 0.1 # size lower bound
Uz<- 10  # size upper bound
Lt <- 0 # touch lower bound
Ut <- 1 # touch upper bound
hz <- (Uz-Lz)/mz # size mesh point breaks
yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
ht <- (Ut-Lt)/mt # touch mesh point breaks
yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points

### compute 4D kernel and 2D iteration matrix



# Function eta to put kernel values in their proper place in A
eta_ij <- function(i, j, mz) {(j-1)*mz+i}

# matrix whose (i,j) entry is eta(ij)
Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)

# code modified from IPM book, Ch 6, SizeQualityExample.R
A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
for(i in 1:mz){
  for(j in 1:mt){
    for(k in 1:mt){
      kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
      A[Eta[,k],Eta[i,j]]=kvals
      Kvals[,k,i,j]=kvals
      
    }}
  cat(i,"\n");
}
A<-hz*ht*A

# # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# 
# Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# 
# A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# dim(A) <- c(mz * mt, mz * mt)
# A <- hz*ht*A

### calculate Lambda, w, and v
out <- domEig(A) # returns lambda and a vector proportional to w
out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
lam.stable <- out$lambda
lam.stable.t <- out2$lambda
w <- Re(matrix(out$w, mz, mt))
w <- w/(hz*ht*sum(w))
v <- Re(matrix(out2$w, mz, mt))
v <- v/sum(v) 

# ---- IPM size touch - aHrH: perturbation analysis ----

### Compute elasticity matrix 
repro.val_aHrH <- matrix(v, mz, mt)
stable.state_aHrH <- matrix(w, mz, mt) 
stable.state_aHrH.vector <- apply(stable.state_aHrH, 1, sum)
v.dot.w <- sum(hz*ht*stable.state_aHrH*repro.val_aHrH)
sens <- outer(repro.val_aHrH,stable.state_aHrH)/v.dot.w
elas <- sens*Kvals/lam.stable

### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
total.elas <- hz*ht*apply(elas,c(3,4),sum) 

### Checks
cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 

### Plots

## stable size distribution

plot(yz
     , (stable.state_aHrH.vector/sum(stable.state_aHrH.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)

# # code following IPM book (I modified this to make it a probability)
# plot(yz
#      , stable.state_aHrH.vector
#      , xlab = "Size x"
#      , ylab = "Frequency"
#      , type = "l"
# )




## stable size and crowding distribution
image(yz
      , yt
      , stable.state_aHrH
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , stable.state_aHrH
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## reproductive value
image(yz
      , yt
      , repro.val_aHrH
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , repro.val_aHrH
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## total elasticity
# image(yt
#       , yz
#       , t(total.elas)
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Crowding"
#       , ylab = "Operculum length, mm"
# )
# contour(yt
#         , yz
#         , t(total.elas)
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

image(yz
      , yt
      , total.elas
      , col = grey(seq(0.5, 1, length=100))
      , xlab = "Operculum length, mm"
      , ylab = "Crowding"
)
contour(yz
        , yt
        , total.elas
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

# # ---- *slow, ~ 34 hrs* IPM size touch - aHrH: lambda CI ----
# # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# 
# # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets.
# 
# ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# boot.lam.st <- function(dataset, sample.index) {
# 
#   ### extract the data used to make this fit
#   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# 
#   ### fit the functions
# 
#   ## survival
#   mod.Surv <- glm(Surv ~ z
#                   , family = binomial
#                   , data = boot.data
#   )
# 
#   ## growth
#   
#   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
#   
#   #z.growth <- boot.data.growth$z
#   mod.Grow_st <- lmer(z1 ~ z * touch_pct
#                       + (1 | subsite)
#                       , data = boot.data.growth
#   )
#   
#   
#   ## recruit size distribution
#   mod.Rcsz <- lm(z ~ 1
#                  , data = rec
#   )
#   
#   ### Store the estimated parameters
#   
#   m.par_st <- c(
#     surv = coef(mod.Surv)
#     , grow.int = summary(mod.Grow_st)$coefficients[1]
#     , grow.z = summary(mod.Grow_st)$coefficients[2]
#     , grow.t = summary(mod.Grow_st)$coefficients[3]
#     , grow.zt = summary(mod.Grow_st)$coefficients[4]
#     , grow.sd = summary(mod.Grow_st)$sigma
#     , rcsz = coef(mod.Rcsz)
#     , rcsz.sd = summary(mod.Rcsz)$sigma
#     , rcst = coef(betafit_t)
#     , p.r = p.r.est
#     , repr = coef(reprodyn.1b)
#     #, U1 = U1
#   )
# 
#   names(m.par_st) <- c(
#     'surv.int'
#     , 'surv.z'
#     , 'grow.int'
#     , 'grow.z'
#     , 'grow.t'
#     , 'grow.zt'
#     , 'grow.sd'
#     , 'rcsz.int'
#     , 'rcsz.sd'
#     , 'rcst.s1'
#     , 'rcst.s2'
#     , 'p.r'
#     , 'repr.int' # model is for volume, NOT operc length
#     , 'repr.t' # model is for volume, NOT operc length
#     , 'repr.v' # v here b/c it's volume, NOT operc length
#     #, 'U1'
#   )
# 
#   ### implement IPM and calculate lambda
#   # following IPM book pg. 162-164 and SizeQualityExample.R
# 
#   ## compute meshpoints
#   mz <- 100 # number of size mesh points
#   mt <- 50  # number of touch mesh points
#   Lz <- 0.1 # size lower bound
#   Uz<- 10  # size upper bound
#   Lt <- 0 # touch lower bound
#   Ut <- 1 # touch upper bound
#   hz <- (Uz-Lz)/mz # size mesh point breaks
#   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
#   ht <- (Ut-Lt)/mt # touch mesh point breaks
#   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
#   ## compute 4D kernel and 2D iteration matrix
#   # Function eta to put kernel values in their proper place in A
#   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
#   # matrix whose (i,j) entry is eta(ij)
#   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
#   # code modified from IPM book, Ch 6, SizeQualityExample.R
#   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
#   for(i in 1:mz){
#     for(j in 1:mt){
#       for(k in 1:mt){
#         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#         A[Eta[,k],Eta[i,j]]=kvals
#         Kvals[,k,i,j]=kvals
# 
#       }}
#     cat(i,"\n");
#   }
#   A<-hz*ht*A
# 
#   out <- domEig(A) # returns lambda and a vector proportional to w
#   lam.boot <- out$lambda
#   cat(lam.boot, "\n")
#   return(lam.boot)
# 
# }
# 
# ### do the bootstrap (code takes ~ X min to run)
# starttime <- Sys.time()
# boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
#                     R = 1000)
# endtime <- Sys.time()
# 
# boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# 
# 
# 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# ~Graveyard~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ aFrU ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# 
# # ---- IPM size touch - aFrU: define kernel components ----
# 
# # IPM book pg. 24-25
# 
# ### survival: s(z)
# 
# s_z_st <- function(z, m.par_st){
#   
#   ## linear predictor
#   
#   # below size ceiling
#   linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
#   # above size ceiling
#   linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
#   
#   # impose ceiling and backtransform from logit
#   p <- ifelse(z <= Uz1
#               , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
#               , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
#   )
#   
#   ## output
#   return( unname(p) )
#   
#   # ## back transform from logit
#   # p <- 1/(1+exp(-linear.p))
# }
# 
# ### growth: G(z', z, t)
# 
# G_z1zt_st <- function(z1, z, t, m.par_st){
#   
#   ## based on growth regression coefficients, define mu and sigma of growth PDF
#   # below size ceiling
#   mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
#   # above size ceiling
#   mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
#   
#   sig <- m.par_st['grow.sd']
#   
#   ## calculate growth PDF
#   # have to loop b/c it's a PDF, not a probability
#   p.den.grow <- c()
#   for(q in 1:length(z1)){
#     p.den.grow[q] <- ifelse(z <= Uz1
#                             , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
#                             , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
#     )
#   }
#   
#   ## output
#   return(p.den.grow)
#   
# }
# ### touch dynamics: T(t', t)
# 
# T_t1t_st <- function(t1, t, m.par_st){
#   
#   p.den.touch <- dbeta(t1
#                        , shape1 = m.par_st['rcst.s1']
#                        , shape2 = m.par_st['rcst.s2']
#   )
#   
#   ## output
#   return(p.den.touch)
# }
# 
# ### prob reproduce: p_bz
# 
# p_bz_st <- function(z, t, m.par_st){
#   ### convert z to v
#   v <- ifelse(z <= Uz1
#               , z * 131.18 - 152.52 # below ceiling
#               , Uz1 * 131.18 - 152.52 # above ceiling
#   )
#   
#   # linear predictor
#   linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
#   
#   # back transform from logit
#   p <- 1/(1+exp(-linear.p))
#   
#   # output
#   return( unname(p) )
# }
# 
# 
# ### fecundity: b_z
# 
# b_z <- function(z){
#   # from Strathmann
#   
#   p.bz <- ifelse(z <= Uz1
#                  , 0.147 * (10 * z)^2.74
#                  , 0.147 * (10 * Uz1)^2.74
#   )
#   
#   return(p.bz)
# }
# 
# ### recruit size dist: C_0(z')
# 
# C_0z1 <- function(z1, m.par_st){
#   mu <- m.par_st['rcsz.int']
#   sig <- m.par_st['rcsz.sd']
#   p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
#   return(p.den.rcsz)
# }
# 
# ### recruit touch dist: C_0(t')
# 
# C_0t1 <- function(t1, m.par_st){
#   
#   p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
#                       , shape1 = 1
#                       , shape2 = 1
#   )
#   
#   # shp1 <- m.par_st['rcst.s1']
#   # shp2 <- m.par_st['rcst.s2']
#   # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
#   
#   
#   return(p.den.rcst)
# }
# 
# 
# ### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)
# 
# P_z1z <- function(z1, z, t1, t, m.par_st){
#   return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
#   )
# }
# 
# ### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')
# 
# F_z1z <- function(z1, z, t1, t, m.par_st){
#   return(
#     unname(
#       s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
#     )
#   )
# }
# 
# 
# ### define K kernel
# k_st <- function(z1, t1,  z, t, m.par_st){
#   P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
# }
# 
# # ---- *slow, ~ 2 min* IPM size touch - aFrU: numerical implementation ----
# 
# # following IPM book pg. 162-164 and SizeQualityExample.R
# 
# ### compute meshpoints
# mz <- 100 # number of size mesh points 
# mt <- 50  # number of touch mesh points
# Lz <- 0.1 # size lower bound
# Uz<- 10  # size upper bound
# Lt <- 0 # touch lower bound
# Ut <- 1 # touch upper bound
# hz <- (Uz-Lz)/mz # size mesh point breaks
# yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
# ht <- (Ut-Lt)/mt # touch mesh point breaks
# yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
# ### compute 4D kernel and 2D iteration matrix
# 
# 
# 
# # Function eta to put kernel values in their proper place in A
# eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
# # matrix whose (i,j) entry is eta(ij)
# Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
# # code modified from IPM book, Ch 6, SizeQualityExample.R
# A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
# for(i in 1:mz){
#   for(j in 1:mt){
#     for(k in 1:mt){
#       kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#       A[Eta[,k],Eta[i,j]]=kvals
#       Kvals[,k,i,j]=kvals
#       
#     }}
#   cat(i,"\n");
# }
# A<-hz*ht*A
# 
# # # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# # 
# # Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# # 
# # A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# # dim(A) <- c(mz * mt, mz * mt)
# # A <- hz*ht*A
# 
# ### calculate Lambda, w, and v
# out <- domEig(A) # returns lambda and a vector proportional to w
# out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
# lam.stable <- out$lambda
# lam.stable.t <- out2$lambda
# w <- Re(matrix(out$w, mz, mt))
# w <- w/(hz*ht*sum(w))
# v <- Re(matrix(out2$w, mz, mt))
# v <- v/sum(v) 
# 
# # ---- IPM size touch - aFrU: perturbation analysis ----
# 
# ### Compute elasticity matrix 
# repro.val_aFrU <- matrix(v, mz, mt)
# stable.state_aFrU <- matrix(w, mz, mt) 
# stable.state_aFrU.vector <- apply(stable.state_aFrU, 1, sum)
# 
# v.dot.w <- sum(hz*ht*stable.state_aFrU*repro.val_aFrU)
# sens <- outer(repro.val_aFrU,stable.state_aFrU)/v.dot.w
# elas <- sens*Kvals/lam.stable
# 
# ### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
# total.elas <- hz*ht*apply(elas,c(3,4),sum) 
# 
# ### Checks
# cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
# cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 
# 
# ### Plots
# 
# ## stable size distribution
# 
# plot(yz
#      , (stable.state_aFrU.vector/sum(stable.state_aFrU.vector))*10
#      , xlab = "Operculum length, mm"
#      , ylab = "Probability"
#      , type = "l"
#      , ylim = c(0, 1)
# )
# 
# # # code following IPM book (I modified this to make it a probability)
# # plot(yz
# #      , stable.state_aFrU.vector
# #      , xlab = "Size x"
# #      , ylab = "Frequency"
# #      , type = "l"
# # )
# 
# # # stable size dist with small ones excluded
# # plot(yz
# #      , (stable.state_aFrU.vector/sum(stable.state_aFrU.vector))*10
# #      , xlab = "Operculum length, mm"
# #      , ylab = "Probability"
# #      , type = "l"
# #      , ylim = c(0, 0.1)
# #      , xlim = c(2.5, 10)
# # )
# 
# 
# ## stable size and crowding distribution
# image(yz
#       , yt
#       , stable.state_aFrU
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , stable.state_aFrU
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## reproductive value
# image(yz
#       , yt
#       , repro.val_aFrU
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , repro.val_aFrU
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## total elasticity
# # image(yt
# #       , yz
# #       , t(total.elas)
# #       , col = grey(seq(0.5, 1, length=100))
# #       , xlab = "Crowding"
# #       , ylab = "Operculum length, mm"
# # )
# # contour(yt
# #         , yz
# #         , t(total.elas)
# #         , add = TRUE
# #         , nlevels = 6
# #         , labcex = 0.8
# # )
# 
# image(yz
#       , yt
#       , total.elas
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , total.elas
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# # # ---- *slow, ~ 34 hrs* IPM size touch - aFrU: lambda CI ----
# # # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# # 
# # # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets. 
# # 
# # ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# # boot.lam.st <- function(dataset, sample.index) {
# #   
# #   ### extract the data used to make this fit
# #   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# #   
# #   ### fit the functions
# #   
# #   ## survival
# #   mod.Surv <- glm(Surv ~ z
# #                   , family = binomial
# #                   , data = boot.data
# #   )
# #   
# #   ## growth
# #   
# #   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
# #   
# #   #z.growth <- boot.data.growth$z
# #   mod.Grow_st <- lm(z1 ~ z * touch_pct
# #                     , data = boot.data.growth
# #                     
# #   )
# #   
# #   ## recruit size distribution
# #   mod.Rcsz <- lm(z ~ 1 
# #                  , data = rec
# #   )
# #   
# #   ### Store the estimated parameters
# #   
# #   m.par_st <- c(
# #     surv = coef(mod.Surv)
# #     , grow = coef(mod.Grow_st)
# #     , grow.sd = summary(mod.Grow_st)$sigma
# #     , rcsz = coef(mod.Rcsz)
# #     , rcsz.sd = summary(mod.Rcsz)$sigma
# #     , rcst = coef(betafit_t)
# #     , p.r = p.r.est
# #     , repr = coef(reprodyn.1b)
# #     #, U1 = U1
# #   )
# #   
# #   names(m.par_st) <- c(
# #     'surv.int'
# #     , 'surv.z'
# #     , 'grow.int'
# #     , 'grow.z'
# #     , 'grow.t'
# #     , 'grow.zt'
# #     , 'grow.sd'
# #     , 'rcsz.int'
# #     , 'rcsz.sd'
# #     , 'rcst.s1'
# #     , 'rcst.s2'
# #     , 'p.r'
# #     , 'repr.int' # model is for volume, NOT operc length
# #     , 'repr.t' # model is for volume, NOT operc length
# #     , 'repr.v' # v here b/c it's volume, NOT operc length
# #     #, 'U1'
# #   )
# #   
# #   ### implement IPM and calculate lambda
# #   # following IPM book pg. 162-164 and SizeQualityExample.R
# #   
# #   ## compute meshpoints
# #   mz <- 100 # number of size mesh points 
# #   mt <- 50  # number of touch mesh points
# #   Lz <- 0.1 # size lower bound
# #   Uz<- 10  # size upper bound
# #   Lt <- 0 # touch lower bound
# #   Ut <- 1 # touch upper bound
# #   hz <- (Uz-Lz)/mz # size mesh point breaks
# #   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
# #   ht <- (Ut-Lt)/mt # touch mesh point breaks
# #   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# #   
# #   ## compute 4D kernel and 2D iteration matrix
# #   # Function eta to put kernel values in their proper place in A
# #   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# #   
# #   # matrix whose (i,j) entry is eta(ij)
# #   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# #   
# #   # code modified from IPM book, Ch 6, SizeQualityExample.R
# #   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
# #   for(i in 1:mz){
# #     for(j in 1:mt){
# #       for(k in 1:mt){
# #         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
# #         A[Eta[,k],Eta[i,j]]=kvals
# #         Kvals[,k,i,j]=kvals
# #         
# #       }}
# #     cat(i,"\n");
# #   }
# #   A<-hz*ht*A
# #   
# #   out <- domEig(A) # returns lambda and a vector proportional to w
# #   lam.boot <- out$lambda
# #   cat(lam.boot, "\n")
# #   return(lam.boot)
# #   
# # }
# # 
# # ### do the bootstrap (code takes ~ X min to run)
# # starttime <- Sys.time()
# # boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
# #                     R = 1000)
# # endtime <- Sys.time()
# # 
# # boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# # 
# # 
# # 
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ aUrM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# # ---- IPM size touch - aUrM: define kernel components ----
# 
# # IPM book pg. 24-25
# 
# ### survival: s(z)
# 
# s_z_st <- function(z, m.par_st){
#   
#   ## linear predictor
#   
#   # below size ceiling
#   linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
#   # above size ceiling
#   linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
#   
#   # impose ceiling and backtransform from logit
#   p <- ifelse(z <= Uz1
#               , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
#               , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
#   )
#   
#   ## output
#   return( unname(p) )
#   
#   # ## back transform from logit
#   # p <- 1/(1+exp(-linear.p))
# }
# 
# ### growth: G(z', z, t)
# 
# G_z1zt_st <- function(z1, z, t, m.par_st){
#   
#   ## based on growth regression coefficients, define mu and sigma of growth PDF
#   # below size ceiling
#   mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
#   # above size ceiling
#   mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
#   
#   sig <- m.par_st['grow.sd']
#   
#   ## calculate growth PDF
#   # have to loop b/c it's a PDF, not a probability
#   p.den.grow <- c()
#   for(q in 1:length(z1)){
#     p.den.grow[q] <- ifelse(z <= Uz1
#                             , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
#                             , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
#     )
#   }
#   
#   ## output
#   return(p.den.grow)
#   
# }
# ### touch dynamics: T(t', t)
# 
# T_t1t_st <- function(t1, t, m.par_st){
#   
#   p.den.touch <- dbeta(t1
#                        , shape1 = 1
#                        , shape2 = 1
#   )
#   
#   ## output
#   return(p.den.touch)
# }
# 
# ### prob reproduce: p_bz
# 
# p_bz_st <- function(z, t, m.par_st){
#   ### convert z to v
#   v <- ifelse(z <= Uz1
#               , z * 131.18 - 152.52 # below ceiling
#               , Uz1 * 131.18 - 152.52 # above ceiling
#   )
#   
#   # linear predictor
#   linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
#   
#   # back transform from logit
#   p <- 1/(1+exp(-linear.p))
#   
#   # output
#   return( unname(p) )
# }
# 
# 
# ### fecundity: b_z
# 
# b_z <- function(z){
#   # from Strathmann
#   
#   p.bz <- ifelse(z <= Uz1
#                  , 0.147 * (10 * z)^2.74
#                  , 0.147 * (10 * Uz1)^2.74
#   )
#   
#   return(p.bz)
# }
# 
# ### recruit size dist: C_0(z')
# 
# C_0z1 <- function(z1, m.par_st){
#   mu <- m.par_st['rcsz.int']
#   sig <- m.par_st['rcsz.sd']
#   p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
#   return(p.den.rcsz)
# }
# 
# ### recruit touch dist: C_0(t')
# 
# C_0t1 <- function(t1, m.par_st){
#   
#   p.den.rcst <- dbeta(t1 # recruits start with low crowd
#                       , shape1 = 300
#                       , shape2 = 300
#   )
#   
#   # shp1 <- m.par_st['rcst.s1']
#   # shp2 <- m.par_st['rcst.s2']
#   # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
#   
#   
#   return(p.den.rcst)
# }
# 
# 
# ### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)
# 
# P_z1z <- function(z1, z, t1, t, m.par_st){
#   return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
#   )
# }
# 
# ### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')
# 
# F_z1z <- function(z1, z, t1, t, m.par_st){
#   return(
#     unname(
#       s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
#     )
#   )
# }
# 
# 
# ### define K kernel
# k_st <- function(z1, t1,  z, t, m.par_st){
#   P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
# }
# 
# # ---- *slow, ~ 2 min* IPM size touch - aUrM: numerical implementation ----
# 
# # following IPM book pg. 162-164 and SizeQualityExample.R
# 
# ### compute meshpoints
# mz <- 100 # number of size mesh points 
# mt <- 50  # number of touch mesh points
# Lz <- 0.1 # size lower bound
# Uz<- 10  # size upper bound
# Lt <- 0 # touch lower bound
# Ut <- 1 # touch upper bound
# hz <- (Uz-Lz)/mz # size mesh point breaks
# yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
# ht <- (Ut-Lt)/mt # touch mesh point breaks
# yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
# ### compute 4D kernel and 2D iteration matrix
# 
# 
# 
# # Function eta to put kernel values in their proper place in A
# eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
# # matrix whose (i,j) entry is eta(ij)
# Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
# # code modified from IPM book, Ch 6, SizeQualityExample.R
# A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
# for(i in 1:mz){
#   for(j in 1:mt){
#     for(k in 1:mt){
#       kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#       A[Eta[,k],Eta[i,j]]=kvals
#       Kvals[,k,i,j]=kvals
#       
#     }}
#   cat(i,"\n");
# }
# A<-hz*ht*A
# 
# # # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# # 
# # Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# # 
# # A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# # dim(A) <- c(mz * mt, mz * mt)
# # A <- hz*ht*A
# 
# ### calculate Lambda, w, and v
# out <- domEig(A) # returns lambda and a vector proportional to w
# out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
# lam.stable <- out$lambda
# lam.stable.t <- out2$lambda
# w <- Re(matrix(out$w, mz, mt))
# w <- w/(hz*ht*sum(w))
# v <- Re(matrix(out2$w, mz, mt))
# v <- v/sum(v) 
# 
# # ---- IPM size touch - aUrM: perturbation analysis ----
# 
# ### Compute elasticity matrix 
# repro.val_aUrM <- matrix(v, mz, mt)
# stable.state_aUrM <- matrix(w, mz, mt) 
# stable.state_aUrM.vector <- apply(stable.state_aUrM, 1, sum)
# v.dot.w <- sum(hz*ht*stable.state_aUrM*repro.val_aUrM)
# sens <- outer(repro.val_aUrM,stable.state_aUrM)/v.dot.w
# elas <- sens*Kvals/lam.stable
# 
# ### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
# total.elas <- hz*ht*apply(elas,c(3,4),sum) 
# 
# ### Checks
# cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
# cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 
# 
# ### Plots
# 
# ## stable size distribution
# 
# plot(yz
#      , (stable.state_aUrM.vector/sum(stable.state_aUrM.vector))*10
#      , xlab = "Operculum length, mm"
#      , ylab = "Frequency"
#      , type = "l"
#      , ylim = c(0, 1)
# )
# 
# # # code following IPM book (I modified this to make it a probability)
# # plot(yz
# #      , stable.state_aUrM.vector
# #      , xlab = "Size x"
# #      , ylab = "Frequency"
# #      , type = "l"
# # )
# 
# 
# 
# 
# ## stable size and crowding distribution
# image(yz
#       , yt
#       , stable.state_aUrM
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , stable.state_aUrM
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## reproductive value
# image(yz
#       , yt
#       , repro.val_aUrM
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , repro.val_aUrM
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## total elasticity
# # image(yt
# #       , yz
# #       , t(total.elas)
# #       , col = grey(seq(0.5, 1, length=100))
# #       , xlab = "Crowding"
# #       , ylab = "Operculum length, mm"
# # )
# # contour(yt
# #         , yz
# #         , t(total.elas)
# #         , add = TRUE
# #         , nlevels = 6
# #         , labcex = 0.8
# # )
# 
# image(yz
#       , yt
#       , total.elas
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , total.elas
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# # # ---- *slow, ~ 34 hrs* IPM size touch - aUrM: lambda CI ----
# # # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# # 
# # # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets. 
# # 
# # ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# # boot.lam.st <- function(dataset, sample.index) {
# #   
# #   ### extract the data used to make this fit
# #   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# #   
# #   ### fit the functions
# #   
# #   ## survival
# #   mod.Surv <- glm(Surv ~ z
# #                   , family = binomial
# #                   , data = boot.data
# #   )
# #   
# #   ## growth
# #   
# #   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
# #   
# #   #z.growth <- boot.data.growth$z
# #   mod.Grow_st <- lm(z1 ~ z * touch_pct
# #                     , data = boot.data.growth
# #                     
# #   )
# #   
# #   ## recruit size distribution
# #   mod.Rcsz <- lm(z ~ 1 
# #                  , data = rec
# #   )
# #   
# #   ### Store the estimated parameters
# #   
# #   m.par_st <- c(
# #     surv = coef(mod.Surv)
# #     , grow = coef(mod.Grow_st)
# #     , grow.sd = summary(mod.Grow_st)$sigma
# #     , rcsz = coef(mod.Rcsz)
# #     , rcsz.sd = summary(mod.Rcsz)$sigma
# #     , rcst = coef(betafit_t)
# #     , p.r = p.r.est
# #     , repr = coef(reprodyn.1b)
# #     #, U1 = U1
# #   )
# #   
# #   names(m.par_st) <- c(
# #     'surv.int'
# #     , 'surv.z'
# #     , 'grow.int'
# #     , 'grow.z'
# #     , 'grow.t'
# #     , 'grow.zt'
# #     , 'grow.sd'
# #     , 'rcsz.int'
# #     , 'rcsz.sd'
# #     , 'rcst.s1'
# #     , 'rcst.s2'
# #     , 'p.r'
# #     , 'repr.int' # model is for volume, NOT operc length
# #     , 'repr.t' # model is for volume, NOT operc length
# #     , 'repr.v' # v here b/c it's volume, NOT operc length
# #     #, 'U1'
# #   )
# #   
# #   ### implement IPM and calculate lambda
# #   # following IPM book pg. 162-164 and SizeQualityExample.R
# #   
# #   ## compute meshpoints
# #   mz <- 100 # number of size mesh points 
# #   mt <- 50  # number of touch mesh points
# #   Lz <- 0.1 # size lower bound
# #   Uz<- 10  # size upper bound
# #   Lt <- 0 # touch lower bound
# #   Ut <- 1 # touch upper bound
# #   hz <- (Uz-Lz)/mz # size mesh point breaks
# #   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
# #   ht <- (Ut-Lt)/mt # touch mesh point breaks
# #   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# #   
# #   ## compute 4D kernel and 2D iteration matrix
# #   # Function eta to put kernel values in their proper place in A
# #   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# #   
# #   # matrix whose (i,j) entry is eta(ij)
# #   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# #   
# #   # code modified from IPM book, Ch 6, SizeQualityExample.R
# #   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
# #   for(i in 1:mz){
# #     for(j in 1:mt){
# #       for(k in 1:mt){
# #         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
# #         A[Eta[,k],Eta[i,j]]=kvals
# #         Kvals[,k,i,j]=kvals
# #         
# #       }}
# #     cat(i,"\n");
# #   }
# #   A<-hz*ht*A
# #   
# #   out <- domEig(A) # returns lambda and a vector proportional to w
# #   lam.boot <- out$lambda
# #   cat(lam.boot, "\n")
# #   return(lam.boot)
# #   
# # }
# # 
# # ### do the bootstrap (code takes ~ X min to run)
# # starttime <- Sys.time()
# # boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
# #                     R = 1000)
# # endtime <- Sys.time()
# # 
# # boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# # 
# # 

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ aMrM ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# # ---- IPM size touch - aMrM: define kernel components ----
# 
# # IPM book pg. 24-25
# 
# ### survival: s(z)
# 
# s_z_st <- function(z, m.par_st){
#   
#   ## linear predictor
#   
#   # below size ceiling
#   linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
#   # above size ceiling
#   linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
#   
#   # impose ceiling and backtransform from logit
#   p <- ifelse(z <= Uz1
#               , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
#               , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
#   )
#   
#   ## output
#   return( unname(p) )
#   
#   # ## back transform from logit
#   # p <- 1/(1+exp(-linear.p))
# }
# 
# ### growth: G(z', z, t)
# 
# G_z1zt_st <- function(z1, z, t, m.par_st){
#   
#   ## based on growth regression coefficients, define mu and sigma of growth PDF
#   # below size ceiling
#   mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
#   # above size ceiling
#   mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
#   
#   sig <- m.par_st['grow.sd']
#   
#   ## calculate growth PDF
#   # have to loop b/c it's a PDF, not a probability
#   p.den.grow <- c()
#   for(q in 1:length(z1)){
#     p.den.grow[q] <- ifelse(z <= Uz1
#                             , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
#                             , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
#     )
#   }
#   
#   ## output
#   return(p.den.grow)
#   
# }
# ### touch dynamics: T(t', t)
# 
# T_t1t_st <- function(t1, t, m.par_st){
#   
#   p.den.touch <- dbeta(t1
#                        , shape1 = 300
#                        , shape2 = 300
#   )
#   
#   ## output
#   return(p.den.touch)
# }
# 
# ### prob reproduce: p_bz
# 
# p_bz_st <- function(z, t, m.par_st){
#   ### convert z to v
#   v <- ifelse(z <= Uz1
#               , z * 131.18 - 152.52 # below ceiling
#               , Uz1 * 131.18 - 152.52 # above ceiling
#   )
#   
#   # linear predictor
#   linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
#   
#   # back transform from logit
#   p <- 1/(1+exp(-linear.p))
#   
#   # output
#   return( unname(p) )
# }
# 
# 
# ### fecundity: b_z
# 
# b_z <- function(z){
#   # from Strathmann
#   
#   p.bz <- ifelse(z <= Uz1
#                  , 0.147 * (10 * z)^2.74
#                  , 0.147 * (10 * Uz1)^2.74
#   )
#   
#   return(p.bz)
# }
# 
# ### recruit size dist: C_0(z')
# 
# C_0z1 <- function(z1, m.par_st){
#   mu <- m.par_st['rcsz.int']
#   sig <- m.par_st['rcsz.sd']
#   p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
#   return(p.den.rcsz)
# }
# 
# ### recruit touch dist: C_0(t')
# 
# C_0t1 <- function(t1, m.par_st){
#   
#   p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
#                       , shape1 = 300
#                       , shape2 = 300
#   )
#   
#   # shp1 <- m.par_st['rcst.s1']
#   # shp2 <- m.par_st['rcst.s2']
#   # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
#   
#   
#   return(p.den.rcst)
# }
# 
# 
# ### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)
# 
# P_z1z <- function(z1, z, t1, t, m.par_st){
#   return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
#   )
# }
# 
# ### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')
# 
# F_z1z <- function(z1, z, t1, t, m.par_st){
#   return(
#     unname(
#       s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
#     )
#   )
# }
# 
# 
# ### define K kernel
# k_st <- function(z1, t1,  z, t, m.par_st){
#   P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
# }
# 
# # ---- *slow, ~ 2 min* IPM size touch - aMrM: numerical implementation ----
# 
# # following IPM book pg. 162-164 and SizeQualityExample.R
# 
# ### compute meshpoints
# mz <- 100 # number of size mesh points 
# mt <- 50  # number of touch mesh points
# Lz <- 0.1 # size lower bound
# Uz<- 10  # size upper bound
# Lt <- 0 # touch lower bound
# Ut <- 1 # touch upper bound
# hz <- (Uz-Lz)/mz # size mesh point breaks
# yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
# ht <- (Ut-Lt)/mt # touch mesh point breaks
# yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
# ### compute 4D kernel and 2D iteration matrix
# 
# 
# 
# # Function eta to put kernel values in their proper place in A
# eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
# # matrix whose (i,j) entry is eta(ij)
# Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
# # code modified from IPM book, Ch 6, SizeQualityExample.R
# A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
# for(i in 1:mz){
#   for(j in 1:mt){
#     for(k in 1:mt){
#       kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#       A[Eta[,k],Eta[i,j]]=kvals
#       Kvals[,k,i,j]=kvals
#       
#     }}
#   cat(i,"\n");
# }
# A<-hz*ht*A
# 
# # # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# # 
# # Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# # 
# # A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# # dim(A) <- c(mz * mt, mz * mt)
# # A <- hz*ht*A
# 
# ### calculate Lambda, w, and v
# out <- domEig(A) # returns lambda and a vector proportional to w
# out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
# lam.stable <- out$lambda
# lam.stable.t <- out2$lambda
# w <- Re(matrix(out$w, mz, mt))
# w <- w/(hz*ht*sum(w))
# v <- Re(matrix(out2$w, mz, mt))
# v <- v/sum(v) 
# 
# # ---- IPM size touch - aMrM: perturbation analysis ----
# 
# ### Compute elasticity matrix 
# repro.val_aMrM <- matrix(v, mz, mt)
# stable.state_aMrM <- matrix(w, mz, mt) 
# stable.state_aMrM.vector <- apply(stable.state_aMrM, 1, sum)
# v.dot.w <- sum(hz*ht*stable.state_aMrM*repro.val_aMrM)
# sens <- outer(repro.val_aMrM,stable.state_aMrM)/v.dot.w
# elas <- sens*Kvals/lam.stable
# 
# ### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
# total.elas <- hz*ht*apply(elas,c(3,4),sum) 
# 
# ### Checks
# cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
# cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 
# 
# ### Plots
# 
# ## stable size distribution
# 
# plot(yz
#      , (stable.state_aMrM.vector/sum(stable.state_aMrM.vector))*10
#      , xlab = "Operculum length, mm"
#      , ylab = "Frequency"
#      , type = "l"
#      , ylim = c(0, 1)
# )
# 
# # # code following IPM book (I modified this to make it a probability)
# # plot(yz
# #      , stable.state_aMrM.vector
# #      , xlab = "Size x"
# #      , ylab = "Frequency"
# #      , type = "l"
# # )
# 
# 
# 
# 
# ## stable size and crowding distribution
# image(yz
#       , yt
#       , stable.state_aMrM
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , stable.state_aMrM
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## reproductive value
# image(yz
#       , yt
#       , repro.val_aMrM
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , repro.val_aMrM
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## total elasticity
# # image(yt
# #       , yz
# #       , t(total.elas)
# #       , col = grey(seq(0.5, 1, length=100))
# #       , xlab = "Crowding"
# #       , ylab = "Operculum length, mm"
# # )
# # contour(yt
# #         , yz
# #         , t(total.elas)
# #         , add = TRUE
# #         , nlevels = 6
# #         , labcex = 0.8
# # )
# 
# image(yz
#       , yt
#       , total.elas
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , total.elas
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# # # ---- *slow, ~ 34 hrs* IPM size touch - aMrM: lambda CI ----
# # # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# # 
# # # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets. 
# # 
# # ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# # boot.lam.st <- function(dataset, sample.index) {
# #   
# #   ### extract the data used to make this fit
# #   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# #   
# #   ### fit the functions
# #   
# #   ## survival
# #   mod.Surv <- glm(Surv ~ z
# #                   , family = binomial
# #                   , data = boot.data
# #   )
# #   
# #   ## growth
# #   
# #   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
# #   
# #   #z.growth <- boot.data.growth$z
# #   mod.Grow_st <- lm(z1 ~ z * touch_pct
# #                     , data = boot.data.growth
# #                     
# #   )
# #   
# #   ## recruit size distribution
# #   mod.Rcsz <- lm(z ~ 1 
# #                  , data = rec
# #   )
# #   
# #   ### Store the estimated parameters
# #   
# #   m.par_st <- c(
# #     surv = coef(mod.Surv)
# #     , grow = coef(mod.Grow_st)
# #     , grow.sd = summary(mod.Grow_st)$sigma
# #     , rcsz = coef(mod.Rcsz)
# #     , rcsz.sd = summary(mod.Rcsz)$sigma
# #     , rcst = coef(betafit_t)
# #     , p.r = p.r.est
# #     , repr = coef(reprodyn.1b)
# #     #, U1 = U1
# #   )
# #   
# #   names(m.par_st) <- c(
# #     'surv.int'
# #     , 'surv.z'
# #     , 'grow.int'
# #     , 'grow.z'
# #     , 'grow.t'
# #     , 'grow.zt'
# #     , 'grow.sd'
# #     , 'rcsz.int'
# #     , 'rcsz.sd'
# #     , 'rcst.s1'
# #     , 'rcst.s2'
# #     , 'p.r'
# #     , 'repr.int' # model is for volume, NOT operc length
# #     , 'repr.t' # model is for volume, NOT operc length
# #     , 'repr.v' # v here b/c it's volume, NOT operc length
# #     #, 'U1'
# #   )
# #   
# #   ### implement IPM and calculate lambda
# #   # following IPM book pg. 162-164 and SizeQualityExample.R
# #   
# #   ## compute meshpoints
# #   mz <- 100 # number of size mesh points 
# #   mt <- 50  # number of touch mesh points
# #   Lz <- 0.1 # size lower bound
# #   Uz<- 10  # size upper bound
# #   Lt <- 0 # touch lower bound
# #   Ut <- 1 # touch upper bound
# #   hz <- (Uz-Lz)/mz # size mesh point breaks
# #   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
# #   ht <- (Ut-Lt)/mt # touch mesh point breaks
# #   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# #   
# #   ## compute 4D kernel and 2D iteration matrix
# #   # Function eta to put kernel values in their proper place in A
# #   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# #   
# #   # matrix whose (i,j) entry is eta(ij)
# #   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# #   
# #   # code modified from IPM book, Ch 6, SizeQualityExample.R
# #   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
# #   for(i in 1:mz){
# #     for(j in 1:mt){
# #       for(k in 1:mt){
# #         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
# #         A[Eta[,k],Eta[i,j]]=kvals
# #         Kvals[,k,i,j]=kvals
# #         
# #       }}
# #     cat(i,"\n");
# #   }
# #   A<-hz*ht*A
# #   
# #   out <- domEig(A) # returns lambda and a vector proportional to w
# #   lam.boot <- out$lambda
# #   cat(lam.boot, "\n")
# #   return(lam.boot)
# #   
# # }
# # 
# # ### do the bootstrap (code takes ~ X min to run)
# # starttime <- Sys.time()
# # boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
# #                     R = 1000)
# # endtime <- Sys.time()
# # 
# # boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# # 
# # 
# # 
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ aMrU ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ----
# # ---- IPM size touch - aMrU: define kernel components ----
# 
# # IPM book pg. 24-25
# 
# ### survival: s(z)
# 
# s_z_st <- function(z, m.par_st){
#   
#   ## linear predictor
#   
#   # below size ceiling
#   linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
#   # above size ceiling
#   linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
#   
#   # impose ceiling and backtransform from logit
#   p <- ifelse(z <= Uz1
#               , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
#               , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
#   )
#   
#   ## output
#   return( unname(p) )
#   
#   # ## back transform from logit
#   # p <- 1/(1+exp(-linear.p))
# }
# 
# ### growth: G(z', z, t)
# 
# G_z1zt_st <- function(z1, z, t, m.par_st){
#   
#   ## based on growth regression coefficients, define mu and sigma of growth PDF
#   # below size ceiling
#   mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
#   # above size ceiling
#   mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t
#   
#   sig <- m.par_st['grow.sd']
#   
#   ## calculate growth PDF
#   # have to loop b/c it's a PDF, not a probability
#   p.den.grow <- c()
#   for(q in 1:length(z1)){
#     p.den.grow[q] <- ifelse(z <= Uz1
#                             , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
#                             , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
#     )
#   }
#   
#   ## output
#   return(p.den.grow)
#   
# }
# ### touch dynamics: T(t', t)
# 
# T_t1t_st <- function(t1, t, m.par_st){
#   
#   p.den.touch <- dbeta(t1
#                        , shape1 = 300
#                        , shape2 = 300
#   )
#   
#   ## output
#   return(p.den.touch)
# }
# 
# ### prob reproduce: p_bz
# 
# p_bz_st <- function(z, t, m.par_st){
#   ### convert z to v
#   v <- ifelse(z <= Uz1
#               , z * 131.18 - 152.52 # below ceiling
#               , Uz1 * 131.18 - 152.52 # above ceiling
#   )
#   
#   # linear predictor
#   linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
#   
#   # back transform from logit
#   p <- 1/(1+exp(-linear.p))
#   
#   # output
#   return( unname(p) )
# }
# 
# 
# ### fecundity: b_z
# 
# b_z <- function(z){
#   # from Strathmann
#   
#   p.bz <- ifelse(z <= Uz1
#                  , 0.147 * (10 * z)^2.74
#                  , 0.147 * (10 * Uz1)^2.74
#   )
#   
#   return(p.bz)
# }
# 
# ### recruit size dist: C_0(z')
# 
# C_0z1 <- function(z1, m.par_st){
#   mu <- m.par_st['rcsz.int']
#   sig <- m.par_st['rcsz.sd']
#   p.den.rcsz <- dnorm(z1, mean = mu, sd = sig)
#   return(p.den.rcsz)
# }
# 
# ### recruit touch dist: C_0(t')
# 
# C_0t1 <- function(t1, m.par_st){
#   
#   p.den.rcst <- dbeta(t1 # recruits start with equal chance across all levels of touch
#                       , shape1 = 1
#                       , shape2 = 1
#   )
#   
#   # shp1 <- m.par_st['rcst.s1']
#   # shp2 <- m.par_st['rcst.s2']
#   # p.den.rcst <- dbeta(t1, shape1 = shp1, shape2 = shp2 )
#   
#   
#   return(p.den.rcst)
# }
# 
# 
# ### define P component of kernel (survival and growth): P(z',z,t) = s(z)*G(z',z, t)*T(t1,t)
# 
# P_z1z <- function(z1, z, t1, t, m.par_st){
#   return( s_z_st(z, m.par_st) * G_z1zt_st(z1, z, t, m.par_st) * T_t1t_st(t1, t, m.par_st)
#   )
# }
# 
# ### define F component of kernel (reproduction): F(z',z,t')=s(z)*Pb(z, t)*b(z)*Pr*Co(z')*Co(t')
# 
# F_z1z <- function(z1, z, t1, t, m.par_st){
#   return(
#     unname(
#       s_z_st(z, m.par_st) * p_bz_st(z, t, m.par_st) * b_z(z) * m.par_st['p.r'] * C_0z1(z1, m.par_st)  * C_0t1(t1, m.par_st)
#     )
#   )
# }
# 
# 
# ### define K kernel
# k_st <- function(z1, t1,  z, t, m.par_st){
#   P_z1z(z1, z, t1, t, m.par_st) + F_z1z(z1, z, t1, t, m.par_st)
# }
# 
# # ---- *slow, ~ 2 min* IPM size touch - aMrU: numerical implementation ----
# 
# # following IPM book pg. 162-164 and SizeQualityExample.R
# 
# ### compute meshpoints
# mz <- 100 # number of size mesh points 
# mt <- 50  # number of touch mesh points
# Lz <- 0.1 # size lower bound
# Uz<- 10  # size upper bound
# Lt <- 0 # touch lower bound
# Ut <- 1 # touch upper bound
# hz <- (Uz-Lz)/mz # size mesh point breaks
# yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
# ht <- (Ut-Lt)/mt # touch mesh point breaks
# yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# 
# ### compute 4D kernel and 2D iteration matrix
# 
# 
# 
# # Function eta to put kernel values in their proper place in A
# eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# 
# # matrix whose (i,j) entry is eta(ij)
# Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# 
# # code modified from IPM book, Ch 6, SizeQualityExample.R
# A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
# for(i in 1:mz){
#   for(j in 1:mt){
#     for(k in 1:mt){
#       kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
#       A[Eta[,k],Eta[i,j]]=kvals
#       Kvals[,k,i,j]=kvals
#       
#     }}
#   cat(i,"\n");
# }
# A<-hz*ht*A
# 
# # # tried using simpler method on IPM book pg. 166 to create A matrix, much slower. nope.
# # 
# # Kd <- expand.grid(z1 = yz, t1 = yt, z = yz, t = yt)
# # 
# # A <- with(Kd, k_st(z1, t1, z, t, m.par_st))
# # dim(A) <- c(mz * mt, mz * mt)
# # A <- hz*ht*A
# 
# ### calculate Lambda, w, and v
# out <- domEig(A) # returns lambda and a vector proportional to w
# out2 <- domEig(t(A)) # returns lambda and a vector proportional to v
# lam.stable <- out$lambda
# lam.stable.t <- out2$lambda
# w <- Re(matrix(out$w, mz, mt))
# w <- w/(hz*ht*sum(w))
# v <- Re(matrix(out2$w, mz, mt))
# v <- v/sum(v) 
# 
# # ---- IPM size touch - aMrU: perturbation analysis ----
# 
# ### Compute elasticity matrix 
# repro.val_aMrU <- matrix(v, mz, mt)
# stable.state_aMrU <- matrix(w, mz, mt) 
# stable.state_aMrU.vector <- apply(stable.state_aMrU, 1, sum)
# v.dot.w <- sum(hz*ht*stable.state_aMrU*repro.val_aMrU)
# sens <- outer(repro.val_aMrU,stable.state_aMrU)/v.dot.w
# elas <- sens*Kvals/lam.stable
# 
# ### Compute matrix of total (=integrated) elasticities for all transitions (x_i,q_j) to anywhere 
# total.elas <- hz*ht*apply(elas,c(3,4),sum) 
# 
# ### Checks
# cat("Forward and transpose iteration: ",lam.stable," should be the same as ",lam.stable.t,"\n");  
# cat("Integrated elasticity =",sum(hz*ht*hz*ht*elas)," and it should = 1","\n") 
# 
# ### Plots
# 
# ## stable size distribution
# 
# plot(yz
#      , (stable.state_aMrU.vector/sum(stable.state_aMrU.vector))*10
#      , xlab = "Operculum length, mm"
#      , ylab = "Frequency"
#      , type = "l"
#      , ylim = c(0, 1)
# )
# 
# # # code following IPM book (I modified this to make it a probability)
# # plot(yz
# #      , stable.state_aMrU.vector
# #      , xlab = "Size x"
# #      , ylab = "Frequency"
# #      , type = "l"
# # )
# 
# 
# 
# 
# ## stable size and crowding distribution
# image(yz
#       , yt
#       , stable.state_aMrU
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , stable.state_aMrU
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## reproductive value
# image(yz
#       , yt
#       , repro.val_aMrU
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , repro.val_aMrU
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# ## total elasticity
# # image(yt
# #       , yz
# #       , t(total.elas)
# #       , col = grey(seq(0.5, 1, length=100))
# #       , xlab = "Crowding"
# #       , ylab = "Operculum length, mm"
# # )
# # contour(yt
# #         , yz
# #         , t(total.elas)
# #         , add = TRUE
# #         , nlevels = 6
# #         , labcex = 0.8
# # )
# 
# image(yz
#       , yt
#       , total.elas
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , total.elas
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )
# 
# # # ---- *slow, ~ 34 hrs* IPM size touch - aMrU: lambda CI ----
# # # from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30
# # 
# # # NOTE: only the survival/growth dataset is boostrapped. other sources of uncertainty (e.g. reproduction) are not quantified b/c their functions draw on separate datasets. 
# # 
# # ### function to compute lambda from a bootstrapped data set in format required by library(boot)
# # boot.lam.st <- function(dataset, sample.index) {
# #   
# #   ### extract the data used to make this fit
# #   boot.data <- dataset[sample.index, ] # don't shuffle using sample index, library(boot) does it.
# #   
# #   ### fit the functions
# #   
# #   ## survival
# #   mod.Surv <- glm(Surv ~ z
# #                   , family = binomial
# #                   , data = boot.data
# #   )
# #   
# #   ## growth
# #   
# #   boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
# #   
# #   #z.growth <- boot.data.growth$z
# #   mod.Grow_st <- lm(z1 ~ z * touch_pct
# #                     , data = boot.data.growth
# #                     
# #   )
# #   
# #   ## recruit size distribution
# #   mod.Rcsz <- lm(z ~ 1 
# #                  , data = rec
# #   )
# #   
# #   ### Store the estimated parameters
# #   
# #   m.par_st <- c(
# #     surv = coef(mod.Surv)
# #     , grow = coef(mod.Grow_st)
# #     , grow.sd = summary(mod.Grow_st)$sigma
# #     , rcsz = coef(mod.Rcsz)
# #     , rcsz.sd = summary(mod.Rcsz)$sigma
# #     , rcst = coef(betafit_t)
# #     , p.r = p.r.est
# #     , repr = coef(reprodyn.1b)
# #     #, U1 = U1
# #   )
# #   
# #   names(m.par_st) <- c(
# #     'surv.int'
# #     , 'surv.z'
# #     , 'grow.int'
# #     , 'grow.z'
# #     , 'grow.t'
# #     , 'grow.zt'
# #     , 'grow.sd'
# #     , 'rcsz.int'
# #     , 'rcsz.sd'
# #     , 'rcst.s1'
# #     , 'rcst.s2'
# #     , 'p.r'
# #     , 'repr.int' # model is for volume, NOT operc length
# #     , 'repr.t' # model is for volume, NOT operc length
# #     , 'repr.v' # v here b/c it's volume, NOT operc length
# #     #, 'U1'
# #   )
# #   
# #   ### implement IPM and calculate lambda
# #   # following IPM book pg. 162-164 and SizeQualityExample.R
# #   
# #   ## compute meshpoints
# #   mz <- 100 # number of size mesh points 
# #   mt <- 50  # number of touch mesh points
# #   Lz <- 0.1 # size lower bound
# #   Uz<- 10  # size upper bound
# #   Lt <- 0 # touch lower bound
# #   Ut <- 1 # touch upper bound
# #   hz <- (Uz-Lz)/mz # size mesh point breaks
# #   yz <- Lz + hz*((1:mz)-0.5) # size actual mesh points
# #   ht <- (Ut-Lt)/mt # touch mesh point breaks
# #   yt <- Lt + ht*((1:mt)-0.5) # touch acutal mesh points
# #   
# #   ## compute 4D kernel and 2D iteration matrix
# #   # Function eta to put kernel values in their proper place in A
# #   eta_ij <- function(i, j, mz) {(j-1)*mz+i}
# #   
# #   # matrix whose (i,j) entry is eta(ij)
# #   Eta <- outer(1:mz, 1:mt, eta_ij, mz = mz)
# #   
# #   # code modified from IPM book, Ch 6, SizeQualityExample.R
# #   A=matrix(0,mz*mt,mz*mt); Kvals=array(0,c(mz,mt,mz,mt));
# #   for(i in 1:mz){
# #     for(j in 1:mt){
# #       for(k in 1:mt){
# #         kvals=k_st(yz,yt[k],yz[i],yt[j],m.par_st)
# #         A[Eta[,k],Eta[i,j]]=kvals
# #         Kvals[,k,i,j]=kvals
# #         
# #       }}
# #     cat(i,"\n");
# #   }
# #   A<-hz*ht*A
# #   
# #   out <- domEig(A) # returns lambda and a vector proportional to w
# #   lam.boot <- out$lambda
# #   cat(lam.boot, "\n")
# #   return(lam.boot)
# #   
# # }
# # 
# # ### do the bootstrap (code takes ~ X min to run)
# # starttime <- Sys.time()
# # boot.out.st <- boot(data = dat, statistic = boot.lam.st, simple = TRUE,
# #                     R = 1000)
# # endtime <- Sys.time()
# # 
# # boot.ci(boot.out.st, type = c('norm', 'basic', 'perc'))
# # 
# # 
# # 