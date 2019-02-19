### Will King
### Ch 3 - population model for B. glandula, IPM
### data from field survey of barnacles after 1 year (Jun 2017 - May 2018)

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

#source('MatrixImage.R') # from Ellner, Childs, and Rees 2016, IPM book

photos <- read.csv('ch3_field-photos.csv', header = T)
rec <- read.csv('ch3_recruits_sizedist.csv', header = T)
totrec <- read.csv('ch3_recruits_total.csv', header = T)




# ---- generate recruits data ----

summary(rec)
rec$operc_length_mm <- round(rec$operc_length_mm, 1)
h <- hist(rec$operc_length_mm
          , freq = F
          , breaks = seq( min(rec$operc_length_mm), max(rec$operc_length_mm), 0.1)
)

rec_sizedist <- h$density * 0.1 # size distribution of recruits (probability)
rec_totalnew <- mean(totrec$number_of_recruits) # total number of recruits
howmany <- round(rec_totalnew * rec_sizedist, 0) # how many of each size recruited

rec_bins <- h$breaks[h$breaks != max(h$breaks)] # size bins of new recruits, excluding final bin
rec_z1 <- data.frame(z1 = rep(NA, sum(howmany))) # empty data frame for z1 recruits
rec_z1$z1 <- rep(rec_bins, howmany)

hist(rec_z1$z1
     , freq = F
     , breaks = seq( min(rec_z1$z1), max(rec_z1$z1), 0.1)
) # checked to see that the calculated distribution matches empirical data 11/14/2018


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

# ### reproduction
# # if survival = 0, reprod and offspring should have NA
# # reproduced:
# # offspring: # embryos = 0.147 * (operculum length in mm^2.74). Based on Strathmann et al. 1981 Oecologia (see Schubart thesis pg. 20)
# 
# # reproduced? 0/1
# dat$Reprod <- ''
# dat$Reprod[dat$Surv == 0] <- NA
# 
# # fecundity metric
# dat$Offspring <- ''
# dat$Offspring[dat$Surv == 0] <- NA
# dat$Offspring[dat$Surv == 1] <- 0.147 * (10 * dat$z[dat$Surv == 1])^2.74


# ---- adjust data: subsite for size, color by touch, surv to factor ----

## remove data for Balanus touching mostly Chthamalus 
#dat <- dat[dat$touch_spp != 'C', ]


## response surface
plot(touch_pct ~ z
     , data = dat
)
## color points by survived or no
# lived
points(touch_pct ~ z
       , data = dat[dat$Surv == 1, ]
       , col = 'blue'
       , pch = 19
)
# died
points(touch_pct ~ z
       , data = dat[dat$Surv == 0, ]
       , col = 'red'
       , pch = 19
)

abline(h = 0.8)
abline(v = 5.5)


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




# ---- explore data ----

### general size structure

hist(dat$z)
hist(dat$z1)

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
  
# ---- analyze: growth ----

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
) # prefer growth.1
BIC(growth.1
    , growth.2
    , growth.3
    , growth.4
    , growth.5
    , growth.6
    , growth.7
) # prefer growth.1


### determine fixed structure
summary(growth.1) # R^2
Anova(growth.1)

### model validation

# resid vs. fitted ### note non-constant variance ###
plot(x = fitted(growth.1),
     y = resid(growth.1),
     xlab = 'fitted values',
     ylab = 'residuals')
lines(lowess(fitted(growth.1), resid(growth.1)), col="blue", lwd=2)


# resid vs. explanatory
plot(x = dat$z[is.na(dat$z1) == F],
     y = resid(growth.1, type = 'deviance')
)

plot(x = dat$touch_pct[is.na(dat$z1) == F],
     y = resid(growth.1, type = 'deviance')
)

qqnorm(resid(growth.1))
qqline(resid(growth.1))


### try to model the non-constant variance ###

# hetergeneity of variance; see Zurr pg. 74 and 90
# variance seems to be greater at smaller body sizes; model fixed variance structure where variance is inversely proportional to body size

dat.growth <- dat[is.na(dat$z1) == F, ]
z.growth <- dat.growth$z

growth.1_lm <- gls(z1 ~ z * touch_pct
                   , data = dat.growth
)
vfix <- varFixed(~(1/z.growth))

growth.1_gls <- gls(z1 ~ z * touch_pct
                    , weights = vfix
                    , data = dat.growth

)

AIC(growth.1_lm, growth.1_gls) # prefers gls
BIC(growth.1_lm, growth.1_gls) # prefers gls

## determine fixed structure of gls
Anova(growth.1_gls) # effect of size, and size*touch

## model validation

# resid vs. fitted
plot(x = fitted(growth.1_gls),
     y = resid(growth.1_gls
               , type = 'normalized'),
     xlab = 'fitted values',
     ylab = 'residuals')
lines(lowess(fitted(growth.1_gls), resid(growth.1_gls)), col="blue", lwd=2)


# resid vs. explanatory
plot(x = dat$z[is.na(dat$z1) == F],
     y = resid(growth.1_gls, type = 'normalized')
)

plot(x = dat$touch_pct[is.na(dat$z1) == F],
     y = resid(growth.1_gls, type = 'normalized')
)

qqnorm(resid(growth.1_gls))
qqline(resid(growth.1_gls))

# ---- plot: growth ----

# model again, for convenience
dat.growth <- dat[is.na(dat$z1) == F, ]
z.growth <- dat.growth$z
vfix <- varFixed(~(1/z.growth))

growth.1_gls <- gls(z1 ~ z * touch_pct
                    , weights = vfix
                    , data = dat.growth
                    
)


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
growthpreds.00 <- predict(growth.1_gls, newgrowth.00, se.fit = T)
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
growthpreds.04 <- predict(growth.1_gls, newgrowth.04, se.fit = T)
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
growthpreds.08 <- predict(growth.1_gls, newgrowth.08, se.fit = T)
# CIs
growthCI.upper.08 <- growthpreds.08$fit + critval*growthpreds.08$se.fit
growthCI.lower.08 <- growthpreds.08$fit - critval*growthpreds.08$se.fit
growthfit.08 <- growthpreds.08$fit

### actually plot

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


mtext('Operculum length (mm), year t'
      , side = 1
      , line = 2
      , cex = 1.2
)
mtext('Operculum length (mm), year t+1'
      , side = 2
      , line = 2.5
      , cex = 1.2
)

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
summary(survival.1)
Anova(survival.1) # only initial size is significant

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


# ---- IPM: prepare suvival and growth regression models ----

### survival

# final model from analysis only has body size anyway, so was prepared in section "analyze: survival", above.
# it is named mod.Surv
# repeated here for ease of use
mod.Surv <- glm(Surv ~ z
                , family = binomial
                , data = dat
)



### growth

# code for modeling non-constant variance repeated here for ease of use. see "analyze: growth"
dat.growth <- dat[is.na(dat$z1) == F, ]
z.growth <- dat.growth$z
vfix <- varFixed(~(1/z.growth))

## model w/ body size only (lumping all crowding); keep it simple to start with
mod.Grow_simple <- gls(z1 ~ z
                       , weights = vfix
                       , data = dat.growth
                       
)

## model w/ bodysize*crowd, based on section "analyze: growth", above.



# ---- IPM: general parameters ----

### ceiling 
# U1 refers to upper limit of size z, beyond which vital rates are the same as U1
U1 = 6

### number of 'bins' for iteration matrix
nBigMatrix <- 500


# ---- IPM_simple: prepare parameter estimates ----

# IPM book pg. 23-24 (regressions), 48 (ceiling)

m.par_simple <- c(
  surv = coef(mod.Surv)
  , grow = coef(mod.Grow_simple)
  , grow.sd = summary(mod.Grow_simple)$sigma
  #, U1 = U1
)

names(m.par_simple) <- c(
  'surv.int'
  , 'surv.z'
  , 'grow.int'
  , 'grow.z'
  , 'grow.sd'
  #, 'U1'
)

# ---- IPM_simple: define kernel components ----

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
  mu.ceiling <- m.par_simple['grow.int'] + m.par_simple['grow.z'] * U1 # above ceiling
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
    

### define P component of kernel (survival and growth): P(z',z) = s(z)*G(z',z)

P_z1z <- function(z1, z, m.par_simple){
  return( s_z_simple(z, m.par_simple) * G_z1z_simple(z1, z, m.par_simple) )
}

# ---- IPM simple: explore kernel components ----

# Growth: see MonocarpGrowthEviction.R and IPM book pg. 45-48

# use a modified growth function here that implements the ceiling manually, for exploration purposes. However, this function is not used for numerical implementation.

### create growth ceiling demo function
G_z1z_simple.ceiling <- function(z1, z, m.par_simple){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  mu <- m.par_simple['grow.int'] + m.par_simple['grow.z'] * z # below ceiling
  mu.ceiling <- m.par_simple['grow.int'] + m.par_simple['grow.z'] * U1 # above ceiling
  sig <- m.par_simple['grow.sd']
  
  ## calculate growth PDF
  
  # make blank vector
  p.den.grow <- c()
  
  #loop through z1
  for( j in 1:length(z1) ){
    p.den.grow[j] <- ifelse(z <= U1
                            , dnorm(z1[j], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[j], mean = mu.ceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
}


### explore growth

zvals <- seq(0, 10, length = 1000)
g1 <- G_z1z_simple.ceiling(zvals, min(zvals), m.par_simple) # small individuals
g2 <- G_z1z_simple.ceiling(zvals, mean(zvals), m.par_simple) # medium individuals
g3 <- G_z1z_simple.ceiling(zvals, max(zvals), m.par_simple) # large individuals

plot(x = zvals
     , y = g1
     , type = 'l'
     , lty = 1
     , lwd = 1
     , xlab = "Final size z'"
     , ylab = 'Probability density'
     , ylim = c(0, 1)
)
lines(x = zvals
      , y = g2
      , lwd = 2
)
lines(x = zvals
      , y = g3
      , lwd = 3
)

### create survival ceiling demo function
s_z_simple.ceiling <- function(z, m.par_simple){
  
  ## linear predictor
  linear.p <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * z # for cases below ceiling
  linear.p.ceiling <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * U1 # for cases > U1
  
  ## back transform from logit
  # code for manually doing ceiling
  p <- ifelse(z <= U1
              , 1/(1+exp(-linear.p)) # below ceiling; use regression coeff of z
              , 1/(1+exp(-linear.p.ceiling)) # above ceiling; use vital rates of U1
  )
  
  
  ## output
  return( unname(p) )
  
  
}



### explore fraction wrongfully evicted

WrongPlace <- function(z, U) {
  fac1 <- s_z_simple.ceiling(z, m.par_simple)
  fac2 <- integrate(function(x) G_z1z_simple.ceiling(x, z, m.par_simple), U, Inf)$value
  return(fac1 * fac2)
}

zvals <- Wvals <- seq(0, 9, length = 900)

for (j in seq_along(zvals)) Wvals[j] <- WrongPlace(zvals[j], 9)

plot(zvals
     , Wvals
     , type = "l"
     , lty = 1
     , lwd = 2
     , col = "black"
     , xlab = "Initial size z"
     , ylab = "Fraction wrongfully evicted"
     , ylim = c(0, 1.1 * max(Wvals))
)

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
  K <- P # IPM kernel only has survival and growth
  
  #matrix.image(K) # look at kernel
  
  return(list(K = K, meshpts = meshpts, P = P))
  
}
#dev.off()

### create IPM kernel

IPM_simple <- mk_K(nBigMatrix, m.par_simple, 0.1, 10)

### dominant eigenvalue, lambda, population growth rate
Re(eigen(IPM_simple$K)$values[1])  # need to put in recruits lol

### eigenvector, w, stable size distribution
meshpts_simple <- IPM_simple$meshpts
w_simple <- Re(eigen(IPM_simple$K)$vectors[,1])
stable.z.dist_simple <- w_simple/sum(w_simple) # scale eigenvector to discrete probability density

plot(stable.z.dist_simple/diff(meshpts_simple)[1]  ~ meshpts_simple
     , type = 'l'
) #from Monocarp Calculations.R
 


# ---- IPM simple: lambda CI ----

# from Monocarp Lambda Bootstrap CI.R, IPM book pg. 30

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
    # non-constant variance
    boot.data.growth <- boot.data[is.na(boot.data$z1) == F, ]
    #z.growth <- boot.data.growth$z
    vfix <- varFixed(~(1/boot.data.growth$z))
    mod.Grow <- gls(z1 ~ z
                    # , weights = vfix
                    , data = boot.data.growth
                    
    )
  
  
  ### Store the estimated parameters
    
    m.par.est <- c(
      surv = coef(mod.Surv)
      , grow = coef(mod.Grow)
      , grow.sd = summary(mod.Grow)$sigma
    )
    
    names(m.par.est) <- c(
      'surv.int'
      , 'surv.z'
      , 'grow.int'
      , 'grow.z'
      , 'grow.sd'
    )  
    
  ### implement IPM and calculate lambda 
  
    IPM.est <- mk_K(nBigMatrix, m.par.est, 0.1, 10)
    
    lam.boot <- Re(eigen(IPM.est$K, only.values = TRUE)$values[1])
    cat(lam.boot, "\n")
  
return(lam.boot)
}

### do the bootstrap (code takes 5.7 min to run)
#starttime <- Sys.time()
boot.out <- boot(data = dat, statistic = boot.lam, simple = TRUE,
                 R = 1000)
#endtime <- Sys.time()

boot.ci(boot.out, type = c('norm', 'basic', 'perc'))

