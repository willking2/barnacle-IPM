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

photos <- read.csv('ch3_field-photos.csv', header = T)
rec <- read.csv('ch3_recruits.csv', header = T)





# ---- generate recruits data ----

summary(rec)
rec$operc_length_mm <- round(rec$operc_length_mm, 1)
h <- hist(rec$operc_length_mm
          , freq = F
          , breaks = seq( min(rec$operc_length_mm), max(rec$operc_length_mm), 0.1)
)

rec_sizedist <- h$density * 0.1 # size distribution of recruits (probability)
rec_totalnew <- 700 # total number of recruits
howmany <- round(rec_totalnew * rec_sizedist, 0) # how many of each size recruited

rec_bins <- h$breaks[h$breaks != max(h$breaks)] # size bins of new recruits, excluding final bin
rec_z1 <- data.frame(z1 = rep(NA, sum(howmany))) # empty data frame for z1 recruits
rec_z1$z1 <- rep(rec_bins, howmany)




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

# ---- make datasets for IPM ----

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


# ---- adjust data: truncate for size ----

## remove data for Balanus touching mostly Chthamalus 
#dat <- dat[dat$touch_spp != 'C', ]


## response surface
plot(touch_pct ~ z
     , data = dat
)
abline(h = 0.8)
abline(v = 5.5)


# limit data to initial operc size =< 5.5 and touch =< 0.8 
dat <- dat[dat$z <= 5.5 & dat$touch_pct <= 0.8, ]


# make column for colors, based on touch
colPal <- colorRampPalette(c('gray80','black'))
# darker means more touch (lightest gray is touch = 0, black is touch = 1)
dat$colors <- colPal(100)[as.numeric(cut(dat$touch_pct,breaks = 100))]

# change survival from character to factor
dat$Surv <- as.factor(dat$Surv)



# ---- explore data ----

### general data structure
str(dat)  
summary(dat)


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

# resid vs. fitted
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



# ---- plot: growth ----

### calculate predictions and CIs

## touch = 0
# make blank area for new data
newgrowth.00 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0
  , z1 = 0
)
# add predictions
newgrowth.00$z1 <- predict(growth.1, newgrowth.00, re.form=NA)
# CIs
confint.00 <- predict(growth.1
                      , newdata = newgrowth.00
                      , interval="confidence"
                      , level = 0.95
)

## touch = 0.4
# make blank area for new data
newgrowth.04 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0.4
  , z1 = 0
)
# add predictions
newgrowth.05$z1 <- predict(growth.1, newgrowth.04, re.form=NA)
# CIs
confint.04 <- predict(growth.1
                      , newdata = newgrowth.04
                      , interval="confidence"
                      , level = 0.95
)

## touch = 0.8
# make blank area for new data
newgrowth.08 <- expand.grid(
  z = seq(0, 6, 0.1)
  , touch_pct = 0.8
  , z1 = 0
)
# add predictions
newgrowth.10$z1 <- predict(growth.1, newgrowth.08, re.form=NA)
# CIs
confint.08 <- predict(growth.1
                      , newdata = newgrowth.08
                      , interval="confidence"
                      , level = 0.95
)

### actually plot

## base plot
plot(z1 ~ z
     , data = dat
     , xlim = c(0, 6)
     , ylim = c(0, 8)
)
abline(a = 0, b = 1, lty = 1, col = 'gray') # 1:1 line of zero growth

## color points by touch
points(z1 ~ z
       , data = dat
       , pch = 19
       , col = dat$colors
)

## prediction lines and CIs
# touch = 0
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0
        y = c(confint.00[ , 2], rev(confint.00[ , 3])),
        col = alpha('gray', 0.5),
        border = NA)
lines(newgrowth.00$z1 ~ seq(0, 6, 0.1)
      , lwd = 2
      , lty = 3
)
# touch = 0.4
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0.4
        y = c(confint.04[ , 2], rev(confint.04[ , 3])),
        col = alpha('gray', 0.5),
        border = NA)
lines(newgrowth.05$z1 ~ seq(0, 6, 0.1)
      , lwd = 2
      , lty = 5
)
# touch = 0.8
polygon(x = c(seq(0, 6, 0.1), rev(seq(0, 6, 0.1))), # touch = 0.8
        y = c(confint.08[ , 2], rev(confint.08[ , 3])),
        col = alpha('gray', 0.5),
        border = NA)
lines(newgrowth.10$z1 ~ seq(0, 6, 0.1)
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
Anova(survival.1)

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






