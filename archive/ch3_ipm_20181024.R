### Will King
### Ch 3 - population model for B. glandula, IPM
### data from field survey of barnacles after 1 year (Jun 2017 - May 2018)

# ---- load data and packages ----

library(lme4)
library(car)

photos <- read.csv('ch3_field-photos.csv', header = T)






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

### reproduction
# if survival = 0, reprod and offspring should have NA
# fill in reprod and offspring using regressions from 2017 field data

# reproduced? 0/1
dat$Reprod <- ''
dat$Reprod[dat$Surv == 0] <- NA

# fecundity metric
dat$Offspring <- ''
dat$Offspring[dat$Surv == 0] <- NA

# ---- adjust data: truncate for size, remove touching mostly Balanus ----

## remove data for Balanus touching mostly Chthamalus 
dat <- dat[dat$touch_spp != 'C', ]


## response surface
plot(touch_pct ~ z
     , data = dat
)
abline(v = 5)

# limit data to initial operc size < 5
dat <- dat[dat$z < 5, ]


# make column for colors, based on touch
colPal <- colorRampPalette(c('gray80','black'))
# darker means more touch (lightest gray is touch = 0, black is touch = 1)
dat$colors <- colPal(100)[as.numeric(cut(dat$touch_pct,breaks = 100))]

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


# # spp competition?
# points(z1 ~ z
#        , data = dat[dat$touch_spp == '', ]
#        , col = 'black'
# )
# points(z1 ~ z
#        , data = dat[dat$touch_spp == 'B', ]
#        , col = 'green'
# )
# # points(z1 ~ z
# #        , data = dat[dat$touch_spp == 'C', ]
# #        , col = 'orange'
# # )

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
plot(Surv ~ z
     , data = dat
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
) # prefe growth.1
BIC(growth.1
    , growth.2
    , growth.3
    , growth.4
    , growth.5
    , growth.6
    , growth.7
) # prefer growth.1


### determine fixed structure
summary(growth.1)
Anova(growth.1)
