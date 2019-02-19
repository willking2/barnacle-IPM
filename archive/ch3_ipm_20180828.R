### Will King
### Ch 3 - population model for B. glandula, IPM
### data from field survey of barnacles after 1 year (Jun 2017 - May 2018)

# ---- load data and packages ----

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

# ---- explore data ----

### general data structure
str(dat)  
summary(dat)

### growth
plot(z1 ~ z
     , data = dat#[dat$subsite == 'RT2' | dat$subsite == 'RT3', ]
     , xlim = c(0, 8)
     , ylim = c(0, 8)
     )
abline(a = 0, b = 1, lty = 2)
abline(lm(z1 ~ z, data = dat))

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


### survival
plot(Surv ~ z
     , data = dat
)

  