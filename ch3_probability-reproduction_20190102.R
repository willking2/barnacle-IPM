### Will King
### for Ch3, using Expt 2.4.1: air warm x density on individual vital rates (B. glandula)
### field (summer 2017) - reproduction and body weight from dissections

# analysis_a: use subsite (categorical) as a fixed effect
#analysis_b: use temperature (grand mean per subsite; cont') as a fixed effect

# * means has code for plots used in ms

# ---- load data/packages ----

library(lme4)
library(lmerTest)
library(MuMIn)
library(car)
library(scales)

rd <- read.csv('expt2.4.1_field-dissection-shell-vs-somagonad_20180306.csv', 
               header = T)


# ---- adjust data part 1 ----

#### Note on reproduction data collection methods
# Barnacles for reprod were collected separately from barnacles in photo analysis. These barnacles were not in any density-manipulated plots; just from in situ adjacent to the photo plots. 

# Did not record identity (spp) of competing barnacles when collecting B. glandula for dissection. 

# Dissection/weighing was done as follows. If not yolky, dissected into shell and whole body (soma and any minimal gonad). Weight for whole body recorded in "Wboat + gonad + soma weight" column. If yes yolky, trisected into shell, soma only, and gonad only. Weight for soma recorded in "Wboat + gonad + soma weight" column, but only has soma weight. Weight for gonad recorded in "Wboat + gonad only column"
####


# remove rows where fly maggot was found
rd <- rd[rd$notes != 'fly maggot found', ]

## calc pct touch
rd$touch_pct <- round(rd$touch_raw/360, digits = 2)


## adjustments for reproduction data

rd$yolky_yn <- as.numeric(rd$yolky_yn) 
rd$yolky_yn[rd$yolky_yn == 1] <- 0 # change y/n to 1/0, 0 means not yolky, 1 means yes yolky
rd$yolky_yn[rd$yolky_yn == 2] <- 1
rd$f.yolky_yn <- as.factor(rd$yolky_yn) ## make new factor variable for analysis


## calc approximate volume
rd$vol <- (rd$basal_diameter_mm^2)*rd$height_mm # basal length^2 * height


# remove one outlier that is obviously wrong (extremely small, but yes reprod)
rd <- rd[-c(which(rd$f.yolky_yn == 1 & rd$vol < 50)), ]


# ---- adjust data part 2 ----

# make dataset for only reproductive barnacles' gonad mass
rd.gonad <- rd[rd$f.yolky_yn == 1, ]
rd.gonad <- rd.gonad[is.na(rd.gonad$gonad_g) == F, ]
# change gonad mass to mg
rd.gonad$gonad_mg <- rd.gonad$gonad_g * 1000


## calculate total weight
rd.gonad$totalweight_g <- rd.gonad$shell_g + rd.gonad$somagonad_g


## calculate ratios
# shell to volume
rd.gonad$shelltovol <- rd.gonad$shell_g / rd.gonad$vol

# gonad to shell
rd.gonad$gonadtoshell <- rd.gonad$gonad_g / rd.gonad$shell_g

# gonad to volume
rd.gonad$gonadtovol  <- rd.gonad$gonad_mg / rd.gonad$vol



# make a copy before limiting dataset to well sampled regions
rd.gonad.full <- rd.gonad

## limit datasets to well sampled regions
rd <- rd[rd$vol < 600, ]
rd.gonad <- rd.gonad[rd.gonad$vol > 200 & rd.gonad$vol < 600, ]



# scale data to appease glmer's warning messages
rd.s <- rd
rd.s[ , names(rd.s) == 'touch_pct'] <- 
  scale(rd.s[ , names(rd.s) == 'touch_pct'])
rd.s[ , names(rd.s) == 'vol'] <- 
  scale(rd.s[ , names(rd.s) == 'vol'])


# scaled gonad dataset
rd.s.gonad <- rd.gonad
rd.s.gonad[ , names(rd.s.gonad) == 'touch_pct'] <- 
  scale(rd.s.gonad[ , names(rd.s.gonad) == 'touch_pct'])
rd.s.gonad[ , names(rd.s.gonad) == 'vol'] <- 
  scale(rd.s.gonad[ , names(rd.s.gonad) == 'vol'])


# scaled gonad dataset [with full size range]
rd.s.gonad.full <- rd.gonad.full
rd.s.gonad.full[ , names(rd.s.gonad.full) == 'touch_pct'] <- 
  scale(rd.s.gonad.full[ , names(rd.s.gonad.full) == 'touch_pct'])
rd.s.gonad.full[ , names(rd.s.gonad.full) == 'vol'] <- 
  scale(rd.s.gonad.full[ , names(rd.s.gonad.full) == 'vol'])


# ---- analysis_a: yolky y/n ~ touch*size ----

# check assumption: #0's ~ #1's, see Zuur (2009) pg. 251
summary(rd.s) # not really... 68:29. # gonna go forward w/ logit link for now
# check assumption: overdispersion?
# not possible here, b/c treating each barnacle as individual bernoulli trial; NOT using proprotions where N > 1.


# make final model after dropping intxn
reprodyn.1b <- glm(f.yolky_yn ~ 
                    touch_pct + 
                    vol,
                  family = binomial,
                  data = rd
)

# make a model w/ only volume, for ch 3 analysis. 1/2/2019
reprodyn.1_sizeonly <- glm(f.yolky_yn ~ 
                              vol,
                            family = binomial,
                            data = rd
)


# ---- graveyard ----

