### Will King
### Ch 3 - population model for B. glandula, IPM

# need to run regressions from ch3_IPM_20190225 main file before able to plot these

# ---- plot both figs in one figure ----

pdf('plots/surv-and-growth.pdf', width = 10, height = 5)

par(mfrow = c(1, 2)
    , cex = 1.2
    , oma = c(4, 3, 1, 0.5)
    # , mar = c(1, 3, 1, 0.5)
)
# ---- plot: survival ----

par(mar = c(0, 0.5, 1, 2)
)

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


mtext('Body size (mm), year t'
      , side = 1
      , line = 2
      , cex = 1.2
)
mtext('Probability of survival'
      , side = 2
      , line = 2
      , cex = 1.2
)


mtext('a)'
      , side = 3
      , line = 0.5
      , cex = 1.2
      , adj = 0
)

# ---- plot: growth ----

par(mar = c(0, 2, 1, 0.5)
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


mtext('Body size (mm), year t'
      , side = 1
      , line = 2
      , cex = 1.2
)
mtext('Body size (mm), year t+1'
      , side = 2
      , line = 1.5
      , cex = 1.2
)

mtext('b)'
      , side = 3
      , line = 0.5
      , cex = 1.2
      , adj = 0
)


dev.off()


