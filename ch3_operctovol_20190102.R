### ch 3 using expt 2.4.1 ch 2 data
### Will King
### comparing operculum length to volume for banracles used in ch 2 mesocosm and field experiments. All barnacles pooled together and chosen randomly for measuring (all temps, meso/field, crowding, size pooled together)

### use this to do regression to predict vol from operc length, for ch 3 probability of reprod
### for operc-to-vol regression, only use data of volume 600mm3 or smaller, b/c that was the range used for the vol-to-prob-reproduction regression in ch 2, which you want to use in ch 3


# ---- load data and packages ----

setwd("~/PhD/PhD_projects/ch 3 IPM/analysis")
ov <- read.csv('expt2.4.1_operctovol_20180514.csv', header = T)
library(scales)


# ----- adjust data part 1 ----

# calculate volume
ov$vol <- ov$basal_diameter_mm^2 * ov$height_mm

# ----- explore: operc against various measures ----

plot(ov$vol ~ ov$operculum_length_mm)
abline(h = 600, lty = 2)

# ---- adjust data: limit to 600 mm3 or below ----

ov2 <- ov[ov$vol <= 600, ]

plot(ov2$vol ~ ov2$operculum_length_mm)

# ---- analyze: vol ~ operc length ----

m1 <- lm(vol ~ operculum_length_mm
         , data = ov2
)
summary(m1)

# ---- plot: vol ~ operc length, for supp fig ----

pdf('plots/operctovol_ch3.pdf', width = 5.5, height = 5.5)
par(cex = 1.2, mar = c(4.5, 5, 0.5, 0.5))

# make blank graph
plot(vol ~ operculum_length_mm
     , data = ov2
     , pch = 2
     , col = alpha('black', 0.55)
     , axes = F
     , xlab = ''
     , ylab = ''
     , xlim = c(0, 6)
     , ylim = c(0, 600)
)

# add regression line

lines(seq(1.15, 5.75, by = 0.01)
      , predict(m1
                , newdata = data.frame(operculum_length_mm = seq(1.15, 5.75, by = 0.01) )
                , type = 'response'
      )
      , lwd = 2
)


# axes and labels

axis(1
     , tck = 0.02
     , at = seq(0, 6, by = 1)
     #, labels = seq(0.5, 6.5)
     , pos = 0
)
axis(2
     , las = 1
     , tck = 0.02
     , at = seq(0, 600, by = 100)
     , pos = 0
)
axis(3
     , tck = 0.02
     , at = seq(0, 6, by = 1)
     , labels = F
     , pos = 600
)
axis(4
     , tck = 0.02
     , at = seq(0, 600, by = 100)
     , labels = F
     , pos = 6
)
mtext(expression(paste('Body volume, ', mm^3))
      , side = 2, line = 2, cex = 1.2)
mtext('Operculum length, mm', side = 1, line = 2, cex = 1.2)


dev.off()
