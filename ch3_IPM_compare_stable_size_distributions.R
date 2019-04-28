### Will King
### Ch 3 - population model for B. glandula, IPM

# need to run models from ch3_IPM_20190225 main file before able to plot these

### compare size distributions of aLrH, and aHrL


# ---- calculate stable size distributions ----

stablesizes_aLrH <- (stable.state_aLrH.vector/sum(stable.state_aLrH.vector))*10
stablesizes_aHrL <- (stable.state_aHrL.vector/sum(stable.state_aHrL.vector))*10

# ---- make dataset of stable size distributions ----

# full distributions
stablesizes_trts <- data.frame(scenario = rep(NA, 200)
                               , size = rep(NA, 200)
                               , density = rep(NA, 200)
)
stablesizes_trts$scenario <- c(rep("aLrH", 100)
                               , rep("aHrL", 100)
)
stablesizes_trts$size <- yz  
stablesizes_trts$density <- c(stablesizes_aLrH, stablesizes_aHrL)

# distribution of barnacles less than 2 mm
stablesizes_trts.small <- stablesizes_trts[stablesizes_trts$size < 2
                                         , ]

# distribution of barnacles greater than 2 mm
stablesizes_trts.big <- stablesizes_trts[stablesizes_trts$size > 2
                                         #& stablesizes_trts$size < 6
                                         , ]

# ---- KS tests ----

# for entire distribution
ks.test(stablesizes_aLrH, stablesizes_aHrL)

# for barnacles below 2 mm
# aLrH vs. aHrL: NS
ks.test(stablesizes_trts.small$density[stablesizes_trts.small$scenario == 'aLrH']
        , stablesizes_trts.small$density[stablesizes_trts.small$scenario == 'aHrL']
)

# for barnacles above 2 mm
# aLrH vs. aHrL: significant
ks.test(stablesizes_trts.big$density[stablesizes_trts.big$scenario == 'aLrH']
        , stablesizes_trts.big$density[stablesizes_trts.big$scenario == 'aHrL']
)
       

# ---- plot as densities: full size range ----

### full plot

#par(fig = c(0, 1, 0, 1)) # use this for doing main/inset figure

pdf('plots/compare_sizedistributions.pdf', width = 5, height = 5)
par(cex = 1.2
    , mar = c(4, 4, 0.5, 0.5)
)

hist(c(rec$z, dat$z1)
     , freq = F
     , ylim = c(0, 1)
     , xlim = c(0, 10)
     , xlab = ''
     , ylab = ''
     , axes = F
     , main = ''
     , col = 'gray80'
     , border = 'white'
)

# # aUrU
# plot(yz
#      , stablesizes_aHrL
#      # , xlab = "Body size (Operculum length, mm)"
#      # , ylab = "Probability density"
#      , type = 'n'
#      # , col = 'gray'
#      , ylim = c(0, 1)
#      , xlab = ''
#      , ylab = ''
#      , axes = F
# )

# line showing portion compared w/ KS test
lines(x = c(2, 2)
      , y = c(0, 0.3)
      #, col = 'gray'
      , lty = 2
)
lines(x = c(2, 10)
      , y = c(0.3, 0.3)
      #, col = 'gray'
      , lty = 2
)

# aLrH
lines(yz
     , stablesizes_aLrH
     , type = "l"
     , ylim = c(0, 1)
     , lwd = 1
)

# aHrL
lines(yz
      , stablesizes_aHrL
      , type = "l"
      , ylim = c(0, 1)
      , lwd = 2
)


# line labels
text(x = 5
     , y = 0.25
     , 'eLrH'
)

lines(x = c(2.6, 4.3)
      , y = c(0.11, 0.23)
)
  
  
text(x = 5
     , y = 0.15
     , 'eHrL'
)
lines(x = c(4, 4.35)
      , y = c(0.1, 0.13)
)


# axes and axes labels

axis(1
     , las = 1
     , pos = 0
     , tck = 0.02
)
axis(2 
     , pos = 0
     , tck = 0.02
     , las = 2
)
axis(3
     , las = 1
     , pos = 1
     , tck = 0.02
     , labels = F
)

axis(4
     , pos = 10
     , tck = 0.02
     , labels = F
)


mtext('Body size (mm)'
      , side = 1
      , line = 2
      , cex = 1.2
)

mtext('Probability density'
      , side = 2
      , line = 2.25
      , cex = 1.2
)

text(x = 7.5
     , y = 0.2
     , 'P = 0.023'
)

lines(x = c(6, 6)
      , y = c(0.15, 0.25)
)

dev.off()


# 
# ### inset: reproductive values
# 
# par(fig = c(0.5, 1, 0.5, 1)
#     , new = T)
# 
# image(yz
#       , yt
#       , repro.val_aHrL
#       , col = grey(seq(0.5, 1, length=100))
#       , xlab = "Operculum length, mm"
#       , ylab = "Crowding"
# )
# contour(yz
#         , yt
#         , repro.val_aHrL
#         , add = TRUE
#         , nlevels = 6
#         , labcex = 0.8
# )

# ---- graveyard ----

# ---- plot as boxplots ----
# # 
# boxplot(density ~ scenario
#         , data = stablesizes_trts.small[stablesizes_trts.small$scenario != 'aUrU',]
# )


# ---- other plotting code ----
# 
# plot(yz
#      , (stable.state_aFrU.vector/sum(stable.state_aFrU.vector))*10
#      , xlab = "Operculum length, mm"
#      , ylab = "Probability density"
#      , type = "l"
#      , ylim = c(0, 0.11)
# )
# # varying adult crowding
# lines(yz
#       , (stable.state_aLrU.vector/sum(stable.state_aLrU.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'blue'
# )
# lines(yz
#       , (stable.state_aHrU.vector/sum(stable.state_aHrU.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'green'
# )
# lines(yz
#       , (stable.state_aMrU.vector/sum(stable.state_aMrU.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'red'
# )
# lines(yz
#       , (stable.state_aUrU.vector/sum(stable.state_aUrU.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'purple'
# )
# 
# # varying recruit crowding
# lines(yz
#       , (stable.state_aUrL.vector/sum(stable.state_aUrL.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'lightblue'
#       , lwd = 2
# )
# lines(yz
#       , (stable.state_aUrH.vector/sum(stable.state_aUrH.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'green'
#       , lwd = 2
# )
# lines(yz
#       , (stable.state_aUrM.vector/sum(stable.state_aUrM.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'orange'
#       , lwd = 2
# )
# # varying adult and recruit crowding
# lines(yz
#       , (stable.state_aLrL.vector/sum(stable.state_aLrL.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'blue'
#       , lwd = 2
# )
# lines(yz
#       , (stable.state_aLrH.vector/sum(stable.state_aLrH.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'purple'
#       , lwd = 2
# )
# lines(yz
#       , (stable.state_aHrH.vector/sum(stable.state_aHrH.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'green'
#       , lwd = 2
# )
# lines(yz
#       , (stable.state_aHrL.vector/sum(stable.state_aHrL.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'red'
#       , lwd = 2
# )
# lines(yz
#       , (stable.state_aMrM.vector/sum(stable.state_aMrM.vector))*10
#       , xlab = "Operculum length, mm"
#       , ylab = "Probability density"
#       , type = "l"
#       , ylim = c(0, 1)
#       , col = 'orange'
#       , lwd = 2
# )
