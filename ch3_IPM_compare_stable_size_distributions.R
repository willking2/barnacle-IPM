### Will King
### Ch 3 - population model for B. glandula, IPM

# need to run models from ch3_IPM_20190225 main file before able to plot these

### compare size distributions of aUrU, aLrH, and aHrL

# ---- calculate stable size distributions ----

stablesizes_aUrU <- (stable.state_aUrU.vector/sum(stable.state_aUrU.vector))*10
stablesizes_aLrH <- (stable.state_aLrH.vector/sum(stable.state_aLrH.vector))*10
stablesizes_aHrL <- (stable.state_aHrL.vector/sum(stable.state_aHrL.vector))*10

# ---- make dataset of stable size distributions ----

# full distributions
stablesizes_trts <- data.frame(scenario = rep(NA, 300)
                               , size = rep(NA, 300)
                               , density = rep(NA, 300)
)
stablesizes_trts$scenario <- c(rep("aUrU", 100)
                               , rep("aLrH", 100)
                               , rep("aHrL", 100)
)
stablesizes_trts$size <- yz  
stablesizes_trts$density <- c(stablesizes_aUrU, stablesizes_aLrH, stablesizes_aHrL)

# distribution of barnacles between 2 and 6 mm
stablesizes_trts.big <- stablesizes_trts[stablesizes_trts$size > 2
                                         #& stablesizes_trts$size < 6
                                         , ]

# ---- KS tests ----

# for entire distribution
ks.test(stablesizes_aUrU, stablesizes_aLrH)
ks.test(stablesizes_aUrU, stablesizes_aHrL)
ks.test(stablesizes_aLrH, stablesizes_aHrL)

# for barnacles between 2 and 6 mm
ks.test(stablesizes_trts.big$density[stablesizes_trts.big$scenario == 'aUrU']
       , stablesizes_trts.big$density[stablesizes_trts.big$scenario == 'aLrH']
)
ks.test(stablesizes_trts.big$density[stablesizes_trts.big$scenario == 'aUrU']
        , stablesizes_trts.big$density[stablesizes_trts.big$scenario == 'aHrL']
)
ks.test(stablesizes_trts.big$density[stablesizes_trts.big$scenario == 'aLrH']
        , stablesizes_trts.big$density[stablesizes_trts.big$scenario == 'aHrL']
)
       
# ---- plot as boxplots ----

boxplot(density ~ scenario
        , data = stablesizes_trts.big
)

# ---- plot as densities ----


# aUrU
plot(yz
     , stablesizes_aUrU
     , xlab = "Body size (Operculum length, mm)"
     , ylab = "Probability density"
     , type = "l"
     #, ylim = c(0, 0.15)
)

# aLrH
lines(yz
     , stablesizes_aLrH
     , xlab = "Body size (Operculum length, mm)"
     , ylab = "Probability density"
     , type = "l"
     #, ylim = c(0, 0.15)
     , col = 'blue'
)

# aHrL
lines(yz
      , stablesizes_aHrL
      , xlab = "Body size (Operculum length, mm)"
      , ylab = "Probability density"
      , type = "l"
      #, ylim = c(0, 0.15)
      , col = 'red'
)

# ---- graveyard ----
# 
# ### other code
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
