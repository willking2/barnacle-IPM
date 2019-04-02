### Will King
### Ch 3 - population model for B. glandula, IPM

# need to run models from ch3_IPM_20190225 main file before able to plot these

### compare size distributions of aUrU, aLrH, and aHrL

# aUrU
plot(yz
     , (stable.state_aUrU.vector/sum(stable.state_aUrU.vector))*10
     , xlab = "Body size (Operculum length, mm)"
     , ylab = "Probability density"
     , type = "l"
     , ylim = c(0, 0.11)
)

# aLrH
lines(yz
     , (stable.state_aLrH.vector/sum(stable.state_aLrH.vector))*10
     , xlab = "Body size (Operculum length, mm)"
     , ylab = "Probability density"
     , type = "l"
     , ylim = c(0, 0.11)
     , col = 'blue'
)

# aHrL
lines(yz
      , (stable.state_aHrL.vector/sum(stable.state_aHrL.vector))*10
      , xlab = "Body size (Operculum length, mm)"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 0.11)
      , col = 'red'
)



### other code

plot(yz
     , (stable.state_aFrU.vector/sum(stable.state_aFrU.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Probability density"
     , type = "l"
     , ylim = c(0, 0.11)
)
# varying adult crowding
lines(yz
      , (stable.state_aLrU.vector/sum(stable.state_aLrU.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'blue'
)
lines(yz
      , (stable.state_aHrU.vector/sum(stable.state_aHrU.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'green'
)
lines(yz
      , (stable.state_aMrU.vector/sum(stable.state_aMrU.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'red'
)
lines(yz
      , (stable.state_aUrU.vector/sum(stable.state_aUrU.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'purple'
)

# varying recruit crowding
lines(yz
      , (stable.state_aUrL.vector/sum(stable.state_aUrL.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'lightblue'
      , lwd = 2
)
lines(yz
      , (stable.state_aUrH.vector/sum(stable.state_aUrH.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'green'
      , lwd = 2
)
lines(yz
      , (stable.state_aUrM.vector/sum(stable.state_aUrM.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'orange'
      , lwd = 2
)
# varying adult and recruit crowding
lines(yz
      , (stable.state_aLrL.vector/sum(stable.state_aLrL.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'blue'
      , lwd = 2
)
lines(yz
      , (stable.state_aLrH.vector/sum(stable.state_aLrH.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'purple'
      , lwd = 2
)
lines(yz
      , (stable.state_aHrH.vector/sum(stable.state_aHrH.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'green'
      , lwd = 2
)
lines(yz
      , (stable.state_aHrL.vector/sum(stable.state_aHrL.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'red'
      , lwd = 2
)
lines(yz
      , (stable.state_aMrM.vector/sum(stable.state_aMrM.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Probability density"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'orange'
      , lwd = 2
)
