plot(yz
     , (stable.state_aFrU.vector/sum(stable.state_aFrU.vector))*10
     , xlab = "Operculum length, mm"
     , ylab = "Frequency"
     , type = "l"
     , ylim = c(0, 1)
)
lines(yz
      , (stable.state_aLrU.vector/sum(stable.state_aLrU.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Frequency"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'blue'
)
lines(yz
      , (stable.state_aHrU.vector/sum(stable.state_aHrU.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Frequency"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'green'
)
lines(yz
      , (stable.state_aMrU.vector/sum(stable.state_aMrU.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Frequency"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'red'
)
lines(yz
      , (stable.state_aUrU.vector/sum(stable.state_aUrU.vector))*10
      , xlab = "Operculum length, mm"
      , ylab = "Frequency"
      , type = "l"
      , ylim = c(0, 1)
      , col = 'purple'
)
