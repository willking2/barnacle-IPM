

pdf('plots/compare_elasticities_break.pdf', width = 8, height = 5)

nl <- 6

par(mfrow = c(2, 2)
    , cex = 1.2
    , oma = c(4, 4, 1.5, 1)
    , mar = c(0, 0.5, 0.2, 1)
)

############# PANEL 1 (TOPLEFT): aLrH top half


image(yz
      , yt[40:50]
      , total.elas_aLrH[ , 40:50]
      , col = grey(seq(0.5, 1, length=100))
      , xlab = ''
      , ylab = ''
      , axes = F
)
contour(yz
        , yt[40:50]
        , total.elas_aLrH[ , 40:50]
        , add = TRUE
        , nlevels = nl
        , labcex = 0.8
)

# axis(1, 
#      pos = 0.78,
#      tck = -0.02, 
#      at = c(0, 2, 4, 6, 8, 10)
#      #, labels = c(0, 2, 4, 6, 8, 10)
#      , labels = F
# )
axis(2
     , pos = 0
     , las = 2
     , tck = -0.02
     , at = seq(0.7, 1.0, by = 0.1)
)
axis(3
     , pos = 1
     , at = c(0, 2, 4, 6, 8, 10)
     , labels = F
     , tck = -0.02
)
axis(4
     , pos = 10
     , labels = F
     , tck = -0.02
     , at = seq(0.7, 1.0, by = 0.1)
)

mtext('Crowding'
      , side = 2
      , line = 2.75
      , cex = 1.2
      , outer = T
)

mtext('a) eLrH'
      , side = 3
      , line = 0.5
      , cex = 1.2
      , adj = 0
)

############# PANEL 2 (TOP RIGHT): aHrL top half


image(yz
      , yt[40:50]
      , total.elas_aHrL[ , 40:50]
      , col = grey(seq(0.5, 1, length=100))
      , xlab = ''
      , ylab = ''
      , axes = F
)
contour(yz
        , yt[40:50]
        , total.elas_aHrL[ , 40:50]
        , add = TRUE
        , nlevels = nl
        , labcex = 0.8
)

# axis(1, 
#      pos = 0.78,
#      tck = -0.02, 
#      at = c(0, 2, 4, 6, 8, 10)
#      #, labels = c(0, 2, 4, 6, 8, 10)
#      , labels = F
# )
axis(2
     , pos = 0
     , las = 2
     , tck = -0.02
     , labels = F
     , at = seq(0.7, 1.0, by = 0.1)
)
axis(3
     , pos = 1
     , at = c(0, 2, 4, 6, 8, 10)
     , labels = F
     , tck = -0.02
)
axis(4
     , pos = 10
     , labels = F
     , tck = -0.02
     , at = seq(0.7, 1.0, by = 0.1)
)

mtext('b) eHrL'
      , side = 3
      , line = 0.5
      , cex = 1.2
      , adj = 0
)


############# PANEL 3 (BOTTOM LEFT): aLrH bottom half

image(yz
      , yt[1:10]
      , total.elas_aLrH[, 1:10]
      , col = grey(seq(0.5, 1, length=100))
      , xlab = ''
      , ylab = ''
      , axes = F
)
contour(yz
        , yt[1:10]
        , total.elas_aLrH[, 1:10]
        , add = TRUE
        , nlevels = nl
        , labcex = 0.8
)

axis(1, 
     pos = 0,
     tck = -0.02, 
     at = c(0, 2, 4, 6, 8, 10),
     labels = c(0, 2, 4, 6, 8, 10)
)
axis(2
     , pos = 0
     , las = 2
     , tck = -0.02
     , at = seq(0, 0.2, 0.1)
)
# axis(3
#      , pos = 0.2
#      , at = c(0, 2, 4, 6, 8, 10)
#      , labels = F
#      , tck = -0.02
# )
axis(4
     , pos = 10
     , labels = F
     , tck = -0.02
     , at = seq(0.0, 0.2, by = 0.1)
)


mtext('Body size (mm)'
      , side = 1
      , line = 2.5
      , cex = 1.2
      # , outer = T
)



############# PANEL 4 (BOTTOM RIGHT): aHrL bottom half


image(yz
      , yt[1:10]
      , total.elas_aHrL[ , 1:10]
      , col = grey(seq(0.5, 1, length=100))
      , xlab = ''
      , ylab = ''
      , axes = F
)
contour(yz
        , yt[1:10]
        , total.elas_aHrL[ , 1:10]
        , add = TRUE
        , nlevels = nl
        , labcex = 0.8
)

axis(1, 
     pos = 0,
     tck = -0.02, 
     at = c(0, 2, 4, 6, 8, 10),
     labels = c(0, 2, 4, 6, 8, 10)
)
axis(2
     , pos = 0
     , las = 2
     , tck = -0.02
     , labels = F
     , at = seq(0.0, 2.0, by = 0.1)
)
# axis(3
#      , pos = 0.2
#      , at = c(0, 2, 4, 6, 8, 10)
#      , labels = F
#      , tck = -0.02
# )
axis(4
     , pos = 10
     , labels = F
     , tck = -0.02
)

mtext('Body size (mm)'
      , side = 1
      , line = 2.5
      , cex = 1.2
)

dev.off()
