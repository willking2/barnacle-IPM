### Will King
### Ch 3 - population model for B. glandula, IPM

# compare total elasticities of aLrH and aHrL

## need to run models and perturbation analyses from main file first

# ---- plot both figs in one figure ----



pdf('plots/compare_elasticities_full.pdf', width = 8, height = 5)

par(mfrow = c(1, 2)
    , cex = 1.2
    , oma = c(4, 4, 1, 0.5)
    , mar = c(0, 0.5, 1, 0.5)
)

### aLrH

image(yz
      , yt
      , total.elas_aLrH
      , col = grey(seq(0.5, 1, length=100))
      , xlab = ''
      , ylab = ''
      , axes = F
)
contour(yz
        , yt
        , total.elas_aLrH
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

## axes and labels
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
)

mtext('Body size (mm)'
      , side = 1
      , line = 2.5
      , cex = 1.2
)

mtext('Crowding'
      , side = 2
      , line = 2.75
      , cex = 1.2)

mtext('a) eLrH'
      , side = 3
      , line = 0.5
      , cex = 1.2
      , adj = 0
)

### aHrL
image(yz
      , yt
      , total.elas_aHrL
      , col = grey(seq(0.5, 1, length=100))
      , xlab = ''
      , ylab = ''
      , axes = F
)
contour(yz
        , yt
        , total.elas_aHrL
        , add = TRUE
        , nlevels = 6
        , labcex = 0.8
)

### axes and labels
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
)

mtext('Body size (mm)'
      , side = 1
      , line = 2.5
      , cex = 1.2
)

mtext('b) eHrL'
      , side = 3
      , line = 0.5
      , cex = 1.2
      , adj = 0
)

dev.off()