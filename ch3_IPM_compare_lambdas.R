### Will King
### Ch 3 - population model for B. glandula, IPM

# ---- load data and packages ----
library(scales)

lambdas <- read.csv('ch3_IPM_compare_lambdas.csv', header = T)
lambdas.model <- lambdas[lambdas$type == 'IPM', ]

# ---- adjust data and calculate terms ----

## order lambdas.model
# display data in order of: field, then uniform, then from lowest to highest lambda

moop <- lambdas.model[-c(1:3), ]
meep <- moop[order(moop$lambda), ]
lambdas.model <- rbind(lambdas.model[1:3, ], meep)
rm(moop, meep)

# make continuous variable for easier plotting
lambdas.model$cont <- seq(1:length(lambdas.model$lambda))

# calculate field mean and se with all data included
lambda.field.all_mean <- mean(lambdas$lambda[lambdas$type == 'field']) # 1.53
lambda.field.all_sd <- sd(lambdas$lambda[lambdas$type == 'field']) # 0.98
lambda.field.all_se <- sd(lambdas$lambda[lambdas$type == 'field'])/sqrt(12)
lambda.field.all_serange <- c(lambda.field.all_mean - lambda.field.all_se # lower
                              , lambda.field.all_mean + lambda.field.all_se # upper
)

# # calculate field mean and se with lambdas outside 1 sd of mean excluded
# # mean is 2.02, SD is 1.41, so range to be included is 0.6 to 3.43
# # to exclude: 3.73 (FH1 H), 0.57 (FH2 M), and 5.28 (RT3 L)
# 
# lambdas.subset <- lambdas[lambdas$lambda != 3.73 &
#                                   lambdas$lambda != 0.57 &
#                                   lambdas$lambda != 5.28
#                                 , ]
# lambda.field.subset_mean <- mean(lambdas.subset$lambda)
# lambda.field.subset_se <- sd(lambdas.subset$lambda)/sqrt(9)
# lambda.field.subset_serange <- c(lambda.field.subset_mean - lambda.field.subset_se # lower
#                               , lambda.field.subset_mean + lambda.field.subset_se # upper
# )


# ---- plot: vertical style ----

pdf('plots/compare_lambdas.pdf', width = 5, height = 7)
par(xpd = NA
    , cex = 1.2
    , mar = c(4, 5, 1, 0.5)
)

## blank plot
plot(cont ~ lambda
     , data = lambdas.model
     , ylim = c(max(cont)+1, min(cont)-2)
     , xlim = c(0.5, 2.0)
     , type = 'n'
     , axes = F
     , xlab = ''
     , ylab = ''
)

## add line of zero population growth
lines(x = c(1, 1)
      , y = c(min(lambdas.model$cont)-2, max(lambdas.model$cont)+1)
      , col = 'gray'
      #, lty = 2
)

# ## add line of aUrU comparison
# lines(x = c(1.17, 1.17)
#       , y = c(3, max(lambdas.model$cont)+1)
#       , col = 'gray'
#       #, lty = 2
# )

# ## add line separating comparisons you want to make
# lines(x = c(0.5, 2.5)
#       , y = c(2.5, 2.5)
#       # , col = 'gray'
#       # , lty = 2
# )

# mean and se of all field lambdas
lines(y = c(min(lambdas.model$cont) - 1, min(lambdas.model$cont) - 1)
      , x = c(lambda.field.all_serange[1], lambda.field.all_serange[2])
)
points(y = min(lambdas.model$cont) - 1
       , x = lambda.field.all_mean
       , pch = 21
       , bg = 'white'
)

# # mean and se of field lambdas < 2
# lines(y = c(min(lambdas.model$cont) - 1, min(lambdas.model$cont) - 1)
#       , x = c(lambda.field.subset_serange[1], lambda.field.subset_serange[2])
# )
# points(y = min(lambdas.model$cont) - 1
#        , x = lambda.field.subset_mean
#        , pch = 21
#        , bg = 'white'
# )

## 95% CIs of model lambdas
for (i in min(lambdas.model$cont):max(lambdas.model$cont)){
  lines(y = c(i, i)
        , x = c(lambdas.model$confint95_lower[lambdas.model$cont == i]
                , lambdas.model$confint95_upper[lambdas.model$cont == i])
  )
}

## model lambdas in black
points(cont ~ lambda
       , data = lambdas.model
       , xlim = c(0.5, 2.0)
       , ylim = c(min(cont)-1, max(cont))
       , col = 'black'
       , pch = 19
)


# axes and axes labels

axis(1
     , las = 1
     , pos = max(lambdas.model$cont) + 1
     , tck = 0.02
)
axis(2 
     , at = seq(from = min(lambdas.model$cont) - 2
                , to = max(lambdas.model$cont) + 1
                , by = 1)
     , pos = 0.5
     , tck = 0.02
     , labels = c('', 'field', as.character(lambdas.model$scenario), '')
     , las = 2
)
axis(3
     , las = 1
     , pos = min(lambdas.model$cont) - 2
     , tck = 0.02
     , labels = F
)

axis(4
     , at = seq(from = min(lambdas.model$cont) - 2
                , to = max(lambdas.model$cont) + 1
                , by = 1)
     , pos = 2
     , tck = 0.02
     , labels = F
)


mtext( expression('Population growth rate ('*lambda*')')
      , side = 1
      , line = 2
      , cex = 1.2
)

mtext('IPMs'
      , side = 2
      , line = 3.6
      , cex = 1.2
)

lines(y = c(min(lambdas.model$cont), 4.7)
      , x = c(0.005, 0.005)
)

lines(y = c(6.3, max(lambdas.model$cont))
      , x = c(0.005, 0.005)
)

dev.off()

# ---- graveyard ----
# # ---- plot: horizontal style ----
# 
# ## blank plot
# plot(lambda ~ cont
#      , data = lambdas.model
#      , ylim = c(0.5, 2.5)
#      , xlim = c(min(cont)-2, max(cont))
#      , type = 'n'
#      , axes = F
#      , xlab = ''
#      , ylab = ''
#      )
# 
# ## add line of zero population growth
# lines(x = c(min(lambdas.model$cont)-3, max(lambdas.model$cont)+1)
#       , y = c(1, 1)
#       , col = 'gray'
#       , lty = 2
# )
# 
# # mean and se of all field lambdas
# lines(x = c(min(lambdas.model$cont) - 2, min(lambdas.model$cont) - 2)
#       , y = c(lambda.field.all_serange[1], lambda.field.all_serange[2])
# )
# points(x = min(lambdas.model$cont) - 2
#        , y = lambda.field.all_mean
#        , pch = 21
#        , bg = 'white'
# )
# 
# # mean and se of field lambdas < 2
# lines(x = c(min(lambdas.model$cont) - 1, min(lambdas.model$cont) - 1)
#       , y = c(lambda.field.subset_serange[1], lambda.field.subset_serange[2])
# )
# points(x = min(lambdas.model$cont) - 1
#        , y = lambda.field.subset_mean
#        , pch = 21
#        , bg = 'white'
# )
# 
# ## 95% CIs of model lambdas
# for (i in min(lambdas.model$cont):max(lambdas.model$cont)){
#   lines(x = c(i, i)
#         , y = c(lambdas.model$confint95_lower[lambdas.model$cont == i]
#                 , lambdas.model$confint95_upper[lambdas.model$cont == i])
#   )
# }
# 
# ## model lambdas in black
# points(lambda ~ cont
#        , data = lambdas.model
#        , ylim = c(0.5, 2.5)
#        , xlim = c(min(cont)-2, max(cont))
#        , col = 'black'
#        , pch = 19
# )
# 
# 
# # axes and axes labels
# axis(1 
#      , at = seq(from = min(lambdas.model$cont) - 3
#                 , to = max(lambdas.model$cont) + 1
#                 , by = 1)
#      , pos = 0.5
#      , tck = 0.02
#      , labels = c('', 'field all', 'field subset', as.character(lambdas.model$scenario), '')
#      , las = 2
# )
# axis(2
#      , las = 1
#      , pos = min(lambdas.model$cont) - 3
#      , tck = 0.02
# )
# axis(3
#      , at = seq(from = min(lambdas.model$cont) - 3
#                 , to = max(lambdas.model$cont) + 1
#                 , by = 1)
#      , pos = 2.5
#      , tck = 0.02
#      , labels = F
# )
# axis(4
#      , las = 1
#      , pos = max(lambdas.model$cont) + 1
#      , tck = 0.02
#      , labels = F
# )
# 
# 
# mtext('Lambda'
#       , side = 2
#       , line = 3
#       , cex = 1.2
# )

# ---- subscripts ----
# ## make vector of nice axes labels
# lamlabels <- c(''
#                , 'field all'
#                , 'field subset'
#                , expression(IPM[size-only])
#                , expression(IPM[aFrF])
#                , expression(IPM[aUrU])
#                , expression(IPM[aLrH])
#                , expression(IPM[aUrH])
#                , expression(IPM[aHrH])
#                , expression(IPM[aLrU])
#                , expression(IPM[aHrU])
#                , expression(IPM[aLrL])
#                , expression(IPM[aUrL])
#                , expression(IPM[aHrL])
#                , ''
# )
