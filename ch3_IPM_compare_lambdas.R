### Will King
### Ch 3 - population model for B. glandula, IPM

# ---- load data and packages ----
library(scales)

lambdas <- read.csv('ch3_IPM_compare_lambdas.csv', header = T)
lambdas.model <- lambdas[lambdas$type == 'IPM', ]

# ---- adjust data and calculate terms ----

# make continuous variable for easier plotting
lambdas.model$cont <- as.numeric(rownames(lambdas.model))

# calculate field mean and se with all data included
lambda.field.all_mean <- mean(lambdas$lambda[lambdas$type == 'field'])
lambda.field.all_se <- sd(lambdas$lambda[lambdas$type == 'field'])/sqrt(12)
lambda.field.all_serange <- c(lambda.field.all_mean - lambda.field.all_se # lower
                              , lambda.field.all_mean + lambda.field.all_se # upper
)

# calculate field mean and se with lambdas 2.00 or greater excluded (4 points excluded)
lambda.field.subset_mean <- mean(lambdas$lambda[lambdas$type == 'field' & lambdas$lambda < 2])
lambda.field.subset_se <- sd(lambdas$lambda[lambdas$type == 'field'& lambdas$lambda < 2])/sqrt(12)
lambda.field.subset_serange <- c(lambda.field.subset_mean - lambda.field.subset_se # lower
                              , lambda.field.subset_mean + lambda.field.subset_se # upper
)

# ---- plot: vertical style ----

## blank plot
plot(cont ~ lambda
     , data = lambdas.model
     , ylim = c(max(cont), min(cont)-2)
     , xlim = c(0.5, 2.5)
     , type = 'n'
     , axes = F
     , xlab = ''
     , ylab = ''
)

## add line of zero population growth
lines(x = c(1, 1)
      , y = c(min(lambdas.model$cont)-3, max(lambdas.model$cont)+1)
      , col = 'gray'
      , lty = 2
)

# mean and se of all field lambdas
lines(y = c(min(lambdas.model$cont) - 2, min(lambdas.model$cont) - 2)
      , x = c(lambda.field.all_serange[1], lambda.field.all_serange[2])
)
points(y = min(lambdas.model$cont) - 2
       , x = lambda.field.all_mean
       , pch = 21
       , bg = 'white'
)

# mean and se of field lambdas < 2
lines(y = c(min(lambdas.model$cont) - 1, min(lambdas.model$cont) - 1)
      , x = c(lambda.field.subset_serange[1], lambda.field.subset_serange[2])
)
points(y = min(lambdas.model$cont) - 1
       , x = lambda.field.subset_mean
       , pch = 21
       , bg = 'white'
)

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
       , xlim = c(0.5, 2.5)
       , ylim = c(min(cont)-2, max(cont))
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
     , at = seq(from = min(lambdas.model$cont) - 3
                , to = max(lambdas.model$cont) + 1
                , by = 1)
     , pos = 0.5
     , tck = 0.02
     , labels = c('', 'field all', 'field subset', as.character(lambdas.model$scenario), '')
     , las = 2
)
axis(3
     , las = 1
     , pos = min(lambdas.model$cont) - 3
     , tck = 0.02
     , labels = F
)

axis(4
     , at = seq(from = min(lambdas.model$cont) - 3
                , to = max(lambdas.model$cont) + 1
                , by = 1)
     , pos = 2.5
     , tck = 0.02
     , labels = F
)


mtext('Lambda'
      , side = 1
      , line = 3.5
      , cex = 1.2
)

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
