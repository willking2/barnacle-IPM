# space to build ceilings into functions
# only doing size ceilings, no touch ceilings (makes things too complicated)

### manual ceiling example, from IPM_simple ####

s_z_simple <- function(z, m.par_simple){
  
  ## linear predictor
  linear.p <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * z # for cases below ceiling
  # linear.p.ceiling <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * U1 # for cases > U1
  
  ## back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  
  
  ## output
  return( unname(p) )
  
  # # code for manually doing ceiling
  # p <- ifelse(z <= U1
  #             , 1/(1+exp(-linear.p)) # below ceiling; use regression coeff of z
  #             , 1/(1+exp(-linear.p.ceiling)) # above ceiling; use vital rates of U1
  # )
}

### ceilings for IPM size-touch~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

# ceiling for size is Uz1 = 6

### IPM size touch: survival with ceiling ####

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  
  # below size ceiling
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z
  # above size ceiling
  linear.p.sizeceiling <- m.par_st['surv.int'] + m.par_st['surv.z'] * Uz1
  
  # impose ceiling and backtransform from logit
  p <- ifelse(z <= Uz1
              , 1/(1+exp(-linear.p)) # below ceiling, use regression coeff of z
              , 1/(1+exp(-linear.p.sizeceiling)) # above ceiling, use coeff of Uz1
  )
  
  ## output
  return( unname(p) )
  
  # ## back transform from logit
  # p <- 1/(1+exp(-linear.p))
}

# test
s_z_st(3, m.par_st)
s_z_st(6, m.par_st)
s_z_st(7, m.par_st)


### IPM size touch: Growth with ceiling ####


G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  # below size ceiling
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t
  # above size ceiling
  mu.sizeceiling <- m.par_st['grow.int'] + m.par_st['grow.z'] * Uz1 + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * Uz1 * t

  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  # have to loop b/c it's a PDF, not a probability
  p.den.grow <- c()
  for(q in 1:length(z1)){
    p.den.grow[q] <- ifelse(z <= Uz1
                            , dnorm(z1[q], mean = mu, sd = sig) # below ceiling
                            , dnorm(z1[q], mean = mu.sizeceiling, sd = sig) # above ceiling
    )
  }
  
  ## output
  return(p.den.grow)
  
}

# test
G_z1zt_st(3, 2, 0.5, m.par_st)
G_z1zt_st(6, 6, 0.5, m.par_st)
G_z1zt_st(6, 7, 0.5, m.par_st)
G_z1zt_st(7, 6, 0.5, m.par_st)
G_z1zt_st(7, 7, 0.5, m.par_st)

wombat <- c(1, 2, 3)
meerkat <- c(2, 3, 4)
G_z1zt_st(c(2, 3, 4), 1, 0.5, m.par_st)
G_z1zt_st(c(6, 6, 6), 6, 0.5, m.par_st)
G_z1zt_st(c(6, 6, 6), 7, 0.5, m.par_st)


### IPM size touch: probability reproduction with ceiling ####

p_bz_st <- function(z, t, m.par_st){
  ### convert z to v
  v <- ifelse(z <= Uz1
              , z * 131.18 - 152.52 # below ceiling
              , Uz1 * 131.18 - 152.52 # above ceiling
  )
  
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  # output
  return( unname(p) )
}

# test
p_bz_st(3, 0.5, m.par_st)
p_bz_st(3, 1, m.par_st)
p_bz_st(6, 0.5, m.par_st)
p_bz_st(7, 0.5, m.par_st)


### archive: functions without ceilings ####

s_z_st <- function(z, m.par_st){
  
  ## linear predictor
  linear.p <- m.par_st['surv.int'] + m.par_st['surv.z'] * z # for cases below ceiling
  # linear.p.ceiling <- m.par_simple['surv.int'] + m.par_simple['surv.z'] * U1 # for cases > U1
  ## Note: ceiling is implemented in mk_K, below. just keeping code here to see it.
  
  ## back transform from logit
  p <- 1/(1+exp(-linear.p))
  
  
  
  ## output
  return( unname(p) )
  
  # # code for manually doing ceiling
  # p <- ifelse(z <= U1
  #             , 1/(1+exp(-linear.p)) # below ceiling; use regression coeff of z
  #             , 1/(1+exp(-linear.p.ceiling)) # above ceiling; use vital rates of U1
  # )
}

G_z1zt_st <- function(z1, z, t, m.par_st){
  
  ## based on growth regression coefficients, define mu and sigma of growth PDF
  mu <- m.par_st['grow.int'] + m.par_st['grow.z'] * z + m.par_st['grow.t'] * t + m.par_st['grow.zt'] * z * t # below ceiling
  #mu.ceiling <- m.par_simple['grow.int'] + m.par_simple['grow.z'] * U1 # above ceiling
  ## Note: ceiling is implemented in mk_K, below. just keeping code here to see it.
  sig <- m.par_st['grow.sd']
  
  ## calculate growth PDF
  
  p.den.grow <- dnorm(z1, mean = mu, sd = sig)
  
  ## output
  return(p.den.grow)
  
  
  
  
  ## graveyard
  
  # p.den.grow <- dnorm(z1, mean = mu, sd = sig)
  
  # #code for prevent shrinking. makes model output 'wonky'. not using.
  # p.den.grow <- ifelse(z1 > z # condition
  #                      , dnorm(z1, mean = mu, sd = sig) # if true. pdf of new size z1
  #                      , 0 # if false. forbid barnacles from shrinking
  # )
  
  # ##code for doing ceiling manually within function. However, do not use, b/c it messes up
  ##the numerical implementation of the kernel. Do the ceiling when forming mk_k.
  # # make blank vector
  # p.den.grow <- c()
  # 
  # loop through z1
  #   for( j in 1:length(z1) ){
  #     p.den.grow[j] <- ifelse(z <= U1
  #                             , dnorm(z1[j], mean = mu, sd = sig) # below ceiling
  #                             , dnorm(z1[j], mean = mu.ceiling, sd = sig) # above ceiling
  #     )
  #   }
  # }
  
}

p_bz_st <- function(z, t, m.par_st){
  # convert z to v
  v <- z * 131.18 - 152.52
  # linear predictor
  linear.p <- m.par_st['repr.int'] + m.par_st['repr.v'] * v + m.par_st['repr.t'] * t
  # back transform from logit
  p <- 1/(1+exp(-linear.p))
  # output
  return( unname(p) )
}

