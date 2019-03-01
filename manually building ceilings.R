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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### for IPM size-touch
# ceilings are Uz1 = 6 and Ut1 = 0.8, respectively

### IPM size touch: survival ####

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


### IPM size touch: Growth ####


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

