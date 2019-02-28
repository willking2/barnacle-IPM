

### manual ceiling example, from IPM_simple

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

### survival: s(z)

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
