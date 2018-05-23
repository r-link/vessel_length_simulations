################################################################################
################################################################################
########
########  Functions for the analysis of the accuracy of vessel length estimates
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# Contents ---------------------------------------------------------------------
# 1. Functions for handling the Christman model
# 2. Internals for the evaluation of similarity between distributions 
# 3. Percent overlap between two distributions

# 1. Functions for handling the Christman model --------------------------------

# Christmans equation for the VL distribution
Fchristman <-function(x, c, k){c * k ^ c * x ^ (c - 1) * exp(-(k * x) ^ c) * (c * (k * x) ^ c - c + 1)}

# normalized to be a true pdf (by numerical integration)
Fchristman1 <-function(x, c, k){
  lower <- ifelse(c < 1, 0, (1 / k) * ((c - 1) / c) ^ (1 / c))                               
  ifelse(x< lower, 0, Fchristman(x, c, k) / 
           integrate(function(x) Fchristman(x, c, k), lower, Inf)$value) }


# 2. Internals for the evaluation of similarity between distributions ---------- 
# a) getter function for true distribution based on the model settings
.truedist <- function(settings) {
  if(settings$distribution == "exponential"){
    return(function(x) dexp(x, rate = 1 / settings$mean_true))}
  if(settings$distribution == "erlang2"){
    return(function(x) dgamma(x, shape = 2, rate = 2 / settings$mean_true))} 
  if(settings$distribution == "weibull"){
    return(function(x) dweibull(x, shape = settings$par2, scale = settings$mean_true / 
                                        gamma(1 + 1 / settings$par2)))} 
  if(settings$distribution == "gamma"){
    return(function(x) dgamma(x, shape = settings$par2, rate = settings$par2 / settings$mean_true))} 
  if(settings$distribution == "lnorm"){
    return(function(x) dlnorm(x, meanlog = log(settings$mean_true) - settings$par2 ^ 2 / 2, 
                                    sdlog = settings$par2))}
}

# b) getter function for predicted distribution based on a model object
.preddist <- function(mod, type) {
  parm <- coef(mod)
  if(type == "exponential"){
    return(function(x) dexp(x, rate = 1 / parm[[1]]))}
  if(type == "erlang2"){
    return(function(x) dgamma(x, shape = 2, rate = 2 / parm[[1]]))} 
  if(type == "weibull"){
    return(function(x) dweibull(x, shape = parm[[2]], scale = parm[[1]] / gamma(1 + 1 / parm[[2]])))} 
  if(type == "gamma"){
    return(function(x) dgamma(x, shape = parm[[2]], rate = parm[[2]] / parm[[1]]))} 
  if(type == "lnorm"){
    return(function(x) dlnorm(x, meanlog = log(parm[[1]]) - parm[[2]] ^ 2 / 2, sdlog = parm[[2]]))} 
  if(type == "cohen"){
    return(function(x)  dgamma(x, shape = 2, rate = 2 / parm[[1]]))} 
  if(type == "christman"){
    return(function(x)  Fchristman1(x, c = parm[[1]], k = parm[[2]]))} 
}

# 3. Percent overlap between two distributions ---------------------------------
percOL <- function(mod, type, settings){
  # true distribution
  fx    <- .truedist(settings)  
  # predicted distribution
  fhatx <- .preddist(mod, type)
  # overlap between distributions: integral over minimum of their density functions
  return(integrate(function(x) 100 * pmin(fhatx(x), fx(x)), lower = 0, upper = Inf)$value)
  }    


