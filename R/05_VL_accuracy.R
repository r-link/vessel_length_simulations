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
# a) getter function for true distribution based on a set of parameters 
#    (i.e., a row of the results table)
.truedist <- function(parm) {
  if(parm$distribution == "exponential"){
    return(function(x, parm) dexp(x, rate = 1 / parm$truemean))}
  if(parm$distribution == "erlang2"){
    return(function(x, parm) dgamma(x, shape = 2, rate = 2 / parm$truemean))} 
  if(parm$distribution == "weibull"){
    return(function(x, parm) dweibull(x, shape = parm$trueshape, scale = parm$truemean / 
                                        gamma(1 + 1 / parm$trueshape)))} 
  if(parm$distribution == "gamma"){
    return(function(x, parm) dgamma(x, shape = parm$trueshape, rate = parm$trueshape / parm$truemean))} 
  if(parm$distribution == "lnorm"){
    return(function(x, parm) dlnorm(x, meanlog = log(parm$truemean) - parm$truesdlog ^ 2 / 2, 
                                    sdlog = parm$truesdlog))}
}

# b) getter function for predicted distribution based on a set of parameters
.preddist <- function(parm) {
  if(parm$model == "exponential"){
    return(function(x, parm) dexp(x, rate = 1 / parm$mean))}
  if(parm$model == "erlang2"){
    return(function(x, parm) dgamma(x, shape = 2, rate = 2 / parm$mean))} 
  if(parm$model == "weibull"){
    return(function(x, parm) dweibull(x, shape = parm$shape, scale = parm$mean / gamma(1 + 1 / parm$shape)))} 
  if(parm$model == "gamma"){
    return(function(x, parm) dgamma(x, shape = parm$shape, rate = parm$shape / parm$mean))} 
  if(parm$model == "lnorm"){
    return(function(x, parm) dlnorm(x, meanlog = log(parm$mean) - parm$sdlog ^ 2 / 2, sdlog = parm$sdlog))} 
  if(parm$model == "cohen"){
    return(function(x, parm)  dgamma(x, shape = 2, rate = 2 / parm$mean))} 
  if(parm$model == "christman"){
    return(function(x, parm)  Fchristman1(x, c = parm$christman_c, k = parm$christman_k))} 
}

# 3. Percent overlap between two distributions ---------------------------------
percOL <- function(parm){
  # true distribution
  fx    <- .truedist(parm)  
  # predicted distribution
  fhatx <- .preddist(parm)
  # overlap between distributions: integral over minimum of their density functions
  return(integrate( function(x) 100 * pmin(fhatx(x), fx(x)), lower = 0, upper = Inf)$value)}    


