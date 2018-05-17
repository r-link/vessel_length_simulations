################################################################################
################################################################################
########
########         Utility functions for the maximum likelihood estimation
########                       of xylem vessel lengths
########
########         Script: Roman Link, 2017         Contact: rlink@gwdg.de
########
################################################################################
################################################################################

# This R script contains convenience functions for the maximum likelihood 
# estimation of vessel length distributions under the assumption that the 
# vessel lengths follow one of the following distributions: 
# - exponential (_exponential_)
# - Erlang(2) (_erlang2_) [basically a gamma distribution with shape = 2]
# - Weibull (_weibull_)
# - Gamma (_gamma_)
# - Lognormal (_lnorm_)

# Provided are functions both for a standard parameterization (with a rate para-
# meter for exponential, Erlang(2)-, Weibull-, and Gamma-Distribution) and a
# reparameterized version using the expected value (_mu) instead of the rate.
# For the lognormal distribution, the standard distribution uses natural log-
# transformed mean and standard deviation (meanlog/sdlog), the reparameterized  
# version replaces meanlog by the expected value on the exponentiated scale.

# The functions are named using the following nomenclature (where _dist_ refers
# to the distribution name given in brackets above):

# f_dist(_mu)   Probability density of uncut vessel lengths
# g_dist(_mu)   Size biased probability density of vessel lenghts of the 
#                    vessels that are situated at the cutting plane
# h_dist(_mu)   Probability density of the lenghts of the cut-off vessel
#                    fragments
# p_dist(_mu)   Expected proportion of open vessels (infusion profile)

################################################################################
########  f(x): Probability density function of distribution of uncut vessels
################################################################################
# Unconditional probability density functions of vessel length 
# (overall vessel length distribution).

# In the standard parameterization, this is equivalent to the standard density
# functions (dexp(), dweibull(), dgamma() etc.), and for Weibull, Gamma and 
# lognormal distribution, the code takes advantage of these functions.

## Exponential distribution
f_exponential    <-function(x, rate) rate * exp(-rate * x)
f_exponential_mu <-function(x, mu)  exp(-x / mu) / mu

## Erlang(2) distribution
f_erlang2    <-function(x, rate)  rate^2 * x * exp(-rate * x)
f_erlang2_mu <-function(x, mu)  (4/mu^2) * x * exp(-2 * x/mu)

## Weibull distribution
f_weibull    <-function(x, shape, rate)  
                 dweibull(x, shape = shape, scale = 1/rate)
f_weibull_mu <-function(x, shape, mu) 
                 dweibull(x, shape = shape, scale = mu / gamma(1 + 1/shape))  
                
## Gamma distribution
f_gamma    <-function(x, shape, rate) dgamma(x, shape = shape, rate = rate)
f_gamma_mu <-function(x, shape, mu) dgamma(x, shape = shape, rate = shape/mu)   

## Lognormal distribution
f_lnorm    <-function(x, meanlog, sdlog) dlnorm(x, meanlog = meanlog, sdlog = sdlog)
f_lnorm_mu <-function(x, mu, sdlog) dlnorm(x, meanlog = (log(mu) - sdlog^2/2), 
                                           sdlog = sdlog)



################################################################################
########    g(x): Size biased length distributions at the cut surface
################################################################################
# Conditional distribution of vessel length given that a vessel crosses the
# injection point (size biased vessel length distribution).
# Calculated as the density of the unbiased distribution multiplied with 
# vessel length, divided by average vessel length.
  
      
## Exponential distribution
g_exponential    <-function(x, rate) rate^2 * x * exp(-rate * x)
g_exponential_mu <-function(x, mu)   (1/mu^2) * x * exp(-x/mu)

## Erlang(2) distribution
g_erlang2     <-function(x, rate) rate^3 * x^2/2 * exp(- rate * x)
g_erlang2_mu  <-function(x, mu)   (4/mu^3) * x^2 * exp(-2 * x/mu)

## Weibull distribution
g_weibull    <-function(x, shape, rate)  
                 x * rate/gamma(1 + 1/shape) * dweibull(x, shape = shape, 
                                                        scale = 1/rate)
g_weibull_mu <-function(x, shape, mu) 
                  x/mu * dweibull(x, shape = shape, scale = mu/gamma(1 + 1/shape)) 
                
## Gamma distribution
g_gamma    <-function(x, shape, rate) x * rate/shape * dgamma(x, shape = shape, 
                                                             rate = rate)
g_gamma_mu <-function(x, shape, mu) x/mu * dgamma(x, shape = shape,
                                                  rate = shape/mu)     
          
## Lognormal distribution
g_lnorm    <-function(x, meanlog, sdlog) x/(exp(meanlog + (sdlog^2)/2))*
               dlnorm(x, meanlog = meanlog, sdlog = sdlog)
g_lnorm_mu <-function(x, mu, sdlog) x/mu*
               dlnorm(x, meanlog = (log(mu) - sdlog^2/2), sdlog = sdlog)


################################################################################
########   h(z): Length distribution of cut-off vessel fragments
################################################################################
# Probability density function of the length distributions of the fraction of 
# each vessel that is injected with silicone.
# Given by the survivor function over the true length distribution with the 
# inverse of the average of the true length distribution as a normalizing 
# constant.

# For Gamma and lognormal distribution, the code takes advantage of the 
# fact that the survivor function can be obtained from the corresponding 
# cumulative distribution function (pgamma() and plnorm()) if specifying
# lower.tail = FALSE.

## Exponential distribution
h_exponential  <-function(z, rate) rate * exp(-rate * z)
h_exponential_mu  <-function(z, mu)  exp(-z/mu)/mu

## Erlang(2) distribution
h_erlang2     <-function(z, rate) (rate/2) * (rate * z + 1) * exp(-rate * z)
h_erlang2_mu  <-function(z, mu)   (1/mu) * (2 * z/mu +1) * exp(- z * 2/mu)

## Weibull distribution             
h_weibull    <-function(z, shape, rate)  (rate/(gamma(1 + 1/shape)))*
                 exp(-(rate * z)^shape)
h_weibull_mu <-function(z, shape, mu)  (1/mu) *
                  exp(-(z * gamma(1 + 1/shape)/mu)^shape)      
                  
## Gamma distribution
h_gamma <-function(z, shape, rate)  (rate/shape)*
             pgamma(z, shape = shape, rate = rate, lower.tail = FALSE)
h_gamma_mu <-function(z, shape, mu)  (1/mu)*
             pgamma(z, shape = shape, rate = shape/mu, lower.tail = FALSE)
             
## Lognormal distribution
h_lnorm <-function(z, meanlog, sdlog) (1/exp(meanlog + (sdlog^2)/2)) * 
             plnorm(z, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE)
h_lnorm_mu <-function(z, mu, sdlog) (1/mu) * 
             plnorm(z, meanlog = (log(mu) - (sdlog^2)/2), sdlog = sdlog, 
                    lower.tail=FALSE)


  
################################################################################
########   p(z): Infusion profiles
################################################################################
# Equations for the expected fraction of open vessels. 
# Calculated as the integral over the survivor function of the true length 
# distribution using the inverse of the true length distribution as a  
# normalizing constant.

# For Weibull and Gamma distribution, the code uses the fact that the upper  
# regularized gamma function can be calculated based on the survivor function
# of the gamma function via pgamma(x, a, lower.tail = FALSE).
# The lognormal distribution uses pnorm() with the standard arguments (mean = 0
# and sd = 1) to obtain the CDF of the standard normal distribution.

## Ezponential distribution
p_exponential    <- function(z, rate) exp(-rate*z)
p_exponential_mu <- function(z, mu)   exp(-z/mu)

## Erlang(2) distribution
p_erlang2     <-function(z, rate)  (rate/2) * (z + 2/rate) * exp(- rate * z)
p_erlang2_mu  <-function(z, mu)  (z/mu + 1) * exp(-2 * z/mu)

## Weibull distribution
p_weibull    <- function(z, shape, rate) pgamma((rate * z)^shape, 1/shape, 
                                                lower.tail = FALSE)
p_weibull_mu <- function(z, shape, mu) pgamma(exp((lgamma(1 + 1/shape)+
                      log( z)-log(mu))*shape), 1/shape, lower.tail = FALSE)
                      
## Gamma distribution
p_gamma <- function(z, shape, rate) { 
   pgamma(rate * z, shape + 1, lower.tail = F) - 
   (rate/shape) * z * pgamma(rate * z, shape)
   }
p_gamma_mu <- function(z, shape, mu) { 
   pgamma(shape*z/mu, shape+1, lower =F) -
   (z/mu) * pgamma(shape * z/mu, shape, lower.tail =F)
   }
  
## Lognormal distribution
p_lnorm <- function(z, meanlog, sdlog) {
   pnorm((meanlog + sdlog^2 - log(z))/(sdlog), 0, 1) - 
   z * pnorm((meanlog - log(z))/(sdlog), 0, 1)/(exp(meanlog + (sdlog^2)/2))
   }
p_lnorm_mu <- function(z, mu, sdlog) {
   pnorm((log(mu) - log(z) + (sdlog^2)/2)/(sdlog) ,0,1) -
   z/mu * pnorm((log(mu) - log(z) -(sdlog^2)/2)/(sdlog) ,0,1)
   }

