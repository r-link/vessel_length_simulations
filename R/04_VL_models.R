################################################################################
################################################################################
########
########           Model fitting for vessel length distributions
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# Contents ---------------------------------------------------------------------
# 1. Internals   - getter functions for starting values and lower bounds for 
#                  bounded estimation
# 2. mle_VL()    - maximum likelihood estimation of vessel length distributions
# 3. fit_mod()   - function operator that either passes its arguments to mle_VL,
#                  or fits simple (non)linear models based on the Cohen or 
#                  Christman methods
# 4. get_coefs() - getter function for the coefficients of the fitted models

# 1. Internals -----------------------------------------------------------------
# a) .get_start() - getter function for starting values
.get_start <- function(z, counts, dist){
  # return empty starting values if dist == "cohen"
  if (dist == "cohen") return(NULL)
  
  # crude estimate of average vessel length: x intercept of line through first 
  # two counts (exluding counts that do not differ from first count)
  which <- c(1, which(counts < counts[1])[1])
  coefs <- coef(lm(counts[which] ~ z[which]))
  mu_start <- - coefs[1]/coefs[2]
  
  # return starting values depending on distribution
  if (dist %in% c("exp", "erlang2")) return(list(mu = mu_start))
  if (dist %in% c("weibull", "gamma")) return(list(mu = mu_start, shape = 1))
  if (dist == "lnorm") return(list(mu = mu_start, sdlog = 1))
  if (dist  == "christman") return(list(k = 1/mu_start, c = 1))
}

# b) .get_lower() - getter function for lower bounds in bounded estimation
.get_lower <- function(dist){
  # return nothing for literature models
  if (dist %in% c("cohen", "christman")) return(NULL)
  
  # return starting values depending on distribution
  if (dist %in% c("exp", "erlang2")) return(list(mu = 0.1))
  if (dist %in% c("weibull", "gamma")) return(list(mu = 0.1, shape = 0.01))
  if (dist == "lnorm") return(list(mu = 0.1, sdlog = 0.01))
}


# 2. mle_VL() ------------------------------------------------------------------
# function for the maximum likelihood estimation of the parameters of vessel 
# length distributions (based on the mle2 function from package bbmle)
mle_VL <-   function(z, counts, dist, subs = FALSE, size = NULL, start, lower, 
                     control = list(eval.max = 999, iter.max = 999)){
  # load bbmle package
  require(bbmle)
  
  # get function for the shape of the infusion profile p(z)
  p_dist <- eval(parse(text = paste0("p_", dist, "_mu")), envir = .GlobalEnv)

  # get arguments of p_dist
  args   <- paste(formalArgs(p_dist)[-1], collapse = ", ")

  # model fitting with exhaustive sampling estimator
  if (!subs){
    # reshape data for exhaustive sampling
    m <- length(counts)
    data <- list(z1 = z[2:m],
                 z0 = z[1:(m-1)],
                 counts1 = counts[2:m] ,
                 counts0 = counts[1:(m-1)],
                 p_dist = p_dist) # has to be included here due to the weird evaluation
                                  # tweaks used in bbmle
    
    # get model formula
    form <- eval(parse(text = paste0("counts1 ~ dbinom(size = counts0,
                                                       prob = p_dist(z1, ", args, ") /
                                                              p_dist(z0, ", args, "))")))
    cat(deparse(form))
    # fit model
    mod <- bbmle::mle2(minuslogl = form,
                       data = data,
                       optimizer = "nlminb",
                       start = start,
                       lower = lower,
                       control = control)
    
  } else {  # model fitting with subsampling estimator
    # get model formula
    form <-eval(parse(text = paste0("counts ~ dbinom(size = size,
                                                      prob = p_dist(z, ", args, "))")))
    # fit model
    mod <- bbmle::mle2(minuslogl = form,
                       data = list(counts = counts, distances = distances,
                                   size = size, p_dist= p_dist), # for p_dist see above
                       optimizer = "nlminb",
                       start = start,
                       lower = lower,
                       control = control)
  }
  return(mod)
}

# 3. fit_mod() -----------------------------------------------------------------
# function operator that either passes its arguments to mle_VL, 
# or fits simple (non)linear models based on the Cohen or Christman 
# methods
# for the Christman method and the mle_VL models, sometimes the models do not 
# converge with the standard starting values provided by .get_start(). To reduce
# the amount of failed models, these models are re-run 10 times with different 
# starting values

fit_mod <- function(type, z, counts, dist, subs, size){
  # get initial starting values
  start <- .get_start(z = z, counts = counts, dist = dist)
  # set counter
  t <- 1
  # set output value if model fails to NA
  mod <- NA
  # fit cohen model
  if (type == "cohen"){
    try(mod <- lm(log(counts[counts != 0]) ~ z[counts != 0]))
  }
  # fit christman model
  if (type == "christman"){
    freq <- counts / size
    repeat {
      # try to fit model
      try(mod <- nls(freq ~ exp(- (k * z) ^ c), start = start))
      # break if model converged (second if clause checks whether it truly converged)
      if (!is.na(mod)) break
      # advance counter
      t <- t + 1
      # break if model did not converge after 10 trials
      if (t > 10) break 
      # perturb starting values
      start[[1]] <- 1 / runif(1, 0.1, 0.5 * max(z))
      start[[2]] <- runif(1, 0.2, 1.2)
    }
  }
  
  # fit all other models
  else{
    lower <- .get_lower(dist = dist)
    repeat {
      try(mod <- mle_VL(z = z, counts = counts, dist = type, subs  = subs, 
                        size = size, start = start, lower = lower))
      # break if model converged (second if clause checks whether it truly converged)
      if (!is.na(mod)) {if (!is.na(summary(mod)@coef[2]) )  break}
      # advance counter
      t <- t + 1
      # break if model did not converge after 10 trials
      if (t > 10) break
      # perturb starting values
      start[[1]] <- runif(1, 0.1, 0.5 * max(z))
      if (length(start) > 1) start[[2]] <- ifelse(dist == "lnorm", runif(1, 0.2, 1.2), runif(1, 0.6, 6))
    }
  }
  return(list(mod = mod, start = start, t = t))
}

# 4. get_coefs() ---------------------------------------------------------------
# function operator that gets the estimated average VL and other coefficients from each model
get_coefs <- function(type, mod){
  out <- c(mean_est = NA, par2_est = NA, par2_est_name = NA)

  # cohen model
  if (type == "cohen"){
    out[1] <- -2/coef(mod)[2]
  }
  
  # christman model
  if (type == "christman"){
    k <- coef(mod)[["k"]] 
    c <- coef(mod)[["c"]]        
    # minimum vessel length   
    lower <- ifelse(c < 1, 0, (1 / k) * ((c - 1) / c) ^ (1 / c))
    # average vessel length (by numerical integration)
    out[1] <- integrate(function(x) x * Fchristman1(x, c, k), lower, Inf)$value # below minimum length, integral is 0
    } 
  
  # all other models
  else{ 
    coefs  <- coefs(mod)
    out[1] <- coefs[1]
    if (length(coefs) > 1){
      out[2] <- coefs[2]
      out[3] <- names(coefs)[2]
      }
  }
  # return results
  return(out)
}



## important: starting value have to be in a list!!!

### experimental
# 
# with(VLP, attach(list(z = distances, counts = counts, dist = "gamma", 
#                  subs = FALSE, size = NULL, start = list(shape = 1, mu = 30), 
#                  lower = list(shape = 0.001, mu = 0.001))))
# 
# mmm <- with(VLP, mle_VL(z = distances, counts = counts, dist = "erlang2", 
#                         subs = FALSE, size = NULL, start = c(mu = 1), 
#                         lower = c(mu = 0.000001)))
# mmm
# 
# mmm <- with(VLP, mle_VL(z = distances, counts = counts, dist = "weibull", 
#                         subs = FALSE, size = NULL,  start = list(shape = 1, mu = 2), 
#                         lower = list(shape = 0.0001, mu = 0.0001)))
# mmm
# 
# settings <- read.csv("settings/monte_carlo_settings.csv") %>% tbl_df
# 
# (settings1 <- settings[46711, ])
# 
# VLP  <- VL_profile(settings = settings1)
# Vls  <- VL_subs(VLP, 50)
# 


