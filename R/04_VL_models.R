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
# 1. mle_VL() - maximum likelihood estimation of vessel length distributions

# 1. mle_VL() ------------------------------------------------------------------
# function for the maximum likelihood estimation of the parameters of vessel 
# length distributions (based on the mle2 function from package bbmle)
mle_VL <-   function(z, counts, dist, subs = FALSE, size = NULL, start, lower){
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
                       lower = lower)
    
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
                       lower = lower)
  }
  return(mod)
}

## important: starting value have to be in a list!!!

### experimental

with(VLP, attach(list(z = distances, counts = counts, dist = "gamma", 
                 subs = FALSE, size = NULL, start = list(shape = 1, mu = 30), 
                 lower = list(shape = 0.001, mu = 0.001))))

mmm <- with(VLP, mle_VL(z = distances, counts = counts, dist = "erlang2", 
                        subs = FALSE, size = NULL, start = c(mu = 1), 
                        lower = c(mu = 0.000001)))
mmm

mmm <- with(VLP, mle_VL(z = distances, counts = counts, dist = "weibull", 
                        subs = FALSE, size = NULL,  start = list(shape = 1, mu = 2), 
                        lower = list(shape = 0.0001, mu = 0.0001)))
mmm

settings <- read.csv("settings/monte_carlo_settings.csv") %>% tbl_df

(settings1 <- settings[46711, ])

VLP  <- VL_profile(settings = settings1)
Vls  <- VL_subs(VLP, 50)



