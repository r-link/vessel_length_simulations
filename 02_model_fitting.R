################################################################################
################################################################################
########
########      Monte Carlo simulations of vessel length measurements 
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# Aim: Monte Carlo simulation of the vessel length estimation process 

# Sections:
# 1. Prerequisites     
# 2. Prepare output files
# 3. Modeling loop

#	1. Prerequisites -------------------------------------------------------------

## load packages
# list of packages 
pkgs <-c("tidyverse", "magrittr", "bbmle", "foreach", "doMC")        

# check for existence of packages and install if necessary
to_install<-pkgs[!(pkgs %in% installed.packages()[, 1])]
if (length(to_install) > 0)  for (pkg in to_install) install.packages(pkg)

# load all required packages 
for (pkg in pkgs) require(pkg, character.only = T)

# list all scripts in the R folder 
scripts <- list.files("R")
# load scripts 
for (i in seq(scripts)) source(file.path("R", scripts[i]))

# load modeling settings
settings <- read.csv("settings/monte_carlo_settings.csv") %>% tbl_df

# keeping track of progress in multicore loops based on the foreach package is
# tricky, as console printouts like cat() and message() do nor work.
# I therefore create a socket to track progress (see following StackOverflow 
# thread for details: https://stackoverflow.com/a/24408477)
# [to make this work, open a UNIX terminal and for instance use netcat 
# with the command nc -l 4000 to connect to the right port before executing
# the next line!]
log.socket <- make.socket(port = 4000)

# the following function is called in the foreach loop
log_fun <- function(text, i) {
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), i)
  cat(msg)
  write.socket(log.socket, msg)
}
# test connection 
log_fun("(%d) <- Does this appear as a number in the console?", i = 100)


#	2. Prepare output files ------------------------------------------------------
# vector for model types
models <-  c("exp", "erlang2","weibull", "gamma", "lnorm","cohen", "christman")

# generate tibble with output columns
out <- tibble(# relevant model output
  model         = models,      # type of fitted model
  converged     = NA,          # logical indicator for model convergence
  mean_est      = NA_real_,    # estimated mean VL
  par1_est      = NA_real_,    # estimated first parameter
  par1_est_name = NA_real_,    # name of first parameter
  par2_est      = NA_real_,    # estimated second parameter (if applicable)
  par2_est_name = NA_real_,    # name of second parameter (if applicable)
  OVL           = NA_real_,    # percent overlap between true and estimated VL distribution
  # additional model documentation
  nvess_true    = NA_real_,    # true number of counted vessels at the cut surface
  meanstart     = NA_real_,    # starting value for mean VL
  par2start     = NA_real_,    # starting value for second parameter (if applicable)
  lmean         = NA_real_,    # lower 95% confidence bound of mean VL (from likelihood profile)
  umean         = NA_real_,    # upper 95% confidence bound of mean VL (from likelihood profile)
  lpar2         = NA_real_,    # lower 95% confidence bound of par2 (from likelihood profile)
  upar2         = NA_real_,    # upper 95% confidence bound of par2 (from likelihood profile)
  lmeanquad     = NA_real_,    # lower 95% confidence bound of mean VL (from quadratic approximation)
  umeanquad     = NA_real_,    # upper 95% confidence bound of mean VL (from quadratic approximation)
  lpar2quad     = NA_real_,    # lower 95% confidence bound of par2 (from quadratic approximation)
  upar2quad     = NA_real_,    # upper 95% confidence bound of par2 (from quadratic approximation)
  t             = NA_real_     # number of tries before convergence
)
out

# prepare empty output object with all column titles
empty <- bind_cols(settings[1, ], out[1, ]) 
empty[1, ] <- NA

# create empty output files
write.csv(empty, "output/output_exhaustive.csv", row.names = FALSE)
write.csv(empty, "output/output_subsamples.csv", row.names = FALSE)


#	3. Modeling loop -----------------------------------------------------------
# The computations are performed in parallel on 4 cores using the multicore 
# features of the parallel package via the foreach::foreach()
# The package DoMC is used as a multicore backend. DoMC uses the fork system
# call, which means the following code will NOT work on Windows!

# set multicore options (computation on 4 cores)
registerDoMC()     # register doMC backend
options(cores = 4) # set number of cores
getDoParWorkers()  # check if number of cores is set to correct value

# outer loop: loop over all modeling settings (done in parallel on 4 cores)
output <- foreach(i = 1:nrow(settings)) %dopar% {
# output <- foreach(i = 401:414) %dopar% {
  # paste number of iteration to the socket created above to track progress 
  log_fun("Iteration No. %d of 50000", i)
  
  # set random seed
  set.seed(settings$seed[i]) # calling set.seed() from within the loop allows for 
                             # reproducible dopar() calls without having to use the
                             # doRNG package
  
  # generate random sample of vessels with the functions from "R/02_VL_samples.R"
  sample <- VL_profile(settings[i, ])
  
  # 3.1) modeling with the exhaustive sampling estimator
  # create temporary object for model results
  temp <- out
  
  # get true number of counted vessels
  temp$nvess_true <- sample$nvess_true
  
  # inner loop: loop over all model types
  for(j in 1:7){
    # get current model
    type <- models[j]
    
    fit   <- NA
    # try fitting model with the fit_mod function from "R/04_VL_models.R"
    fit <- fit_mod(type   = type, 
                   z      = sample$distances, 
                   counts = sample$counts, 
                   subs   = FALSE, 
                   size   = sample$nvess_true  # not available in real data [only used in Christman method]
    )
    
    # extract model object
    mod <- fit$mod
    
    # get starting values, convergence indicator and number of trials
    if (!(type %in% c("cohen", "christman") )) {
      if (length(fit$start) > 0) temp[j, "meanstart"] <- fit$start[[1]]
      if (length(fit$start) > 1) temp[j, "par2start"] <- fit$start[[2]]
    }
    temp[j, "converged"] <- !is.null(mod) 
    if (inherits(mod, "mle2")) temp[j, "converged"] <- !(inherits(mod, "mle2") & is.na(summary(mod)@coef[2]))
    temp[j, "t"]         <- fit$t
    
    # get coefficients, percent overlap, confidence intervals and coverage for converged models
    if (!is.null(mod)){
      # define empty object for coefficients
      coefs <- rep(NA, 5)
      # try to extract coefficients from model object with the get_coefs function from "R/04_VL_models.R"
      try(coefs <- get_coefs(type = type, mod = mod))
      temp[j, "mean_est"]      <- coefs[1]
      temp[j, "par1_est"]      <- coefs[2]
      temp[j, "par1_est_name"] <- coefs[3]
      temp[j, "par2_est"]      <- coefs[4]
      temp[j, "par2_est_name"] <- coefs[5]
      
      # get percent overlap
      try(temp[j, "OVL"] <- percOL(mod = mod, type = type, settings = settings[i, ]))
      
      # get confidence intervals and coverage using the get_CI function from "R/05_VL_accuracy.R"
      if (!(type %in% c("cohen", "christman"))) {
        try(temp[j,c ("lmean", "umean", "lpar2", "upar2",
                      "lmeanquad", "umeanquad", "lpar2quad", "upar2quad")] <- get_CI(mod)
        )
      }
    }
  }
  
  # export results
  res <- bind_cols(settings[rep(i, 7), ], temp)
  write.table(res,  "output/output_exhaustive.csv", sep = ",", dec = ".",
              row.names = F, col.names = F, append = T)
  
  
  # 3.2) modeling with the subsampling estimator
  if (settings[i, "nvess"] >= 500){
    # generate random subsample with the functions from "R/03_VL_subsamples.R"
    subsample <- VL_subs(VL_profile = sample, nsub = settings$nvess[i] / 10)
    
    # create new temporary object for results
    temp1 <- out
    
    # get true number of counted vessels
    temp1$nvess_true <- subsample$nsub
    
    # second inner loop for random subsampling
    for(j in 1:7){
      # get current model
      type <- models[j]
      
      # try fitting model with the fit_mod function from "R/04_VL_models.R"
      fit <- fit_mod(type   = type, 
                     z      = subsample$distances, 
                     counts = subsample$subs, 
                     subs   = TRUE, 
                     size   = subsample$nsub
      )
      
      # extract model object
      mod <- fit$mod
      
      # get starting values, convergence indicator and number of trials
      if (!(type %in% c("cohen", "christman") )) {
        if (length(fit$start) > 0) temp1[j, "meanstart"] <- fit$start[[1]]
        if (length(fit$start) > 1) temp1[j, "par2start"] <- fit$start[[2]]
      }
      temp1[j, "converged"] <- !is.null(mod) 
      if (inherits(mod, "mle2")) temp1[j, "converged"] <- !(inherits(mod, "mle2") & is.na(summary(mod)@coef[2]))
      temp1[j, "t"]         <- fit$t
      
      # get coefficients, percent overlap, confidence intervals and coverage for converged models
      if (!is.null(mod)){
        # define empty object for coefficients
        coefs <- rep(NA, 5)
        # try to extract coefficients from model object with the get_coefs function from "R/04_VL_models.R"
        try(coefs <- get_coefs(type = type, mod = mod))
        temp1[j, "mean_est"]      <- coefs[1]
        temp1[j, "par1_est"]      <- coefs[2]
        temp1[j, "par1_est_name"] <- coefs[3]
        temp1[j, "par2_est"]      <- coefs[4]
        temp1[j, "par2_est_name"] <- coefs[5]
        
        # get percent overlap
        try(temp1[j, "OVL"] <- percOL(mod = mod, type = type, settings = settings[i, ]))
        
        # get confidence intervals and coverage using the get_CI function from "R/05_VL_accuracy.R"
        if (!(type %in% c("cohen", "christman"))) {
          try(temp1[j,c ("lmean", "umean", "lpar2", "upar2",
                        "lmeanquad", "umeanquad", "lpar2quad", "upar2quad")] <- get_CI(mod)
          )
        }
      }
    }
    # export results
    res1 <- bind_cols(settings[rep(i, 7), ], temp1)
    write.table(res1,  "output/output_subsamples.csv", sep = ",", dec = ".",
                row.names = F, col.names = F, append = T)
    
  }
  # return nothing inside the foreach call (I found the runtime for foreach loops that return
  # large objects to grow exponentially with the number of iterations, with the tendency to jam
  # all memory after a while. As long as the issue with objects growing in memory is not fixed,
  # storing the outcome of each iteration on HD seems the better option to me)
  return(NA)
}
