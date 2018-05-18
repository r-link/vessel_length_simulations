################################################################################
################################################################################
########
########       Functions for the sampling of vessel length distributions
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# Contents ---------------------------------------------------------------------
# 1. rdist()      - function for random samples from the distribution in the
#                   corresponding line of the sampling settings
# 2. VL_profile() - takes the arguments from a line of the sampling settings to
#                   generate a full sample of the vessel length distribution 
#                   based on exhaustive sampling

# 1. rdist() -------------------------------------------------------------------
# this function takes the sampling settings as an argument and produces a random
# sample of size n of the desired distribution with the specified parameter 
# combination
rdist  <- function(n, settings){
  # for the Erlang(2) cases, sampling is done from the rgamma function
  if(settings$distribution == "erlang2") settings$distribution <- "gamma"
  
  # get correct call 
  # (ifelse statement removes second argument for exponential distribution)
  call <- with(settings, 
               paste0("r", distribution, "(", "n = ", n, ", ",
                      par1_name, " = ", par1,
                      ifelse(distribution != "exp", 
                             paste0(", ", par2_name, " = ", par2), 
                             ""), ")"))
  # evaluate call in parent frame
  eval(parse(text = call), envir = parent.frame())
}

# 2. VL_profile() --------------------------------------------------------------

