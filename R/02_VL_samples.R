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
# 1. rdist()      - random samples from the vessel length distribution
# 2. qdist()      - quantile function for the vessel length distribution
# 3. spacing()    - function that computes cuts with exponential spacing 
# 4. VL_profile() - constructor function for the "VL_profile" class

# 1. rdist() -------------------------------------------------------------------
# random sample from the true vessel length distribution:
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

# 2. qdist() -------------------------------------------------------------------
# quantile function of the true vessel length distribution:
# this function takes the sampling settings and a target probability as 
# arguments and returns the corresponding quantiles of the desired distribution
# with the specified parameter combination
# combination
qdist  <- function(p, settings){
  # for the Erlang(2) cases, sampling is done from the rgamma function
  if(settings$distribution == "erlang2") settings$distribution <- "gamma"
  
  # get correct call 
  # (ifelse statement removes second argument for exponential distribution)
  call <- with(settings, 
               paste0("q", distribution, "(", "p = ", p, ", ",
                      par1_name, " = ", par1,
                      ifelse(distribution != "exp", 
                             paste0(", ", par2_name, " = ", par2), 
                             ""), ")"))
  # evaluate call in parent frame
  eval(parse(text = call), envir = parent.frame())
}

# 3. spacing() -----------------------------------------------------------------
# function for exponential cut spacing (corrected from Sperry et al., 2005)
spacing <- function(lmin, lmax, n) lmin * (lmax / lmin) ^ ((1:n - 1) / (n - 1))

# 4. VL_profile() --------------------------------------------------------------
# constructor function for S3 objects of the class "VL_profile"
# this function takes the arguments from a line of the sampling settings to
# generate a full sample of the vessel length distribution based on exhaustive 
# sampling
VL_profile <- function(settings){
  # length before injection point (chosen to exclude 1 in a million vessels)
  L <- qdist(p = 0.999999, settings = settings)
  
  # number of vessels to simulate in order to obtain an average of nvess vessels
  # at each cutting plane
  n <- with(settings, ceiling(2 * nvess * L / truemean)) 
  
  # random sample of vessel lengths
  VL_true  <- rdist(n = n, settings = settings)
  
  # left end positions from uniform distribution
  left_end <- runif(n, min = -L, L)
  
  # right end positions
  right_end <- VL_true + left_end
  
  # subset of vessels that are situated at injection point
  VL_biased <- VL_true[left_end <= 0 & right_end >= 0]
  
  # sample of cut-off vessels
  VL_cut <- right_end[left_end <= 0 & right_end >= 0]
  
  # maximum distance for the cut positions from the injection point: 
  # 98% quantile of the observed vessel lengths (some of the parameter 
  # combinations result in extremely right skewed VL distributions, so
  # using the maximum leads to convergence issues)
  lcutoff <- quantile(VL_cut, 0.98)
  
  # cutting positions (exponential spacing from a minimum sample length of 
  # 0.3 cm to the cutoff value)
  distances <- spacing(lmin = 0.3, lmax = lcutoff, n = settings$ncuts)
  
  # counts of filled vessels in function of the distance from the cutting plane
  counts <- sapply(distances, function(x) sum(VL_cut >= x))
  
  # store all created objects in a list
  out <- list(counts    = counts,
              distances = distances,
              VL_true   = VL_true,
              VL_biased = VL_biased,
              VL_cut    = VL_cut,
              left_end  = left_end,
              right_end = right_end,
              L         = L,
              n         = n,
              lcutoff   = lcutoff,
              settings  = settings
              )
  
  # return object with the results
  structure(out, class = "VL_profile")
}
