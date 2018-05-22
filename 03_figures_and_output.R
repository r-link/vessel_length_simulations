################################################################################
################################################################################
########
########      Monte Carlo simulations of vessel length measurements 
########           based on the subsampling estimator
########
########                Author: Roman Link (rlink@gwdg.de)
########
################################################################################
################################################################################

# Aim: Monte Carlo simulation of the vessel length estimation process based on 
#      the assumption that only a subsample of the filled and empty vessels in 
#      each cross-section are counted

# Sections:
# 1. Prerequisites     
# 2. Preparation of sampling settings
# 3. 

#	1. Prerequisites -------------------------------------------------------------

## load packages
# list of packages 
pkgs <-c("tidyverse", "magrittr", "bbmle")        

# check for existence of packages and install if necessary
to_install<-pkgs[!(pkgs %in% installed.packages()[, 1])]
if (length(to_install) > 0)  for (pkg in to_install) install.packages(pkg)

# load all required packages 
for (pkg in pkgs) require(pkg, character.only = T)

# load scripts for model fitting and for extracting results
# list all scripts in the R folder 
# grep is used to remove rkward autosaves (filename contains "~")
scripts <- grep("~", list.files("R"), value = T, invert = T)
for (i in seq(scripts)) source(file.path("R", scripts[i]))


################################################################################
########	Preparation of sampling settings
################################################################################
